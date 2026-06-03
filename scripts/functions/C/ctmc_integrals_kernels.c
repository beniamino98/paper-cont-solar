#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <string.h>

static SEXP getListElt(SEXP list, const char *name) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) {
    return R_NilValue;
  }
  for (R_len_t i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

static void copy_Rmat2x2_rowmajor(SEXP M, double out[4]) {
  const double *r = REAL(M);
  out[0] = r[0];
  out[1] = r[2];
  out[2] = r[1];
  out[3] = r[3];
}

static void write_rowmajor_to_Rmat2x2(const double in[4], SEXP M) {
  double *r = REAL(M);
  r[0] = in[0];
  r[1] = in[2];
  r[2] = in[1];
  r[3] = in[3];
}

static void mat2x2_set_I(double A[4]) {
  A[0] = 1.0;
  A[1] = 0.0;
  A[2] = 0.0;
  A[3] = 1.0;
}

static void mat2x2_mul(double A[4], const double B[4]) {
  double c00 = A[0] * B[0] + A[1] * B[2];
  double c01 = A[0] * B[1] + A[1] * B[3];
  double c10 = A[2] * B[0] + A[3] * B[2];
  double c11 = A[2] * B[1] + A[3] * B[3];
  A[0] = c00;
  A[1] = c01;
  A[2] = c10;
  A[3] = c11;
}

static int is_finite_2x2(const double A[4]) {
  return R_FINITE(A[0]) && R_FINITE(A[1]) && R_FINITE(A[2]) && R_FINITE(A[3]);
}

static void validate_generator_2state(const double Q[4]) {
  const double tol = 1e-8;
  if (!is_finite_2x2(Q)) {
    error("Q contains non-finite values.");
  }
  if (fabs(Q[0] + Q[1]) > tol || fabs(Q[2] + Q[3]) > tol) {
    error("Rows of Q must sum to 0.");
  }
  if (Q[1] < -tol || Q[2] < -tol) {
    error("Off-diagonal entries of Q must be non-negative.");
  }
}

static void exp_kQ_2state(double k, const double Q[4], double out[4]) {
  if (!R_FINITE(k) || k < -1e-12) {
    error("k must be a finite non-negative scalar horizon.");
  }

  k = fmax(k, 0.0);
  validate_generator_2state(Q);

  const double q12 = fmax(Q[1], 0.0);
  const double q21 = fmax(Q[2], 0.0);
  const double lambda = q12 + q21;

  if (lambda < 1e-14) {
    out[0] = 1.0 + k * Q[0];
    out[1] =       k * Q[1];
    out[2] =       k * Q[2];
    out[3] = 1.0 + k * Q[3];
    return;
  }

  const double scale = (1.0 - exp(-lambda * k)) / lambda;
  out[0] = 1.0 + scale * Q[0];
  out[1] =       scale * Q[1];
  out[2] =       scale * Q[2];
  out[3] = 1.0 + scale * Q[3];
}

static void check_matrix_2x2(SEXP M, const char *name) {
  if (!isReal(M) || !isMatrix(M)) {
    error("%s must be a numeric 2x2 matrix.", name);
  }
  SEXP dim = getAttrib(M, R_DimSymbol);
  if (dim == R_NilValue || INTEGER(dim)[0] != 2 || INTEGER(dim)[1] != 2) {
    error("%s must be a numeric 2x2 matrix.", name);
  }
}

SEXP ctmc_matrix_exponential_call(SEXP k_sexp, SEXP Q_sexp) {
  if (!isReal(k_sexp) || length(k_sexp) != 1) {
    error("k must be a numeric scalar.");
  }
  check_matrix_2x2(Q_sexp, "Q");

  double Q[4];
  double P[4];
  copy_Rmat2x2_rowmajor(Q_sexp, Q);
  exp_kQ_2state(REAL(k_sexp)[0], Q, P);

  SEXP out = PROTECT(allocMatrix(REALSXP, 2, 2));
  write_rowmajor_to_Rmat2x2(P, out);
  UNPROTECT(1);
  return out;
}

SEXP ctmc_get_month_index_call(SEXP tau_sexp, SEXP bounds) {
  if (!isReal(tau_sexp)) {
    error("tau must be a numeric vector.");
  }
  if (!isNewList(bounds)) {
    error("bounds must be a data.frame/list.");
  }

  SEXP n_raw = getListElt(bounds, "n");
  SEXP N_raw = getListElt(bounds, "N");
  SEXP month_raw = getListElt(bounds, "Month");
  if (n_raw == R_NilValue || N_raw == R_NilValue || month_raw == R_NilValue) {
    error("bounds must have columns 'n', 'N', and 'Month'.");
  }

  SEXP n_col = PROTECT(coerceVector(n_raw, REALSXP));
  SEXP N_col = PROTECT(coerceVector(N_raw, REALSXP));
  SEXP month_col = PROTECT(coerceVector(month_raw, INTSXP));

  const int M = length(n_col);
  const int n_tau = length(tau_sexp);
  const double *n = REAL(n_col);
  const double *N = REAL(N_col);
  const int *month = INTEGER(month_col);
  const double *tau = REAL(tau_sexp);

  SEXP out = PROTECT(allocVector(REALSXP, n_tau));
  double *res = REAL(out);

  for (int k = 0; k < n_tau; k++) {
    int pos = -1;
    for (int j = 0; j < M; j++) {
      if (n[j] <= tau[k] && tau[k] < N[j]) {
        pos = j;
      }
    }
    if (pos >= 0) {
      res[k] = (double) month[pos];
    } else if (tau[k] >= N[M - 1]) {
      res[k] = (double) month[M - 1];
    } else if (tau[k] <= n[0]) {
      res[k] = (double) month[0];
    } else {
      error("tau is outside the supplied CTMC bounds.");
    }
  }

  UNPROTECT(4);
  return out;
}

SEXP ctmc_phi_call(SEXP a_sexp, SEXP b_sexp, SEXP bounds) {
  if (!isReal(a_sexp) || !isReal(b_sexp)) {
    error("a and b must be numeric vectors.");
  }
  if (length(a_sexp) != length(b_sexp)) {
    error("a and b must have the same length after R-side recycling.");
  }
  if (!isNewList(bounds)) {
    error("bounds must be a data.frame/list.");
  }

  SEXP n_raw = getListElt(bounds, "n");
  SEXP N_raw = getListElt(bounds, "N");
  SEXP Q_list = getListElt(bounds, "Q");
  if (n_raw == R_NilValue || N_raw == R_NilValue || Q_list == R_NilValue) {
    error("bounds must have columns 'n', 'N', and 'Q'.");
  }

  SEXP n_col = PROTECT(coerceVector(n_raw, REALSXP));
  SEXP N_col = PROTECT(coerceVector(N_raw, REALSXP));

  const int M = length(n_col);
  const int n_pairs = length(a_sexp);
  const double *n = REAL(n_col);
  const double *N = REAL(N_col);
  const double *a = REAL(a_sexp);
  const double *b = REAL(b_sexp);

  SEXP out = PROTECT(allocVector(VECSXP, n_pairs));

  for (int k = 0; k < n_pairs; k++) {
    if (!R_FINITE(a[k]) || !R_FINITE(b[k])) {
      error("a and b must be finite.");
    }
    if (b[k] < a[k] - 1e-12) {
      error("Phi_C expects a <= b.");
    }

    SEXP P_out = PROTECT(allocMatrix(REALSXP, 2, 2));
    double P[4];
    mat2x2_set_I(P);

    if (fabs(b[k] - a[k]) > 1e-14) {
      int has_overlap = 0;
      for (int j = 0; j < M; j++) {
        if (n[j] < b[k] && N[j] > a[k]) {
          const double left = fmax(a[k], n[j]);
          const double right = fmin(b[k], N[j]);
          if (right - left > 1e-14) {
            SEXP Q_sexp = VECTOR_ELT(Q_list, j);
            check_matrix_2x2(Q_sexp, "Q[[j]]");
            double Q[4];
            double E[4];
            copy_Rmat2x2_rowmajor(Q_sexp, Q);
            exp_kQ_2state(right - left, Q, E);
            mat2x2_mul(P, E);
            has_overlap = 1;
          }
        }
      }
      if (!has_overlap) {
        error("No CTMC interval overlaps [a, b).");
      }
    }

    write_rowmajor_to_Rmat2x2(P, P_out);
    SET_VECTOR_ELT(out, k, P_out);
    UNPROTECT(1);
  }

  UNPROTECT(3);
  return out;
}

static const R_CallMethodDef CallEntries[] = {
  {"ctmc_matrix_exponential_call", (DL_FUNC) &ctmc_matrix_exponential_call, 2},
  {"ctmc_get_month_index_call",    (DL_FUNC) &ctmc_get_month_index_call,    2},
  {"ctmc_phi_call",                (DL_FUNC) &ctmc_phi_call,                3},
  {NULL, NULL, 0}
};

void R_init_ctmc_integrals_kernels(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
