// functions/bounds_kernels.c
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <string.h>

/* -------- helpers -------- */

static SEXP getListElt(SEXP list, const char *name) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (R_len_t i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(list, i);
  }
  return R_NilValue;
}

/* Copy 2x2 R (column-major) -> row-major [a11 a12; a21 a22] */
static inline void copy_Rmat2x2_rowmajor(SEXP M, double out[4]) {
  const double *r = REAL(M); /* [a11, a21, a12, a22] */
out[0] = r[0]; out[1] = r[2];
out[2] = r[1]; out[3] = r[3];
}

static inline void mat2x2_mul(double A[4], const double B[4]) {
  double c00 = A[0]*B[0] + A[1]*B[2];
  double c01 = A[0]*B[1] + A[1]*B[3];
  double c10 = A[2]*B[0] + A[3]*B[2];
  double c11 = A[2]*B[1] + A[3]*B[3];
  A[0]=c00; A[1]=c01; A[2]=c10; A[3]=c11;
}

static inline void mat2x2_set_I(double A[4]) {
  A[0]=1.0; A[1]=0.0; A[2]=0.0; A[3]=1.0;
}

/* exp(k Q) with lambda = Q12 + Q21: I + (1 - exp(lambda k))*(Q/lambda)
 Fallback to I + k Q when |lambda| is tiny */
static inline void exp_kQ(double k, const double Q[4], double out[4]) {
  const double lam = Q[1] + Q[2];  /* row-major: Q12=Q[1], Q21=Q[2] */
if (fabs(lam) < 1e-14) {
  /* I + k Q */
  out[0] = 1.0 + k*Q[0];
  out[1] =       k*Q[1];
  out[2] =       k*Q[2];
  out[3] = 1.0 + k*Q[3];
  return;
}
const double s = (1.0 - exp(lam * k)) / lam;
out[0] = 1.0 + s*Q[0];
out[1] =       s*Q[1];
out[2] =       s*Q[2];
out[3] = 1.0 + s*Q[3];
}

/* ---------- 1) get_monthly_index (vectorized) ---------- */
/* tau_vec: numeric, bounds: tibble/list with n_idx (double) and Month (int/double) */
SEXP get_monthly_index_vec_call(SEXP tau_vec, SEXP bounds) {
  if (!isReal(tau_vec)) error("tau must be a numeric vector");
  if (!isNewList(bounds)) error("bounds must be a list/data.frame");
  
  SEXP n_idx = getListElt(bounds, "n_idx");
  SEXP Month = getListElt(bounds, "Month");
  if (n_idx == R_NilValue || Month == R_NilValue)
    error("bounds must have columns 'n_idx' and 'Month'");
  
  const int M = length(n_idx);
  const double *NIDX = REAL(n_idx);
  const double *MONTHd = REAL(coerceVector(Month, REALSXP)); /* Month can be integer or double */

const int N = length(tau_vec);
const double *TAU = REAL(tau_vec);

SEXP out = PROTECT(allocVector(INTSXP, N));
int *OUT = INTEGER(out);

for (int k = 0; k < N; ++k) {
  double tau = TAU[k];
  int pos = -1;
  /* max i such that n_idx[i] <= tau */
  for (int i = 0; i < M; ++i) {
    if (NIDX[i] <= tau) pos = i;
    else break; /* n_idx should be increasing; early exit */
  }
  if (pos < 0) { OUT[k] = (int)MONTHd[0]; }            /* before first row → first Month */
  else          { OUT[k] = (int)MONTHd[pos]; }
}

UNPROTECT(1);
return out;
}

/* ---------- 2) Phi (vectorized over (a,b)) ---------- */
/* a_vec, b_vec numeric (same length). Returns a LIST of 2x2 double matrices. */
SEXP Phi_vec_call(SEXP a_vec, SEXP b_vec, SEXP bounds) {
  if (!isReal(a_vec) || !isReal(b_vec)) error("a and b must be numeric vectors");
  const int Npairs = length(a_vec);
  if (length(b_vec) != Npairs) error("a and b must have same length");
  
  if (!isNewList(bounds)) error("bounds must be a list/data.frame");
  
  SEXP n_col   = getListElt(bounds, "n");
  SEXP N_col   = getListElt(bounds, "N");
  SEXP Q_list  = getListElt(bounds, "Q");
  SEXP Qp_list = getListElt(bounds, "Q_prod");
  
  if (n_col == R_NilValue || N_col == R_NilValue || Q_list == R_NilValue || Qp_list == R_NilValue)
    error("bounds must have columns 'n', 'N', 'Q', 'Q_prod'");
  
  const int M = length(n_col);
  const double *n_v = REAL(n_col);
  const double *N_v = REAL(N_col);
  
  SEXP out = PROTECT(allocVector(VECSXP, Npairs));
  
  for (int k = 0; k < Npairs; ++k) {
    const double a = REAL(a_vec)[k];
    const double b = REAL(b_vec)[k];
    
    /* find rows where N >= a && n <= b */
    int i_first = -1, i_last = -1;
    for (int i = 0; i < M; ++i) {
      if (N_v[i] >= a && n_v[i] <= b) {
        if (i_first < 0) i_first = i;
        i_last = i;
      } else if (N_v[i] > b) {
        break; /* since N likely increasing */
      }
    }
    
    /* allocate result 2x2 (column-major for R) */
    SEXP Mout = PROTECT(allocMatrix(REALSXP, 2, 2));
    double *Rmat = REAL(Mout);
    
    if (i_first < 0) {
      /* no overlap → identity */
      Rmat[0]=1.0; Rmat[1]=0.0; Rmat[2]=0.0; Rmat[3]=1.0;
      SET_VECTOR_ELT(out, k, Mout);
      UNPROTECT(1);
      continue;
    }
    
    double Prod[4];
    mat2x2_set_I(Prod);
    
    if (i_first == i_last) {
      /* single row: exp((b-a) * Q[first]) */
      SEXP Qmat = VECTOR_ELT(Q_list, i_first);
      if (!isReal(Qmat) || !isMatrix(Qmat)) error("Q[[i]] must be 2x2 real");
      SEXP dm = getAttrib(Qmat, R_DimSymbol);
      if (INTEGER(dm)[0] != 2 || INTEGER(dm)[1] != 2) error("Q[[i]] must be 2x2");
      double Qrm[4]; copy_Rmat2x2_rowmajor(Qmat, Qrm);
      
      double Ek[4];
      exp_kQ(b - a, Qrm, Ek);
      /* copy Ek to Prod */
      Prod[0]=Ek[0]; Prod[1]=Ek[1]; Prod[2]=Ek[2]; Prod[3]=Ek[3];
    } else {
      /* first partial */
      {
        SEXP Qmat = VECTOR_ELT(Q_list, i_first);
        if (!isReal(Qmat) || !isMatrix(Qmat)) error("Q[[i_first]] must be 2x2");
        SEXP dm = getAttrib(Qmat, R_DimSymbol);
        if (INTEGER(dm)[0] != 2 || INTEGER(dm)[1] != 2) error("Q[[i_first]] must be 2x2");
        double Qrm[4]; copy_Rmat2x2_rowmajor(Qmat, Qrm);
        double Ek[4];  exp_kQ(N_v[i_first] - a, Qrm, Ek);
        Prod[0]=Ek[0]; Prod[1]=Ek[1]; Prod[2]=Ek[2]; Prod[3]=Ek[3];
      }
      /* middle full blocks: Q_prod[[j]] for j=i_first+1..i_last-1 */
      for (int j = i_first + 1; j <= i_last - 1; ++j) {
        SEXP Pj = VECTOR_ELT(Qp_list, j);
        if (!isReal(Pj) || !isMatrix(Pj)) error("Q_prod[[j]] must be 2x2 real");
        SEXP dmj = getAttrib(Pj, R_DimSymbol);
        if (INTEGER(dmj)[0]!=2 || INTEGER(dmj)[1]!=2) error("Q_prod[[j]] must be 2x2");
        double Pjrm[4]; copy_Rmat2x2_rowmajor(Pj, Pjrm);
        mat2x2_mul(Prod, Pjrm);
      }
      /* last partial */
      {
        SEXP QmatL = VECTOR_ELT(Q_list, i_last);
        if (!isReal(QmatL) || !isMatrix(QmatL)) error("Q[[i_last]] must be 2x2");
        SEXP dml = getAttrib(QmatL, R_DimSymbol);
        if (INTEGER(dml)[0] != 2 || INTEGER(dml)[1] != 2) error("Q[[i_last]] must be 2x2");
        double Qlrm[4]; copy_Rmat2x2_rowmajor(QmatL, Qlrm);
        double EkL[4];  exp_kQ(b - N_v[i_last - 1], Qlrm, EkL);
        mat2x2_mul(Prod, EkL);
      }
    }
    
    /* write row-major Prod back to R (column-major) */
    Rmat[0] = Prod[0]; Rmat[1] = Prod[2];
    Rmat[2] = Prod[1]; Rmat[3] = Prod[3];
    
    SET_VECTOR_ELT(out, k, Mout);
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return out;
}

/* ---- registration ---- */
static const R_CallMethodDef CallEntries[] = {
  {"get_monthly_index_vec_call", (DL_FUNC) &get_monthly_index_vec_call, 2},
  {"Phi_vec_call",               (DL_FUNC) &Phi_vec_call,               3},
  {NULL, NULL, 0}
};

void R_init_bounds_kernels(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}