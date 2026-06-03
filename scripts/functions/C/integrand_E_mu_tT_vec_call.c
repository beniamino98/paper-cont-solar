// functions/mu_integrals.c
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <string.h>  // strcmp

/* Utility: get named element from a list (not strictly needed here but handy) */
static SEXP getListElt(SEXP list, const char *name) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (R_len_t i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

/* Copy a 2x2 R matrix (column-major) into row-major [a11 a12; a21 a22] */
static inline void copy_Rmat2x2_rowmajor(SEXP M, double out[4]) {
  const double *r = REAL(M); // r = [a11, a21, a12, a22]
  out[0] = r[0]; out[1] = r[2];
  out[2] = r[1]; out[3] = r[3];
}

/* (1x2) * (2x2) -> (1x2) */
static inline void vec2T_mul_mat2x2(const double v[2], const double A[4], double out[2]) {
  out[0] = v[0]*A[0] + v[1]*A[2];
  out[1] = v[0]*A[1] + v[1]*A[3];
}

/*========================
 1) e_mu_B_t_vec_call (vectorized)
   Args (in order):
     s_vec                   numeric vector (length N)
     t0                      scalar double
     bounds                  list (passed to R closures)
     p0                      numeric length 2
     mu_list                 list of 12 numeric length 2
     Phi_fun                 closure: Phi(t_from, t_to, bounds) -> 2x2 matrix
     get_monthly_index_fun   closure: get_monthly_index(s, bounds) -> integer 1..12

   Returns numeric vector length N: e_mu(s_i)
========================*/
SEXP e_mu_B_t_vec_call(SEXP s_vec,
                            SEXP t0,
                            SEXP bounds,
                            SEXP p0,
                            SEXP mu_list,
                            SEXP Phi_fun,
                            SEXP get_monthly_index_fun) {
  if (!isReal(s_vec)) error("s must be a numeric vector");
  const int N = length(s_vec);
  const double *S = REAL(s_vec);

  if (length(p0) != 2) error("p0 must have length 2");
  const double p0v[2] = { REAL(p0)[0], REAL(p0)[1] };

  if (!isNewList(mu_list)) error("mu must be a list");
  if (!isNewList(bounds))  error("bounds must be a list");

  const double t0d = REAL(t0)[0];

  SEXP out = PROTECT(allocVector(REALSXP, N));
  double *OUT = REAL(out);

  for (int k = 0; k < N; ++k) {
    double s = S[k];

    /* month(s) */
    int m_s;
    {
      SEXP xs = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang3(get_monthly_index_fun, xs, bounds));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      m_s = asInteger(res);
      UNPROTECT(3);
    }
    if (m_s < 1 || m_s > length(mu_list)) error("get_monthly_index out of range");

    /* mu_s (length 2) */
    SEXP mu_s_sexp = VECTOR_ELT(mu_list, m_s - 1);
    if (!isReal(mu_s_sexp) || length(mu_s_sexp) < 2) error("mu[[m]] must be numeric length 2");
    const double mu_s[2] = { REAL(mu_s_sexp)[0], REAL(mu_s_sexp)[1] };

    /* Phi(t0, s, bounds) */
    double P_t0s[4];
    {
      SEXP a = PROTECT(ScalarReal(t0d));
      SEXP b = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang4(Phi_fun, a, b, bounds));
      SEXP M = PROTECT(eval(call, R_GlobalEnv));
      if (!isReal(M) || !isMatrix(M)) error("Phi(t0,s) must return 2x2 real matrix");
      SEXP dm = getAttrib(M, R_DimSymbol);
      if (INTEGER(dm)[0]!=2 || INTEGER(dm)[1]!=2) error("Phi(t0,s): not 2x2");
      copy_Rmat2x2_rowmajor(M, P_t0s);
      UNPROTECT(4);
    }

    /* v = p0' %*% Phi(t0,s)  (1x2) */
    double v[2];
    vec2T_mul_mat2x2(p0v, P_t0s, v);

    /* e_mu(s) = v dot mu_s */
    OUT[k] = v[0]*mu_s[0] + v[1]*mu_s[1];
  }

  UNPROTECT(1);
  return out;
}

/*========================
 2) integrand_E_mu_tT (vectorized)
   Args (in order):
     s_vec, bounds, p0, mu_list, theta, sigma_bar_fun, t0, T_, Phi_fun, get_monthly_index_fun

   Returns numeric vector length N:
     exp(-theta*(T - s)) * sigma_bar(s) * e_mu(s)
========================*/
SEXP integrand_E_mu_tT_vec_call(SEXP s_vec,
                                SEXP bounds,
                                SEXP p0,
                                SEXP mu_list,
                                SEXP theta,
                                SEXP sigma_bar_fun,
                                SEXP t0,
                                SEXP T_,
                                SEXP Phi_fun,
                                SEXP get_monthly_index_fun) {
  if (!isReal(s_vec)) error("s must be a numeric vector");
  const int N = length(s_vec);
  const double *S = REAL(s_vec);

  if (length(p0) != 2) error("p0 must have length 2");
  const double p0v[2] = { REAL(p0)[0], REAL(p0)[1] };

  const double t0d   = REAL(t0)[0];
  const double Td    = REAL(T_)[0];
  const double Theta = REAL(theta)[0];

  if (!isNewList(mu_list)) error("mu must be a list");
  if (!isNewList(bounds))  error("bounds must be a list");

  SEXP out = PROTECT(allocVector(REALSXP, N));
  double *OUT = REAL(out);

  for (int k = 0; k < N; ++k) {
    double s = S[k];

    /* month(s) */
    int m_s;
    {
      SEXP xs = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang3(get_monthly_index_fun, xs, bounds));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      m_s = asInteger(res);
      UNPROTECT(3);
    }
    if (m_s < 1 || m_s > length(mu_list)) error("get_monthly_index out of range");

    /* mu_s */
    SEXP mu_s_sexp = VECTOR_ELT(mu_list, m_s - 1);
    if (!isReal(mu_s_sexp) || length(mu_s_sexp) < 2) error("mu[[m]] must be numeric length 2");
    const double mu_s[2] = { REAL(mu_s_sexp)[0], REAL(mu_s_sexp)[1] };

    /* Phi(t0, s, bounds) */
    double P_t0s[4];
    {
      SEXP a = PROTECT(ScalarReal(t0d));
      SEXP b = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang4(Phi_fun, a, b, bounds));
      SEXP M = PROTECT(eval(call, R_GlobalEnv));
      if (!isReal(M) || !isMatrix(M)) error("Phi(t0,s) must return 2x2 real matrix");
      SEXP dm = getAttrib(M, R_DimSymbol);
      if (INTEGER(dm)[0]!=2 || INTEGER(dm)[1]!=2) error("Phi(t0,s): not 2x2");
      copy_Rmat2x2_rowmajor(M, P_t0s);
      UNPROTECT(4);
    }

    /* v = p0' %*% Phi(t0,s)  (1x2) */
    double v[2];
    vec2T_mul_mat2x2(p0v, P_t0s, v);

    /* e_mu(s) */
    double e_mu = v[0]*mu_s[0] + v[1]*mu_s[1];

    /* sigma_bar(s) */
    double sig_s;
    {
      SEXP xs = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang2(sigma_bar_fun, xs));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      sig_s = asReal(res);
      UNPROTECT(3);
    }

    OUT[k] = exp(-Theta * (Td - s)) * sig_s * e_mu;
  }

  UNPROTECT(1);
  return out;
}

/* registration */
static const R_CallMethodDef CallEntries[] = {
  {"e_mu_B_t_vec_call",     (DL_FUNC) & e_mu_B_t_vec_call,     8},
  {"integrand_E_mu_tT_vec_call", (DL_FUNC) &integrand_E_mu_tT_vec_call, 11},
  {NULL, NULL, 0}
};

void R_init_mu_integrals(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}