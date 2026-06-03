// integrand.c
/* --- helpers --- */
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <string.h>  // for strcmp


static SEXP getListElt(SEXP list, const char *name) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (R_len_t i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

/* R 2x2 (column-major) -> C row-major [a11 a12; a21 a22] */
static inline void copy_Rmat2x2_rowmajor(SEXP M, double out[4]) {
  const double *r = REAL(M); // [a11, a21, a12, a22]
  out[0] = r[0]; out[1] = r[2];
  out[2] = r[1]; out[3] = r[3];
}

/* (1x2) * (2x2) -> (1x2) */
static inline void vec2T_mul_mat2x2(const double v[2], const double A[4], double out[2]) {
  out[0] = v[0]*A[0] + v[1]*A[2];
  out[1] = v[0]*A[1] + v[1]*A[3];
}

/* (2x2) * (2x1) -> (2x1) */
static inline void mat2x2_mul_vec2(const double A[4], const double x[2], double out[2]) {
  out[0] = A[0]*x[0] + A[1]*x[1];
  out[1] = A[2]*x[0] + A[3]*x[1];
}

/* elementwise scale (1x2) by diag(d0,d1): v <- v .* d */
static inline void vec2T_scale_diag(double v[2], const double d[2]) {
  v[0] *= d[0];
  v[1] *= d[1];
}

/* .Call entry */
SEXP integrand_V_mu_tT_cond_vec_call(SEXP su,
                                     SEXP bounds,
                                     SEXP p0,
                                     SEXP mu_list,
                                     SEXP theta,
                                     SEXP sigma_bar_fun,
                                     SEXP i_arg,
                                     SEXP p_T,
                                     SEXP Phi_fun,
                                     SEXP get_month_index_fun) {
  if (!isReal(su) || !isMatrix(su)) error("su must be a real 2xN matrix");
  SEXP dim = getAttrib(su, R_DimSymbol);
  if (INTEGER(dim)[0] != 2) error("su must have 2 rows: (s; u)");
  const int N = INTEGER(dim)[1];
  const double *SU = REAL(su);
  
  if (!isNewList(bounds)) error("bounds must be a list");
  SEXP n_vec   = getListElt(bounds, "n");
  SEXP tau_vec = getListElt(bounds, "tau");
  if (n_vec == R_NilValue || tau_vec == R_NilValue)
    error("bounds must contain named elements 'n' and 'tau'");
  const double t0  = REAL(n_vec)[0];
  const double T_  = REAL(tau_vec)[0];
  
  if (length(p0) != 2) error("p0 must have length 2");
  const double p0v[2] = { REAL(p0)[0], REAL(p0)[1] };
  
  if (!isReal(theta) || length(theta) < 1) error("theta must be scalar");
  const double Theta = REAL(theta)[0];
  
  if (!isInteger(i_arg) && !isReal(i_arg)) error("i must be scalar");
  const int    i_val = asInteger(i_arg);
  const double ei[2] = { (double)i_val, 1.0 - (double)i_val };
  
  if (!isReal(p_T) || length(p_T) < 1) error("p_T must be scalar");
  const double pT = REAL(p_T)[0];
  
  SEXP ans = PROTECT(allocVector(REALSXP, N));
  double *OUT = REAL(ans);
  
  for (int k = 0; k < N; ++k) {
    double s = SU[0 + 2*k];
    double u = SU[1 + 2*k];
    if (s > u) { double tmp = s; s = u; u = tmp; }
    
    /* month(s), month(u) */
    int m_s, m_u;
    {
      SEXP xs = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang3(get_month_index_fun, xs, bounds));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      m_s = asInteger(res);
      UNPROTECT(3);
    }
    {
      SEXP xu = PROTECT(ScalarReal(u));
      SEXP call = PROTECT(lang3(get_month_index_fun, xu, bounds));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      m_u = asInteger(res);
      UNPROTECT(3);
    }
    if (m_s < 1 || m_s > length(mu_list) || m_u < 1 || m_u > length(mu_list))
      error("get_month_index out of range");
    
    /* mu_s, mu_u */
    SEXP mu_s_sexp = VECTOR_ELT(mu_list, m_s - 1);
    SEXP mu_u_sexp = VECTOR_ELT(mu_list, m_u - 1);
    if (!isReal(mu_s_sexp) || length(mu_s_sexp) < 2 ||
        !isReal(mu_u_sexp) || length(mu_u_sexp) < 2)
        error("mu_list[[m]] must be numeric length 2");
    const double mu_s[2] = { REAL(mu_s_sexp)[0], REAL(mu_s_sexp)[1] };
    const double mu_u[2] = { REAL(mu_u_sexp)[0], REAL(mu_u_sexp)[1] };
    
    /* Phi(t0,s), Phi(s,u), Phi(u,T) -> row-major copies */
    double P_t0s[4], P_su[4], P_uT[4];
    {
      SEXP a = PROTECT(ScalarReal(t0));
      SEXP b = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang4(Phi_fun, a, b, bounds));
      SEXP M = PROTECT(eval(call, R_GlobalEnv));
      if (!isReal(M) || !isMatrix(M) || INTEGER(getAttrib(M, R_DimSymbol))[0]!=2 ||
          INTEGER(getAttrib(M, R_DimSymbol))[1]!=2)
        error("Phi(t0,s) must be 2x2");
      copy_Rmat2x2_rowmajor(M, P_t0s);
      UNPROTECT(4);
    }
    {
      SEXP a = PROTECT(ScalarReal(s));
      SEXP b = PROTECT(ScalarReal(u));
      SEXP call = PROTECT(lang4(Phi_fun, a, b, bounds));
      SEXP M = PROTECT(eval(call, R_GlobalEnv));
      if (!isReal(M) || !isMatrix(M) || INTEGER(getAttrib(M, R_DimSymbol))[0]!=2 ||
          INTEGER(getAttrib(M, R_DimSymbol))[1]!=2)
        error("Phi(s,u) must be 2x2");
      copy_Rmat2x2_rowmajor(M, P_su);
      UNPROTECT(4);
    }
    {
      SEXP a = PROTECT(ScalarReal(u));
      SEXP b = PROTECT(ScalarReal(T_));
      SEXP call = PROTECT(lang4(Phi_fun, a, b, bounds));
      SEXP M = PROTECT(eval(call, R_GlobalEnv));
      if (!isReal(M) || !isMatrix(M) || INTEGER(getAttrib(M, R_DimSymbol))[0]!=2 ||
          INTEGER(getAttrib(M, R_DimSymbol))[1]!=2)
        error("Phi(u,T) must be 2x2");
      copy_Rmat2x2_rowmajor(M, P_uT);
      UNPROTECT(4);
    }
    
    /* w = Phi(u,T) %*% ei (2x1) */
    double w[2];
    mat2x2_mul_vec2(P_uT, ei, w);
    
    /* v = p0' %*% Phi(t0,s)  -> (1x2) */
    double v[2];
    vec2T_mul_mat2x2(p0v, P_t0s, v);
    
    /* v = v %*% diag(mu_s)  (elementwise scale by mu_s) */
    vec2T_scale_diag(v, mu_s);
    
    /* v = v %*% Phi(s,u) -> (1x2) */
    double vP[2];
    vec2T_mul_mat2x2(v, P_su, vP);
    
    /* vP = vP %*% diag(w)  (elementwise scale by w) */
    vec2T_scale_diag(vP, w);
    
    /* now multiply by mu_u (1x2 · 2x1) */
    double numer = vP[0]*mu_u[0] + vP[1]*mu_u[1];
    double C_su  = numer / pT;
    
    /* sigma_bar(s), sigma_bar(u) */
    double sig_s, sig_u;
    {
      SEXP xs = PROTECT(ScalarReal(s));
      SEXP call = PROTECT(lang2(sigma_bar_fun, xs));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      sig_s = asReal(res);
      UNPROTECT(3);
    }
    {
      SEXP xu = PROTECT(ScalarReal(u));
      SEXP call = PROTECT(lang2(sigma_bar_fun, xu));
      SEXP res  = PROTECT(eval(call, R_GlobalEnv));
      sig_u = asReal(res);
      UNPROTECT(3);
    }
    
    /* kernel */
    double expo = exp(-Theta * (2.0*T_ - s - u));
    OUT[k] = expo * sig_s * sig_u * C_su;
  }
  
  UNPROTECT(1);
  return ans;
}

/* registration */
static const R_CallMethodDef CallEntries[] = {
  {"integrand_V_mu_tT_cond_vec_call", (DL_FUNC) &integrand_V_mu_tT_cond_vec_call, 10},
  {NULL, NULL, 0}
};

void R_init_integrand(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}