#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void check_real_vector(SEXP x, const char *name) {
  if (!isReal(x)) {
    error("%s must be a numeric vector.", name);
  }
}

static void check_real_matrix(SEXP x, const char *name) {
  if (!isReal(x) || !isMatrix(x)) {
    error("%s must be a numeric matrix.", name);
  }
}

static double normal_density(double x, double mean, double sd) {
  const double z = (x - mean) / sd;
  return exp(-0.5 * z * z) / (sd * sqrt(2.0 * M_PI));
}

SEXP ctmc_density_propagate_call(
    SEXP y_grid_sexp,
    SEXP Y0_sexp,
    SEXP p0_sexp,
    SEXP ds_sexp,
    SEXP P_step_sexp,
    SEXP Ybar0_sexp,
    SEXP Ybar1_sexp,
    SEXP sigbar1_sexp,
    SEXP mu_step_sexp,
    SEXP sig_step_sexp,
    SEXP theta_sexp,
    SEXP normalize_sexp) {

  check_real_vector(y_grid_sexp, "y_grid");
  check_real_vector(p0_sexp, "p0");
  check_real_vector(ds_sexp, "ds");
  check_real_matrix(P_step_sexp, "P_step");
  check_real_vector(Ybar0_sexp, "Ybar0");
  check_real_vector(Ybar1_sexp, "Ybar1");
  check_real_vector(sigbar1_sexp, "sigbar1");
  check_real_matrix(mu_step_sexp, "mu_step");
  check_real_matrix(sig_step_sexp, "sig_step");

  if (!isReal(Y0_sexp) || length(Y0_sexp) != 1) {
    error("Y0 must be a numeric scalar.");
  }
  if (!isReal(theta_sexp) || length(theta_sexp) != 1) {
    error("theta must be a numeric scalar.");
  }
  if (!isLogical(normalize_sexp) || length(normalize_sexp) != 1) {
    error("normalize must be a logical scalar.");
  }

  const R_xlen_t N = XLENGTH(y_grid_sexp);
  const R_xlen_t K = XLENGTH(ds_sexp);
  if (N < 2) {
    error("y_grid must have at least two points.");
  }
  if (XLENGTH(p0_sexp) != 2) {
    error("p0 must have length two.");
  }
  if (XLENGTH(Ybar0_sexp) != K || XLENGTH(Ybar1_sexp) != K ||
      XLENGTH(sigbar1_sexp) != K) {
    error("Ybar0, Ybar1, and sigbar1 must have one value per time step.");
  }

  SEXP dim_P = getAttrib(P_step_sexp, R_DimSymbol);
  SEXP dim_mu = getAttrib(mu_step_sexp, R_DimSymbol);
  SEXP dim_sig = getAttrib(sig_step_sexp, R_DimSymbol);
  if (INTEGER(dim_P)[0] != K || INTEGER(dim_P)[1] != 4) {
    error("P_step must be a K x 4 matrix.");
  }
  if (INTEGER(dim_mu)[0] != K || INTEGER(dim_mu)[1] != 2) {
    error("mu_step must be a K x 2 matrix.");
  }
  if (INTEGER(dim_sig)[0] != K || INTEGER(dim_sig)[1] != 2) {
    error("sig_step must be a K x 2 matrix.");
  }

  const double *y = REAL(y_grid_sexp);
  const double Y0 = REAL(Y0_sexp)[0];
  const double *p0 = REAL(p0_sexp);
  const double *ds = REAL(ds_sexp);
  const double *P_step = REAL(P_step_sexp);
  const double *Ybar0 = REAL(Ybar0_sexp);
  const double *Ybar1 = REAL(Ybar1_sexp);
  const double *sigbar1 = REAL(sigbar1_sexp);
  const double *mu_step = REAL(mu_step_sexp);
  const double *sig_step = REAL(sig_step_sexp);
  const double theta = REAL(theta_sexp)[0];
  const int do_normalize = LOGICAL(normalize_sexp)[0] == TRUE;
  const double dy = y[1] - y[0];

  if (!R_FINITE(Y0) || !R_FINITE(theta) || !R_FINITE(dy) || dy <= 0.0) {
    error("Invalid Y0, theta, or y_grid spacing.");
  }
  if (!R_FINITE(p0[0]) || !R_FINITE(p0[1])) {
    error("p0 contains non-finite values.");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, N, 2));
  double *f_out = REAL(out);
  double *f0 = (double *) R_alloc(N, sizeof(double));
  double *f1 = (double *) R_alloc(N, sizeof(double));
  double *new0 = (double *) R_alloc(N, sizeof(double));
  double *new1 = (double *) R_alloc(N, sizeof(double));

  R_xlen_t idx0 = 0;
  double best = fabs(y[0] - Y0);
  for (R_xlen_t j = 1; j < N; j++) {
    const double d = fabs(y[j] - Y0);
    if (d < best) {
      best = d;
      idx0 = j;
    }
  }

  for (R_xlen_t j = 0; j < N; j++) {
    f0[j] = 0.0;
    f1[j] = 0.0;
  }
  f0[idx0] = p0[0] / dy;
  f1[idx0] = p0[1] / dy;

  for (R_xlen_t k = 0; k < K; k++) {
    const double dsk = ds[k];
    if (!R_FINITE(dsk) || dsk <= 0.0) {
      error("All time steps must be positive and finite.");
    }

    const double p11 = P_step[k];
    const double p12 = P_step[k + K];
    const double p21 = P_step[k + 2 * K];
    const double p22 = P_step[k + 3 * K];
    const double ybar0 = Ybar0[k];
    const double ybar1 = Ybar1[k];
    const double sb1 = sigbar1[k];
    const double phi = exp(-theta * dsk);
    const double sqrt_ds = sqrt(dsk);

    if (!R_FINITE(p11) || !R_FINITE(p12) || !R_FINITE(p21) || !R_FINITE(p22) ||
        !R_FINITE(ybar0) || !R_FINITE(ybar1) || !R_FINITE(sb1)) {
      error("Non-finite value in precomputed time-step inputs.");
    }

    for (R_xlen_t l = 0; l < N; l++) {
      new0[l] = 0.0;
      new1[l] = 0.0;
    }

    for (int state = 0; state < 2; state++) {
      const double mu_i = mu_step[k + K * state];
      const double sig_i = sig_step[k + K * state];
      const double sd_next = sb1 * sig_i * sqrt_ds;
      const double drift_shift = sb1 * mu_i * dsk;

      if (!R_FINITE(mu_i) || !R_FINITE(sig_i) ||
          !R_FINITE(sd_next) || sd_next <= 0.0) {
        error("Invalid CTMC density Gaussian kernel standard deviation.");
      }

      for (R_xlen_t j = 0; j < N; j++) {
        const double g_prev = (state == 0) ?
          (p11 * f0[j] + p21 * f1[j]) :
          (p12 * f0[j] + p22 * f1[j]);
        if (g_prev == 0.0) {
          continue;
        }

        const double mean_next = ybar1 + phi * (y[j] - ybar0) + drift_shift;
        const double weight = g_prev * dy;

        if (state == 0) {
          for (R_xlen_t l = 0; l < N; l++) {
            new0[l] += normal_density(y[l], mean_next, sd_next) * weight;
          }
        } else {
          for (R_xlen_t l = 0; l < N; l++) {
            new1[l] += normal_density(y[l], mean_next, sd_next) * weight;
          }
        }
      }
    }

    if (do_normalize) {
      double mass = 0.0;
      for (R_xlen_t l = 0; l < N; l++) {
        mass += new0[l] + new1[l];
      }
      mass *= dy;
      if (R_FINITE(mass) && mass > 0.0) {
        const double inv_mass = 1.0 / mass;
        for (R_xlen_t l = 0; l < N; l++) {
          new0[l] *= inv_mass;
          new1[l] *= inv_mass;
        }
      }
    }

    for (R_xlen_t l = 0; l < N; l++) {
      f0[l] = new0[l];
      f1[l] = new1[l];
    }
  }

  for (R_xlen_t l = 0; l < N; l++) {
    f_out[l] = f0[l];
    f_out[l + N] = f1[l];
  }

  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef CallEntries[] = {
  {"ctmc_density_propagate_call", (DL_FUNC) &ctmc_density_propagate_call, 12},
  {NULL, NULL, 0}
};

void R_init_ctmc_density_kernel(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
