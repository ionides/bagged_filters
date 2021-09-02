library(spatPomp)
#' Constrained motion spatPomp simulator
#'
#' Generate a class \sQuote{spatPomp} object representing a \code{U}-dimensional
#' Brownian motion constrained to sum to zero across units, with a nonlinear
#' drift that applies if this constraint is broken.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @param delta_t Process simulations are performed every \code{delta_t} time units
#' whereas observations occur every one time unit
#' @importFrom utils data
#' @return An object of class \sQuote{spatPomp} representing a simulation from a \code{U}-dimensional
#' Constrained motion
#' @examples
#' co <- con(U=4, N=20)
#' # See all the model specifications of the object
#' spy(co)
#' @export

con <- function(U=10,N=20,delta_t=0.2){
# debug
# U=10 ; N=20; delta_t=0.2
  obs_names <- paste0("U",1:U)
  con_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)
  con_unitnames <- unique(con_data[["unit"]])
  con_unitnames_level <- paste("U",sort(as.numeric(stringr::str_remove(con_unitnames, "U"))),sep='')

  con_unit_statenames <- c("X")
  con_statenames <- paste0(con_unit_statenames,1:U)

  con_IVPnames <- paste0(con_statenames,"_0")
  con_RPnames <- c("sigma","tau")
  con_paramnames <- c(con_RPnames,con_IVPnames)

  con_rprocess <- spatPomp_Csnippet("
    double dW[U];
    double dWsum = 0;
    double Xsum = 0;
    int u;

    for (u = 0 ; u < U ; u++) {
      dW[u] = rnorm(0,sigma*sqrt(dt));
      dWsum += dW[u];
      Xsum += X[u];
    }
    for (u = 0 ; u < U ; u++) {
      dW[u] -= dWsum/U;
    }
    for (u = 0 ; u < U ; u++) {
      X[u] += dW[u] + Xsum*dt;
    }
  ", unit_statenames = c("X"))


  con_skel2 <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_vfnames = c("X"),
    code = "
      double Xsum = 0;
      int u;
      for (u = 0 ; u < U ; u++) {
        DX[u] = 0;
      }
    "
  )

  con_skel <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_vfnames = c("X"),
    code = "
      double Xsum = 0;
      int u;  
      for (u = 0 ; u < U ; u++) {
        Xsum += X[u];
      }
      for (u = 0 ; u < U ; u++) {
        DX[u] = Xsum;
      }
    "
  )

  con_rinit <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_ivpnames = c("X"),
    code = "
      for (int u = 0; u < U; u++) {
        X[u]=X_0[u];
      }
    "
  )

  con_dmeasure <- Csnippet("
    const double *X = &X1;
    const double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    lik=0;
    for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau,1);
    if(!give_log) lik = exp(lik) + tol;
  ")

  con_eunit_measure <- Csnippet("
    ey = X;
  ")

  con_munit_measure <- Csnippet("
    M_tau = sqrt(vc);
  ")

  con_vunit_measure <- Csnippet("
    vc = tau*tau;
  ")

  con_rmeasure <- Csnippet("
    const double *X = &X1;
    double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau+tol);
  ")

  con_dunit_measure <- Csnippet("
    //double tol = 1.0e-18;
    lik = dnorm(Y,X,tau,1);
    if(!give_log) lik = exp(lik);
  ")

  con_runit_measure <- Csnippet("
    double tol = pow(1.0e-18,U);
    double Y;
    Y = rnorm(X,tau+tol);
  ")

  con_spatPomp <- spatPomp(con_data %>% dplyr::arrange(time, factor(.data$unit, levels = con_unitnames_level)),
                 times="time",
                 t0=0,
                 units="unit",
                 unit_statenames = con_unit_statenames,
                 rprocess=euler(con_rprocess,delta.t = delta_t),
                 skeleton=vectorfield(con_skel),
                 paramnames=con_paramnames,
                 rmeasure=con_rmeasure,
                 dmeasure=con_dmeasure,
                 eunit_measure=con_eunit_measure,
                 munit_measure=con_munit_measure,
                 vunit_measure=con_vunit_measure,
                 dunit_measure=con_dunit_measure,
                 runit_measure=con_runit_measure,
                 partrans = parameter_trans(log = c("sigma", "tau")),
                 rinit=con_rinit
    )

  ## We need a parameter vector. For now, we initialize the process at zero.
  test_ivps <- rep(0,U)
  names(test_ivps) <- con_IVPnames
  test_params <- c(sigma=1, tau=1, test_ivps)
  simulate(con_spatPomp,params=test_params)
}

