##' Joint estimation of initial incidence and reproduction number
##'
##'
##' @title Jointly estimate I0 and Rt
##' @param window length of the window used for estimation. If you
##' have the full incidence time-series, only data in this window are
##' used for estimatiom.
##'
##' @param window_back length of time for which incidence should be
##' estimated.
##' @param incid Incidence used to fit the model. Numeric vector of
##' length \code{window}
##' @param si_distr Discretised serial interval distribution
##' @param ... Arguments passed to `rstan::fit` (e.g. iter, chains).
##' @return an object of class stanfit
##' @author Sangeeta Bhatia
##' @export
jointly_estimate <- function(window, window_back, incid, si_distr, ...) {

  omega <- rev(si_distr)
  log_incid_init <- log(mean(incid))
  standata <- list(
    window = window, window_back = window_back, incid = incid,
    omega = omega, log_incid_init = log_incid_init,
    si_trunc = length(omega) - 1
  )
  fit <- rstan::sampling(stanmodels$rti0_bayesian, data = standata, ...)
  fit
}
