functions {
  real force_of_infection_obs(vector incid_init, vector incid_obs,
                              vector omega, real rt, int back,
                              int si_trunc) {

    int T = rows(incid_init) + rows(incid_obs);
    vector[T] incid_full;
    real lambda;
    // At this point you have the incidence from days 1 through to
    // `back` which here is 100.
    // Stack this on top of the observed deaths to get a time series
    // from 1 to 110.
    incid_full[1:back] = incid_init;
    incid_full[(back + 1):rows(incid_full)] = incid_obs;
    lambda = dot_product(incid_full[(T - si_trunc):T], omega);
    return(lambda);
  }

  vector force_of_infection_unobs(real log_i0, vector omega, real rt,
                                  int back, int si_trunc) {
    
    vector[back] incid_init;
    int tmp;
    incid_init = rep_vector(0, back);
    incid_init[1] = exp(log_i0);
    for (index in 2:back) {
      if ((index - si_trunc) > 1) tmp = index - si_trunc;
      else tmp = 1;
      incid_init[index] = rt * dot_product(incid_init[tmp:index],
                                           omega[(si_trunc + 1 - (index - tmp)):(si_trunc + 1)]);
    }
    return(incid_init);
  }
}
data {
  int window; // The window size used for model calibration
  // Size of the window over which to estimate incidence
  int window_back;
  // Incidence in the window. This will be integer of course but here
  // we use a real data type because estimated incidence is on log-scale
  // and then exponentiated. So that it is not guranteed to be an integer.
  vector[window] incid;
  int si_trunc;
  // flipped discretised SI distribution, +1 because of probability
  // mass at day 0
  vector[si_trunc + 1] omega; 
  //set to log of the mean of the incidence in the calibration window
  real log_incid_init;

}
parameters {
  real <lower = log_incid_init - log(10) * 100/6, upper = log_incid_init - log(0.5) * 100/6> log_i0;
  real <lower = 0.1, upper = 10> rt_est;
}
model {
  vector[window_back] incid_est;
  real lambda;
  incid_est = force_of_infection_unobs(log_i0, omega, rt_est,
                                       window_back, si_trunc);
  // Get the daily force of infection
  for (index in 1:window) {
    lambda = force_of_infection_obs(incid_est, incid[1:index], 
                                    omega,  rt_est, window_back, si_trunc);
    target += -rt_est * lambda + incid[index] * log(rt_est * lambda);
    //target += 1;
  }

}
generated quantities {
  vector[window_back + window + 7] incid_est;
  // Incidence prior to the window
  incid_est = force_of_infection_unobs(log_i0, omega, rt_est,
                                       window_back + window + 7,
                                       si_trunc);

}
