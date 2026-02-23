// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> N;
  real<lower=0,upper=1> p_flu;
  real<lower=0,upper=1> p_fever;
  real<lower=0,upper=1> p_fever_flu;
}

parameters {
  real<lower=0,upper=1> p_flu_fever;
  real<lower=0> sigma;
}

model {
  p_flu_fever ~ normal((p_fever_flu * p_flu)/(p_fever_flu), sigma) T[0,1];
}
