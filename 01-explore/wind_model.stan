data {
  int<lower=1> n_sites;               // Total number of sites
  int<lower=0, upper=1> observed_detections[n_sites]; // Observed detections
  matrix[n_sites, n_sites] wind_effect;  // Wind kernel matrix
}

parameters {
  vector<lower=0, upper=1>[n_sites] z;  // Latent true presence/absence
  real<lower=0, upper=1> p_true_pos;    // True positive detection rate
  real<lower=0, upper=1> p_false_pos;   // False positive detection rate
}

model {
  // Priors
  p_true_pos ~ beta(2, 2);   // Prior for true positive rate
  p_false_pos ~ beta(2, 2);  // Prior for false positive rate
  z ~ uniform(0, 1);         // Weak prior for true presence

  // Observation model for all sites
  for (i in 1:n_sites) {
    real detection_prob = p_true_pos * z[i] + 
                          p_false_pos * (1 - z[i]) + 
                          dot_product(wind_effect[i, ], z);  // Wind effect term

    // Observed detection model
    observed_detections[i] ~ bernoulli_logit(detection_prob);
  }
}
