data {
    int<lower=1> n; // number of instances
    int<lower=2> m; // number of algorithms
     
    // Matrix with all the rankings, one per row
    int ranks [n,m];
    
    real weights[n];
    
    // Parameters for Dirichlet prior.
    vector[m] alpha; 
}


transformed data {
  // The implementation of the probability of the PL model uses the order, rather
  // than the rank
  int order [n,m];
  for (s in 1:n){
    for (i in 1:m){
      order[s, ranks[s, i]]=i;
    }
  }
}

parameters {
    // Vector of ratings for each team.
    // The simplex constrains the ratings to sum to 1 
    simplex[m] ratings;
}
 
transformed parameters{
  real loglik;
  real rest;
  
  loglik=0;
  for (s in 1:n){
    for (i in 1:(m-1)){
      rest=0;
      for (j in i:m){
        rest = rest + ratings[order[s, j]];
      }
      loglik = loglik + log(weights[s] * ratings[order[s, i]] / rest);
    }
  }
}
 
model {
    ratings ~ dirichlet(alpha); // Dirichlet prior
    target += loglik;
}
