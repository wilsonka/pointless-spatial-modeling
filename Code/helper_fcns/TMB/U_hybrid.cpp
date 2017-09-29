
#include <TMB.hpp>
#include <vector>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR( y_i );  // counts for observation i
  DATA_VECTOR( e_i );  // expected counts for observation i
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(A);
  
  // Parameters
  PARAMETER_VECTOR( theta );
  
  // Random effects
  PARAMETER_VECTOR( S_j );
  PARAMETER( alpha );

  // Objective function
  int n_i = y_i.size();
  vector<Type> jnll_comp(3);
  jnll_comp.setZero(); 
  
  // Probability of thetas
  // Flat prior on thetas
  
  // Probability of intercept
  jnll_comp(2) -= dnorm( alpha, Type(0), Type(10), true);
  
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(2*theta(0)) * ( exp(4*theta(1))*M0 + Type(2.0)*exp(2*theta(1))*M1 + M2);
  jnll_comp(1) += GMRF(Q) ( S_j );  
  
  // Probability of data conditional on random effects
  vector<Type> S_j_star = A * exp(S_j);
  for( int i=0; i<n_i; i++){
  jnll_comp(0) -= dpois( y_i(i), e_i(i) * exp(alpha) * S_j_star(i), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  return jnll;
}

