// State-Space Leslie Model

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_VECTOR(CPUE);
  DATA_VECTOR(CAT);
  DATA_VECTOR(WEEK);  
  DATA_IVECTOR(YEAR);
  DATA_IVECTOR(START);
  DATA_IVECTOR(END);
  DATA_INTEGER(Y);
  DATA_SCALAR(M);
  DATA_SCALAR(x_lo);
  DATA_SCALAR(x_up);
  DATA_VECTOR(nodes);
  DATA_VECTOR(wt);
  DATA_SCALAR(g);
  DATA_INTEGER(model);     // 0: HS, 1: MR, 2: BHS, 3: BH, 4: RI
  
  // PARAMETER //
  PARAMETER(log_a);
  PARAMETER(log_b);
  PARAMETER(tilde_n0);
  PARAMETER(log_tilde_q);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER(log_eta);
  PARAMETER_VECTOR(log_q);  
  PARAMETER_VECTOR(n0);
  
  // PARAMETER TRANSFORMATION //
  Type a = exp(log_a);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);
  Type eta = exp(log_eta);
  int T = CPUE.size();
  int ng = wt.size();
  
  vector<Type> q(Y);
  vector<Type> n_last(Y);
  vector<Type> N_S(Y);
  vector<Type> N_E(Y);
  vector<Type> pred_n0(Y);
  vector<Type> U(T);
  vector<Type> n(T);
  vector<Type> N(T);
  vector<Type> pred_n(T);
     
  Type nll = 0.0;
  
  // GL integration
  
  Type xm_x = Type(0.5)*(x_up+x_lo);
  Type xr_x = Type(0.5)*(x_up-x_lo);
  vector<Type> dx_gi_x = xr_x*nodes;
  vector<Type> xp_x = xm_x+dx_gi_x;
  
  // DeLury Model
  
  for (int i=0;i<T;i++){
    if (START(i)==1) {
      if (i==0) nll += -dnorm(log_q(YEAR(i)),log_tilde_q,eta,true); 
      if (i > 0) nll += -dnorm(log_q(YEAR(i)),log_q(YEAR(i)-1),eta,true);
      n(i) = n0(YEAR(i))-M*Type(7.0)*WEEK(i);
      N(i) = exp(n(i));
      N_S(YEAR(i)) = exp(n0(YEAR(i)));
      U(i) = CAT(i)/N(i);
    }
    if (START(i)==0) {
      n(i) = n(i-1)-M*Type(7.0)*(WEEK(i)-WEEK(i-1))+log(Type(1.0)-U(i-1));
      N(i) = exp(n(i));
      U(i) = CAT(i)/N(i);
    }
    U(i) = CppAD::CondExpLe(Type(1.0)-U(i),Type(0.0),Type(0.99),U(i));
    if (END(i)==1){
      n_last(YEAR(i)) = n(i)-M*Type(7.0)*(Type(26.0)-WEEK(i))+log(Type(1.0)-U(i));
      N_E(YEAR(i)) = exp(n_last(YEAR(i)));
    }

    pred_n(i) = log_q(YEAR(i))+n(i);
    nll += -dnorm(log(CPUE(i)),pred_n(i),tau,true);
  }
  
  // Density-Dependent Model
  
  vector<Type> S = exp(n_last);
  Type max_S = max(S);
  Type min_S = min(S);
  Type pred;
  
  Type b = CppAD::CondExpLe(exp(log_b), min_S, min_S, exp(log_b));
  b = CppAD::CondExpGe(exp(log_b), max_S, max_S, exp(log_b));
  
  nll += -dnorm(n0(0),tilde_n0,sigma,true);
  for (int i=1;i<Y;i++){ 
    if (model==0) pred_n0(i) = CppAD::CondExpGe(S(i-1), b, log_a + log(b), log_a+log(S(i-1)));
    if (model==1) pred_n0(i) = log_a+log(S(i-1)+sqrt(pow(b,Type(2.0))+pow(g,Type(2.0))/Type(4.0))-sqrt(pow(S(i-1)-b,Type(2.0))+pow(g,Type(2.0))/Type(4.0)));
    if (model==2) {
      pred = 0.0;
      for (int j=0;j<ng;j++){
        pred += xr_x*wt(j)*a*b*pow(S(i-1)/b,1-pow(S(i-1)/b,1/xp_x(j)));
      }
      pred /= 2.0*xr_x;
      pred_n0(i) = CppAD::CondExpGe(S(i-1), b, log_a + log(b), log(pred));
    }
    if (model==3) pred_n0(i) = log_a + log(S(i-1)) - log(Type(1.0)+S(i-1)/b);
    if (model==4) pred_n0(i) = log_a + log(S(i-1)) - S(i-1)/b;
    nll += -dnorm(n0(i), pred_n0(i), sigma, true);
  }
  
  q = exp(log_q);
  
  Type n0_new = 0.0;
  if (model==0) n0_new = CppAD::CondExpGe(S(Y-1), b, log_a + log(b), log_a+log(S(Y-1)));
  if (model==1) n0_new = log_a+log(S(Y-1)+sqrt(pow(b,Type(2.0))+pow(g,Type(2.0))/Type(4.0))-sqrt(pow(S(Y-1)-b,Type(2.0))+pow(g,Type(2.0))/Type(4.0)));
  if (model==2) {
    pred=0.0;
    for (int j=0;j<ng;j++){
        pred += xr_x*wt(j)*a*b*pow(S(Y-1)/b,1-pow(S(Y-1)/b,1/xp_x(j)));
    }
    pred /= 2.0*xr_x;
    n0_new = CppAD::CondExpGe(S(Y-1), b, log_a + log(b), log(pred));
  }
  if (model==3) n0_new = log_a + log(S(Y-1)) - log(Type(1.0)+S(Y-1)/b);
  if (model==4) n0_new = log_a + log(S(Y-1)) - S(Y-1)/b;
  
  // ADREPORTS
  
  REPORT(a);
  REPORT(b);
  REPORT(max_S);
  REPORT(min_S);
          
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(eta);
  ADREPORT(q);
  ADREPORT(U);
  ADREPORT(N_S);
  ADREPORT(N_E);
  ADREPORT(pred_n);
  ADREPORT(N);
  ADREPORT(n0_new);
    
  return nll;
}

