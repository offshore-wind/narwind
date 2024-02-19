#ifndef BIOENERGETICS_H
#define BIOENERGETICS_H

#include <RcppEigen.h>
#include <random>
#include <cmath>
#include <cstdio>
#include <vector>
#include "spline.h"
#include <algorithm>    // std::all_of
#include <array>        // std::array

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

/// Truncated Normal RNG ----------------------------
// Functions taken from the RcppTN package
// https://github.com/olmjo/RcppTN/tree/master

inline bool CheckSimple(const double low, ///< lower bound of distribution
                        const double high ///< upper bound of distribution
) {
  // Init Values Used in Inequality of Interest
  double val1 = (2 * sqrt(exp(1.0))) / (low + sqrt(pow(low, 2.0) + 4));
  double val2 = exp((pow(low, 2.0) - low * sqrt(pow(low, 2.0) + 4)) / (4)) ;
  //
  
  // Test if Simple is Preferred
  if (high > low + val1 * val2) {
    return true ;
  } else {
    return false ;
  }
}


/// Draw using algorithm 1.

/// 
/// Naive Accept-Reject algorithm.
/// 
inline double UseAlg1(const double low, ///< lower bound of distribution
                      const double high ///< upper bound of distribution
) {
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Init Draw Storage
  double z = 0.0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    z = Rf_rnorm(0.0, 1.0) ;
    
    if (z <= high && z >= low) {
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}

/// Draw using algorithm 2.

/// 
///  Accept-Reject Algorithm
///

inline double UseAlg2(const double low ///< lower bound of distribution
) {
  // Init Values
  const double alphastar = (low +
                            sqrt(pow(low, 2.0) + 4.0)
  ) / (2.0) ;
  const double alpha = alphastar ;
  double e = 0 ;
  double z = 0 ;
  double rho = 0 ;
  double u = 0 ;
  //
  
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    e = Rf_rexp(1.0) ;
    z = low + e / alpha ;
    
    rho = exp(-pow(alpha - z, 2.0) / 2) ;
    u = Rf_runif(0, 1) ;
    if (u <= rho) {
      // Keep Successes
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}

/// Draw using algorithm 3.

/// 
/// Accept-Reject Algorithm
/// 

inline double UseAlg3(const double low, ///< lower bound of distribution
                      const double high ///< upper bound of distribution
) {
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Declare Qtys
  double rho = 0 ;
  double z = 0 ;
  double u = 0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    z = Rf_runif(low, high) ;
    if (0 < low) {
      rho = exp((pow(low, 2.0) - pow(z, 2.0)) / 2) ;
    } else if (high < 0) {
      rho = exp((pow(high, 2.0) - pow(z, 2.0)) / 2) ;
    } else if (0 < high && low < 0) {
      rho = exp(- pow(z, 2.0) / 2) ;
    }
    
    u = Rf_runif(0, 1) ;
    if (u <= rho) {
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}


/// Draw from an arbitrary truncated normal distribution.

///
/// See Robert (1995): <br />
/// Reference Type: Journal Article <br />
/// Author: Robert, Christian P. <br />
/// Primary Title: Simulation of truncated normal variables <br />
/// Journal Name: Statistics and Computing <br />
/// Cover Date: 1995-06-01 <br />
/// Publisher: Springer Netherlands <br />
/// Issn: 0960-3174 <br />
/// Subject: Mathematics and Statistics <br />
// Start Page: 121 <br />
// End Page: 125 <br />
/// Volume: 5 <br />
/// Issue: 2 <br />
/// Url: http://dx.doi.org/10.1007/BF00143942 <br />
/// Doi: 10.1007/BF00143942 <br />
///
// [[Rcpp::export]]
double rtnorm(const double mean,
              const double sd,
              const double low,
              const double high
) {
  // Namespace
  using namespace Rcpp ;
  //

  // Init Useful Values
  double draw = 0;
  int type = 0 ;
  int valid = 0 ; // used only when switching to a simplified version
  // of Alg 2 within Type 4 instead of the less
  // efficient Alg 3
  //

  // Set Current Distributional Parameters
  const double c_mean = mean ;
  double c_sd = sd ;
  const double c_low = low ;
  const double c_high = high ;
  double c_stdlow = (c_low - c_mean) / c_sd ;
  double c_stdhigh = (c_high - c_mean) / c_sd ; // bounds are standardized
  //

  // Map Conceptual Cases to Algorithm Cases
  // Case 1 (Simple Deterministic AR)
  // mu \in [low, high]
  if (0 <= c_stdhigh &&
      0 >= c_stdlow
  ) {
    type = 1 ;
  }

  // Case 2 (Robert 2009 AR)
  // mu < low, high = Inf
  if (0 < c_stdlow &&
      c_stdhigh == INFINITY
  ) {
    type = 2 ;
  }

  // Case 3 (Robert 2009 AR)
  // high < mu, low = -Inf
  if (0 > c_stdhigh &&
      c_stdlow == -INFINITY
  ) {
    type = 3 ;
  }

  // Case 4 (Robert 2009 AR)
  // mu -\in [low, high] & (abs(low) =\= Inf =\= high)
  if ((0 > c_stdhigh || 0 < c_stdlow) &&
      !(c_stdhigh == INFINITY || c_stdlow == -INFINITY)
  ) {
    type = 4 ;
  }

  ////////////
  // Type 1 //
  ////////////
  if (type == 1) {
    draw = UseAlg1(c_stdlow, c_stdhigh) ;
  }

  ////////////
  // Type 3 //
  ////////////
  if (type == 3) {
    c_stdlow = -1 * c_stdhigh ;
    c_stdhigh = INFINITY ;
    c_sd = -1 * c_sd ; // hack to get two negative signs to cancel out

    // Use Algorithm #2 Post-Adjustments
    type = 2 ;
  }

  ////////////
  // Type 2 //
  ////////////
  if (type == 2) {
    draw = UseAlg2(c_stdlow) ;
  }

  ////////////
  // Type 4 //
  ////////////
  if (type == 4) {
    if (CheckSimple(c_stdlow, c_stdhigh)) {
      while (valid == 0) {
        draw = UseAlg2(c_stdlow) ;
        // use the simple
        // algorithm if it is more
        // efficient
        if (draw <= c_stdhigh) {
          valid = 1 ;
        }
      }
    } else {
      draw = UseAlg3(c_stdlow, c_stdhigh) ; // use the complex
      // algorithm if the simple
      // is less efficient
    }
  }



  // Returns
  return  c_mean + c_sd * draw ;
  //
}

// [[Rcpp::export]]
Rcpp::NumericMatrix transpose_c(Rcpp::NumericMatrix m, int k){

  // int nc = k - (j + 1);
  Rcpp::NumericMatrix out(1, k);
  for(int i = 0; i < k ; i++){
    out(0, i) = m(i + 1, 0);
  }
  return(out);
}

// [[Rcpp::export]]
int multinomial(Rcpp::NumericVector probs) {
  int k = probs.size();
  Rcpp::IntegerVector ans(k);
  R::rmultinom(1, probs.begin(), k, ans.begin());
  int pos;
  for (int j=0; j<ans.size(); j++){
    if(ans[j]==1) pos= j;
  }
  return(pos);
}

// [[Rcpp::export]]
Rcpp::NumericVector random_int(int n, int lwr = 0, int uppr = 3){
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist(lwr, uppr); // distribution in range [1, 6]
  Rcpp::NumericVector out(n);
  for(int i = 0 ; i < n ; i++) out(i) = dist(rng);
  return(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector random_multivariate_normal(const Eigen::MatrixXd mu, const Eigen::MatrixXd Sigma){

  Rcpp::NumericVector out(2);
  int P = mu.rows(), i = 0;
  Eigen::MatrixXd y(Eigen::MatrixXd(P, 1).setZero()); // Create a column matrix of zeroes
  Eigen::MatrixXd z(Eigen::MatrixXd(P, 1).setZero());

  for(i = 0 ; i < P ; i++) z(i, 0) = Rf_rnorm(0, 1); // To get original R api, use Rf_*

  y = mu + Sigma.llt().matrixL() * z;

  out[0] = y(0,0);
  out[1] = y(1,0);
  return out;
}

// // [[Rcpp::export]]
// Rcpp::NumericVector csample_num( NumericVector x,
//                            int size,
//                            bool replace,
//                            NumericVector prob = Rcpp::NumericVector::create()
// ) {
//   Rcpp::NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
//   return ret;
// }


// [[Rcpp::export]]
Rcpp::NumericVector prob_migration(int n, std::string destination, int cohortID){
  
  // Migration is condition-dependent in NARW.
  // 
  // Gowan et al. (2019) show that:
  // (1) All individuals have the potential to travel to the SEUS, which is an important wintering ground used by all cohorts
  // (2) Migration probabilities are higher for juveniles than for non-breeding adults
  // (3) Resting females have close to 0 probability of migrating to the SEUS - this is called skipped breeding partial migration
  // (4) This probability is 1 for calving females
  // (5) This probability is more variable for non-breeders
  // 
  // Crowe et al. (2021) show that 40% of the population migrates to the GSL
  
  Rcpp::NumericVector out(n);
  
  for(int i = 0; i<n; i++){
    
    if(destination == "GSL"){
      
      out(i) = R::rbinom(1, 0.4);
        // multinomial(Rcpp::NumericVector::create(0.6, 0.4));
      
    } else if(destination == "SEUS"){
      
      if(cohortID <= 2){ // Juveniles (males and females)
        // Higher migration probabilities for juveniles than adults -- Gowan et al. (2019)
        out(i) = R::rbinom(1, 0.25);
        
      } else if (cohortID == 3){ // Adult males
        // Males more likely to visit SEUS than non-calving females -- Gowan et al. (2019)
        out(i) = R::rbinom(1,0.1);
        
      } else if (cohortID == 4){ // Adult females (pregnant)
        // Lack of use by females in years preceding calving -- Gowan et al. (2019)
        out(i) = 0;
        
      } else if (cohortID == 5){ // Adult females (lactating)
        out(i) = 1;
        
      } else if (cohortID == 6){ // Adult females (resting)
        // Lack of use by females in years following calving -- Gowan et al. (2019)
        out(i) = 0;
      }
      
    }
  } // End for loop
  return(out);
}

// // [[Rcpp::export]]
// int southward_migration(int cohortID){
//   
//   double p = 0;
//   
//   if(cohortID <= 2){ // Juveniles (males and females)
//     // Higher migration probabilities for juveniles than adults -- Gowan et al. (2019)
//     p = 0.25;
//   } else if (cohortID == 3){ // Adult males
//     // Males more likely to visit SEUS than non-calving females -- Gowan et al. (2019)
//     p = 0.1;
//   } else if (cohortID == 4){ // Adult females (pregnant)
//     // Lack of use by females in years preceding calving -- Gowan et al. (2019)
//     p = 0;
//   } else if (cohortID == 5){ // Adult females (lactating)
//     p = 1;
//   } else if (cohortID == 6){ // Adult females (resting)
//     // Lack of use by females in years following calving -- Gowan et al. (2019)
//     p = 0;
//   }
// 
//   return(R::rbinom(1,p));
// }

// // [[Rcpp::export]]
// double response_threshold(Rcpp::NumericVector db){
//   std::random_device rd;     // Only used once to initialize (seed) engine
//   std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
//   std::uniform_int_distribution<int> uniform(0,9999);
//   int d = uniform(rd);
//   return(db[d]);
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector testseed(int seed){
//   std::random_device rd;     // Only used once to initialize (seed) engine
//   std::mt19937 mt;    // Random-number engine used (Mersenne-Twister in this case)
//   mt.seed(seed);
//   std::uniform_int_distribution<int> uniform(0,9999);
//   Rcpp::NumericVector db(366);
//   for(int i = 0; i <366; ++i){
//     db[i] = uniform(mt);
//   }
//   return(db);
// }

// [[Rcpp::export]]
double response_threshold(Rcpp::NumericVector dose, int day, int simduration, int seed){
  std::random_device rd; // Only used once to initialize (seed) engine
  std::mt19937 mt; // Random-number engine used (Mersenne-Twister in this case)
  mt.seed(seed); // Set the seed
  std::uniform_int_distribution<int> uniform(0,9999);
  Rcpp::NumericVector db(simduration);
  for(int i = 0; i <simduration; ++i){
    db[i] = uniform(mt);
  }
  double thresh = dose[db(day)];
  return(thresh);
}

// // [[Rcpp::export]]
// std::vector<double> dose_range(double lwr, double uppr, double n){
//   std::vector<double> v(n);
//   for (int i = 0; i < n; i++) v[i] = lwr + i * ((uppr - lwr)/(n-1));
//   return v;
// }

// [[Rcpp::export]]
double starvation_mortality(double bc,
                            Eigen::MatrixXd mortality,
                            double starvation_death,
                            double starvation_onset)
{
  double out;
  int n = mortality.rows();
  std::vector<double> p_x(n);
  std::vector<double> p_y(n);
  
  for(int i = 0; i < n; i++){
    p_x[i] = mortality(i,0);
    p_y[i] = mortality(i,1);
  }
  
  tk::spline s(p_x, p_y); // Fit cubic spline
  
  if(bc > starvation_onset){
    out = 0.0L;
  } else if (bc < starvation_death){
    out = 1.0L;
  } else {
    out = s(bc);
  }

  return out;
}
                        
// [[Rcpp::export]]
Rcpp::NumericVector starvation_mortality_vec(Rcpp::NumericVector bc,
                      Eigen::MatrixXd mortality,
                      double starvation_death,
                      double starvation_onset)
{
  int v = bc.size();
  Rcpp::NumericVector out(v);
  int n = mortality.rows();
  std::vector<double> p_x(n);
  std::vector<double> p_y(n);
  
  for(int i = 0; i < n; i++){
    p_x[i] = mortality(i,0);
    p_y[i] = mortality(i,1);
  }
  
  tk::spline s(p_x, p_y); // Fit cubic spline  
  
  for(int j = 0; j < v; j++){
    if(bc(j) > starvation_onset){
      out(j) = 0.0L;
    } else if (bc(j) < starvation_death){
      out(j) = 1.0L;
    } else {
      out(j) = s(bc(j));
    }
  }

  return out;
}

// // [[Rcpp::export]]
// double prob_response(std::vector<double> x, // Range of dose - from 80 to 200 dB
//                      Eigen::MatrixXd p, // Dose-response functions (5,000 realizations from expert elicitation)
//                      int id, // Which dose-response function to use
//                      double z) // Dose at which response is to be evaluated
// {
//   double out;
//   p = p.col(id); // Extract relevant column
//   std::vector<double> pr(p.data(), p.data() + p.rows()); // Convert to vector<double>
//   tk::spline s(x,pr); // Fit cubic spline
//   out = s(z); // Retrieve value
//   if(out > 1) out = 1;
//   if(out < 0) out = 0;
//   return out; 
// }

// //' Random deviate from a truncated Normal distribution
// //' @name rtnorm
// //' @param location Location parameter
// //' @param scale Scale parameter
// //' @param L Lower bound
// //' @param U Upper bound
// // [[Rcpp::export]]
// double rtnorm(double location,
//                double scale,
//                double L,
//                double U) {
// 
//    double tot = R::pnorm(L, location, scale, TRUE, FALSE) + R::runif(0,1) * (R::pnorm(U, location, scale, TRUE, FALSE) - R::pnorm(L, location, scale, TRUE, FALSE));
//    double out = location + scale * R::qnorm(tot, 0, 1, TRUE, FALSE);
//    return out;
// }
// 
// //' Random deviate from a truncated Normal distribution
// //' @name rtnorm_vec
// //' @param location Location parameter
// //' @param scale Scale parameter
// //' @param L Lower bound
// //' @param U Upper bound
// // [[Rcpp::export]]
// Rcpp::NumericVector rtnorm_vec(int n,
//                                 double location,
//                                 double scale,
//                                 double L,
//                                 double U) {
// 
//    Rcpp::NumericVector out(n);
//    Rcpp::NumericVector tot(n);
// 
//    for(int i = 0; i<n; i++){
//      tot[i] = R::pnorm(L, location, scale, TRUE, FALSE) + R::runif(0,1) * (R::pnorm(U, location, scale, TRUE, FALSE) - R::pnorm(L, location, scale, TRUE, FALSE));
//      out[i] = location + scale * R::qnorm(tot[i], 0, 1, TRUE, FALSE);
//    }
// 
//    return out;
//  }

// [[Rcpp::export]]
double is_female(int cohort){
  if(cohort == 0){
    return R::rbinom(1, 0.46);
  } else if(cohort == 1 | cohort == 3){
    return 0;
  } else {
    return 1;
  }
}


//' Initialize age
//' @name start_age
 //' @description Performs a random draw from a uniform distribution to initialize the age of simulated animals
 //' @param cohort Integer indicating which cohort the animal belongs to
 // [[Rcpp::export]]
 
 double start_age(int cohort){
   if(cohort == 0){
     return 0;
   } else if(cohort >= 1 && cohort <= 2){
     return R::runif(1, 8);
   } else {
     return R::runif(9, 68); // Longevity of 69 years (-1 for 1 year simulation)
   }
 }


// [[Rcpp::export]]

Rcpp::NumericVector start_age_vec(Rcpp::NumericVector cohort){
  
  int n = cohort.size();
  Rcpp::NumericVector a(n);
  
  for(int i = 0; i < n; i++) {
    if(cohort[i] == 0){
      a[i] = 0;
    } else if(cohort[i] >= 1 && cohort[i] <= 2){
      a[i] = R::runif(1, 8);
    } else {
      a[i] = R::runif(9, 68);
    }
  }
  return(a);
}

// //' Initialize entanglement load
//  //' @description Performs a random draw from a discrete uniform distribution to initialize the entanglement state of simulated animals
//  // [[Rcpp::export]]
//  int start_entangled(double p_entangled = 0){                     
// 
//    
// 
//    int is_entangled = R::rbinom(1, p_entangled); // Incidence of entanglement
//    return is_entangled;
//  }

// // [[Rcpp::export]]
// double survivorship(double age,
//                     double a1 = 0.1,
//                     double a2 = 0,
//                     double a3 = 0.01,
//                     double b1 = 60,
//                     double b3 = 8,
//                     double longevity = 69){
// 
//   // From Barlow and Boveng (1991)
// 
//   double lj, lc, ls, out;
// 
//   // Exponentially decreasing risk due to juvenile mortality factors
//   lj = std::exp((-a1/b1) * (1-std::exp(-b1*age/longevity)));
//   
//   // Constant risk due to other factors
//   lc = std::exp(-a2*age/longevity);
//   
//   // Exponentially increasing risk due to senescent mortality factors
//   ls = std::exp((a3/b3) * (1-exp(b3*age/longevity)));
//   
//   out = lj * ls * lc;
//   
//   return(out);
// 
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector survivorship_vec(Rcpp::NumericVector age,
//                                      double a1 = 0.1,
//                                      double a2 = 0,
//                                      double a3 = 0.01,
//                                      double b1 = 60,
//                                      double b3 = 8,
//                                      double longevity = 69){
//   
//   // From Barlow and Boveng (1991)
//   
//   int n = age.size();
//   Rcpp::NumericVector out(n);
//   double lj, lc, ls;
//   
//   for(int i = 0; i<n; i++){
//     // Exponentially decreasing risk due to juvenile mortality factors
//     lj = std::exp((-a1/b1) * (1-std::exp(-b1*age(i)/longevity)));
//     
//     // Constant risk due to other factors
//     lc = std::exp(-a2*age(i)/longevity);
//     
//     // Exponentially increasing risk due to senescent mortality factors
//     ls = std::exp((a3/b3) * (1-exp(b3*age(i)/longevity)));
//     
//     out(i) = lj * ls * lc;
//   }
// 
//   return(out);
//   
// }


//' Entanglement event
//' @name entanglement_event
//' @param p_head Probability that the entanglement involves the anterior region of the body (mouth, head, rostrum)
// [[Rcpp::export]]
 
 Rcpp::NumericVector entanglement_event(double p_entangled = 0,
                                        double p_head = 0.732, 
                                        Rcpp::NumericVector p_severity = Rcpp::NumericVector::create(0.819, 0.140, 0.041)
                                        ){       
   
   // p(head) weighted mean of entries in spreadsheet of model parameters,
   
   Rcpp::NumericVector out (7); // Store results
   
   // Is the animal entangled?
   // Annual entanglement rate = 0.259 -- Knowlton et al. (2012)
   // So probability of becoming entangled on each day is 0.000821 (see SI in Pirotta et al. 2023)
   int is_entangled = R::rbinom(1, p_entangled);
   
   if(is_entangled){
     
     out(0) = is_entangled; // First value

     // How severe is the entanglement (0 = minor, 1 = moderate, 2 = severe)
     out(1) = multinomial(p_severity); // Second value
     
     // Does the entanglement involve the anterior region of the body?
     out(2) = R::rbinom(1, p_head);  // Third value
     
     if(out(2) == 1){
       // Reduction in the mouth gape
       if(out(1)==0){
         out(3) = R::runif(0,0.15);  // Fourth value
       } else if(out(1)==1){
         out(3) = R::runif(0.15,0.3);
       } else if(out(1)==2){
         out(3) = R::runif(0.3,0.5);
       }
     } 
     
     // Duration of entanglement (days)
     // From Pirotta et al. (2023 Oikos)
     // See entgl_durations() function
     // All long-term events involve moderate or severe injuries (Knowlton et al. 2016 Cons Biol)
     if(out(1) == 0){
       out(4) = Rf_rnbinom_mu(0.7363271, 111.3911746); // Fifth value
     } else if(out(1) == 1){
       out(4) = Rf_rnbinom_mu(0.7493824, 124.2184062);
     } else {
       out(4) = Rf_rnbinom_mu(0.6527522, 211.8187153);
     }
     
   }
   
   return out;
 }

//' Mouth-width-to-length ratio
//' @name start_mouth
 //' @description Returns the ratio between the width of the animal's mouth (defined as body width measured at 10% of the body length from the snout) and its body length, as per Miller et al. (2012)
 //' @param cohort Population cohort to which the animal belongs
 //' @param age Age of the animal (years)
 //' @note Assumes that juveniles follow the ratio defined by older calves and adult males that of resting females
 //' @return Estimate mouth-to-width ratio
 // [[Rcpp::export]]
 
 double start_mouth(int cohort, double age){             
   double r = 0.0;
   if(age < 30/365) { // calf (1st month of suckling)
     r = R::rnorm(0.125, 0.00430);
   } else if (age < 9){ // calf (3rd month of suckling) + Juveniles
     r = R::rnorm(0.144, 0.00686);
   } else {
     if(cohort == 3 | cohort == 6){
       r = R::rnorm(0.149, 0.00240);
     } else if(cohort == 4){
       r = R::rnorm(0.157, 0.00610);
     } else if (cohort == 5){
       r = R::rnorm(0.149, 0.00634);
     }
   }
   return r;
 }

//' Initialize body condition
//' @name start_bcondition
//' @description Performs a random draw from a beta distribution to initialize 
//' the body condition of simulated animals, taken as the ratio of fat mass to total mass.
//' @param age Age in years
//' @param shape1 First shape parameter of the beta distribution
//' @param shape2 Second shape parameter of the beta distribution
// [[Rcpp::export]]
long double start_bcondition(double cohort){
  
  // Mean and SD values obtained from the <init_bc> function based on Christiansen et al.
  // Lower bound = starvation threshold
  // Upper bound = maximum allowable BC (see <find_maxBC> function)
  
  long double bc;
  double lower = 0.05;
  double upper = 0.665289; // find_maxBC()
  
  // Calves
  if(cohort == 0){
    
    bc = rtnorm(0.35884621, 0.07788818, lower, upper);
    
    // Juveniles (males and females)
  } else if(cohort == 1 | cohort == 2){
    
    // bc = rtnorm(0.28574784, 0.02027885, 0.05, 0.665289);
    bc = rtnorm(0.2904040, 0.0681619, lower, upper);
    
    // Adult males, resting females, pregnant females
  } else if(cohort == 3 | cohort == 4 | cohort == 6) {
    
    // bc = rtnorm(0.25940715, 0.01519702, 0.05, 0.665289);
    bc = rtnorm(0.27617455, 0.06969136, lower, upper);
    
    // Lactating females, which start simulation as late pregnant
  } else {
    
    bc = rtnorm(0.40950094, 0.05721548, lower, upper);
    
  }
  return bc;
  
}

//' Initialize body condition 
//' @name start_bcondition_vec
//' @description Performs a random draw from a beta distribution to initialize
//' the body condition of simulated animals, taken as the ration of fat mass to total mass.
//' @param age Age in years
//' @param shape1 First shape parameter of the beta distribution
//' @param shape2 Second shape parameter of the beta distribution
// [[Rcpp::export]]
Rcpp::NumericVector start_bcondition_vec(Rcpp::NumericVector cohort, int month = 10){

  int n = cohort.size();
  Rcpp::NumericVector bc(n);
  double lower = 0.05;
  double upper = 0.665289; // find_maxBC()

  for(int i = 0; i < n; i++) {

    // Calves
    if(cohort[i] == 0){ 
      
      bc[i] = rtnorm(0.35884621, 0.07788818, lower, upper);
      
      // Juveniles (males and females)
    } else if(cohort[i] == 1 | cohort[i] == 2){
      
      // bc[i] = rtnorm(0.28574784, 0.02027885, 0.05, 0.665289);
      bc[i] = rtnorm(0.2904040, 0.0681619, lower, upper);
      
      // Adult males, resting females, pregnant females
    } else if(cohort[i] == 3 | cohort[i] == 4 | cohort[i] == 6) {
      
      // bc[i] = rtnorm(0.25940715, 0.01519702, 0.05, 0.665289);
      bc[i] = rtnorm(0.27617455, 0.06969136, lower, upper);
      
      // Lactating females, which start simulation as late pregnant
    } else if(cohort[i] == 5){ 
      
      bc[i] = rtnorm(0.40950094, 0.05721548, lower, upper);
      
    }
    
  }

  return bc;
}

// ' Parameters of the length-at-age relationship
// ' @description Generates vectors of coefficients from the Gompertz growth
//  curve describing the change in length as a function of age relationship, as described in Fortune et al. (2021)
// ' @return Matrix of coefficients
// [[Rcpp::export]]
Eigen::MatrixXd agL(double age, int n = 1, bool sd = false){
  
  double a, b, c; // Parameters of the Gompertz growth curves        
  Eigen::MatrixXd out(n,3);
  
  for(int i = 0; i < n; i++){
    
    if(age <= 0.79){
      if(sd){
        a = R::rnorm(1067.19, 19.67);
        b = R::rnorm(-0.93, 0.08);
        c = R::rnorm(-3.11, 0.28);
      } else {
        a = 1067.19;
        b = -0.93;
        c = -3.11;
      }
    } else {
      if(sd){
        a = R::rnorm(1362.755, 22.88); 
        b = R::rnorm(-0.37, 0.03); 
        c = R::rnorm(-0.18, 0.03);
      } else {
        a = 1362.755;
        b = -0.37;
        c = -0.18;
      }
    }
    
    out(i,0) = a;
    out(i,1) = b;
    out(i,2) = c;
    
  }
  return(out);
}

// [[Rcpp::export]]
Eigen::MatrixXd agL_vec(Rcpp::NumericVector age, bool sd = false){
  
  int n = age.size();
  Eigen::MatrixXd out(n,3);
  
  for(int i = 0; i < n; i++){
    
    double a, b, c; // Parameters of the Gompertz growth curves  
    
    if(age[i] <= 0.79){
      if(sd){
        a = R::rnorm(1067.19, 19.67);
        b = R::rnorm(-0.93, 0.08);
        c = R::rnorm(-3.11, 0.28);
      } else {
        a = 1067.19;
        b = -0.93;
        c = -3.11;
      }
    } else {
      if(sd){
        a = R::rnorm(1362.755, 22.88); 
        b = R::rnorm(-0.37, 0.03); 
        c = R::rnorm(-0.18, 0.03);
      } else {
        a = 1362.755;
        b = -0.37;
        c = -0.18;
      }
    }
    
    out(i,0) = a;
    out(i,1) = b;
    out(i,2) = c;
    
  }
  return(out);
}

//' Age to length conversion
//' @name age2length
 //' @description Calculates body length from age according to the length-at-age relationship described
 //' in Fortune et al. (2021)
 //' @param age Age of the animal (years)
 //' @param gompertz Parameters of the Gompertz curve, as returned by agL()
 //' @return Estimated length (m)
 // [[Rcpp::export]]
 double age2length(double age, Eigen::MatrixXd gompertz){             
   return gompertz(0,0) * exp(gompertz(0,1)*exp(gompertz(0,2)*age)) / 100; // cm to m
 }

// // [[Rcpp::export]]
// Rcpp::NumericVector age2length_vec(Rcpp::NumericVector age, Eigen::MatrixXd gompertz){        
//   int n = age.size();
//   int nr = gompertz.rows();
//   if(n != nr) Rcpp::stop("Arguments have different lengths");
//   Rcpp::NumericVector a(n);
//   
//   for(int i = 0; i < n; i++) {
//     a[i] = gompertz(i,0) * exp(gompertz(i,1)*exp(gompertz(i,2)*age[i])) / 100; // cm to m
//   }
//   return a;
// }

// [[Rcpp::export]]
Rcpp::NumericVector age2length_vec(Rcpp::NumericVector age){        
  int n = age.size();
  Rcpp::NumericVector a(n);
  Eigen::MatrixXd gompertz = agL_vec(age);
  for(int i = 0; i < n; i++) {
    a[i] = gompertz(i,0) * exp(gompertz(i,1)*exp(gompertz(i,2)*age[i])) / 100; // cm to m
  }
  return a;
}

// // [[Rcpp::export]]
// Rcpp::NumericVector age2length_vec(Rcpp::NumericVector age){        
//   int n = age.size();
//   Rcpp::NumericVector a(n);
// 
//   for(int i = 0; i < n; i++) {
//     if(age[i] <= 0.79){
//       a[i] = 1067.19 * exp(-0.93*exp(-3.11*age[i])) / 100; // cm to m
//     } else {
//       a[i] = 1362.755 * exp(-0.37*exp(-0.18*age[i])) / 100; // cm to m
//     }
//     
//   }
//   return a;
// }

// // [[Rcpp::export]]
// Eigen::MatrixXd create_mat(){
//   return Eigen::MatrixXd(2,1);
// }

// ' Parameters of the mass-at-length relationship
// ' @description Generates pairs of coefficients (intercept and slope) from the logarithmic mass-at-length relationship
//   described in Fortune et al. (2021)
// ' @return Matrix of coefficients
// [[Rcpp::export]]
Eigen::MatrixXd mL(int n = 1, bool sd = false){
  
  Eigen::MatrixXd out = Eigen::MatrixXd(n,2);
  
  if(sd){
    
    // Multivariate parameters of the mass-at-length function
    // See growth_curves.R script
    Eigen::MatrixXd mu = Eigen::MatrixXd(2,1);
    Eigen::MatrixXd sigma = Eigen::MatrixXd(2,2);
    
    // Means
    mu(0,0) = -4.834189;
    mu(1,0) = 2.984353;
    
    // Variance-covariance matrix
    sigma(0,0) = 0.21128304;
    sigma(0,1) = -0.07515154;
    sigma(1,0) = -0.07515154;
    sigma(1,1) = 0.02686353;
    
    for (int i = 0; i < n; i++) {
      Rcpp::NumericVector deviates = random_multivariate_normal(mu, sigma);
      out(i,0) = deviates[0];
      out(i,1) = deviates[1];
    }
  } else {
    for (int i = 0; i < n; i++) {
      // Atlantic and Pacific right whales
      // out(i,0) = -5.096379;  
      // out(i,1) = 3.079408;
      // Atlantic only
      out(i,0) = -4.834189;
      out(i,1) = 2.984353;
    }
  }
  return(out);
}

// ' Length to mass conversion
// ' @description Calculates body mass from body length according to the logarithmic
// ' mass-at-length relationship of Fortune et al. (2021)
// ' @param L Input body length (m)
// ' @param param Parameters of the mass-at-length relationship, as returned by mL
// ' @return Estimated mass (kg)
// [[Rcpp::export]]
double length2mass(double L,
                   Eigen::MatrixXd param,
                   double lean = 0.5554227){
  
  // lean determined by find_lean()
  
  double a = param(0,0);
  double b = param(0,1);
  double L_cm = L * 100;
  return lean * std::pow(10, a)*std::pow(L_cm,b); // m to cm
}

// // [[Rcpp::export]]
// Rcpp::NumericVector length2mass_vec(Rcpp::NumericVector L,
//                        Eigen::MatrixXd param,
//                        double lean = 0.5435686){
//   
//   int n = L.size();
//   int nr = param.rows();
//   if(n != nr) Rcpp::stop("Arguments have different lengths");
//   Rcpp::NumericVector mass(n);
//   
//   for(int i = 0; i < n; i++) {
//     mass[i] = lean * std::pow(10, param(i,0))*std::pow(L[i] * 100,param(i,1)); // m to cm
//   }
//   return(mass);
// }

// [[Rcpp::export]]
Rcpp::NumericVector length2mass_vec(Rcpp::NumericVector L, double lean = 0.5435686){
  
  // lean determined by find_lean(45000)
  int n = L.size();
  Rcpp::NumericVector mass(n);
  
  for(int i = 0; i < n; i++) {
    mass[i] = lean * std::pow(10, -4.834189)*std::pow(L[i] * 100, 2.984353); // m to cm
  }
  return(mass);
}
// double length2mass(double L, 
//                    double age,
//                    bool lean){
//   
//   double a, b;
//   double L_cm = L * 100;        
//   double perc = 1;
//   if(lean) perc = 0.73; 
//   
//   // Eigen::MatrixXd mv = random_multivariate_normal(mu, sigma);
//   // a = mv(0,0);
//   // b = mv(1,0);
//   
//   if(age <= 0.79){
//     a = R::rnorm(-5.091821, 0.2578327);
//     b = R::rnorm(3.077823, 0.08325852);
//   } else {
//     a = R::rnorm(-5.096379, 0.2592405);
//     b = R::rnorm(3.079408, 0.08360103);
//   }
//   
//   return perc * pow(10, a + b*std::log10(L*100));
// }


// double length2mass(double L,
//                    Eigen::MatrixXd mu,
//                    Eigen::MatrixXd sigma,
//                    // double age,
//                    bool lean){
// 
//   double a, b;
//   double L_cm = L * 100;
//   double perc = 1;
//   if(lean) perc = 0.73;
// 
//   Eigen::MatrixXd mv = random_multivariate_normal(mu, sigma);
//   a = mv(0,0);
//   b = mv(1,0);
// 
//   // if(age <= 0.79){
//   //   a = R::rnorm(-5.091821, 0.2578327);
//   //   b = R::rnorm(3.077823, 0.08325852);
//   // } else {
//   //   a = R::rnorm(-5.096379, 0.2592405);
//   //   b = R::rnorm(3.079408, 0.08360103);
//   // }
// 
//  return std::pow(10, a)*std::pow(L_cm,b); // m to cm
// 
//   // return perc * pow(10, a + b*std::log10(L*100));
// }

// // [[Rcpp::export]]
// int separation_duration(){
//   // The returned value is automatically converted to an integer by defining a function of type <int>
//   return rtnorm(0, 6.152344, 1, 23);
// } 

// [[Rcpp::export]]
Rcpp::NumericVector increment_cohort(Rcpp::NumericVector alive,
                                     Rcpp::NumericVector cohort, 
                                     Rcpp::NumericVector age, 
                                     Rcpp::NumericVector female,
                                     Rcpp::NumericVector bc,
                                     Rcpp::NumericVector min_bc,
                                     Rcpp::NumericVector reprod,
                                     double abort){
  int nal = alive.size();
  int nc = cohort.size();
  int na = age.size();
  int nf = female.size();
  int nbc = bc.size();
  int nmin = min_bc.size();
  int nrep = reprod.size();
  
  // Check that all vectors have the same size
  std::array<int,7> sizes = {nal, nc, na, nf, nbc, nmin, nrep};
  int allequal = std::equal(sizes.begin() + 1, sizes.end(), sizes.begin());
  if(allequal == 0) Rcpp::stop("Arguments have different lengths");
  
  Rcpp::NumericVector coh(nc);
  
  for(int i = 0; i < nc; i++) {
    
    if(alive[i] == 0){
      
      coh[i] = cohort[i];
      
    } else {
      
      // ----------------------
      // Growth and maturity
      // ----------------------
      
      if(cohort[i] == 0 & age[i] >= 1 & female[i] == 0){              // Calves (m) -> Juveniles (m)
        coh[i] = 1;
      } else if (cohort[i] == 0 & age[i] >= 1 & female[i] == 1){      // Calves (f) -> Juveniles (f)
        coh[i] = 2;
      } else if (cohort[i] == 1 & age[i] >= 9 & female[i] == 0){      // Juveniles (m) -> Adults (m)
        coh[i] = 3;
      } else if (cohort[i] == 2 & age[i] >= 9 & female[i] == 1){      // Juveniles (f) -> Adults, resting (f)
        coh[i] = 6;
        
        // ----------------------
        // Transitions between reproductive states
        // ----------------------
        
      } else if (cohort[i] == 5) {                                    // Lactating (f) -> Resting (f)
        coh[i] = 6;
      } else if(cohort[i] == 6){  
        if(bc[i] > min_bc[i] & reprod[i] == 1){                                        // Resting (f) -> Pregnant (f)
          coh[i] = 4; 
        } else {
          coh[i] = 6;                                                 // Resting (f) -> Resting (f)
        }
      } else if(cohort[i] == 4){                                      // Pregnant (f) -> Lactating (f)
        int has_aborted = R::rbinom(1, abort);
        if(has_aborted){
          coh[i] = 6;
        } else {
          coh[i] = 5;
        }
      } else {
        coh[i] = cohort[i];
      }
      
    }
    
  }
    
  return coh;
}

//' Incidence of foraging behavior
//' @name feeding_threshold
//' @description Determines whether the prey concentration encountered by an animal
//' is sufficient to support foraging
//' @param min_prey Minimum prey density threshold that triggers foraging 
//' @param D Prey concentration
// [[Rcpp::export]]
double feeding_threshold(double min_prey, double prey_concentration){
   return std::round(1/(1 + exp(min_prey - prey_concentration)));
 }

//' Incidence of foraging behavior (vectorized)
//' @name feeding_threshold_vec
//' @description Determines whether the prey concentration encountered by an animal
//' is sufficient to support foraging
//' @param min_prey Minimum prey density threshold that triggers foraging
//' @param D Prey concentration
// [[Rcpp::export]]
Rcpp::NumericVector feeding_threshold_vec(double min_prey, 
                                           Rcpp::NumericVector prey_concentration){
   int n = prey_concentration.size();
   Rcpp::NumericVector out(n);
     for(int i = 0; i < n; i++) {
       out(i) = std::round(1/(1 + exp(min_prey - prey_concentration(i))));
     };
   return out;
 }


// //' Swim speed during foraging
// //' @description Provides an estimate of swim speed during foraging based on body length
// //' @param L Body length (cm)
// //' @return Estimated swim speed (body length/s)
// // [[Rcpp::export]]
//  
//  double swimspeed(double L){
//    return 0.09 * L;
//  }

// //' Variation in feeding/suckling effort with body condition
// //' @description Defines how feeding/suckling effort varies in response to the animal's body condition
// //' @param eta Parameter controlling the non-linearity of the response
// //' @param beta Target body condition (d.u.)
// //' @param M Total body mass (kg)
// //' @param R Reserve (fat) mass (kg)
// //' @note An increase in feeding effort during periods of compromised health is simulated as a 
// //' behavioral mechanism to compensate for periods of disturbance, should adequate resources be available. 
// //' This is achieved by using a sigmoidal decreasing function (bounded by 0 and 1) of body condition, which 
// //' is designed to equal 0.5 at the animal’s target body condition and ensures that fat reserves do not grow out
// //' of bounds under favorable conditions. While the possibility of early weaning resulting from good body 
// //' condition of both mother and calf is not considered explicitly, suckling effort is made to vary as a 
// //' function of the calf’s reserve mass, in the same way that foraging effort does in juveniles/adults.
// //' @return A scalar on feeding/suckling effort
// // [[Rcpp::export]]
//  
//  double scale_effort(double eta, double beta, double M, double R){
//    return 1/(1 + exp(-eta * ((beta * M / R) - 1)));
//  }

// //' Variation in feeding/suckling effort with body condition
// //'  @name scale_effort
//  //' @description Calculates feeding/suckling effort given an animal's body condition
//  //' @param A Horizontal asymptote when x tends to -Inf
//  //' @param D Horizontal asymptote when x tends to Inf
//  //' @param B Steepness of the transition between the two asymptotes
//  //' @param C Location parameter
//  //' @param S Curve asymmetry (the curve us symmetric when S = 1)
//  //' @note An increase in feeding effort during periods of compromised health is simulated as a 
//  //' behavioral mechanism to compensate for periods of disturbance, should adequate resources be available. 
//  //' This is achieved by using a 5-parameter logistic function (bounded by 0 and 1) of body condition, which 
//  //' is designed to reduce to 0 at the animal’s target body condition and ensures that fat reserves do not grow out
//  //' of bounds under favorable conditions. While the possibility of early weaning resulting from good body 
//  //' condition of both mother and calf is not considered explicitly, suckling effort is made to vary as a 
//  //' function of the calf’s reserve mass, in the same way that foraging effort does in juveniles/adults.
//  //' @return A scalar on feeding/suckling effort
//  // [[Rcpp::export]]
//  
//  double scale_effort(double x, double A, double D, double B, double C, double S){
//    return A + (D-A) / std::pow(1 + std::exp(B*(C-x)), S);
//  }
// 
// // [[Rcpp::export]]
// 
// Rcpp::NumericVector scale_effort_vec(Rcpp::NumericVector x, 
//                                      double A, double D, double B, double C, double S){
//   
//   int n = x.size();
//   Rcpp::NumericVector out(n);
//   for(int i = 0; i<n; i++){
//     out(i) = A + (D-A) / std::pow(1 + std::exp(B*(C-x(i))), S);
//   }
//   return out;
// }


// [[Rcpp::export]]
 double feeding_effort(double eta, double rho, double bc){
   return 1 / (1 + std::exp(-eta*(rho*(1/bc)-1)));
 }

// [[Rcpp::export]]
Rcpp::NumericVector feeding_effort_vec(double eta, double rho, Rcpp::NumericVector bc){
  int n = bc.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i<n; i++){
    out(i) = 1 / (1 + std::exp(-eta*(rho*(1/bc(i))-1)));
  }
  return out;
}

//' Convert degrees to radians
//' @name deg2radians
//' @param angle Input angle in degrees
//' @return Equivalent angle in radians
// [[Rcpp::export]]
 double deg2radians(double angle){
   return angle * M_PI / 180;
 }

//' Area of the gape
//' @name gape_size
//' @description Provides an estimate of the area of the mouth gape
//' @param L Total body length (m)
//' @param omega Width of the mouth (m)
//' @param alpha Angle between the tip of the baleen plates and outer edge of the baleen racks (rad)
//' @note The gape is assumed to be a disk sector defined by the height of the tallest baleen plate 
//' and the width of the mouth. The central angle can be estimated based on 3D models 
//' of right whales and bowhead whales developed from measurements taken during whaling,
//' necropsies, and aerial photogrammetry studies
//' @return Estimated gape area
// [[Rcpp::export]]
double gape_size(double L, double omega, double alpha){ 
   double a = deg2radians(alpha);
   return a * (std::pow(omega,2)/4 + std::pow((0.2077*L - 1.095),2))/2;
 }

// //' Active feeding during foraging
// //' @description Determines the proportion of time spent actively feeding while in foraging mode
// //' @param F Filtration rate (m^3/s)
// //' @param tau_clear Time required to clear the forestomach (s)
// //' @param stomach_size Forestomach capacity (m^3)
// //' @return Proportion of active feeding (\%)
// // [[Rcpp::export]]
//  
//  double active_feed(double F, double tau_clear, double stomach_size){ 
//    return 1 / (1 + (F * tau_clear / stomach_size));
//  }

// //' Energy content of the prey
// //' @name Econtent_cop
//  //' @description Calculates the amount of energy that a whale can obtain from an average individual prey
//  //' @param m Average mass of the prey (g/cop)
//  //' @param rho Energy density of the prey (J/g) 
//  //' @param E_digest Digestive efficiency (fecal and urinary, \%)
//  //' @param E_hif Metabolizing efficiency (1-heat increment of feeding, \%)
//  // [[Rcpp::export]]
//  
//  double Econtent_cop(double m, double rho, double E_digest, double E_hif){ 
//    return m * rho * E_digest * E_hif;
//  }

// //' Filtration rate
// //' @name filtration_rate
//  //' @description Calculates the volume of water filtered per unit time
//  //' @param A Area of the mouth gape (m^2)
//  //' @param lambda_gape Percent reduction in the gape during an entanglement event (\%)
//  //' @param V Swimming speed while filtering (m/s)
//  //' @param E_capt Capture efficiency (\%)
//  //' @param lambda_capt Percent reduction in capture efficiency during an entanglement event (\%)
//  //' @return Estimated filtration rate (m^3/s)
//  // [[Rcpp::export]]
//  
//  double filtration_rate(double A, double lambda_gape, double V, double E_capt, double lambda_capt){ 
//    return A * (1 - lambda_gape) * V * E_capt * (1 - lambda_capt);
//  }


// //' Milk ingestion rate
// //' @name milk_ingestion
//  //' @description Calculates the maximum volume of milk that a calf can ingest per unit time
//  //' @param E_milk Milk transfer/assimilation efficiency (\%)
//  //' @param M Scalar on milk provisioning by the mother (d.u.)
//  //' @param E_gland Mammary gland efficiency (\%)
//  //' @param mu_female Total mass of mammary glands (kg)
//  //' @param delta_female Milk production rate by the mother (m^2)
//  //' @param D Density of milk (kg/m^3)
//  // [[Rcpp::export]]
//  
//  double milk_ingestion(double E_milk, double M, double E_gland, double mu_female, double delta_female, double D){ 
//    return E_milk * M * E_gland * (mu_female * delta_female/D);
//  }

//' Milk assimilation efficiency
//' @name milk_assimilation
 //' @description Calculates the efficiency with which maternal milk is transferred from mother to calf
 //' @param t Time (d)
 //' @param T_lac Duration of the lactation period (i.e., age at weaning) (d)
 //' @param a Age at which milk consumption starts to decrease (d)
 //' @param zeta Non-linearity between milk assimilation and calf age
 // [[Rcpp::export]]
 
 double milk_assimilation(double t, int T_lac, double a, double zeta){ 
   double out;
   if(t <= a){
     out = 1;
   } else if (t >= T_lac){
     out = 0;
   } else {
     out = (1 - (t - a)/(T_lac - a))/(1 - (zeta * (t - a)/(T_lac - a)));
   }
   return out;
 }

// [[Rcpp::export]]

Rcpp::NumericVector milk_assimilation_vec(Rcpp::NumericVector t,
                             int T_lac,
                             double a,
                             double zeta){ 
  int n = t.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i<n; i++){
    if(t(i) <= a){
      out(i) = 1;
    } else if (t(i) >= T_lac){
      out(i) = 0;
    } else {
      out(i) = (1 - (t(i) - a)/(T_lac - a))/(1 - (zeta * (t(i) - a)/(T_lac - a)));
    }
  }
  return out;
}

//' Milk provisioning
//' @name milk_supply
 //' @description Defines how milk provisioning varies in response to the female's body condition
 //' @param kappa Starvation threshold (d.u.)
 //' @param target_condition Target body condition (d.u.)
 //' @param M Total body mass (kg)
 //' @param R Reserve (fat) mass (kg)
 //' @param zeta Non-linearity between milk supply and the body condition of the mother
 //' @note Milk production is only possible when the female’s reserve (i.e., blubber) mass 
 //' is sufficient to offset the costs of maintenance and lactation. Milk supply therefore
 //' varies as a function of the female’s body condition such that it equals 1 if the female
 //' is at or above her target body condition, and ceases altogether when the female is at
 //' or below the starvation threshold. This allows us to emulate the abandonment of a calf
 //' in response to poor female health.
 //' @return A scalar on milk supply
 // [[Rcpp::export]]
 
 double milk_supply(double kappa, double target_condition, double M, double R, double zeta){ 
   double out;
   double bc_female = R/M;
   if(bc_female <= kappa){
     out = 0;
   } else if (bc_female >= target_condition){
     out = 1;
   } else {
     out = ((1-zeta)*(R - kappa * M))/(M * (target_condition - kappa) - zeta * (R - kappa * M));
   }
   return out;
 } 

// [[Rcpp::export]]

Rcpp::NumericVector milk_supply_vec(double kappa,
                                    double target_condition, 
                                    Rcpp::NumericVector M, 
                                    Rcpp::NumericVector R, 
                                    double zeta){ 
  int nm = M.size();
  int nr = R.size();
  if (nr!=nm) Rcpp::stop("Arguments have different lengths");
  
  Rcpp::NumericVector out (nm);

  for(int i = 0; i<nm; i++){
    if(R(i)/M(i) <= kappa){
      out(i) = 0;
    } else if (R(i)/M(i) >= target_condition){
      out(i) = 1;
    } else {
      out(i) = ((1-zeta)*(R(i) - kappa * M(i)))/(M(i) * (target_condition - kappa) - zeta * (R(i) - kappa * M(i)));
    }
  }
  
  return out;
} 

//' Total mass of mammary glands
//' @name mammary_mass
 //' @description Predicts mammary mass from total mass in females
 //' @param M Total body mass (kg)
 //' @note Oftedal (1997)
 // [[Rcpp::export]]
 
 double mammary_mass(double M){ 
   return std::pow(10, 0.902 * std::log10(M) - 1.965);
 }  

//' Milk production rate
//' @name milk_production
 //' @description Predicts the amount of milk produced by unit time from mammary mass
 //' @param m Mass of mammary glands (kg)
 //' @return Milk yield/production rate (kg/s)
 //' @note Equation taken from Hanwell & Peaker ()1997)
 // [[Rcpp::export]]
 
 double milk_production(double m){ 
   // return 0.0835 * std::pow(m, 1.965);
   return (1.67 * std::pow(m, 0.95));
 } 

//' Resting metabolic rate
//' @name RMR
 //' @description Predicts the resting metabolic rate of an animal from its mass 
 //' using the allometric relationship proposed proposed by Williams & Maresh (2015)
 //' @param M Total body mass (kg)
 //' @param phi Scalar constant for immature animals
 //' @note The RMR is scaled up in immature animals (Fortune et al., 2013; Rechsteiner et al., 2013)
 //' in order to account for the elevated metabolic demand associated with active growth (Lavigne et al., 1986)
 // [[Rcpp::export]]
 
 double RMR(double M){ 
  
   return(581 * std::pow(M, 0.68) / 1000);
   
   // double out;
   // 
   // if(function_ID == 0) {  // Kleiber
   // 
   //   out = phi * 3.771 * std::pow(M, 0.75) * 86400 / 1000000; // Watts (J/s) to MJ/day
   //   
   // } else if(function_ID == 1){ // Williams & Maresh (2015)
   //   
   //   out = phi * 581 * std::pow(M, 0.68) / 1000;
   //   
   // } else { // George et al. (2021)
   //   
   //   out = phi * 3.693 * std::pow(M, 0.667) * 86400 / 1000000;
   //   
   // }
   // 
   // return out;
 } 

//' Costs of locomotion
 // [[Rcpp::export]]
 double locomotor_costs(double mass,
                        double distance,
                        double strokerate_foraging,
                        double strokerate,
                        double glide_foraging,
                        double glide,
                        double t_feed,
                        double t_activ,
                        double scalar){

   // Check that vector arguments are of the same size
   // if (t_activ.size()!=scalars.size()) Rcpp::stop("Arguments have different lengths");

   // Total number of strokes per day
   double tot_strokes = scalar * (strokerate_foraging * t_feed * (1 - glide_foraging) + strokerate * t_activ * (1 - glide));

   // Cost per stroke
   // 0.4 J per kg per m (Williams et al. 1999)
   // Distance traveled in m
   // Average stroke rate during routine diving (0.11375 Hz) 
   double J_stroke = 0.4 * distance/(0.11375*86400);

   // Total cost per day (MJ)
   return J_stroke * mass * tot_strokes/1000000;
 }

// ' Costs of locomotion
// ' @name locomotor_costs
// ' @description Predicts total locomotor costs from the mass-specific, stroke-based allometric
// ' relationships proposed by Williams et al. (2017).
// ' @param M Total body mass (kg)
// ' @param Sroutine Routine stroke frequency
// ' @param delta Time spent in activity a (s)
// ' @param phi Locomotory cost scalars (d.u.)
// ' @note We only consider stroke rates associated with routine swimming. Williams et al. (2017) do also
// ' present equations for performance (maximum aerobic) swimming, which likely capture startle/flight responses.
// ' To our knowledge, no data on the swimming kinematics of right whales exhibiting this behavior exist at present.
// ' The equations used were derived from data on odontocetes up to 3,000 kg in weight,
// ' including killer whales (Orcinus orca). It is assumed that these relationships hold when
// ' extrapolating to the larger body mass range exhibited by North Atlantic right whales.
// ' @return Total daily locomotor costs (kJ)

// double locomotor_costs(double mass, 
//                        double strokerate_foraging,
//                        double strokerate, 
//                        double glide_foraging,
//                        double glide,
//                        double t_feed,
//                        double t_activ, 
//                        double scalar){ 
//   
//   // Check that vector arguments are of the same size
//   // if (t_activ.size()!=scalars.size()) Rcpp::stop("Arguments have different lengths");
//   
//   // Total number of strokes per day
//   double tot_strokes = scalar * (strokerate_foraging * t_feed * (1 - glide_foraging) + strokerate * t_activ * (1 - glide));
//   
//   // Cost per stroke (Williams et al. 2017) (J / kg / stroke)
//   double J_stroke = (1.46 + 0.0005 * mass);
//   // double J_stroke = 3.38;
//   
//   // Total cost per day (MJ)
//   return J_stroke * mass * tot_strokes/1000000;
//   
//   // int n = t_activ.size();
//   // Rcpp::NumericVector scalar (n);
//   // 
//   // for(int i = 0; i < n; ++i){
//   //   scalar[i] += delta[i] * phi[i];
//   // }
//   
//   // return sum(scalar) * M * strokerate * 
// } 



// // [[Rcpp::export]]
// double testmin(double n){
//   
//   Rcpp::NumericVector distances(n);
//   for(int i = 0; i<n;i++){
//     distances(i) = R::rnorm(150, 20);
//     std::cout << distances(i) << std::endl;
//   }
//   
//   return(Rcpp::which_min(distances));
//   
// }

// //' Energetic cost of fetal growth during pregnancy
//  //' @param M_muscle Mass of muscles in fetus (kg)
//  //' @param M_viscera Mass of viscera in fetus (kg)
//  //' @param M_bones Mass of bones in fetus (kg)
//  //' @param M_blubber Mass of blubber in fetus (kg)
//  //' @param rho_lipid Energy density of lipids (kJ/kg)
//  //' @param rho_protein Energy density of protein (kJ/kg)
//  //' @param P_lip_muscle Lipid concentration in muscle (\%)
//  //' @param P_pro_muscle Protein concentration in muscle (\%)
//  //' @param P_lip_viscera Lipid concentration in viscera (\%)
//  //' @param P_pro_viscera Protein concentration in viscera (\%)
//  //' @param P_lip_bones Lipid concentration in bones (\%)
//  //' @param P_pro_bones Protein concentration in bones (\%)
//  //' @param P_lip_blubber Lipid concentration in blubber (\%)
//  //' @param P_pro_blubber Protein concentration in blubber (\%)
//  // [[Rcpp::export]]
//  
//  double fetal_growth(double M_muscle, double M_viscera, 
//                      double M_bones, double M_blubber,
//                      double rho_lipid, double rho_protein, 
//                      double P_lip_muscle, double P_pro_muscle,
//                      double P_lip_viscera, double P_pro_viscera, 
//                      double P_lip_bones, double P_pro_bones,
//                      double P_lip_blubber, double P_pro_blubber){ 
//    
//    double fg_m, fg_v, fg_b, fg_bl;
//    fg_m = M_muscle * (rho_lipid * P_lip_muscle + rho_protein * P_pro_muscle);
//    fg_v = M_viscera * (rho_lipid * P_lip_viscera + rho_protein * P_pro_viscera);
//    fg_b = M_bones * (rho_lipid * P_lip_bones + rho_protein * P_pro_bones);
//    fg_bl = M_blubber * (rho_lipid * P_lip_blubber + rho_protein * P_pro_blubber);
//    
//    return fg_m + fg_v + fg_b + fg_bl;
//    
//  }  

//' Energetic cost of placental maintenance during pregnancy (vectorised)
//' @name placental_maintenance
//' @param G Energetic cost of fetal growth (kJ)
// [[Rcpp::export]]
 Rcpp::NumericVector placental_maintenance_vec(Rcpp::NumericVector G){ 
   int n = G.size();
   Rcpp::NumericVector out(n);
   for(int i = 0; i < n; i++){
     out(i) = (G(i)/0.807)*(1-0.807);
   }
   return(out);
 } 

//' Energetic cost of placental maintenance during pregnancy
//' @name placental_maintenance
//' @param G Energetic cost of fetal growth (kJ)
// [[Rcpp::export]]
 
 double placental_maintenance(double G){ 
   return (G/0.807)*(1-0.807);
 }  

//' Heat increment of gestation
//' @name heat_gestation
//' @param birth_mass Birth mass of the fetus (kg)
//' @param delta_m Daily growth rate of the fetus (kg/day)
// [[Rcpp::export]]
 
 double heat_gestation(double birth_mass, double delta_m){ 
   return 18409.6 * std::pow(birth_mass, 1.2) * (delta_m/birth_mass);
 }  

//' Heat increment of gestation (vectorised)
//' @name heat_gestation
//' @param birth_mass Birth mass of the fetus (kg)
//' @param delta_m Daily growth rate of the fetus (kg/day)
// [[Rcpp::export]]
 
 Rcpp::NumericVector heat_gestation_vec(double birth_mass,
                                        Rcpp::NumericVector delta_m){ 
   int n = delta_m.size();
   Rcpp::NumericVector out(n);
   for(int i = 0; i < n; i++){
    out(i) = 18409.6 * std::pow(birth_mass, 1.2) * (delta_m(i)/birth_mass);
   }
   return(out);
 }  


//' Fetal tissue mass
//' @name fetal_tissue_mass
//' @param P_b Proportion of the body volume comprised of tissue b
//' @param L Length of the fetus (m)
//' @note This relationship only applies to muscles, bones, and viscera
// [[Rcpp::export]]
 
 double fetal_tissue_mass(double P_b, double L){ 
   return 1000 * P_b * std::exp(-4.115 + 3.016 * std::log(L));
 }  

//' Fetal tissue mass (vectorised)
//' @name fetal_tissue_mass
//' @param P_b Proportion of the body volume comprised of tissue b
//' @param L Length of the fetus (m)
//' @note This relationship only applies to muscles, bones, and viscera
// [[Rcpp::export]]
 
 Rcpp::NumericVector fetal_tissue_mass_vec(double P_b, 
                                           Rcpp::NumericVector L){ 
   int n = L.size();
   Rcpp::NumericVector out(n);
   for(int i = 0; i < n; i++){
     out(i) = 1000 * P_b * std::exp(-4.115 + 3.016 * std::log(L(i)));
   }
   return(out);
 }  


//' Fetal blubber mass
//' @name fetal_blubber_mass
//' @param L Length of the fetus (m)
//' @param M_muscle Mass of muscles in the fetus (kg)
//' @param M_viscera Mass of viscera in the fetus (kg)
//' @param M_bones Mass of bone tissues in the fetus (kg)
//' @param D_blubber Average blubber density (kg/m^3)
//' @param D_muscle Average muscle density (kg/m^3)
//' @param D_viscera Average density of viscera (kg/m^3)
//' @param D_bones Average bone density (kg/m^3)
//' @note The original equation from Christiansen et al. (2022) (DOI: 10.3354/meps14009) includes an additional term
//' designed to account for the calf's body condition at birth. However, Christiansen et al. rely on a metric of body condition (BC)
//' that differs from, and is not readily comparable to, ours. Here, we assume that BC = 0, which corresponds to an animal of average
//' body condition. 
// [[Rcpp::export]]
 
 double fetal_blubber_mass(double L,
                           double BC,
                           double M_muscle, 
                           double M_viscera, 
                           double M_bones,
                           double D_blubber,
                           double D_muscle, 
                           double D_viscera, 
                           double D_bones){ 
   return D_blubber * (std::exp(-4.115 + 3.016 * std::log(L)) * (1+BC) - (M_muscle/D_muscle) - (M_viscera/D_viscera) - (M_bones/D_bones));
 }  

//' Fetal blubber mass (vectorised)
//' @name fetal_blubber_mass
//' @param L Length of the fetus (m)
//' @param M_muscle Mass of muscles in the fetus (kg)
//' @param M_viscera Mass of viscera in the fetus (kg)
//' @param M_bones Mass of bone tissues in the fetus (kg)
//' @param D_blubber Average blubber density (kg/m^3)
//' @param D_muscle Average muscle density (kg/m^3)
//' @param D_viscera Average density of viscera (kg/m^3)
//' @param D_bones Average bone density (kg/m^3)
//' @note The original equation from Christiansen et al. (2022) (DOI: 10.3354/meps14009) includes an additional term
//' designed to account for the calf's body condition at birth. However, Christiansen et al. rely on a metric of body condition (BC)
//' that differs from, and is not readily comparable to, ours. Here, we assume that BC = 0, which corresponds to an animal of average
//' body condition. 
// [[Rcpp::export]]
 
 Rcpp::NumericVector fetal_blubber_mass_vec(Rcpp::NumericVector L,
                                            double BC,
                                            Rcpp::NumericVector M_muscle, 
                                            Rcpp::NumericVector M_viscera, 
                                            Rcpp::NumericVector M_bones,
                                            double D_blubber,
                                            double D_muscle, 
                                            double D_viscera, 
                                            double D_bones){ 
   
   int n = L.size();
   int mn = M_muscle.size();
   int mv = M_viscera.size();
   int mb = M_bones.size();
   
   // Check that all vectors have the same size
   std::array<int,4> sizes = {n, mn, mv, mb};
   int allequal = std::equal(sizes.begin() + 1, sizes.end(), sizes.begin());
   if(allequal == 0) Rcpp::stop("Arguments have different lengths");
   
   Rcpp::NumericVector out(n);
   for(int i = 0; i < n; i++){
     out(i) = D_blubber * (std::exp(-4.115 + 3.016 * std::log(L(i))) * 
       (1+BC) - (M_muscle(i)/D_muscle) - (M_viscera(i)/D_viscera) - (M_bones(i)/D_bones));
   }
   
   return(out);
 }  

//' Fetal mass 
//' @name fetal_mass
//' @param days_to_birth Number of days until birth, assuming a 365-day gestation period (d)
//' @param mother_length Body length of the mother (m)
//' @param bbc Body condition, as defined by Christiansen et al. (2022). Defaults to 0 for an individual of average condition.
//' @param body_density Average body density (kg/m^3)
//' @note In this parameterization, birth corresponds to t=0 and conception corresponds to t=-365
// [[Rcpp::export]]
 
 double fetal_mass(int days_to_birth, double mother_length, int bbc = 0, double body_density = 805.07){
   return (std::exp(-4.115 + 3.016 * std::log((std::exp(-1.050376 + 0.007685 * days_to_birth) + (0.021165/365)*days_to_birth)*mother_length))*(1+bbc))*body_density;
 }

//' Fetal mass 
//' @name fetal_mass (vectorised)
//' @param days_to_birth Number of days until birth, assuming a 365-day gestation period (d)
//' @param mother_length Body length of the mother (m)
//' @param bbc Body condition, as defined by Christiansen et al. (2022). Defaults to 0 for an individual of average condition.
//' @param body_density Average body density (kg/m^3)
//' @note In this parameterization, birth corresponds to t=0 and conception corresponds to t=-365
// [[Rcpp::export]]
 
 Rcpp::NumericVector fetal_mass_vec(Rcpp::NumericVector days_to_birth, 
                                    Rcpp::NumericVector mother_length, 
                                    int bbc = 0, 
                                    double body_density = 805.07){
   int n = days_to_birth.size();
   int m = mother_length.size();
   double ML;
   
   if(m > 1){
     // Check that all vectors have the same size
     std::array<int,2> sizes = {n, m};
     int allequal = std::equal(sizes.begin() + 1, sizes.end(), sizes.begin());
     if(allequal == 0) Rcpp::stop("Arguments have different lengths");
   }
   
   Rcpp::NumericVector out(n);
   for(int i = 0; i < n; i++){
     if(m > 1){
       ML = mother_length(i);
     } else {
       ML = mother_length(0);
     }
     out(i) = (std::exp(-4.115 + 3.016 * std::log((std::exp(-1.050376 + 0.007685 * days_to_birth(i)) + 
       (0.021165/365)*days_to_birth(i))*ML))*(1+bbc))*body_density;
   }
   return(out);
 }

//' Fetal length 
//' @name fetal_length
//' @param days_to_birth Number of days until birth, assuming a 365-day gestation period (d)
//' @param mother_length Body length of the mother (m)
//' @note In this parameterization, birth corresponds to t=0 and conception corresponds to t=-365
// [[Rcpp::export]]
 
 double fetal_length(int days_to_birth, double mother_length){
   return (std::exp(-1.050376 + 0.007685 * days_to_birth) + (0.021165/365) * days_to_birth) * mother_length;
 }

//' Fetal length (vectorised)
//' @name fetal_length
//' @param days_to_birth Number of days until birth, assuming a 365-day gestation period (d)
//' @param mother_length Body length of the mother (m)
//' @note In this parameterization, birth corresponds to t=0 and conception corresponds to t=-365
// [[Rcpp::export]]
 
Rcpp::NumericVector fetal_length_vec(Rcpp::NumericVector days_to_birth, 
                                     Rcpp::NumericVector mother_length){
  int n = days_to_birth.size();
  int m = mother_length.size();
  double ML;
  
  // Check that all vectors have the same size
  if(m > 1){
  std::array<int,2> sizes = {n, m};
  int allequal = std::equal(sizes.begin() + 1, sizes.end(), sizes.begin());
  if(allequal == 0) Rcpp::stop("Arguments have different lengths");
  }
  
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; i++){
    if(m > 1){
      ML = mother_length(i);
    } else {
      ML = mother_length(0);
    }
    out(i) = (std::exp(-1.050376 + 0.007685 * days_to_birth(i)) + (0.021165/365) * days_to_birth(i)) * ML;
  }
  return(out);
}


// // [[Rcpp::export]]
// double fetal_growth(double L0, double L1, 
//              double Prop_m, double Prop_v, double Prop_b, 
//              double Dens_m, double Dens_v, double Dens_bo, double Dens_bl,
//              double Lip_m, double Lip_v, double Lip_bo, double Lip_bl,
//              double Pro_m, double Pro_v, double Pro_bo, double Pro_bl,
//              double EL, double EP){
//   
//   // Mass of muscles in fetus -- (kg)
//   double mass_muscle = fetal_tissue_mass(Prop_m, L0);
//   double mass_muscle_next = fetal_tissue_mass(Prop_m, L1);
//   
//   // Mass of viscera in fetus -- (kg)
//   double mass_viscera = fetal_tissue_mass(Prop_v, L0);
//   double mass_viscera_next = fetal_tissue_mass(Prop_v, L1);
//   
//   // Mass of bones in fetus -- (kg)
//   double mass_bones = fetal_tissue_mass(Prop_b, L0);
//   double mass_bones_next = fetal_tissue_mass(Prop_b, L1);
//   
//   // Mass of blubber in fetus -- (kg)
//   double mass_blubber = fetal_blubber_mass(L0, mass_muscle, mass_viscera, mass_bones, Dens_bl, Dens_m, Dens_v, Dens_bo);
//   double mass_blubber_next = fetal_blubber_mass(L1, mass_muscle_next, mass_viscera_next, mass_bones_next, Dens_bl, Dens_m, Dens_v, Dens_bo);
//   
//   double FG0 = mass_muscle * (EL * Lip_m[current_animal] + EP * Pro_m[current_animal]) + 
//     mass_viscera * (EL * Lip_v[current_animal] + EP * Pro_v[current_animal]) +
//     mass_bones * (EL * Lip_bo[current_animal] + EP * Pro_bo) +
//     mass_blubber * (EL * Lip_bl[current_animal] +  EP * Pro_bl[current_animal]);
//   
//   double FG1 = mass_muscle_next * (EL * Lip_m[current_animal] + EP * Pro_m[current_animal]) + 
//     mass_viscera_next * (EL * Lip_v[current_animal] + EP * Pro_v[current_animal]) +
//     mass_bones_next * (EL * Lip_bo[current_animal] + EP * Pro_bo) +
//     mass_blubber_next * (EL * Lip_bl[current_animal] +  EP * Pro_bl[current_animal]);
//   
//   return(FG1 - FG0);
// }


// 
// //' Energetic cost of growth
// //' @name growth_cost
// //' @param delta_m Body mass growth increment (kg/day)
// //' @param prop_blubber Proportion of the body that is blubber (\%)
// //' @param prop_water Proportion of lean body mass that is water (\%)
// //' @param P_lipid_blubber Proportion of blubber that is lipid (\%)
// //' @param rho_lipid Energy density of lipids (kJ/kg)
// //' @param rho_protein Energy density of protein (kJ/kg)
// //' @param D_lipid Efficiency of deposition of lipids (\%)
// //' @param D_protein Efficiency of deposition of protein (\%)
// // [[Rcpp::export]]
//  
//  double growth_cost_old(double delta_m, double prop_blubber, double prop_water, double P_lipid_blubber,
//                     double rho_lipid, double rho_protein, double D_lipid, double D_protein){
//    
//    return delta_m * ((prop_blubber * P_lipid_blubber * rho_lipid * D_lipid) + 
//                      ((1 - prop_blubber) * (1 - prop_water) * rho_protein * D_protein));
//  }

// [[Rcpp::export]]
double growth_cost(double leanmass_increment,
                   double EDens_lipids, 
                   double EDens_protein,
                   double lipid_in_muscle,
                   double lipid_in_viscera,
                   double lipid_in_bones,
                   double protein_in_muscle,
                   double protein_in_viscera,
                   double protein_in_bones,
                   double prop_muscle,
                   double prop_viscera,
                   double prop_bones){
  
  // double lipid_in_blubber,
  // double protein_in_blubber,
  // double prop_blubber,
  
  // double cost_blubber = EDens_lipids * lipid_in_blubber + EDens_protein * protein_in_blubber;
  double cost_muscle = (EDens_lipids * lipid_in_muscle) + (EDens_protein * protein_in_muscle);
  double cost_viscera = (EDens_lipids * lipid_in_viscera )+ (EDens_protein * protein_in_viscera);
  double cost_bones = (EDens_lipids * lipid_in_bones) + (EDens_protein * protein_in_bones);
  
  // return mass_increment * (cost_blubber * prop_blubber + cost_muscle * prop_muscle + cost_viscera * prop_viscera + cost_bones * prop_bones);
  return leanmass_increment * (cost_muscle * prop_muscle + cost_viscera * prop_viscera + cost_bones * prop_bones);
}

// // [[Rcpp::export]]
// 
// double fatdeposition_cost(double fatmass_increment,
//                           double energy_density_lipids, 
//                           double energy_density_protein,
//                           double lipid_in_blubber,
//                           double protein_in_blubber){
//   
//   double cost_blubber = (energy_density_lipids * lipid_in_blubber) + (energy_density_protein * protein_in_blubber);
//   return fatmass_increment * cost_blubber;
// }


// §newind[year,"alive", 1] <- 1
// newind[year,"cohort", 1] <- 0
// newind[year,"female", 1] <- rbinom(n = 1, size = 1, prob = 0.5) # 1:1 sex ratio at birth
// newind[year, "age", 1] <- 0
// 
// l.params <- agL(0)
//   newind[year,"length", 1] <- age2length(0, l.params)
//   newind[year:dim(arr)[1],"length_a", 1] <- l.params[,1]
// newind[year:dim(arr)[1],"length_b", 1] <- l.params[,2]
// newind[year:dim(arr)[1],"length_c", 1] <- l.params[,3]
// 
// m.params <- mL(1)
//   mass <- length2mass(newind[year,"length", 1], m.params, FALSE)
//   newind[year:dim(arr)[1],"mass_a", 1] <- m.params[, 1]
// newind[year:dim(arr)[1],"mass_b", 1] <- m.params[, 2]
// newind[year,"tot_mass", 1] <- mass
// 
// newind[year,"bc",] <- start_bcondition(0)
//   
//   newind[year,"lean_mass", 1] <- mass - (newind[year,"bc",] * mass)
//   
//   newind[year,"trest",] <- 0
// newind[year,"t2calf",] <- 0
// newind[year,"min_bc",] <- 0
// newind[year,"birth",] <- 0
// newind[year,"p_surv",] <- 1

// [[Rcpp::export]]
Rcpp::NumericMatrix add_calf(int n, 
                             Rcpp::StringVector attr, 
                             Rcpp::NumericVector sex,
                             Rcpp::NumericVector nonreprod){
  
  int nattr = attr.size();
  Rcpp::NumericMatrix out(nattr,n);
  
  // Alive vs. dead
  Rcpp::NumericMatrix::Row alive = out.row(0);
  alive = alive + 1;
  
  // Sex (male, female)
  out(2,Rcpp::_) = sex;

  // Body length
  Rcpp::NumericVector L = age2length_vec(out.row(3));
  out(4,Rcpp::_) = L;
  
  // Lean_mass
  Rcpp::NumericVector M = length2mass_vec(L);
  out(6,Rcpp::_) = M;
  
  // Body condition
  Rcpp::NumericVector BC = start_bcondition_vec(out(3,Rcpp::_));
  out(7,Rcpp::_) = BC;
  
  // Total mass
  out(5,Rcpp::_) = M / (1-BC);
  
  // Probability of survival
  Rcpp::NumericMatrix::Row psurv = out.row(8);
  psurv = psurv + 1;
  
  // Sterility
  out(12,Rcpp::_) = nonreprod;
  // out(12,Rcpp::_) = Rcpp::rbinom(n, 1, 1-nonreprod);
  
  Rcpp::rownames(out) = attr;
  
  return(out);
}

// 
// Rcpp::NumericMatrix add_calf(int n, Rcpp::StringVector attr){
//   int nattr = attr.size();
//   Rcpp::NumericMatrix out(nattr,n);
//   
//   Rcpp::NumericMatrix::Row alive = out.row(0);
//   alive = alive + 1;
//   
//   out(2,Rcpp::_) = Rcpp::rbinom(n, 1, 0.5); // Sex
//   // Rcpp::NumericMatrix::Row female = out.row(2);
//   // female = Rcpp::rbinom(n, 1, 0.5);
//   
//   Eigen::MatrixXd lparams = agL_vec(n);
//   // Rcpp::NumericMatrix::Row length = out.row(4);
//   Rcpp::NumericVector L = age2length_vec(out.row(3));
//   out(4,Rcpp::_) = L;
//   
//   Eigen::MatrixXd mparams = mL(n);
//   Rcpp::NumericVector M = length2mass_vec(L);
//   out(8,Rcpp::_) = M; // Total mass
//   
//   Rcpp::NumericVector BC = start_bcondition_vec(out(3,Rcpp::_));
//   out(10,Rcpp::_) = BC; // Body condition
//   
//   Rcpp::NumericMatrix::Row psurv = out.row(13);
//   psurv = psurv + 1;
//   
//   for(int i = 0; i < n; i++){
//     for(int j = 0; j < 3; j++){
//       out(5+j, i) = lparams(i,j); // age to length parameters
//     }
//     out(9,i) = M(i) - (BC(i) * M(i)); // lean mass
//     for(int k = 0; k < 2; k++){
//       out(11+k, i) = mparams(i,k); // length to mass parameters
//     }
//   }
//   
//   return(out);
// }


// //' Tissue deposition
//  //' @description Returns the mass obtained after deposition of any surplus energy as lipid and lean tissue
//  //' @param M Start mass (kg)
//  //' @param E_net Net energy balance (kJ)
//  //' @param add_protein Percent of protein synthesis during anabolism (\%)
//  //' @param add_lipid Percent of lipid synthesis during anabolism (\%)
//  //' @param D_lipid Efficiency of deposition of lipids during anabolism (\%)
//  //' @param D_protein Efficiency of deposition of proteins during anabolism (\%)
//  //' @param break_lipid Percent of lipid breakdown during catabolism (\%)
//  //' @param break_protein Percent of protein breakdown during catabolism (\%)
//  //' @param C_lipid Efficiency of breakdown of lipids during catabolism (\%)
//  //' @param C_protein Efficiency of breakdown of protein during catabiolism (\%)
//  //' @param rho_lipid Energy density of lipids (kJ/kg)
//  //' @param rho_protein Energy density of proteins (kJ/kg)
//  // [[Rcpp::export]]
//  
//  double new_mass(double M, double E_net,
//                  double add_protein, double add_lipid,
//                  double D_lipid, double D_protein,
//                  double break_lipid, double break_protein,
//                  double C_lipid, double C_protein,
//                  double rho_lipid, double rho_protein){
//    
//    double fat_growth = E_net/rho_lipid;
//    double lean_growth = E_net/rho_protein;
//    
//    // Fat growth
//    if(E_net > 0){
//      fat_growth *= add_lipid * D_lipid;
//    } else if(E_net < 0){
//      fat_growth *= break_lipid * C_lipid;
//    }
//    
//    // Lean tissue growth
//    if(E_net > 0){
//      lean_growth *= add_protein * D_protein;
//    } else if(E_net < 0){
//      lean_growth *= break_protein * C_protein;
//    }
//    
//    return M + fat_growth + lean_growth;
//  }
// // [[Rcpp::export]]
// 
// double findminval(double num1, double num2){
//   if (num1 < num2){
//     return num1;
//   } else {
//     return num2;
//   }
// }

// // [[Rcpp::export]]
// double pbirth(float now, 
//               float enter,
//               float timespan = 60){
//   double out = (now - enter)/timespan;
//   if(out > 1) out = 1;
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector pbirth_vec(Rcpp::NumericVector now, 
//                                float enter,
//                                float timespan = 60){
//   int n = now.size();
//   Rcpp::NumericVector out(n);
//   for (int i = 0; i<n; i++){
//     out(i)= (now(i) - enter)/timespan;
//     if(out(i) > 1) out(i) = 1;
//   }
//   return out;
// }

// [[Rcpp::export]]
double pleave(float now, 
              float enter,
              float cohortID,
              float factor,
              Rcpp::NumericMatrix resid){
  
  int days_in_SEUS = now - enter;
  double resid_m = resid(cohortID - 1,0);
  double resid_sd = resid(cohortID - 1,1);
  double q = std::round(days_in_SEUS*factor);
  double out = R::pnorm5(q, resid_m, resid_sd, 1, 0);

  return out;
}

  
// // [[Rcpp::export]]
// Rcpp::NumericVector seq_cpp(double start,
//                             double end,
//                             int npts){
// 
//   Rcpp::NumericVector out(npts);
//   double interval = (end-start)/(npts + 1);
//   for (int i = 0; i<npts; i++){
//     out(i) = start + (i+1)*interval; 
//   }
//   return out;
// }

// // [[Rcpp::export]]
// double starvation_mortality(double bc,
//                             double starve_start, 
//                             double starve_threshold,
//                             double starve_scalar){
//   
//   double out;
//   
//   if(bc < starve_threshold){
//     out = 1.0L;
//   } else if(bc > starve_start) {
//     out = 0.0L;
//   } else {
//     out = starve_scalar *(starve_start * (1/bc) - 1);
//   }
//   
//   if(out < 0) out = 0;
//   if(out > 1) out = 1  ;
//   
//   return out;
// }  
//   
// // [[Rcpp::export]]
// Rcpp::NumericVector starvation_mortality_vec(Rcpp::NumericVector bc,
//                                              double starve_start, 
//                                              double starve_threshold,
//                                              double starve_scalar){
//   
//   int n = bc.size();
//   Rcpp::NumericVector out(n);
//   for (int i = 0; i<n; i++){
//     if(bc(i) < starve_threshold){
//       out(i) = 1.0L;
//     } else if(bc(i) > starve_start) {
//       out(i) = 0.0L;
//     } else {
//       out(i) = starve_scalar *(starve_start * (1/bc(i)) - 1);
//     }
//     if(out(i) < 0) out(i) = 0;
//     if(out(i) > 1) out(i) = 1;
//   }
// 
// return out;
//   }  

// [[Rcpp::export]]
double entanglement_effect(double prob_survival,
                           int severity,
                           int ndays,
                           Rcpp::NumericVector entgl_hit){
  
  double p_surv3mo = 0.0L;
  double p_surv3mo_cloglog;
  double surv_drop = entgl_hit(severity);
  
  // If the entanglement effect estimated under PCoMS is zero, then
  // return 0; otherwise, calculate the modified survival probability
  if(surv_drop != 0){
    
    // If the daily probability of a mortality event (starving) is p,
    // then the daily probability of survival is (1-p). 
    // The probability of survival over 3 months (taken to be 90 days) is then given 
    // by (1-p)^90. Incidentally, the probability of starving during the 3-month period is
    // is p_3mo = 1-((1-p)^90). 
    // It follows that p = 1-((1-p_3mo)^(1/90))
    
    p_surv3mo = std::pow(prob_survival, ndays);
    
    // Convert complementary log-log scale
    p_surv3mo_cloglog = std::log(-std::log(1.0L - p_surv3mo));
    
    if(std::isinf(p_surv3mo_cloglog)) throw std::invalid_argument("Argument out of bounds; reduce --p or change --ndays");
    
    // Add effect of entanglement estimated from PCoMS (Pirotta et al. 2023)
    p_surv3mo_cloglog += entgl_hit(severity);
    // if(severity == 0) p_surv3mo_cloglog += minor;
    // if(severity == 1) p_surv3mo_cloglog += moderate;
    // if(severity == 2) p_surv3mo_cloglog += severe;
    
    // Apply inverse cloglog
    p_surv3mo = 1.0L - std::exp(-std::exp(p_surv3mo_cloglog));
    
    // Convert back to daily probability
    // double out = std::pow(p_surv3mo, 1.0L/ndays);
    
    return p_surv3mo;
    
  } else {
    
    return prob_survival;
    
  }
}


// [[Rcpp::export]]
double survival(double age,
                int regime,
                int female){
  
  // Linden (2023). Population size estimation of North Atlantic right whales from 1990-2022
  // https://www.fisheries.noaa.gov/s3/2023-10/TM314-508-0.pdf
  
  // This function gives the following survival estimates:
  // Calves (0.5 yrs): 0.920195
  // Yearlings (1.5 yrs): 0.9340726
  // Year 2 (2.5 yrs): 0.9456794
  // Year 3 (3.5 yrs): 0.9553405
  // Year 4 (4.5 years):  0.9633499
  // Adult males: 0.9668179
  // Adult females: 0.9521659
  
  // These are comparable to Pace (2021):
  // https://repository.library.noaa.gov/view/noaa/50567
  // Calves: 0.924234
  // Yearlings: 0.937877
  // Year 2: 0.949159
  // Year 3: 0.958477
  // Year 4: 0.966181
  // Adult males: 0.972534
  // Adult females: 0.956194
  
  // And to Reed et al. (2023)
  // https://www.frontiersin.org/articles/10.3389/fmars.2022.994481/full
  // Pre-breeders: 0.96 (range: 0.92-0.98)
  // Breeders: 0.97 (range: 0.94-0.99) 
  // Calves: 0.79 (range: 0.50-0.92)
  
  // logit-linear coefficients for survival probability -- Table 1
  double b_age = 0.206; // Linear effect of age from 0–5
  double b_female = -0.381; // Effect of being an adult female
  double b_regime = -0.718; // Regime effect for years after 2010
  double b0 = 3.060; // Intercept
  
  int is_adult = 0;
  if(age > 9.5) is_adult += 1;
  
  double age_5 = std::min(age, 5.0); // Daytime foraging
  
  double y = b0 + b_female * female * is_adult + regime * b_regime + b_age * age_5;
  double out = std::exp(y)/(1+exp(y));
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector clamp(Rcpp::NumericVector v, double threshold){
  int n = v.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++){
    if(v(i) > threshold){
      out(i) = threshold;
    } else{
      out(i) = v(i);
    }
  }
  return out;
}


#endif