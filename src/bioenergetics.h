#ifndef BIOENERGETICS_H
#define BIOENERGETICS_H

#include <RcppEigen.h>
#include <random>
#include <cmath>

#include <cstdio>
#include <vector>
#include "spline.h"

// [[Rcpp::plugins(cpp11)]]

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

// [[Rcpp::export]]
double response_threshold(Rcpp::NumericVector db){
  std::random_device rd;     // Only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uniform(0,4999);
  int d = uniform(rd);
  return(db[d]);
}

// // [[Rcpp::export]]
// std::vector<double> dose_range(double lwr, double uppr, double n){
//   std::vector<double> v(n);
//   for (int i = 0; i < n; i++) v[i] = lwr + i * ((uppr - lwr)/(n-1));
//   return v;
// }

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

//' Random deviate from a truncated Normal distribution
//' 
//' @param location Location parameter
//' @param scale Scale parameter
//' @param L Lower bound
//' @param U Upper bound
// [[Rcpp::export]]
double rtnorm(double location,
               double scale,
               double L,
               double U) {
   
   double tot = R::pnorm(L, location, scale, TRUE, FALSE) + R::runif(0,1) * (R::pnorm(U, location, scale, TRUE, FALSE) - R::pnorm(L, location, scale, TRUE, FALSE));
   double out = location + scale * R::qnorm(tot, 0, 1, TRUE, FALSE);
   return out;
}

//' Initialize age
//' @name start_age
 //' @description Performs a random draw from a uniform distribution to initialize the age of simulated animals
 //' @param cohort Integer indicating which cohort the animal belongs to
 // [[Rcpp::export]]
 
 double start_age(int cohort){
   if(cohort == 0){
     return 1/365;
   } else if(cohort >= 1 && cohort <= 2){
     return R::runif(1, 9);
   } else {
     return R::runif(9, 69);
   }
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

//' Entanglement event
//' @name entanglement_event
 //' @param p_anterior Probability that the entanglement involves the anterior region of the body (mouth, head, rostrum)
 // [[Rcpp::export]]
 
 Rcpp::NumericVector entanglement_event(double p_entangled = 0,
                                        double p_anterior = 0.732, // weighted mean of entries in spreadsheet of model parameters,
                                        Rcpp::NumericVector p_severity = Rcpp::NumericVector::create(0.761, 0.155, 0.084)){       
   
   Rcpp::NumericVector out (4); // Store results
   
   // Is the animal entangled?
   // Annual entanglement rate = 0.259 -- Knowlton et al. (2012)
   // So probability of becoming entangled on each day is 0.000821 (see SI in Pirotta et al. 2023)
   int is_entangled = R::rbinom(1, p_entangled);
   
   if(is_entangled){
     
     out(0) = is_entangled;

     // Does the entanglement involve the anterior region of the body?
     out(1) = R::rbinom(1, p_anterior); 
     
     // How severe is the entanglement (0 = minor, 1 = moderate, 2 = severe)
     out(2) = multinomial(p_severity);
     
     // Duration  of entanglement (days)
     // From Pirotta et al. (2023 Oikos)
     // See entgl_durations() function
     // All long-term events involve moderate or severe injuries (Knowlton et al. 2016 Cons Biol)
     if(out(2) == 0){
       out(3) = Rf_rnbinom_mu(0.7363271, 111.3911746);
     } else if(out(2) == 1){
       out(3) = Rf_rnbinom_mu(0.7493824, 124.2184062);
     } else {
       out(3) = Rf_rnbinom_mu(0.6527522, 211.8187153);
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
//' @name start_percfat
 //' @description Performs a random draw from a beta distribution to initialize 
 //' the body condition of simulated animals, taken as the ration of fat mass to total mass.
 //' @param age Age in years
 //' @param shape1 First shape parameter of the beta distribution
 //' @param shape2 Second shape parameter of the beta distribution
 // [[Rcpp::export]]
 
 long double start_percfat(double cohort, int shape1 = 6, int shape2 = 20){
   long double bc;
   if(cohort == 0){
     bc = 0.060L;
   } else {
     bc = R::rbeta(shape1, shape2);
   }
   return bc;
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

// [[Rcpp::export]]
Eigen::MatrixXd create_mat(){
  return Eigen::MatrixXd(2,1);
}

// ' Length to mass conversion
// ' @description Calculates body mass from body length according to the logarithmic
// ' mass-at-length relationship of Fortune et al. (2021)
// ' @param L Input body length (m)
// ' @param param Parameters of the mass-at-length relationship, as returned by mL
// ' @param lean If TRUE, returns lean mass only.
// ' @return Estimated mass (kg)
// [[Rcpp::export]]
double length2mass(double L,
                   Eigen::MatrixXd param,
                   bool lean){
  
  double a = param(0,0);
  double b = param(0,1);
  double L_cm = L * 100;
  double perc = 1;
  if(lean) perc = 0.73;
  return perc * std::pow(10, a)*std::pow(L_cm,b); // m to cm
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


// ' Parameters of the length-at-age relationship
// ' @description Generates vectors of coefficients from the Gompertz growth
//  curve describing the change in length as a function of age relationship, as described in Fortune et al. (2021)
// ' @return Matrix of coefficients
// [[Rcpp::export]]
Eigen::MatrixXd agL(double age, int n = 1){
  
  double a, b, c; // Parameters of the Gompertz growth curves        
  Eigen::MatrixXd out(n,3);
  
  for(int i = 0; i < n; i++){
    
    if(age <= 0.79){
      a = R::rnorm(1067.19, 19.67);
      b = R::rnorm(-0.93, 0.08);
      c = R::rnorm(-3.11, 0.28);
    } else {
      a = R::rnorm(1362.755, 22.88); 
      b = R::rnorm(-0.37, 0.03); 
      c = R::rnorm(-0.18, 0.03); 
    }
    
    out(i,0) = a;
    out(i,1) = b;
    out(i,2) = c;
    
  }
  return(out);
}

// ' Parameters of the mass-at-length relationship
// ' @description Generates pairs of coefficients (intercept and slope) from the logarithmic mass-at-length relationship
//   described in Fortune et al. (2021)
// ' @return Matrix of coefficients
// [[Rcpp::export]]
Eigen::MatrixXd mL(int n = 1){
  
  Eigen::MatrixXd out = Eigen::MatrixXd(n,2);
  
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
  return(out);
}

// // [[Rcpp::export]]
// int separation_duration(){
//   // The returned value is automatically converted to an integer by defining a function of type <int>
//   return rtnorm(0, 6.152344, 1, 23);
// } 

//' Incidence of foraging behavior
//' @name feeding_threshold
 //' @description Determines whether the prey concentration encountered by an animal
 //' is sufficient to support foraging
 //' @param min_prey Minimum prey density threshold that triggers foraging (\ifelse{html}{\out{copepods/m<sup>3</sup>}}{\eqn{copepods/m^3})
 //' @param D Prey concentration (\ifelse{html}{\out{copepods/m<sup>3</sup>}}{\eqn{copepods/m^3})
 // [[Rcpp::export]]
 
 double feeding_threshold(double min_prey, double D){
   return std::round(1/(1 + exp(min_prey - D)));
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

//' Variation in feeding/suckling effort with body condition
//'  @name scale_effort
 //' @description Calculates feeding/suckling effort given an animal's body condition
 //' @param A Horizontal asymptote when x tends to -Inf
 //' @param D Horizontal asymptote when x tends to Inf
 //' @param B Steepness of the transition between the two asymptotes
 //' @param C Location parameter
 //' @param S Curve asymmetry (the curve us symmetric when S = 1)
 //' @note An increase in feeding effort during periods of compromised health is simulated as a 
 //' behavioral mechanism to compensate for periods of disturbance, should adequate resources be available. 
 //' This is achieved by using a 5-parameter logistic function (bounded by 0 and 1) of body condition, which 
 //' is designed to reduce to 0 at the animal’s target body condition and ensures that fat reserves do not grow out
 //' of bounds under favorable conditions. While the possibility of early weaning resulting from good body 
 //' condition of both mother and calf is not considered explicitly, suckling effort is made to vary as a 
 //' function of the calf’s reserve mass, in the same way that foraging effort does in juveniles/adults.
 //' @return A scalar on feeding/suckling effort
 // [[Rcpp::export]]
 
 double scale_effort(double x, double A, double D, double B, double C, double S){
   return A + (D-A) / std::pow(1 + std::exp(B*(C-x)), S);
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
 //' @return Estimated gape area (\ifelse{html}{\out{m<sup>2</sup>}}{\eqn{m^2})
 // [[Rcpp::export]]
 
 double gape_size(double L, double omega, double alpha){ 
   double a = deg2radians(alpha);
   return a * (std::pow(omega,2)/4 + std::pow((0.2077*L - 1.095),2))/2;
 }

// //' Active feeding during foraging
// //' @description Determines the proportion of time spent actively feeding while in foraging mode
// //' @param F Filtration rate (\ifelse{html}{\out{m<sup>3</sup>/s}}{\eqn{m^3/s})
// //' @param tau_clear Time required to clear the forestomach (s)
// //' @param stomach_size Forestomach capacity (\ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3})
// //' @return Proportion of active feeding (\%)
// // [[Rcpp::export]]
//  
//  double active_feed(double F, double tau_clear, double stomach_size){ 
//    return 1 / (1 + (F * tau_clear / stomach_size));
//  }

//' Energy content of the prey
//' @name Econtent_cop
 //' @description Calculates the amount of energy that a whale can obtain from an average individual prey
 //' @param m Average mass of the prey (g/cop)
 //' @param rho Energy density of the prey (J/g) 
 //' @param E_digest Digestive efficiency (fecal and urinary, \%)
 //' @param E_hif Metabolizing efficiency (1-heat increment of feeding, \%)
 // [[Rcpp::export]]
 
 double Econtent_cop(double m, double rho, double E_digest, double E_hif){ 
   return m * rho * E_digest * E_hif;
 }

//' Filtration rate
//' @name filtration_rate
 //' @description Calculates the volume of water filtered per unit time
 //' @param A Area of the mouth gape (\ifelse{html}{\out{m<sup>2</sup>}}{\eqn{m^2})
 //' @param lambda_gape Percent reduction in the gape during an entanglement event (\%)
 //' @param V Swimming speed while filtering (m/s)
 //' @param E_capt Capture efficiency (\%)
 //' @param lambda_capt Percent reduction in capture efficiency during an entanglement event (\%)
 //' @return Estimated filtration rate (\ifelse{html}{\out{m<sup>3</sup>/s}}{\eqn{m^3/s})
 // [[Rcpp::export]]
 
 double filtration_rate(double A, double lambda_gape, double V, double E_capt, double lambda_capt){ 
   return A * (1 - lambda_gape) * V * E_capt * (1 - lambda_capt);
 }


//' Milk ingestion rate
//' @name milk_ingestion
 //' @description Calculates the maximum volume of milk that a calf can ingest per unit time
 //' @param E_milk Milk transfer/assimilation efficiency (\%)
 //' @param M Scalar on milk provisioning by the mother (d.u.)
 //' @param E_gland Mammary gland efficiency (\%)
 //' @param mu_female Total mass of mammary glands (kg)
 //' @param delta_female Milk production rate by the mother (\ifelse{html}{\out{m<sup>2</sup>}}{\eqn{m^2}})
 //' @param D Density of milk (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 // [[Rcpp::export]]
 
 double milk_ingestion(double E_milk, double M, double E_gland, double mu_female, double delta_female, double D){ 
   return E_milk * M * E_gland * (mu_female * delta_female/D);
 }

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

//' Total mass of mammary glands
//' @name mammary_mass
 //' @description Predicts mammary mass from total mass in females
 //' @param M Total body mass (kg)
 // [[Rcpp::export]]
 
 double mammary_mass(double M){ 
   return std::pow(10, 0.902 * std::log10(M) - 1.965);
 }  

//' Milk production rate
//' @name milk_production
 //' @description Predicts the amount of milk produced by unit time from mammary mass
 //' @param m Mass of mammary glands (kg)
 //' @return Milk yield/production rate (kg/s)
 // [[Rcpp::export]]
 
 double milk_production(double m){ 
   // return 0.0835 * std::pow(m, 1.965);
   return (1.67 * std::pow(m, 0.95))/86400;
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
 
 double RMR(int function_ID, double M, double phi){ 
  
   double out;
   
   if(function_ID == 0) {  // Kleiber

     out = phi * 3.771 * std::pow(M, 0.75) * 86400 / 1000000; // Watts (J/s) to MJ/day
     
   } else if(function_ID == 1){ // Gavrilchuk et al. (2021)
     
     out = phi * 581 * std::pow(M, 0.68) / 1000;
     
   } else { // George et al. (2021)
     
     out = phi * 3.693 * std::pow(M, 0.667) * 86400 / 1000000;
     
   }
   return out;
 } 

//' Costs of locomotion
//' @name locomotor_costs
 //' @description Predicts total locomotor costs from the mass-specific, stroke-based allometric
 //' relationships proposed by Williams et al. (2017).
 //' @param M Total body mass (kg)
 //' @param Sroutine Routine stroke frequency
 //' @param delta Time spent in activity a (s)
 //' @param phi Locomotory cost scalars (d.u.)
 //' @note We only consider stroke rates associated with routine swimming. Williams et al. (2017) do also 
 //' present equations for performance (maximum aerobic) swimming, which likely capture startle/flight responses.
 //' To our knowledge, no data on the swimming kinematics of right whales exhibiting this behavior exist at present.
 //' The equations used were derived from data on odontocetes up to 3,000 kg in weight,
 //' including killer whales (Orcinus orca). It is assumed that these relationships hold when
 //' extrapolating to the larger body mass range exhibited by North Atlantic right whales.
 //' @return Total daily locomotor costs (kJ)
 // [[Rcpp::export]]
 double locomotor_costs(double M, 
                        double Sroutine, 
                        Rcpp::NumericVector delta, 
                        Rcpp::NumericVector phi){ 
   
   if (delta.size()!=phi.size()) Rcpp::stop("Arguments have different lengths");
   
   int n = delta.size();
   Rcpp::NumericVector scalar (n);
   
   for(int i = 0; i < n; ++i){
     scalar[i] += delta[i] * phi[i];
   }
   
   return sum(scalar) * M * Sroutine * (1.46 + 0.0005 * M);
 } 

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

//' Fetal tissue mass
//' @name fetal_tissue_mass
 //' @param P_b Proportion of the body volume comprised of tissue b
 //' @param L Length of the fetus (m)
 //' @note This relationship only applies to muscles, bones, and viscera
 // [[Rcpp::export]]
 
 double fetal_tissue_mass(double P_b, double L){ 
   return 1000 * P_b * std::exp(-4.115 + 3.016 * std::log(L));
 }  

//' Fetal blubber mass
//' @name fetal_blubber_mass
 //' @param L Length of the fetus (m)
 //' @param M_muscle Mass of muscles in the fetus (kg)
 //' @param M_viscera Mass of viscera in the fetus (kg)
 //' @param M_bones Mass of bone tissues in the fetus (kg)
 //' @param D_blubber Average blubber density (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 //' @param D_muscle Average muscle density (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 //' @param D_viscera Average density of viscera (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 //' @param D_bones Average bone density (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 //' @note The original equation from Christiansen et al. (2022) (DOI: 10.3354/meps14009) includes an additional term
 //' designed to account for the calf's body condition at birth. However, Christiansen et al. rely on a metric of body condition (BC)
 //' that differs from, and is not readily comparable to, ours. Here, we assume that BC = 0, which corresponds to an animal of average
 //' body condition. 
 // [[Rcpp::export]]
 
 double fetal_blubber_mass(double L,
                           double M_muscle, double M_viscera, double M_bones,
                           double D_blubber, double D_muscle, double D_viscera, double D_bones){ 
   return D_blubber * (std::exp(-4.115 + 3.016 * std::log(L)) - 
                       (M_muscle/D_muscle) - (M_viscera/D_viscera) - (M_bones/D_bones));
 }  

//' Fetal mass 
//' @name fetal_mass
 //' @param days_to_birth Number of days until birth, assuming a 365-day gestation period (d)
 //' @param mother_length Body length of the mother (m)
 //' @param bbc Body condition, as defined by Christiansen et al. (2022). Defaults to 0 for an individual of average condition.
 //' @param body_density Average body density (\ifelse{html}{\out{kg/m<sup>3</sup>}}{\eqn{kg/m^3})
 //' @note In this parameterization, birth corresponds to t=0 and conception corresponds to t=-365
 // [[Rcpp::export]]
 
 double fetal_mass(int days_to_birth, double mother_length, int bbc = 0, double body_density = 805.07){
   return (std::exp(-4.115 + 3.016 * std::log((std::exp(-1.050376 + 0.007685 * days_to_birth) + (0.021165/365)*days_to_birth)*mother_length))*(1+bbc))*body_density;
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

//' Energetic cost of growth
//' @name growth_cost
 //' @param delta_m Body mass growth increment (kg/day)
 //' @param prop_blubber Proportion of the body that is blubber (\%)
 //' @param prop_water Proportion of lean body mass that is water (\%)
 //' @param P_lipid_blubber Proportion of blubber that is lipid (\%)
 //' @param rho_lipid Energy density of lipids (kJ/kg)
 //' @param rho_protein Energy density of protein (kJ/kg)
 //' @param D_lipid Efficiency of deposition of lipids (\%)
 //' @param D_protein Efficiency of deposition of protein (\%)
 // [[Rcpp::export]]
 
 double growth_cost(double delta_m, double prop_blubber, double prop_water, double P_lipid_blubber,
                    double rho_lipid, double rho_protein, double D_lipid, double D_protein){
   
   return delta_m * ((prop_blubber * P_lipid_blubber * rho_lipid * D_lipid) + 
                     ((1 - prop_blubber) * (1 - prop_water) * rho_protein * D_protein));
 }

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

#endif