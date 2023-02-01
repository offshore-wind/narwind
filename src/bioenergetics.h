#ifndef BIOENERGETICS_H
#define BIOENERGETICS_H

#include <RcppEigen.h>
#include <random>
#include <cmath>

// [[Rcpp::plugins(cpp11)]]


//' Random deviate from a truncated Normal distribution
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
 
//' Initialize entanglement load
//' @description Performs a random draw from a discrete uniform distribution to initialize the entanglement state of simulated animals
// [[Rcpp::export]]
 
 int start_entangled(){                     
   std::random_device rd;   // Used to obtain a seed for the random number engine
   std::mt19937 gen(rd());  // Pseudo-random number generator
   std::uniform_int_distribution<> distrib(0, 2);
   int rand = distrib(gen);
   return rand;
 }

//' Mouth-width-to-length ratio
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
 //' @description Performs a random draw from a beta distribution to initialize 
 //' the body condition of simulated animals, taken as the ration of fat mass to total mass.
 //' @param shape1 First shape parameter of the beta distribution
 //' @param shape2 Second shape parameter of the beta distribution
 // [[Rcpp::export]]
 
 double start_percfat(int shape1 = 6, int shape2 = 20){
   return R::rbeta(shape1, shape2);
 }

//' Age to length conversion
//' @description Calculates body length from age according to the length-at-age relationship described
//' in Fortune et al. (2021)
//' @param age Age of the animal (years)
//' @return Estimated length (m)
// [[Rcpp::export]]
 
 double age2length(double age){             
   double a, b, c;// Parameters of the Gompertz growth curves                          
   if(age <= 0.79){
     a = R::rnorm(1067.353, 20.479);
     b = R::rnorm(-0.923, 0.058);
     c = R::rnorm(-3.075, 0.315);
   } else {
     a = R::rnorm(1360.675, 19.501); 
     b = R::rnorm(-0.361, 0.023); 
     c = R::rnorm(-0.166, 0.026); 
   }
   return a * exp(b*exp(c*age)) / 100; // cm to m
 }
 
//' Length to mass conversion
//' @description Calculates body mass from body length according to the logarithmic 
//' mass-at-length relationship of Fortune et al. (2021)
//' @param L Input body length (m)
//' @param age Age of the animal (years)
//' @return Estimated mass (kg)
// [[Rcpp::export]]
 
 double length2mass(double L, double age){
   double a, b;                           
   if(age <= 0.79){
     a = R::rnorm(-5.091821, 0.2578327);
     b = R::rnorm(3.077823, 0.08325852);
   } else {
     a = R::rnorm(-5.096379, 0.2592405);
     b = R::rnorm(3.079408, 0.08360103);
   }
   return pow(10, a + b*log10(L*100)); // m to cm
 }

 
//' Incidence of foraging behavior
//' @description Determines whether the prey concentration encountered by an animal
//' is sufficient to support foraging
//' @param gamma Minimum prey density threshold that triggers foraging (\ifelse{html}{\out{copepods/m<sup>3</sup>}}{\eqn{copepods/m^3})
//' @param D Prey concentration (\ifelse{html}{\out{copepods/m<sup>3</sup>}}{\eqn{copepods/m^3})
// [[Rcpp::export]]
 
 double feeding_threshold(double gamma, double D){
   return std::round(1/(1 + exp(gamma - D)));
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
//'
//' @param angle Input angle in degrees
//' @return Equivalent angle in radians
// [[Rcpp::export]]
 
 double deg2radians(double angle){
   return angle * M_PI / 180;
 }
 
//' Area of the gape
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
//' @description Calculates the amount of energy that a whale can obtain from an average individual prey
//' @param m Average mass of the prey (kg/cop)
//' @param rho Energy density of the prey (J/kg) 
//' @param E_digest Digestive efficiency (fecal and urinary, \%)
//' @param E_hif Metabolizing efficiency (1-heat increment of feeding, \%)
// [[Rcpp::export]]
 
 double Econtent_cop(double m, double rho, double E_digest, double E_hif){ 
   return m * rho * E_digest * E_hif;
 }
 
//' Filtration rate
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
//' @description Predicts mammary mass from total mass in females
//' @param M Total body mass (kg)
// [[Rcpp::export]]
 
 double mammary_mass(double M){ 
   return std::pow(10, 0.902 * std::log10(M) - 1.965);
 }  

//' Milk production rate
//' @description Predicts the amount of milk produced by unit time from mammary mass
//' @param m Mass of mammary glands (kg)
//' @return Milk production rate (\ifelse{html}{\out{m<sup>3</sup>/kg/s)}}{\eqn{m^3/kg/s})
// [[Rcpp::export]]
 
 double milk_production(double m){ 
   return 0.0835 * std::pow(m, 1.965);
 } 

//' Resting metabolic rate
//' @description Predicts the resting metabolic rate of an animal from its mass 
//' using the allometric relationship proposed proposed by Williams & Maresh (2015)
//' @param M Total body mass (kg)
//' @param phi Scalar constant for immature animals
//' @note The RMR is scaled up in immature animals (Fortune et al., 2013; Rechsteiner et al., 2013)
//' in order to account for the elevated metabolic demand associated with active growth (Lavigne et al., 1986)
// [[Rcpp::export]]
 
 double RMR(double M, double phi){ 
   return phi * 581 * std::pow(M, 0.68);
 } 
 
//' Costs of locomotion
//' @description Predicts total locomotor costs from mass-specific, stroke-based allometric
//' relationships proposed by Williams et al. (2017), which distinguish between routine and
//' maximum swimming gaits.
//' @param M Total body mass (kg)
//' @param Sroutine Routine stroke frequency
//' @param Smax Maximum (aerobic) stroke frequency
//' @param delta_a Proportion of time spent in activity a
//' @param phi_a Scaling factor for the routine stroke frequency during activity a
//' @param psi_a Scaling factor for the maximum stroke frequency during activity a
//' @note These equations were derived from data on odontocetes up to 3,000 kg in weight,
//' including killer whales (Orcinus orca). It is assumed that these relationships hold when
//' extrapolating to the larger body mass range exhibited by North Atlantic right whales.
//' @return Total daily locomotor costs (kJ)
// [[Rcpp::export]]
 
 double locomotor_costs(double M, double Sroutine, double Smax, 
                       Rcpp::NumericVector delta_a, 
                       Rcpp::NumericVector phi_a, 
                       Rcpp::NumericVector psi_a){ 
   
   if ((delta_a.size()!=phi_a.size()) || (psi_a.size() ^ delta_a.size())) 
     Rcpp::stop("Arguments have different lengths");
   
   int n = delta_a.size();
   double t_routine = 0, t_max = 0;
   
   for(int i = 0; i < n; ++i){
     t_routine += delta_a[i] * phi_a[i];
     t_max += delta_a[i] * psi_a[i];
   }
   
   return M * (t_routine * Sroutine * (1.46 + 0.0005 * M) + t_max * Smax * (5.17 + 0.0002 * M));
 } 
 
//' Energetic cost of fetal growth during pregnancy
//' @param M_muscle Mass of muscles in fetus (kg)
//' @param M_viscera Mass of viscera in fetus (kg)
//' @param M_bones Mass of bones in fetus (kg)
//' @param M_blubber Mass of blubber in fetus (kg)
//' @param rho_lipid Energy density of lipids (kJ/kg)
//' @param rho_protein Energy density of protein (kJ/kg)
//' @param P_lip_muscle Lipid concentration in muscle (\%)
//' @param P_pro_muscle Protein concentration in muscle (\%)
//' @param P_lip_viscera Lipid concentration in viscera (\%)
//' @param P_pro_viscera Protein concentration in viscera (\%)
//' @param P_lip_bones Lipid concentration in bones (\%)
//' @param P_pro_bones Protein concentration in bones (\%)
//' @param P_lip_blubber Lipid concentration in blubber (\%)
//' @param P_pro_blubber Protein concentration in blubber (\%)
// [[Rcpp::export]]
 
 double fetal_growth(double M_muscle, double M_viscera, 
                    double M_bones, double M_blubber,
                    double rho_lipid, double rho_protein, 
                    double P_lip_muscle, double P_pro_muscle,
                    double P_lip_viscera, double P_pro_viscera, 
                    double P_lip_bones, double P_pro_bones,
                    double P_lip_blubber, double P_pro_blubber){ 
   
   double fg_m, fg_v, fg_b, fg_bl;
   fg_m = M_muscle * (rho_lipid * P_lip_muscle + rho_protein * P_pro_muscle);
   fg_v = M_viscera * (rho_lipid * P_lip_viscera + rho_protein * P_pro_viscera);
   fg_b = M_bones * (rho_lipid * P_lip_bones + rho_protein * P_pro_bones);
   fg_bl = M_blubber * (rho_lipid * P_lip_blubber + rho_protein * P_pro_blubber);
   
   return fg_m + fg_v + fg_b + fg_bl;
     
 }  

//' Energetic cost of placental maintenance during pregnancy
//' @param G Energetic cost of fetal growth (kJ)
// [[Rcpp::export]]
 
 double placental_maintenance(double G){ 
   return (G/0.807)*(1-0.807);
 }  
 
//' Heat increment of gestation
//' @param birth_mass Birth mass of the fetus (kg)
//' @param delta_m Daily growth rate of the fetus (kg/day)
// [[Rcpp::export]]
 
 double heat_gestation(double birth_mass, double delta_m){ 
   return 18409.6 * std::pow(birth_mass, 1.2) * (delta_m/birth_mass);
 }  

//' Fetal tissue mass
//' @param P_b Proportion of the body volume comprised of tissue b
//' @param L_birth Length of the fetus at birth (m)
//' @note This relationship only applies to muscles, bones, and viscera
// [[Rcpp::export]]
 
 double fetal_tissue_mass(double P_b, double L_birth){ 
   return 1000 * P_b * std::exp(-4.115 + 3.016 * std::log10(L_birth));
 }  

//' Fetal blubber mass
//' @param L_birth Length at birth (m)
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
 
 double fetal_blubber_mass(double L_birth,
                          double M_muscle, double M_viscera, double M_bones,
                          double D_blubber, double D_muscle, double D_viscera, double D_bones){ 
   return D_blubber * (std::exp(-4.115 + 3.016 * std::log10(L_birth)) - 
                       (M_muscle/D_muscle) - (M_viscera/D_viscera) - (M_bones/D_bones));
 }  

//' Energetic cost of gestation
//' @param cohortID Integer indicating the population cohort to which an animal belongs
//' @param G Energetic cost of fetal growth during pregnancy (kJ/day)
//' @param P Energetic cost of placental maintenance during pregnancy (kJ/day)
//' @param Q Heat increment of gestation (kJ/day)
// [[Rcpp::export]]
 
 double gestation_cost(int cohortID, double G, double P, double Q){ 
   if(cohortID == 4){
     return G + P + Q;
   } else 
     return 0.0f;
 }

//' Energetic cost of growth
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

//' Daily growth increment
//' @param m Current body mass (kg)
//' @param m_previous Body mass in the previous day (kg)
// [[Rcpp::export]]
 
 double growth_increment(double m, double m_previous){
   return (m - m_previous)/365;
 }

//' Tissue deposition
//' @description Returns the mass obtained after deposition of any surplus energy as lipid and lean tissue
//' @param M Start mass (kg)
//' @param E_net Net energy balance (kJ)
//' @param add_protein Percent of protein synthesis during anabolism (\%)
//' @param add_lipid Percent of lipid synthesis during anabolism (\%)
//' @param D_lipid Efficiency of deposition of lipids during anabolism (\%)
 //' @param D_protein Efficiency of deposition of proteins during anabolism (\%)
//' @param break_lipid Percent of lipid breakdown during catabolism (\%)
//' @param break_protein Percent of protein breakdown during catabolism (\%)
//' @param C_lipid Efficiency of breakdown of lipids during catabolism (\%)
//' @param C_protein Efficiency of breakdown of protein during catabiolism (\%)
//' @param rho_lipid Energy density of lipids (kJ/kg)
//' @param rho_protein Energy density of proteins (kJ/kg)
// [[Rcpp::export]]
 
 double new_mass(double M, double E_net,
                double add_protein, double add_lipid,
                double D_lipid, double D_protein,
                double break_lipid, double break_protein,
                double C_lipid, double C_protein,
                double rho_lipid, double rho_protein){
   
   double fat_growth = E_net/rho_lipid;
   double lean_growth = E_net/rho_protein;
   
   // Fat growth
   if(E_net > 0){
     fat_growth *= add_lipid * D_lipid;
   } else if(E_net < 0){
     fat_growth *= break_lipid * C_lipid;
   }
   
   // Lean tissue growth
   if(E_net > 0){
     lean_growth *= add_protein * D_protein;
   } else if(E_net < 0){
     lean_growth *= break_protein * C_protein;
   }
   
   return M + fat_growth + lean_growth;
 }

#endif