// R and 'Eigen' integration using 'Rcpp'. 
// 'Eigen' is a C++ template library 
// for linear algebra: matrices, vectors, 
// numerical solvers and related algorithms.

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppEigen.h>
#include <Rcpp.h> // For lists
#include <iostream>
#include <cstdlib>
#include "environments.h"
#include "movements.h"
#include "geodesic.h"
#include "bioenergetics.h"
#include "progressbar.hpp"
#include <random>
using namespace std;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]


// Structures (also called structs) are a way to group several 
// related variables into one place. 
// Each variable in the structure is known as a member of the structure.
// Unlike an array, a structure can contain many different data types (int, string, bool, etc.).

// The public keyword is an access specifier. 
// Access specifiers define how the members (attributes and methods) 
// of a class can be accessed. In the example above, the members are public 
// - which means that they can be accessed and modified from outside the code.

/**
 * Basic state for an animal
 */
// Using const to declare data members constant, i.e., must have some value before the object is constructed
struct Animal {
  double x, y;
  int cohortID;
  int alive;
  double age;
  float bc;
  Eigen::MatrixXd lengthatage;
  float length;
  Eigen::MatrixXd massatlength;
  float mass;
  float fatmass;
  float mouth_ratio;
  float mouth_angle;
  float mouth_width;
  Rcpp::NumericVector entangled;
  
  // Member initializer list - age initialized by default at birth
  Animal() : x(0), y(0), cohortID(0), alive(1), age(0), bc(0.06), 
  lengthatage(agL(age)), length(age2length(age, lengthatage)),
  massatlength(mL()), mass(length2mass(length, massatlength, false)), fatmass(bc*mass), 
  mouth_ratio(start_mouth(cohortID, age)), mouth_angle(76.7), 
  mouth_width(length*mouth_ratio), entangled() { }
  
  // x: Easting
  // y: Northing
  // cohortID: Unique identifier for the cohort to which an animal belongs.
  // age: Age of the animal in decimal years
  // bc: Body condition, expressed as a relative % of fat
  // lengthatage: Coefficients of the Gompertz length-at-age relationship
  // length: Total body length in cm
  // massatlength: Coefficients of the logarithmic mass-at-length relationship
  // mass: Total body mass in kg
  // fatmass: Total blubber mass in kg
  // mouth_ratio: Ratio of mouth width to body length (d.u.)
  // mouth_angle: Angle between the tip of the baleen plates and the outer edges of the baleen rack (rad)
  // mouth_width: Width of the mouth, i.e. body width at 10% of the body length from the snout (m)
  // entangled: Integer indicating whether the animal is entangled and 
  //   where the attachment site of the ropes is on the animal's body. The variable can take
  //   the following values: 0 (not entangled), (1) entangled (anterior), (2) entangled (posterior).
  // feeding: Integer indicating whether the animal is in a foraging state (1) or not (0)
  
  Animal(const double & x0, 
         const double & y0,
         int cohort) : 
    x(x0), y(y0), 
    cohortID(cohort),
    alive(1),
    age(start_age(cohort)),
    bc(start_percfat(age)),
    lengthatage(agL(age)), 
    length(age2length(age, lengthatage)),
    massatlength(mL()),
    mass(length2mass(length, massatlength, false)),
    fatmass(bc*mass),
    mouth_ratio(start_mouth(cohort, age)), 
    mouth_angle(76.7), 
    mouth_width(length*mouth_ratio),
    entangled(entanglement_event()){ }
};

/**
 //' @title Animal with capacity for attraction to other latent animals
 //'
 //' @description Initialize an animal with latent members, and couple the observed position to one of the latent animals
 //' 
 //' @param init_latent latent animals to be attracted to
 //' @param init_animal index of latent animal to use as starting location
 */

struct CouplingAnimal : public Animal {
  
  std::vector<Animal> latent;
  
  // Member initializer list
  CouplingAnimal() { }
  
  CouplingAnimal(
    const std::vector<Animal> & init_latent,
    std::size_t init_animal
  ) : Animal(init_latent[init_animal].x, 
  init_latent[init_animal].y, 
  init_latent[init_animal].cohortID),
  latent(init_latent) { }
};

/**
 * Simulate a collection of animals moving in a projected coordinate space
 *
 * @tparam AnimalType
 * @tparam EnvironmentType
 * @tparam MovementType
 * @tparam CoordSize the number of coordinates
 * @param animals reference to a collection of pre-initialized animals
 * @param environments reference to a collection of environments for all times
 * @param m object used to update animal's position
 * @return array of simulated locations with indices [coord, animal, time]
 */
// template<typename AnimalType, 
//          typename EnvironmentContainer,
//          typename MovementRule, 
//          int CoordSize = 2>
// Rcpp::NumericVector movesim(
//     std::vector<AnimalType> & animals,
//     EnvironmentContainer & environments,
//     MovementRule & m,
//     Eigen::MatrixXd support,
//     Eigen::VectorXd limits, 
//     Eigen::VectorXd resolution,
//     double stepsize,
//     bool progress) {
// 
//     // Number of animals and time points to simulate
//     std::size_t n = animals.size();
//     std::size_t t = environments.size();
//     
//     // Initialize and label output
//     Rcpp::NumericVector out(n*t*CoordSize);
//     out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
//     out.attr("dimnames") = Rcpp::List::create(
//         Rcpp::CharacterVector::create("easting", "northing"), R_NilValue,
//         R_NilValue
//     );
// 
//     // Loop over time points and environments: One environment per time point
//     Rcpp::NumericVector::iterator data_out = out.begin();
//     
//     auto env_end = environments.end();
//     auto animals_end = animals.end();
//     
//     // Initialize progress bar
//     progressbar bar(t);
//     
//     // 365 steps per animal
//     for(auto env = environments.begin(); env != env_end; ++env) {
//       
//       if(progress) bar.update();
//       
//         for(auto animal = animals.begin(); animal != animals_end; ++animal) {
//           
//             // Save current coordinates
//             *(data_out++) = animal->x;
//             *(data_out++) = animal->y;
//             
//             // Update coordinates
//             m.update(*animal, *env);
//         }
//     }
// 
//     return out;
// }

template<typename AnimalType, 
         typename EnvironmentContainer,
         typename MovementRule, 
         int CoordSize = 4,
         int nparam_E = 6,
         int nparam_animal = 15,
         int nparam_fetus = 15,
         int nparam_stressors = 9,
         int nparam_activity = 9,
         int nparam_feeding = 15,
         int nparam_nursing = 12,
         int nparam_costs = 10,
         int nparam_costs_calves = 5>

Rcpp::List movesim(
    int cohortID,
    std::vector<AnimalType> & animals,
    std::vector<AnimalType> & calves,
    EnvironmentContainer & environments,
    std::vector<Environment> & layers,
    Rcpp::NumericVector doseresp,
    std::vector<std::size_t> densitySeq,
    MovementRule & m,
    Eigen::MatrixXd support,
    Eigen::VectorXd limits, 
    Eigen::VectorXd resolution,
    double stepsize,
    bool progress) {
  
  // Number of animals and time points to simulate
  std::size_t n = animals.size();  // nsim
  std::size_t t = environments.size(); // 365 days
  
  // Initialize list in which results will be stored
  Rcpp::List results;
  
  // Initialize and label array in which x,y coordinates will be stored
  Rcpp::NumericVector out(n*t*CoordSize);
  out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
  out.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("easting", "northing", "region", "d"), R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which energy budget will be stored
  Rcpp::NumericVector attrib_E(n*t*nparam_E);
  attrib_E.attr("dim") = Rcpp::Dimension(nparam_E,n,t);
  attrib_E.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_tot", "E_in", "E_out", "E_tot_calves", "E_in_calves", "E_out_calves"), R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which bioenergetic parameters will be stored
  Rcpp::NumericVector attrib(n*t*nparam_animal);
  attrib.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  attrib.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c", 
                                  "mass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", "mouth_w"),
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for dependent calves
  Rcpp::NumericVector attrib_calves(n*t*nparam_animal);
  attrib_calves.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  attrib_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c", 
                                  "mass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", "mouth_w"), 
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for fetuses
  Rcpp::NumericVector attrib_fetus(n*t*nparam_fetus);
  attrib_fetus.attr("dim") = Rcpp::Dimension(nparam_fetus,n,t);
  attrib_fetus.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("fetus_l", "fetus_m", "birth_l", "birth_m", "muscle_m", "viscera_m", "bones_m", "blubber_m",
                                  "muscle_lip", "muscle_pro", "visc_lip", "visc_pro", "bone_lip", "blubber_lip", "blubber_pro"), 
                                  R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which stressor data will be stored
  Rcpp::NumericVector attrib_stressors(n*t*nparam_stressors);
  attrib_stressors.attr("dim") = Rcpp::Dimension(nparam_stressors,n,t);
  attrib_stressors.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("is_entgl", "entgl_head", "severity", "entgl_d",
                                  "entgl_start", "entgl_end", "strike", "noise", "dB_thresh"),
    R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which activity budgets will be stored
  Rcpp::NumericVector attrib_activity(n*t*nparam_activity);
  attrib_activity.attr("dim") = Rcpp::Dimension(nparam_activity,n,t);
  attrib_activity.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("d_travel", "speed", "t_travel", "t_feed", "t_nurse", "t_rest", "n_zero", "t_sum", "t_remain"),
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_feeding(n*t*nparam_feeding);
  attrib_feeding.attr("dim") = Rcpp::Dimension(nparam_feeding,n,t);
  attrib_feeding.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("feed", "preyconc", "minprey", "gape", "speed", "captEff", "impedance", "daylight",
                                  "feed_effort", "targetBC", "cop_mass", "cop_kJ", "digestEff", "metabEff", "E_cop"), 
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in calves
  Rcpp::NumericVector attrib_nursing(n*t*nparam_nursing);
  attrib_nursing.attr("dim") = Rcpp::Dimension(nparam_nursing,n,t);
  attrib_nursing.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("assim", "provision", "mamm_M", "milk_rate", "Dmilk", "t_suckling",
                                  "targetBC", "nursing", "milk_lip", "milk_pro", "ED_lip", "ED_pro"),
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy costs in adults/juveniles
  Rcpp::NumericVector attrib_costs(n*t*nparam_costs);
  attrib_costs.attr("dim") = Rcpp::Dimension(nparam_costs,n,t);
  attrib_costs.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("rmr", "LC", "stroke", "E_growth", "E_gest", "fgrowth",
                                  "placenta", "hic", "E_lac", "delta_m"),
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy costs in calves
  Rcpp::NumericVector attrib_costs_calves(n*t*nparam_costs_calves);
  attrib_costs_calves.attr("dim") = Rcpp::Dimension(nparam_costs_calves,n,t);
  attrib_costs_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("rmr", "LC", "stroke", "E_growth", "delta_m"),
                                  R_NilValue, R_NilValue
  );
  
  // "fetal_mass", "fm_muscle", "fm_viscera", "fm_bones", "fm_blubber",
  // "muscle_lip", "muscle_pro", "visc_lip", "visc_pro", "bones_lip", "bones_pro",
  // "blubber_lip", "blubber_pro", ")
  
  // Loop over time points and environments: One environment per time point
  Rcpp::NumericVector::iterator data_out = out.begin();
  Rcpp::NumericVector::iterator attrib_out = attrib.begin();
  Rcpp::NumericVector::iterator attrib_E_out = attrib_E.begin();
  Rcpp::NumericVector::iterator attrib_calves_out = attrib_calves.begin();
  Rcpp::NumericVector::iterator attrib_fetus_out = attrib_fetus.begin();
  Rcpp::NumericVector::iterator attrib_stressors_out = attrib_stressors.begin();
  Rcpp::NumericVector::iterator attrib_activity_out = attrib_activity.begin();
  Rcpp::NumericVector::iterator attrib_feeding_out = attrib_feeding.begin();
  Rcpp::NumericVector::iterator attrib_nursing_out = attrib_nursing.begin();
  Rcpp::NumericVector::iterator attrib_costs_out = attrib_costs.begin();
  Rcpp::NumericVector::iterator attrib_costs_calves_out = attrib_costs_calves.begin();
  
  // std::advance(attrib_out, nparam); // Skip update of parameters at first time point as already initialized
  
  auto env_end = environments.end();
  auto animals_end = animals.end();
  
  // Inititalize progress bar
  progressbar bar(t);
  
  // std::cout << cohortID << std::endl;
  
  // double age, length, mtot, mfat, bc;
  // std::default_random_engine generator;
  // std::uniform_int_distribution<int> binary_distribution(0,1);
  
  // ------------------------------------------------------------
  // DECLARE VARIABLES
  // ------------------------------------------------------------
  
  int current_animal, current_month, current_day, current_region;
  double current_x, current_y;
  Rcpp::NumericVector feeding_nursing_effort (2); // feeding, nursing
  Rcpp::NumericVector target_bc (2); // adjuv, calves
  
  // Stressors +++
  
  int vessel_strike = 0;
  int resp_noise = 0;
  double gear_risk;
  int start_entanglement = 0;
  int end_entanglement = 0;
  double strike_risk;
  // double p_response;
  double dB;
  
  // Activity budgets +++
  
  double daylight_hours;
  double feeding_time;
  double travel_time;
  double nursing_time;
  double resting_time;
  double travel_dist;
  int n_zeroes;
  double t_sum;
  double t_remain = 0;
  double tmax;
 
  // Foraging +++
  
  double D_cop;
  double min_calanus;
  bool is_feeding;
  double mouth_gape;
  double swim_speed_foraging;
  double swim_speed;
  double capture_efficiency;
  double gape_reduction = 0;
  double cop_mass;
  double cop_kJ;
  double digestive_efficiency;
  double metabolizing_efficiency;
  double E_calanus;
  double E_in;

  // Nursing +++
  
  double milk_assim = 0; 
  double milk_provisioning = 0;
  double mammary_M = 0;
  double milk_production_rate = 0;
  double D_milk = 0;
  double t_suckling = 0;
  double milk_lipids = 0; 
  double milk_protein = 0; 
  double E_in_calves = 0;
  
  // Energetic costs
  
  double resting_metab, resting_metab_calves; 
  double phi_rmr = 1;
  double stroke_rate;
  double scalar_LC;
  
  if(cohortID == 4){
    scalar_LC = 1.4; // Pregnant females
  } else if(cohortID == 5){
    scalar_LC = 1.7; // Lactating females
  } else if(cohortID == 3) {
    scalar_LC = 0.96; // Adult males
  } else if(cohortID == 1 | cohortID == 2) {
    scalar_LC = 1.1;
  } else if(cohortID == 0) {
    scalar_LC = 0.98;
  } else {
    scalar_LC = 1;
  }
  
  double LC_tot, LC_tot_calves;
  
  Rcpp::NumericVector birth_length (n);
  Rcpp::NumericVector birth_mass (n);
  
  if(cohortID == 4){
    
    // Length and mass at birth (only relevant for pregnant females) -- (m)
    Eigen::MatrixXd birth_agL = agL(0, n);
    Eigen::MatrixXd birth_mL = mL(n);
    for (int i = 0; i < n; i++) {
      birth_length(i) = age2length(0, birth_agL.row(i));
      birth_mass(i) = length2mass(birth_length(i), birth_mL.row(i), false);
    }
  }

  double mass_muscle = 0;
  double mass_viscera = 0;
  double mass_bones = 0;
  double mass_blubber = 0;
  double fetal_growth_cost = 0;
  double placental_cost = 0;
  double fetus_growth = 0;
  double fetus_m_next;
  double HIC = 0;
  double E_gest = 0;
  
  double E_lac = 0;
  
  double lean_mass_current, lean_mass_current_calves;
  double length_nextday, length_nextday_calves;
  double lean_mass_nextday, lean_mass_nextday_calves;
  double mass_increment, mass_increment_calves;
  double E_growth, E_growth_calves;
  double E_out, E_out_calves;
  
  // Range of dose values for the dose-response function
  // int dose_lwr = 85;
  // int dose_uppr = 200;
  double response_dB;
  // std::vector<double> dose = dose_range(dose_lwr, dose_uppr, 1000);
  
  // ------------------------------------------------------------
  // CONSTANTS
  // ------------------------------------------------------------
  
  // Mammary gland efficiency -- (%)
  // As per Brody 1968
  const double mammary_efficiency = 0.9;
  
  // Starvation threshold expressed as relative blubber mass -- (d.u.)
  // As per Pirotta et al. (2018) - to test extreme conditions of leanness
  const double starvation = 0.05;
  
  // Non-linearity between milk supply and the body condition of the mother
  const int zeta = -2;
  
  // Age at which milk consumption starts to decrease -- (days)
  // Fortune et al. (2020)
  const int milk_decrease = 288;
  
  // Non-linearity between milk assimilation and calf age -- (d.u.)
  // Hin et al. (2019)
  const double eta_milk = 0.9;
  
  // Duration of the lactation -- (days)
  const int T_lac = 365;
  
  // Energy density of lipids -- (kJ / kg - converted to MJ)
  const int ED_lipids = 39300/1000;
  
  // Energy density of protein -- (kJ / kg - converted to MJ)
  // mean(c(23600,18000))
  const int ED_protein = 20800/1000;
  
  // Proportion of the body volume comprised of muscle -- (d.u.)
  const double P_muscle = 0.2820;
  
  // Proportion of the body volume comprised of viscera -- (d.u.)
  const double P_viscera = 0.1020;
  
  // Proportion of the body volume comprised of bones -- (d.u.)
  const double P_bones = 0.1250;
  
  // Density of blubber -- (kg / m3) 
  const int blubber_density = 700;
  
  // Density of muscle -- (kg / m3) 
  const int muscle_density = 960;
  
  // Density of bones -- (kg / m3) 
  const int bone_density = 720;
  
  // Density of viscera -- (kg / m3) 
  const int visceral_density = 930;
  
  // Proportion of protein in bones -- (d.u.)
  const double bones_protein = 0.248; // From Christiansen et al. (2022)  
  
  // Proportion of lean mass that is water -- (d.u.)
  const double lean_water = 0.672;
  
  // Efficiency of deposition of lipids -- (d.u.)
  const double deposition_lipids = 0.75;
  
  // Efficiency of deposition of protein -- (d.u.)
  const double deposition_protein = 0.4975;
  
  // Energetic cost of entanglement -- (MJ)
  const double entanglement_cost = 248.8222913;
  
  // ------------------------------------------------------------
  // Fetal development
  // ------------------------------------------------------------
  
  // Proportion of lipids in muscle -- (d.u.)
  Rcpp::NumericVector muscle_lipids = Rcpp::rnorm(n, 0.114, 0.08); // From Christiansen et al. (2022)   
  
  // Proportion of protein in muscle -- (d.u.)
  Rcpp::NumericVector muscle_protein = Rcpp::rnorm(n, 0.221, 0.023); // From Christiansen et al. (2022)   
  
  // Proportion of lipids in viscera -- (d.u.)
  Rcpp::NumericVector visceral_lipids = Rcpp::rnorm(n, 0.758, 0.096); // From Christiansen et al. (2022)   
  
  // Proportion of protein in viscera -- (d.u.)
  Rcpp::NumericVector visceral_protein = Rcpp::rnorm(n, 0.037, 0.017); // From Christiansen et al. (2022)   
  
  // Proportion of lipids in bones -- (d.u.)
  Rcpp::NumericVector bones_lipids = Rcpp::rnorm(n, 0.758, 0.096); // From Christiansen et al. (2022)   
  
  // Proportion of lipids in blubber -- (d.u.)
  Rcpp::NumericVector blubber_lipids = Rcpp::rnorm(n, 0.73056, 0.058482);
  
  // Proportion of protein in blubber -- (d.u.)
  Rcpp::NumericVector blubber_protein = Rcpp::rnorm(n, 0.1020, 0.039); // From Christiansen et al. (2022)  
  
  Rcpp::NumericVector fetus_l (n);
  Rcpp::NumericVector fetus_m (n);
  
  // ------------------------------------------------------------
  // START SIMULATION
  // ------------------------------------------------------------
  
  // 365 steps per animal
  for(auto env = environments.begin(); env != env_end; ++env) {
    
    current_day = env - environments.begin() + 1;
    
    if(progress) bar.update();
    
    // env is a vector iterator, so env - environments.begin() returns the associated index
    current_month = densitySeq[env - environments.begin()];
    
    // std::cout << current_month << " " << environments[env - environments.begin()].id << std::endl;
    
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {
      
      current_x = animal->x;
      current_y = animal->y;
      
      current_animal = animal - animals.begin();
      current_region = layers[current_month](current_x, current_y, 'R');
      
      // std::cout << current_x << "-" << current_y << "-" << current_region << std::endl;
      
      // Save current coordinates (fills by column)
      
      *(data_out++) = animal->x;
      *(data_out++) = animal->y;
      *(data_out++) = current_region;
      
      // ------------------------------------------------------------
      // VESSEL STRIKES
      // ------------------------------------------------------------
      
      if(animal->alive){
      
      // Vessel strikes +++
      strike_risk = layers[current_month](current_x, current_y, 'V');
      
      // Incidence of a vessel strike -- (d.u.)
      vessel_strike = R::rbinom(1, strike_risk);
      // vessel_strike = 0; // For model development only
      }
      
      if(vessel_strike | animal->bc < starvation){
        animal->alive = 0;
      }
      
      if(animal->alive){
        
        // ------------------------------------------------------------
        // COORDINATES
        // ------------------------------------------------------------
        
        // animal->x = coords_xy[0];
        // animal->y = coords_xy[1];
        
        // // Save current coordinates (fills by column)
        // 
        // *(data_out++) = animal->x;
        // *(data_out++) = animal->y;
        
        bool calc_geo = false;
        int gd;
        int *gdptr;
        gdptr = &gd;
        
        
        // Only activate geodesic calculations in some areas
        if((current_x >= 685 & current_x <= 1420 & current_y >= 1070 & current_y <= 1865) ||
           (current_x >= 520 & current_x <= 730 & current_y >= 825 & current_y <= 960)) calc_geo = true;
        
        if(!calc_geo){
          
          // Update the animal's state
          m.update(*animal, *env, TRUE, current_region);
          
        } else {
          
          while(calc_geo == true){
            
            // Update the animal's state
            m.update(*animal, *env, TRUE, current_region);
            
            double next_x = animal->x;
            double next_y = animal->y;
            
            // Calculate geodesic distance
            *gdptr = geoD(support, current_x, current_y, next_x, next_y, limits, resolution);
            
            if(*gdptr >= 0 & *gdptr < stepsize){
              
              calc_geo = false;
              
            } else {
              
              animal->x = current_x;
              animal->y = current_y;
            }
            
          } // End while calc_geo
          
        } // End if(calc_geo)
        
        travel_dist = std::sqrt(std::pow(current_x - animal->x, 2) +  std::pow(current_y - animal->y, 2));
        *(data_out++) = travel_dist;
        
        // ------------------------------------------------------------
        // NOISE
        // ------------------------------------------------------------
        
        resp_noise = 0; // Reset daily
        
        // Pile-driving noise +++
        dB = layers[current_month](current_x, current_y, 'N');
        
        // Response threshold (dB)
        response_dB = response_threshold(doseresp);
        
        // Behavioral response to pile-driving noise exposure -- (d.u.)
        if(dB > response_dB) resp_noise = 1;

        // if(dB <= dose_lwr){
        //   p_response = 0;
        // } else if(dB >= dose_uppr){
        //   p_response = 1;
        // } else {
        //   p_response = prob_response(dose, doseresp, animal->doseID, dB);
        // }
        // resp_noise = R::rbinom(1, p_response);
        
        // ------------------------------------------------------------
        // FISHING GEAR
        // ------------------------------------------------------------
        
        // An animal can only become entangled if not carrying gear already 
        if(animal->entangled(0) == 0){
          gear_risk = layers[current_month](current_x, current_y, 'F');
          animal->entangled = entanglement_event(gear_risk);
          if(animal->entangled(0) == 1){
            start_entanglement = current_day;
            end_entanglement = current_day + animal->entangled(3);
          } else {
            start_entanglement = 0;
            end_entanglement = 0;
          }
        }

        // Terminate entanglement event
        if(current_day == end_entanglement){
          animal->entangled = entanglement_event(0);
          end_entanglement = 0;
          start_entanglement = 0;
        }
        
        // ------------------------------------------------------------
        // PREY
        // ------------------------------------------------------------
        
        // Prey concentration in cell at time t -- (copepods/cubic m)
        D_cop = layers[current_month](current_x, current_y, 'P');
        
        // Minimum prey concentration (Calanus finmarchicus stage V) -- (copepods/cubic m)
        // Meta: From year 1987, by species
        min_calanus = rtnorm(1784.979, 718.77, 0, INFINITY);
        
        // Determine whether the prey concentration is sufficient to support foraging -- (d.u.)
        if(resp_noise){
          is_feeding = 0;
        } else {
          is_feeding = feeding_threshold(min_calanus, D_cop);
        }
        
        // ------------------------------------------------------------
        // ACTIVITY BUDGET
        // ------------------------------------------------------------
        
        // Total time available for foraging activities --(hr converted to s)
        daylight_hours = layers[current_month](current_x, current_y, 'L');
        
        // Based largely on Cusano et al. (2019)
        
        if(current_region == 9){ // SEUS
          
          // On the calving grounds, mothers and calves prioritize 
          // rest, travel and nursing in this order
          
          // Time spent traveling per day -- (hr)
          // Calculated based on distance covered that day and mean swimming speed
          // 4.6 m/s taken as maximum average speed based on satellite tag record described in Mate (1997)
          // Minimum chosen to ensure that travel_time does not exceed 24 hours
          swim_speed = rtnorm(0.7480997, 0.5403338, (travel_dist*1000)/(24*3600), 4.6) * 3600; // m per hr
          travel_time = travel_dist * 1000 / swim_speed; // m
          
          // Time spent resting + nursing per day -- (hr)
          // Meta: By age class, combine calves and mother/calf pairs, bounds of (0,24)
          if(cohortID == 5){ // Lactating females
            resting_time = rtnorm(13.76, 0.76984, 0, 24 - travel_time);
            nursing_time = rtnorm(1.6342, 4.6767, 0, 24 - (resting_time + travel_time)); 
          } else {
            resting_time = rtnorm(5.487199, 0.8499895, 0, 24 - travel_time);
            nursing_time = 0;
          }
          
          feeding_time = 0;
          
        } else { // Feeding grounds and elsewhere
          
          if (current_region == 7 | current_region == 8) { // Migratory corridor (MIDA and SCOS)
            
            // During migration, travel is the priority
            
            swim_speed = rtnorm(2.251177, 1.269262, (travel_dist*1000)/(24*3600), 4.6) * 3600; // m per hour
            travel_time = travel_dist * 1000 / swim_speed;
            
          } else {
            
            swim_speed = rtnorm(0.8368928, 0.3751697, (travel_dist*1000)/(24*3600), 4.6) * 3600; // m per hour
            travel_time = travel_dist * 1000 / swim_speed;
            
          }
            
            tmax = std::min(24 - travel_time, daylight_hours);
            
            feeding_time = is_feeding * rtnorm(12.768, 2.9166, 0, tmax);
            resting_time = rtnorm(5.487199, 0.8499895, 0, 24 - (feeding_time + travel_time));
            
            if(cohortID == 5){ 
              nursing_time = 24 - (feeding_time + travel_time + resting_time);
            } else {
              nursing_time = 0;
            }
        }

          // if(cohortID == 5){ 
          //   resting_time = rtnorm(13.76, 0.76984, 0, 24 - travel_time);
          // } else { 
          //   resting_time = rtnorm(0.6449982, 0.1963043, 0, 24 - travel_time);
          //   }
          // 
          // tmax = std::min(24 - (travel_time + resting_time), daylight_hours);
          // feeding_time = is_feeding * rtnorm(12.768, 2.9166, 0, tmax);
          // 
          // if(cohortID == 5){ 
          //   nursing_time = rtnorm(1.6342, 4.6767, 0, 24 - (feeding_time + resting_time + travel_time)); 
          // } else { 
          //   nursing_time = 0;
          //   }
          
          // resting_time = 24 - (travel_time + nursing_time);
          
          
        // } else { // Feeding grounds and elsewhere
        //   
        //  
        //   
        //   tmax = std::min(24 - travel_time, daylight_hours);
        //   
        //   feeding_time = is_feeding * rtnorm(12.768, 2.9166, 0, tmax);
        //   resting_time = rtnorm(5.487199, 0.8499895, 0, 24 - (feeding_time + travel_time));
        //   
        //   if(cohortID == 5){ 
        //     nursing_time = 24 - (feeding_time + travel_time + resting_time);
        //   } else {
        //     nursing_time = 0;
        //   }
        // }

        t_sum = travel_time + resting_time + nursing_time + feeding_time;
        std::vector<double> t = {travel_time, resting_time, nursing_time, feeding_time};
        n_zeroes = std::count(t.begin(), t.end(), 0);
        
        if(t_sum < 24){
          t_remain = (24 - t_sum)/(4-n_zeroes);
        }
        
        // int can_feed = 0;
        // if(D_cop > min_calanus) can_feed = 1;
        // std::cout << current_region << " - " << resp_noise << " - " << response_dB << " - " << dB << " - " << daylight_hours << " - " << is_feeding << " - " << can_feed << " - " << min_calanus << " - " << D_cop << " - " << travel_time << " - " << 
        //   feeding_time << " - " << resting_time << " - " << nursing_time << " - " << t_remain << std::endl;
        
        if(travel_time > 0) travel_time += t_remain;
        if(resting_time > 0) resting_time += t_remain;
        if(nursing_time > 0) nursing_time += t_remain;
        if(feeding_time > 0) feeding_time += t_remain;
        
        // Convert to seconds
        travel_time = travel_time * 3600;
        feeding_time = feeding_time * 3600;
        resting_time = resting_time * 3600;
        nursing_time = nursing_time * 3600;
        
        // ------------------------------------------------------------
        // ENERGY INTAKE - FORAGING
        // ------------------------------------------------------------
        
        // Mouth gape area -- (m^2)
        mouth_gape = gape_size(animal->length, animal->mouth_width, animal->mouth_angle);
        
        // Swim speed during foraging -- (m/s)
        // Meta: No temporal filter, all records
        swim_speed_foraging = R::rnorm(0.8994045, 0.1830964);
        
        // Capture efficiency -- (%)
        // Meta: No temporal filter, all records
        // draw("beta", 0.91575, 0.0392)
        capture_efficiency = R::rbeta(45.06241, 4.145791);
        
        // Percent reduction in the gape during an entanglement event
        // Meta: Parameter unknown
        if(animal->entangled(0) == 1 & animal->entangled(1) == 1) gape_reduction = R::runif(0,0.5);
        
        // Foraging/nursing response to body condition -- (d.u.)
        // Parameters are estimated by:
        // (1) running the meta-analysis function on target relative blubber mass (see separate R script)
        // (2) drawing the resulting curve and recording the maximum value from the associated distribution
        // (3) using the optim_feeding function with bounds of c(0, max) to find the parameters of the logistic function
        
        // See meta_analysis file for calculations
        
        if(animal->cohortID == 4){  // Pregnant females
          
          target_bc[0] = 0.6439795;
          feeding_nursing_effort[0] = scale_effort(animal->bc, 1, 0, 27.36939, 0.4258282 , 0.5918155);
          
        } else if(animal->cohortID == 5){ // Lactating mothers and their calves
          
          target_bc[0] = 0.5271131;
          target_bc[1] = 0.6439795;
          
          // Feeding effort is zero for lactating mothers (Miller et al. 2012)
          feeding_nursing_effort[1] = scale_effort(calves[current_animal].bc, 1, 0, 27.36939, 0.4258282 , 0.5918155);
          
        } else { // All other cohorts
          
          target_bc[0] = 0.5271131;
          feeding_nursing_effort[0] = scale_effort(animal->bc, 1, 0, 24.18117, 0.2219739, 1.284047);
          
        }
        
        // Average mass of the prey -- (mg converted to g)
        // Meta: From year 2000, all records
        cop_mass = rtnorm(0.38224, 0.095553, 0, INFINITY)/1000;
        
        // Energy density of Calanus stage V -- (kJ / g converted to MJ)
        // Meta: No temporal filter, by species, combine all records for C. finmarchicus
        cop_kJ = R::rnorm(24.054, 3.3191)/1000;
        
        // Digestive efficiency -- (%)
        // Meta: No temporal filter, all records
        digestive_efficiency = R::rnorm(0.88504, 0.0045112);
        
        // Metabolizing efficiency (1 - heat increment of feeding) -- (%)
        // Meta: No temporal filter, by age_class, combining juveniles and calves
        
        if(animal->cohortID <= 2){ // Juveniles
          metabolizing_efficiency = rtnorm(0.740, 0.046, 0, 1);
        } else { // Adults
          metabolizing_efficiency = rtnorm(0.875, 0.036, 0, 1);
        }
        
        // Prey energy content -- (MJ / copepod)
        E_calanus = Econtent_cop(cop_mass, cop_kJ, digestive_efficiency, metabolizing_efficiency);
        
        // ENERGY INTAKE (MJ)
        E_in = (1 - resp_noise) * is_feeding * D_cop * mouth_gape * (1 - gape_reduction) *
          swim_speed_foraging * capture_efficiency * feeding_time * feeding_nursing_effort[0] *
          E_calanus;
        
        // std::cout << current_x << ";" << current_y << " - " << current_month << " - " << D_cop << " - " << min_calanus << " " << is_feeding << " " << E_in<< std::endl;
        
        // ------------------------------------------------------------
        // ENERGY INTAKE - SUCKLING
        // ------------------------------------------------------------
        
        // // Incidence of mother-calf separation outside of SEUS -- (d.u.)
        // // Reported to be between 10 and 40% of sightings as per Hamilton et al. (2022) - take mean = 25%
        // bool is_separated = 0;
        // int separation_days = 0;
        // bool previously_separated = 0;
        // 
        // if(!previously_separated & animal->y > -12){
        //   is_separated = R::rbinom(1, 0.25);
        // }
        // 
        // // Separation duration -- (days)
        // // Using half-normal with mean of 5.9 and max of 23, as per Hamilton et al. (2022)
        // // Coded as truncated Normal centered on zero
        // if(is_separated){
        //   separation_days = separation_duration();
        //   end_separation = current_day + separation_days;
        //   previously_separated = 1;
        // }
        // 
        // // Terminate separation event if needed
        // if(current_day <= end_separation){
        //   is_separated = 0;
        // }
        
        // if(!is_separated){
        
        // Milk assimilation -- (d.u.)
        // Varies with calf age
        milk_assim = milk_assimilation(current_day, T_lac, milk_decrease, eta_milk);
        
        // Milk provisioning -- (d.u.)
        // Varies with mother's body condition
        milk_provisioning = milk_supply(starvation, target_bc[0], animal->mass, animal->fatmass, zeta);
        
        // Total mass of mammary glands -- (kg)
        mammary_M = mammary_mass(animal->mass);
        
        // Milk production rate by the mother -- (cubic m per kg per sec)
        milk_production_rate = milk_production(mammary_M);
        
        // Density of milk -- (kg/m3)
        // Meta: No temporal filter, all records
        D_milk = R::rnorm(1.02, 0.01) * 1000; // L to cubic m
        
        // Time spent nursing/suckling per day -- (hours to sec)
        // Meta:: No filter on time spent nursing, time spent suckling during nursing bouts filtered by species (value for E. australis)
        t_suckling = rtnorm(0, 3.688672, 0, 24) *3600;
        
        // Proportion of lipids in milk -- (%)
        // Meta: All records – as only value for closest relative (Eubalaena mysticetus) deemed low
        // Set maximum to 60% as reported in the literature (White 1953)
        milk_lipids = rtnorm(0.365680, 0.11932, 0, 0.51);
        
        // Proportion of protein in milk -- (%)
        milk_protein = rtnorm(0.13650, 0.031562, 0, 1);
        
        if((milk_lipids + milk_protein) > 1){
          std::cout << "Inconsistent milk composition" << std::endl;
          break;
        }
        
        // ENERGY INTAKE (MJ)
        E_in_calves = (1 - resp_noise) * milk_assim * milk_provisioning * mammary_efficiency *
          mammary_M * milk_production_rate * t_suckling * feeding_nursing_effort[1] * D_milk *
          (milk_lipids * ED_lipids + milk_protein * ED_protein);
        
        // } # Separation conditional
        
        
        // ------------------------------------------------------------
        // ENERGY COSTS
        // ------------------------------------------------------------
        
        // Metabolic scalar for RMR based on age class -- (d.u.)
        if(cohortID == 0) phi_rmr = 1.4; // Calves
        if(cohortID > 0 & cohortID <= 2) phi_rmr = 1.2; // Juveniles
        
        // Resting metabolic rate -- (MJ)
        // Calculated based on lean mass as blubber is more or less metabolically inert (George et al. 2021)
        resting_metab = RMR(2, animal->mass - animal->fatmass, phi_rmr);
        resting_metab_calves = RMR(2, calves[current_animal].mass - calves[current_animal].fatmass, phi_rmr);
        
        // Stroke rate during routine swimming -- (d.u.)
        stroke_rate = rtnorm(0.1505755, 0.01366726, 0, INFINITY);
        
        // Locomotory costs (J converted to MJ)
        Rcpp::NumericVector delta = {feeding_time, (24*3600) - feeding_time}; // Time spent in different activities
        Rcpp::NumericVector phi = {2, scalar_LC}; // Scalars on locomotory costs
        
        LC_tot = locomotor_costs(animal->mass, stroke_rate, delta, phi)/1000000;
        LC_tot_calves = locomotor_costs(calves[current_animal].mass, stroke_rate, delta, phi)/1000000;
        
        // FETAL DEVELOPMENT
        
        if(cohortID == 4){
        
        // Mass of muscles in fetus -- (kg)
        mass_muscle = fetal_tissue_mass(P_muscle, fetus_l[current_animal]);

        // Mass of viscera in fetus -- (kg)
        mass_viscera = fetal_tissue_mass(P_viscera, fetus_l[current_animal]);
        
        // Mass of bones in fetus -- (kg)
        mass_bones = fetal_tissue_mass(P_bones, fetus_l[current_animal]);

        // Mass of blubber in fetus -- (kg)
        mass_blubber = fetal_blubber_mass(fetus_l[current_animal], 
                                          mass_muscle, mass_viscera, 
                                          mass_bones, blubber_density, 
                                          muscle_density, 
                                          visceral_density, bone_density);
        
        // Energetic cost of fetal growth during pregnancy -- (MJ) [G]
        fetal_growth_cost = mass_muscle * (ED_lipids * muscle_lipids[current_animal] + ED_protein * muscle_protein[current_animal]) + 
          mass_viscera * (ED_lipids * visceral_lipids[current_animal] + ED_protein * visceral_protein[current_animal]) +
          mass_bones * (ED_lipids * bones_lipids[current_animal] + ED_protein * bones_protein) +
          mass_blubber * (ED_lipids * blubber_lipids[current_animal] +  ED_protein * blubber_protein[current_animal]);
        
        // Energetic cost of placental maintenance during pregnancy -- (kJ converted to MJ) [P]
        placental_cost = placental_maintenance(fetal_growth_cost)/1000;
        
        // Daily growth rate of the fetus (kg)
        fetus_m_next = fetal_mass(current_day - 365, animal->length, 0, 805.07);
        fetus_growth =  fetus_m_next - fetus_m[current_animal];
        
        // Heat increment of gestation -- (kJ converted to MJ) [Q]
        HIC = heat_gestation(birth_mass[current_animal], fetus_growth)/1000;
        
        // Total cost of gestation -- (kJ)
        E_gest = fetal_growth_cost + placental_cost + HIC;
        
        } // End if cohortID == 4
        
        // Cost of growth (adults and juveniles) -- (MJ)
        lean_mass_current = length2mass(animal->length, animal->massatlength, true);
        length_nextday = age2length(animal->age + 1/365, animal->lengthatage);
        lean_mass_nextday = length2mass(length_nextday, animal->massatlength, true);
        mass_increment = (lean_mass_nextday - lean_mass_current) / 365;
        
        E_growth = growth_cost(mass_increment, animal->bc, lean_water, blubber_lipids[current_animal],
                               ED_lipids, ED_protein, deposition_lipids, deposition_protein);
        
        
        // Cost of growth (calves) -- (MJ)
        lean_mass_current_calves = length2mass(calves[current_animal].length, calves[current_animal].massatlength, true);
        length_nextday_calves = age2length(calves[current_animal].age + 1/365, calves[current_animal].lengthatage);
        lean_mass_nextday_calves = length2mass(length_nextday_calves, calves[current_animal].massatlength, true);
        mass_increment_calves = (lean_mass_nextday_calves - lean_mass_current_calves) / 365;

        E_growth_calves = growth_cost(mass_increment_calves, calves[current_animal].bc, lean_water, blubber_lipids[current_animal],
                               ED_lipids, ED_protein, deposition_lipids, deposition_protein);
        
        // Cost of lactation - estimated from the energy requirements of the calf
        if(cohortID == 5) E_lac = LC_tot_calves + E_growth_calves;
        
        // ENERGY EXPENDITURE (Adults and juveniles -- MJ)
        E_out = (resting_metab + LC_tot + E_gest + E_lac + E_growth);
        if(animal->entangled(0) == 1) E_out = E_out + entanglement_cost;
        
        // ENERGY EXPENDITURE (Adults and juveniles -- MJ)
        E_out_calves = (resting_metab_calves + LC_tot_calves + E_growth_calves);
        if(animal->entangled(0) == 1) E_out_calves = E_out_calves + entanglement_cost;
        
        
      } // End is_alive
      
      // ------------------------------------------------------------
      // STORE VALUES in output matrices
      // ------------------------------------------------------------
      
      // Convert back to relevant units
      travel_time = travel_time / 3600;
      feeding_time = feeding_time / 3600;
      resting_time = resting_time / 3600;
      nursing_time = nursing_time / 3600;
      
      // +++ Animal traits +++
      
      *(attrib_out++) = animal->cohortID;                // Unique cohort identifier
      *(attrib_out++) = animal->alive;                   // Alive or dead
      *(attrib_out++) = animal->age;                     // Age
      *(attrib_out++) = animal->bc;                      // Body condition
      *(attrib_out++) = animal->length;                  // Total body length 
      *(attrib_out++) = animal->lengthatage(0,0);        // Coefficient of the Gompertz growth curve
      *(attrib_out++) = animal->lengthatage(0,1);        // Coefficient of the Gompertz growth curve
      *(attrib_out++) = animal->lengthatage(0,2);        // Coefficient of the Gompertz growth curve
      *(attrib_out++) = animal->mass;                    // Total body mass (including blubber)
      *(attrib_out++) = animal->fatmass;                 // Blubber mass
      *(attrib_out++) = animal->massatlength(0,0);       // Intercept of the logarithmic mass-at-length relationship
      *(attrib_out++) = animal->massatlength(0,1);       // Slope of the logarithmic mass-at-length relationship
      *(attrib_out++) = animal->mouth_ratio;             // Ratio of mouth width to total body length
      *(attrib_out++) = animal->mouth_angle;             // Upper angle of the mouth
      *(attrib_out++) = animal->mouth_width;             // Width of the mouth
      
      // +++ Stressors +++
      
      *(attrib_stressors_out++) = animal->entangled(0);  // Is the animal entangled?
      *(attrib_stressors_out++) = animal->entangled(1);  // Does the entanglement involve the anterior part of the body?
      *(attrib_stressors_out++) = animal->entangled(2);  // Entanglement severity
      *(attrib_stressors_out++) = animal->entangled(3);  // Entanglement duration
      *(attrib_stressors_out++) = start_entanglement;    // Day when entanglement event started
      *(attrib_stressors_out++) = end_entanglement;      // Day when entanglement event ended
      *(attrib_stressors_out++) = vessel_strike;         // Incidence of vessel strike
      *(attrib_stressors_out++) = resp_noise;            // Incidence of behavioral response to noise exposure
      *(attrib_stressors_out++) = response_dB;           // Behavioral response threshold
      
      // +++ Energy budget +++
      
      *(attrib_E_out++) = E_in - E_out;                  // Net energy (adults and juveniles)
      *(attrib_E_out++) = E_in;                          // Energy gain (adults and juveniles)
      *(attrib_E_out++) = E_out;                         // Energy expenditure (adults and juveniles)
      *(attrib_E_out++) = E_in_calves - E_out_calves;    // Net energy (calves)
      *(attrib_E_out++) = E_in_calves;                   // Energy gain (calves)
      *(attrib_E_out++) = E_out_calves;                  // Energy expenditure (calves)
      
      // +++ Energy intake +++
      
      *(attrib_feeding_out++) = is_feeding;                          // Incidence of foraging
      *(attrib_feeding_out++) = D_cop;                               // Prey concentration
      *(attrib_feeding_out++) = min_calanus;                         // Minimum prey concentration needed for foraging
      *(attrib_feeding_out++) = mouth_gape;                          // Energy gain (adults and juveniles)
      *(attrib_feeding_out++) = swim_speed_foraging;                 // Swimming speed while foraging
      *(attrib_feeding_out++) = capture_efficiency;                  // Capture efficiency
      *(attrib_feeding_out++) = gape_reduction;                      // % reduction in the gape when entangled
      *(attrib_feeding_out++) = daylight_hours;                      // Number of daylight hours for given location/season
      *(attrib_feeding_out++) = feeding_nursing_effort[0];           // Feeding/nursing effort
      *(attrib_feeding_out++) = target_bc[0];                        // Target body condition
      *(attrib_feeding_out++) = cop_mass;                            // Copepod mass
      *(attrib_feeding_out++) = cop_kJ;                              // Energy density of Calanus finmarchicus stage V
      *(attrib_feeding_out++) = digestive_efficiency;                // Digestive efficiency
      *(attrib_feeding_out++) = metabolizing_efficiency;             // Metabolizing efficiency
      *(attrib_feeding_out++) = E_calanus;                           // Prey energy content
      
      // +++ Energy costs +++
      
      *(attrib_costs_out++) = resting_metab;           // Resting metabolic rate
      *(attrib_costs_out++) = LC_tot;                  // Locomotory costs
      *(attrib_costs_out++) = stroke_rate;             // Stroke rate
      *(attrib_costs_out++) = E_growth;                // Energetic cost associated with growth
      *(attrib_costs_out++) = E_gest;                  // Energetic cost associated with gestation
      *(attrib_costs_out++) = fetal_growth_cost;       // Energetic cost associated with fetal growth
      *(attrib_costs_out++) = placental_cost;          // Energetic cost associated with placental maintenance
      *(attrib_costs_out++) = HIC;                     // Heat increment of gestation
      *(attrib_costs_out++) = E_lac;                   // Energetic cost associated with lactation
      *(attrib_costs_out++) = mass_increment;          // Daily change in mass
      
      // +++ Activity budgets +++
      
      *(attrib_activity_out++) = travel_dist;          // Distance traveled
      *(attrib_activity_out++) = swim_speed;           // Swimming speed
      *(attrib_activity_out++) = travel_time;          // Time spent traveling
      *(attrib_activity_out++) = feeding_time;         // Time spent feeding
      *(attrib_activity_out++) = nursing_time;         // Time spent nursing
      *(attrib_activity_out++) = resting_time;         // Time spent resting
      *(attrib_activity_out++) = n_zeroes;             
      *(attrib_activity_out++) = t_sum;                
      *(attrib_activity_out++) = t_remain;             
      
      // +++ Fetus development +++
      
      if(cohortID == 4){
        
        *(attrib_fetus_out++) = fetus_l[current_animal];
        *(attrib_fetus_out++) = fetus_m[current_animal];
        *(attrib_fetus_out++) = birth_length[current_animal];
        *(attrib_fetus_out++) = birth_mass[current_animal];
        *(attrib_fetus_out++) = mass_muscle;
        *(attrib_fetus_out++) = mass_viscera;
        *(attrib_fetus_out++) = mass_bones;
        *(attrib_fetus_out++) = mass_blubber;
        *(attrib_fetus_out++) = muscle_lipids[current_animal];
        *(attrib_fetus_out++) = muscle_protein[current_animal];
        *(attrib_fetus_out++) = visceral_lipids[current_animal];
        *(attrib_fetus_out++) = visceral_protein[current_animal];
        *(attrib_fetus_out++) = bones_lipids[current_animal];
        *(attrib_fetus_out++) = blubber_lipids[current_animal];
        *(attrib_fetus_out++) = blubber_protein[current_animal];

      }
      
      // +++ Calf traits +++
      
      if(cohortID == 5){
        
        *(attrib_calves_out++) = calves[current_animal].cohortID;
        *(attrib_calves_out++) = calves[current_animal].alive;
        *(attrib_calves_out++) = calves[current_animal].age;
        *(attrib_calves_out++) = calves[current_animal].bc;
        *(attrib_calves_out++) = calves[current_animal].length;
        *(attrib_calves_out++) = calves[current_animal].lengthatage(0,0);
        *(attrib_calves_out++) = calves[current_animal].lengthatage(0,1);
        *(attrib_calves_out++) = calves[current_animal].lengthatage(0,2);
        *(attrib_calves_out++) = calves[current_animal].mass;
        *(attrib_calves_out++) = calves[current_animal].fatmass;
        *(attrib_calves_out++) = calves[current_animal].massatlength(0,0);
        *(attrib_calves_out++) = calves[current_animal].massatlength(0,1);
        *(attrib_calves_out++) = calves[current_animal].mouth_ratio;
        *(attrib_calves_out++) = calves[current_animal].mouth_angle;
        *(attrib_calves_out++) = calves[current_animal].mouth_width;

        *(attrib_nursing_out++) = milk_assim;
        *(attrib_nursing_out++) = milk_provisioning;
        *(attrib_nursing_out++) = mammary_M;
        *(attrib_nursing_out++) = milk_production_rate;
        *(attrib_nursing_out++) = D_milk;
        *(attrib_nursing_out++) = t_suckling;
        *(attrib_nursing_out++) = target_bc[1];
        *(attrib_nursing_out++) = feeding_nursing_effort[1];
        *(attrib_nursing_out++) = milk_lipids;
        *(attrib_nursing_out++) = milk_protein;
        *(attrib_nursing_out++) = ED_lipids;
        *(attrib_nursing_out++) = ED_protein;
        
        *(attrib_costs_calves_out++) = resting_metab_calves;           // Resting metabolic rate
        *(attrib_costs_calves_out++) = LC_tot_calves;                  // Locomotory costs
        *(attrib_costs_calves_out++) = E_growth_calves;                // Energetic cost associated with growth
        *(attrib_costs_calves_out++) = mass_increment_calves;          // Daily change in mass
        
      }
      
    } // End loop over animals
  } // End loop over time points (days)
  
  results["locs"] = out;
  results["attrib"] = attrib;
  results["attrib_calves"] = attrib_calves;
  results["attrib_fetus"] = attrib_fetus;
  results["stressors"] = attrib_stressors;
  results["activity"] = attrib_activity;
  results["in.kj"] = attrib_feeding;
  results["in.kj_calves"] = attrib_nursing;
  results["out.kj"] = attrib_costs;
  results["out.kj_calves"] = attrib_costs_calves;
  results["E"] = attrib_E;
  
  return results;
}


// //' Retrieve nested list elements
// //' 
// //' @param l Input list
// //' @param loc
// //' @param month
// //' @return List element
// // [[Rcpp::export]]
// Eigen::MatrixXd accessLOL(Rcpp::List l, std::string loc, int month) {
//   Rcpp::List louter = l[loc];
//   return louter[month-1];
// }

// //' Convert integer to character string
// // [[Rcpp::export]]
// std::string int2char(int i) {
//   std::string c = std::to_string(i);
//   return(c);
// }
// 

// void my_func(double threshold){
//   bool obj = false;
//   int i = 0;
//   while(!obj){
//     i++;
//     std::cout << i << "\n";
//     double a = R::runif(0,1);
//     if(a > threshold) obj = true;
//   }
// }


//' Simulate right whale movements
 //' 
 //' @param densities array of density maps, stored in column-major order with 
 //' indices [xind, yind, map] for projected spaces (xind is index for a projected lon coordinate, 
 //' yind is index for a projected lat)
 //' @param densitySeq 1-indexed vector specifying which density maps should be
 //' used to simulate movement at each timepoint; number of timepoints to
 //' simulate comes from length of vector densitySeq
 //' @param latentDensitySeq 1-indexed vector specifying which density maps should be used to simulate movement for each latent animal
 //' @param limits vector (xmin, xmax, ymin, ymax) of spatial coordinate extents for spatial densities
 //' @param resolution vector (xres, yres) of spatial step sizes for spatial densities
 //' @param M Number of proposals used in the importance sampler for movement (Michelot, 2019)
 //' @param stepsize Rradius of the proposal circle for movement (Michelot, 2019)
 //' @param xinit matrix of initial x coordinates for latent animals (nlatent x n)
 //' @param yinit matrix of initial y coordinates for latent animals (nlatent x n)
 //' @return List of coordinates
 // [[Rcpp::export]]
 
 Rcpp::List MovementSimulator(
     int cohortID,
     Eigen::MatrixXd support,
     std::vector<Eigen::MatrixXd> densities,
     std::vector<std::size_t> densitySeq,
     std::vector<std::size_t> latentDensitySeq,
     std::vector<Eigen::MatrixXd> prey,
     std::vector<Eigen::MatrixXd> fishing,
     std::vector<Eigen::MatrixXd> vessels,
     std::vector<Eigen::MatrixXd> noise,
     Rcpp::NumericVector doseresp,
     std::vector<Eigen::MatrixXd> daylight,
     Eigen::MatrixXd regions,
     Eigen::VectorXd limits, 
     Eigen::VectorXd resolution,
     std::size_t M,
     double stepsize, 
     Eigen::MatrixXd xinit, 
     Eigen::MatrixXd yinit,
     bool progress
 ) {
   
   // ------------------------------------------------------------------------
   // 0-index variables
   // ------------------------------------------------------------------------
   
   // Coming from R, we assume latentDensitySeq is 1-indexed
   // .begin = Returns an iterator pointing to the first element in the sequence
   // .end = Returns an iterator pointing to the last element in the sequence
   for(auto d = latentDensitySeq.begin(); d != latentDensitySeq.end(); ++d)
     --(*d);
   
   // Ditto for densitySeq
   for(auto d = densitySeq.begin(); d != densitySeq.end(); ++d)
     --(*d);
   
   // ------------------------------------------------------------------------
   // Initialize environments for latent animals
   // ------------------------------------------------------------------------
   // emplace_back: inserts a new element at the end of a vector
   
   std::vector<Environment> latent_environments;
   auto dend  = latentDensitySeq.end();
   
   for(auto d = latentDensitySeq.begin(); d != dend; ++d)
     latent_environments.emplace_back(densities[*d], prey[*d], fishing[*d], vessels[*d], noise[*d], daylight[*d], regions,
                                      limits, resolution, *d);
   
   // Initialize prey base environments
   // std::vector<PreyEnv> prey_environments;
   // auto pend  = latentDensitySeq.end();
   // for(auto p = latentDensitySeq.begin(); p != pend; ++p)
   //   prey_environments.emplace_back(prey[*p], limits, resolution, *p);
   
   // ------------------------------------------------------------------------
   // Initialize environments with attraction to latent animals
   // ------------------------------------------------------------------------
   
   std::vector<LatentAttractor> environments;
   auto dattrend = densitySeq.end();
   for(auto d = densitySeq.begin(); d != dattrend; ++d)
     environments.emplace_back(*d);
   
   // ------------------------------------------------------------------------
   // Initialize population
   // ------------------------------------------------------------------------
   
   std::vector<CouplingAnimal> animals;
   std::vector<CouplingAnimal> calves;
   
   double *xiter = xinit.data();
   double *yiter = yinit.data();
   
   // Loop over animals
   for(std::size_t j = 0; j < xinit.cols(); ++j) {
     std::vector<Animal> latent_animals;
     std::vector<Animal> latent_calves;
     
     // Loop over months
     for(std::size_t i = 0; i < xinit.rows(); ++i) {
       latent_animals.emplace_back(Animal(*(xiter++), *(yiter++), cohortID));
       latent_calves.emplace_back(Animal());
     }
     animals.push_back(CouplingAnimal(latent_animals, densitySeq[0]));
     calves.push_back(CouplingAnimal());
   }
   
   // ------------------------------------------------------------------------
   // Initialize latent movement rule
   // ------------------------------------------------------------------------
   
   // typedef gives a type a new name
   typedef ReweightedRandomWalk<Animal, Environment> LatentMvmt;
   LatentMvmt rm(M, stepsize);
   
   // Initialize observed movement rule
   CoupledRandomWalk<
     CouplingAnimal, LatentAttractor, LatentMvmt, std::vector<Environment>
   > crm(rm, latent_environments, M, stepsize);
   
   // ------------------------------------------------------------------------
   // Run simulations
   // ------------------------------------------------------------------------
   
   return movesim(cohortID, animals, calves, environments, 
                  latent_environments, doseresp, densitySeq, crm, support, 
                  limits, resolution, stepsize, progress);
   
 }

/**
 * Evaluate the environment at the given x and y coordinates
 *
 * @param density density map, stored in column-major order with
 *   indices [xind, yind] for projected spaces (xind is index for a
 *   projected lon coordinate, yind is index for a projected lat)
 * @param limits vector (xmin, xmax, ymin, ymax) of spatial coordinate extents
 *   for spatial densities
 * @param resolution vector (xres, yres) of spatial step sizes for spatial
 *   densities
 * @param x x coordinate at which to query the density map
 * @param y y coordinate at which to query the density map
 * @return the density map value at location (x,y)
 */
// [[Rcpp::export]]
double evalEnvironment(Eigen::MatrixXd density, 
                       Eigen::MatrixXd prey, 
                       Eigen::MatrixXd fishing, 
                       Eigen::MatrixXd vessels, 
                       Eigen::MatrixXd noise, 
                       Eigen::MatrixXd daylight, 
                       Eigen::MatrixXd regions, 
                       Eigen::VectorXd limits, 
                       Eigen::VectorXd resolution,
                       double x, double y, char layer) {
  Environment e(density, prey, fishing, vessels, noise, daylight, regions, limits, resolution, 0);
  return e(x, y, layer);
}

// Version with regions
// Eigen::VectorXd evalEnvironment(Eigen::MatrixXd density, 
//                                 Eigen::MatrixXd regions,
//                                 Eigen::VectorXd limits, 
//                                 Eigen::VectorXd resolution,
//                                 double x, double y) {
//   Environment e(density, regions, limits, resolution, 0);
//   return e(x, y);
// }
// 