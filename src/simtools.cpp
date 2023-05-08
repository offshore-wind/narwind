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
  int seus;
  int gsl;
  int alive;
  long double age;
  long double bc;
  Eigen::MatrixXd lengthatage;
  long double length;
  Eigen::MatrixXd massatlength;
  long double tot_mass;
  long double fat_mass;
  long double lean_mass;
  double mouth_ratio;
  double mouth_angle;
  double mouth_width;
  Rcpp::NumericVector entangled;
  int strike;
  int respnoise;
  int abort;

  
  // Member initializer list - age initialized by default at birth
  Animal() : 
    x(0), y(0), 
    cohortID(0), 
    seus(0),
    gsl(0), 
    alive(cohortID >=1), 
    age(0),
    bc(start_bcondition(cohortID)),
    lengthatage(agL(age)), 
    length(age2length(age, lengthatage)),
    massatlength(mL()), 
    tot_mass(length2mass(length, massatlength, false)),
    fat_mass(bc*tot_mass),
    lean_mass(tot_mass - fat_mass), 
    mouth_ratio(start_mouth(cohortID, age)), 
    mouth_angle(76.7), 
    mouth_width(length*mouth_ratio), 
    entangled(Rcpp::NumericVector::create(0,0,0,0,0,0)),
    strike(0), respnoise(0),
    abort(0){ }
  
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
         int cohort,
         int to_seus,
         int to_gsl) : 
    x(x0), y(y0), 
    cohortID(cohort),
    seus(to_seus),
    gsl(to_gsl),
    alive(cohortID >=1),
    age(start_age(cohort)),
    bc(start_bcondition(cohort)),
    lengthatage(agL(age)), 
    length(age2length(age, lengthatage)),
    massatlength(mL()),
    tot_mass(length2mass(length, massatlength, false)),
    fat_mass(bc*tot_mass),
    lean_mass(tot_mass - fat_mass),
    mouth_ratio(start_mouth(cohort, age)), 
    mouth_angle(76.7), 
    mouth_width(length*mouth_ratio),
    entangled(Rcpp::NumericVector::create(0,0,0,0,0,0)),
    strike(0), respnoise(0),
    abort(0){ }

  
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
  init_latent[init_animal].cohortID,
  init_latent[init_animal].seus,
  init_latent[init_animal].gsl),
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
         int CoordSize = 3,
         int nparam_E = 8,
         int nparam_E_calves = 4,
         int nparam_animal = 19,
         int nparam_fetus = 17,
         int nparam_stressors = 12,
         int nparam_activity = 9,
         int nparam_feeding = 15,
         int nparam_nursing = 12,
         int nparam_costs = 11,
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
    bool stressors,
    bool growth,
    Rcpp::NumericVector cumgest,
    bool progress) {
  
  // Generate hash table of geodesic distances 
  // phmap::flat_hash_map<std::int64_t, int> geodesic_table;
  // for (int i = 0; i < hashID.rows(); i++) {
  //   geodesic_table[hashID(i,0)] = hashID(i,1);
  // }

  // Number of animals and time points to simulate
  std::size_t n = animals.size();  // nsim
  std::size_t t = environments.size(); // 365 days
  
  // Initialize list in which results will be stored
  Rcpp::List results;
  
  // Initialize and label array in which x,y coordinates will be stored
  Rcpp::NumericVector out(n*t*CoordSize);
  out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
  out.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("easting", "northing", "region"), R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which bioenergetic parameters will be stored
  Rcpp::NumericVector attrib(n*t*nparam_animal);
  attrib.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  attrib.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c", 
                                  "mass", "leanmass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", 
                                  "mouth_w", "gsl", "seus", "abort"),
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for dependent calves
  Rcpp::NumericVector attrib_calves(n*t*(nparam_animal-3));
  attrib_calves.attr("dim") = Rcpp::Dimension(nparam_animal-3,n,t);
  attrib_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort_calf", "alive_calf", "age_calf", "bc_calf", "length_calf", 
                                  "La_calf", "Lb_calf", "Lc_calf", 
                                  "mass_calf", "leanmass_calf", "fatmass_calf", 
                                  "ma_calf", "mb_calf", "mouth_r_calf", 
                                  "mouth_a_calf", "mouth_w_calf"), 
                                  R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which energy budget will be stored
  Rcpp::NumericVector attrib_E(n*t*nparam_E);
  attrib_E.attr("dim") = Rcpp::Dimension(nparam_E,n,t);
  attrib_E.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_tot", "E_in", "E_out", "delta_fat", "DE_lip", "ED_lip", "lip_anab", "lip_catab"),
    R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which energy budget will be stored
  Rcpp::NumericVector attrib_E_calves(n*t*nparam_E_calves);
  attrib_E_calves.attr("dim") = Rcpp::Dimension(nparam_E_calves,n,t);
  attrib_E_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_tot_calf", "E_in_calf", "E_out_calf", "delta_fat_calf"),
    R_NilValue, R_NilValue
  );
  

  // Create a separate array for fetuses
  Rcpp::NumericVector attrib_fetus(n*t*nparam_fetus);
  attrib_fetus.attr("dim") = Rcpp::Dimension(nparam_fetus,n,t);
  attrib_fetus.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("fetus_l", "fetus_m", "delta_fetus_m", "delta_fetus_l", "birth_l", "birth_m", "muscle_m", "viscera_m",
                                  "bones_m", "blubber_m", "muscle_lip", "muscle_pro", "visc_lip", "visc_pro", "bone_lip",
                                  "blubber_lip", "blubber_pro"), 
                                  R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which stressor data will be stored
  Rcpp::NumericVector attrib_stressors(n*t*nparam_stressors);
  attrib_stressors.attr("dim") = Rcpp::Dimension(nparam_stressors,n,t);
  attrib_stressors.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("gear_risk", "is_entgl", "entgl_head", "severity", "entgl_d",
                                  "entgl_start", "entgl_end", "strike_risk", "strike",
                                  "noise_resp", "noise_lvl", "dB_thresh"),
                                  R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which activity budgets will be stored
  Rcpp::NumericVector attrib_activity(n*t*nparam_activity);
  attrib_activity.attr("dim") = Rcpp::Dimension(nparam_activity,n,t);
  attrib_activity.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("d_travel", "swimspeed", "t_travel", "t_feed", "t_nurse", "t_rest", "n_zero", "t_sum", "t_remain"),
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_feeding(n*t*nparam_feeding);
  attrib_feeding.attr("dim") = Rcpp::Dimension(nparam_feeding,n,t);
  attrib_feeding.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("feed", "preyconc", "minprey", "gape", "feedspeed", "captEff", "impedance", "daylight",
                                  "feed_effort", "targetBC", "cop_mass", "cop_kJ", "digestEff", "metabEff", "E_cop"), 
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in calves
  Rcpp::NumericVector attrib_nursing(n*t*nparam_nursing);
  attrib_nursing.attr("dim") = Rcpp::Dimension(nparam_nursing,n,t);
  attrib_nursing.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("assim", "provision", "mamm_M", "milk_rate", "Dmilk", "t_suckling",
                                  "targetBC_calf", "nursing", "milk_lip", "milk_pro", "EDlip", "EDpro"),
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy costs in adults/juveniles
  Rcpp::NumericVector attrib_costs(n*t*nparam_costs);
  attrib_costs.attr("dim") = Rcpp::Dimension(nparam_costs,n,t);
  attrib_costs.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("rmr", "LC", "scalar_LC", "stroke", "E_growth", "E_gest", "fgrowth",
                                  "placenta", "hic", "E_lac", "delta_m"),
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy costs in calves
  Rcpp::NumericVector attrib_costs_calves(n*t*nparam_costs_calves);
  attrib_costs_calves.attr("dim") = Rcpp::Dimension(nparam_costs_calves,n,t);
  attrib_costs_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("rmr_calf", "LC_calf","scalar_LC", "E_growth_calf", "delta_m_calf"),
    R_NilValue, R_NilValue
  );
  

  Rcpp::NumericVector::iterator data_out = out.begin();
  Rcpp::NumericVector::iterator attrib_out = attrib.begin();
  Rcpp::NumericVector::iterator attrib_E_out = attrib_E.begin();
  Rcpp::NumericVector::iterator attrib_E_calves_out = attrib_E_calves.begin();
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
  
  double gear_risk;
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
  double stroke_rate_foraging, stroke_rate, percent_glide_foraging, percent_glide;
  double scalar_LC;
  
  // if(cohortID == 4){
  //   scalar_LC = 1.4; // Pregnant females
  // } else if(cohortID == 5){
  //   scalar_LC = 1.7; // Lactating females
  // } else if(cohortID == 3) {
  //   scalar_LC = 0.96; // Adult males
  // } else if(cohortID == 1 | cohortID == 2) {
  //   scalar_LC = 1.1;
  // } else if(cohortID == 0) {
  //   scalar_LC = 0.98;
  // } else {
  //   scalar_LC = 1;
  // }
  
  double LC_tot, LC_tot_calves;
  
  Rcpp::NumericVector birth_length (n);
  Rcpp::NumericVector birth_mass (n);
  
  // if(cohortID == 4){
  //   
  //   // Length and mass at birth (only relevant for pregnant females) -- (m)
  //   Eigen::MatrixXd birth_agL = agL(0, n);
  //   Eigen::MatrixXd birth_mL = mL(n);
  //   for (int i = 0; i < n; i++) {
  //     birth_length[i] = age2length(0, birth_agL.row(i));
  //     birth_mass[i] = length2mass(birth_length[i], birth_mL.row(i), false);
  //   }
  // }
  
  if(cohortID == 4){
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {
      int i = animal - animals.begin();
      birth_length[i] = fetal_length(0, animal->length);
      birth_mass[i] = fetal_mass(0, animal->length);
    }
  }

  
  double fetus_l, fetus_m;
  double bc_gest;
  
  double mass_muscle, mass_muscle_next;
  double mass_viscera, mass_viscera_next;
  double mass_bones, mass_bones_next;
  double mass_blubber, mass_blubber_next;
  double FG0, FG1;
  double fetal_growth_cost = 0;
  double placental_cost = 0;
  double mass_increment_fetus = 0;
  double length_increment_fetus = 0;
  double fetus_m_next, fetus_l_next;
  double HIC = 0;
  double E_gest = 0;
  
  double E_lac = 0;
  
  double length_nextday, length_nextday_calves;
  double lean_mass_nextday, lean_mass_nextday_calves;
  double mass_increment, mass_increment_calves;
  double E_growth, E_growth_calves;
  double E_out, E_out_calves;
  int E_balance, E_balance_calves;
  
  // Range of dose values for the dose-response function
  // int dose_lwr = 85;
  // int dose_uppr = 200;
  double response_dB = 0;
  // std::vector<double> dose = dose_range(dose_lwr, dose_uppr, 1000);
  
  // Percent of protein synthesis during anabolism / breakdown during catabolism
  Rcpp::NumericVector lipid_anacat(2);
  Rcpp::NumericVector lipid_anacat_calves(2);
  lipid_anacat[0] = 0.8; // Lipid anabolism -- (d.u.)
  lipid_anacat_calves[0] = 0.8;
  
  double delta_blubber, delta_blubber_calves;
  
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
  const double ED_lipids = 39300.0L/1000;
  
  // Energy density of protein -- (kJ / kg - converted to MJ)
  // mean(c(23600,18000))
  const double ED_protein = 20800.0L/1000;
  
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
  
  // Rcpp::NumericVector fetus_l (n);
  // Rcpp::NumericVector fetus_m (n);
  
  // ------------------------------------------------------------
  // START SIMULATION
  // ------------------------------------------------------------
  
  // 365 steps per animal
  for(auto env = environments.begin(); env != env_end; ++env) {

    current_day = env - environments.begin();
    
    if(progress) bar.update();
    
    // env is a vector iterator, so env - environments.begin() returns the associated index
    current_month = densitySeq[env - environments.begin()];
    
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {

      current_x = animal->x;
      current_y = animal->y;
      
      current_animal = animal - animals.begin();
      current_region = layers[current_month](current_x, current_y, 'R');

      // std::cout << animal->seus << " - " << animal->gsl << std::endl;
      
      // Coordinates coordinates (fills by column)
      
      *(data_out++) = animal->x;
      *(data_out++) = animal->y;
      *(data_out++) = current_region;
      
      // ------------------------------------------------------------
      // TARGET BC
      // ------------------------------------------------------------
      
      if(animal->cohortID == 4){  // Pregnant females
        
        target_bc[0] = 0.6439795;
        target_bc[1] = 0;
        
      } else if(animal->cohortID == 5){ // Lactating mothers and their calves
        
        target_bc[0] = 0.5271131;
        target_bc[1] = 0.6439795;
        
      } else { // All other cohorts
        
        target_bc[0] = 0.5271131;
        target_bc[1] = 0;
        
      }
      
      // ------------------------------------------------------------
      // VESSEL STRIKES
      // ------------------------------------------------------------
      
      if(current_day > 0){
        
        if(stressors){
          
          if(animal->alive){
            
            // Vessel strikes +++
            strike_risk = layers[current_month](current_x, current_y, 'V');
            
            // Incidence of a vessel strike -- (d.u.)
            animal->strike = R::rbinom(1, strike_risk);
            
          } // End if alive
        } // End if stressors

        // ------------------------------------------------------------
        // MORTALITY
        // ------------------------------------------------------------
        
        if(animal->strike){
          animal->alive = 0;
          calves[current_animal].alive = 0;
        }
        
      } // End if current_day > 0

      if(animal->alive){
        
        if(current_day > 0){
        
        // ------------------------------------------------------------
        // COORDINATES
        // ------------------------------------------------------------
        
        // animal->x = coords_xy[0];
        // animal->y = coords_xy[1];
        
        // // Save current coordinates (fills by column)
        // 
        // *(data_out++) = animal->x;
        // *(data_out++) = animal->y;
        
        // bool calc_geo = false;
        // int gd;
        // int *gdptr;
        // gdptr = &gd;
        
        m.update(*animal, *env, TRUE, current_region, support, limits, resolution);
          
        // // Only activate geodesic calculations in some areas
        // if((current_x >= 685 & current_x <= 1420 & current_y >= 1070 & current_y <= 1865) ||
        //    (current_x >= 520 & current_x <= 730 & current_y >= 825 & current_y <= 960)) calc_geo = true;
        // 
        // if(!calc_geo){
        //   
        //   // Update the animal's state
        //   m.update(*animal, *env, TRUE, current_region, support, limits, resolution);
        //   
        // } else {
        //   
        //   while(calc_geo == true){
        //     
        //     // Update the animal's state
        //     m.update(*animal, *env, TRUE, current_region, support, limits, resolution);
        //     
        //     double next_x = animal->x;
        //     double next_y = animal->y;
        //     
        //     // Calculate geodesic distance
        //     *gdptr = geoD(support, current_x, current_y, next_x, next_y, limits, resolution);
        //     
        //     if(*gdptr >= 0 & *gdptr < (2*stepsize)){
        //       
        //       calc_geo = false;
        //       
        //     } else {
        //       
        //       animal->x = current_x;
        //       animal->y = current_y;
        //     }
        //     
        //   } // End while calc_geo
        //   
        // } // End if(calc_geo)
        
        if(cohortID == 5 & current_region == 9){
          
          if(animal->y > current_y) calves[current_animal].alive = 1;

        }
        
        travel_dist = std::sqrt(std::pow(current_x - animal->x, 2) +  std::pow(current_y - animal->y, 2));
        
        // Percent lipid breakdown during catabolism
        // See Beltran et al. (2017)
        lipid_anacat[1] = R::runif(0.53, 0.6);
        
        // ------------------------------------------------------------
        // NOISE
        // ------------------------------------------------------------
        
        if(stressors){
          
          animal->respnoise = 0; // Reset daily
          
          // Pile-driving noise +++
          dB = layers[current_month](current_x, current_y, 'N');
          
          // Response threshold (dB)
          response_dB = response_threshold(doseresp);
          
          // Behavioral response to pile-driving noise exposure -- (d.u.)
          if(dB > response_dB) animal->respnoise = 1;
          
          // if(dB <= dose_lwr){
          //   p_response = 0;
          // } else if(dB >= dose_uppr){
          //   p_response = 1;
          // } else {
          //   p_response = prob_response(dose, doseresp, animal->doseID, dB);
          // }
          // animal->respnoise = R::rbinom(1, p_response);
          
          // ------------------------------------------------------------
          // FISHING GEAR
          // ------------------------------------------------------------
          
          // An animal can only become entangled if not carrying gear already 
        
          if(animal->entangled(0) == 0){
            
            gear_risk = layers[current_month](current_x, current_y, 'F');
            animal->entangled = entanglement_event(gear_risk);
            
            if(animal->entangled(0) == 1){
              animal->entangled(4) = current_day;
              animal->entangled(5) = current_day + animal->entangled(3);
            }
            
          } else {

            // Terminate entanglement event
            if(current_day == animal->entangled(5)){
              animal->entangled = entanglement_event(0); // Reset entanglement state
            }
            
          } // End if entangled

        } // End if stressors
        
        // ------------------------------------------------------------
        // PREY
        // ------------------------------------------------------------
        
        // Prey concentration in cell at time t -- (copepods/cubic m)
        D_cop = layers[current_month](current_x, current_y, 'P');
        
        // Minimum prey concentration (Calanus finmarchicus stage V) -- (copepods/cubic m)
        // Meta: From year 1987, by species
        min_calanus = rtnorm(1784.979, 718.77, 0, INFINITY);
        
        // Determine whether the prey concentration is sufficient to support foraging -- (d.u.)
        if(animal->respnoise){
          is_feeding = 0;
        } else {
          is_feeding = feeding_threshold(min_calanus, D_cop);
        }
        
        // ------------------------------------------------------------
        // ACTIVITY BUDGET
        // ------------------------------------------------------------
        
        // Total time available for foraging activities --(hr converted to s)
        daylight_hours = layers[current_month](current_x, current_y, 'L');
        
        // Swimming speed (during travel)
        // 4.6 m/s taken as maximum average speed based on satellite tag record described in Mate (1997)
        // Minimum chosen to ensure that travel_time does not exceed 24 hours
        swim_speed = rtnorm(0.8119258, 0.4530960, (travel_dist*1000)/(24*3600), 4.6) * 3600; // m per hr
        
        // Time spent traveling per day -- (hr)
        // Calculated based on distance covered that day and swimming speed
        travel_time = travel_dist * 1000 / swim_speed;
        
        // Based largely on Cusano et al. (2019)
        if(current_region == 9){ // SEUS
          
          feeding_time = 0;
          
          // Time spent resting + nursing per day -- (hr)
          // Meta: By age class, combine calves and mother/calf pairs, bounds of (0,24)
          if(cohortID == 5){ // Lactating females
            resting_time = rtnorm(13.76, 0.76984, 0, 24 - travel_time);
            nursing_time = rtnorm(1.6342, 4.6767, 0, 24 - (resting_time + travel_time)); 
          } else {
            resting_time = rtnorm(5.487199, 0.8499895, 0, 24 - travel_time);
            nursing_time = 0;
          }
          
        } else { // Feeding grounds and elsewhere
          
          tmax = std::min(24 - travel_time, daylight_hours);
          feeding_time = is_feeding * rtnorm(12.768, 2.9166, 0, tmax);
          resting_time = rtnorm(5.487199, 0.8499895, 0, 24 - (feeding_time + travel_time));
          
          if(cohortID == 5){ 
            nursing_time = 24 - (feeding_time + travel_time + resting_time);
          } else {
            nursing_time = 0;
          }
        }
        
        // Allocate any remaining time equally to each behavior
        t_sum = travel_time + resting_time + nursing_time + feeding_time;
        std::vector<double> t = {travel_time, resting_time, nursing_time, feeding_time};
        n_zeroes = std::count(t.begin(), t.end(), 0);
        
        if(t_sum < 24){
          t_remain = (24 - t_sum)/(4-n_zeroes);
        }
        
        // int can_feed = 0;
        // if(D_cop > min_calanus) can_feed = 1;
        // std::cout << current_region << " - " << animal->respnoise << " - " << response_dB << " - " << dB << " - " << daylight_hours << " - " << is_feeding << " - " << can_feed << " - " << min_calanus << " - " << D_cop << " - " << travel_time << " - " << 
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
          
          // Values for 5-parameter logistic curve obtained through optimization using target BC
          // passed to optim_feeding
          feeding_nursing_effort[0] = scale_effort(animal->bc, 1, 0, 27.36939, 0.4258282 , 0.5918155);
          
        } else if(animal->cohortID == 5){ // Lactating mothers and their calves
          
          // Lactating mothers fast in the initial months and only resume feeding midway through lactation (Miller et al. 2012)
          feeding_nursing_effort[0] = scale_effort(animal->bc, 1, 0, 24.18117, 0.2219739, 1.284047);
          feeding_nursing_effort[1] = scale_effort(calves[current_animal].bc, 1, 0, 27.36939, 0.4258282 , 0.5918155);
          
        } else { // All other cohorts
          
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
        E_in = (1 - animal->respnoise) * is_feeding * D_cop * mouth_gape * (1 - gape_reduction) *
          swim_speed_foraging * capture_efficiency * feeding_time * feeding_nursing_effort[0] *
          E_calanus;
        
        // std::cout << current_x << ";" << current_y << " - " << current_month << " - " << D_cop << " - " << min_calanus << " " << is_feeding << " " << E_in<< std::endl;
        
        // ------------------------------------------------------------
        // ENERGY INTAKE - SUCKLING
        // ------------------------------------------------------------
        
        if(cohortID == 5 && calves[current_animal].alive == 1){
          
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
          milk_provisioning = milk_supply(starvation, target_bc[0], animal->tot_mass, animal->fat_mass, zeta);
          
          // Total mass of mammary glands -- (kg)
          mammary_M = mammary_mass(animal->tot_mass);
          
          // Milk production rate by the mother -- (kg per sec)
          milk_production_rate = milk_production(mammary_M);
          
          // Density of milk -- (kg/m3)
          // Meta: No temporal filter, all records
          D_milk = R::rnorm(1.02, 0.01) * 1000; // L to cubic m
          
          // Time spent nursing/suckling per day -- (hours to sec)
          // Meta:: No filter on time spent nursing, time spent suckling during nursing bouts filtered by species (value for E. australis)
          t_suckling = rtnorm(0, 3.688672, 0, 24) * 3600;
          
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
          E_in_calves = (1 - animal->respnoise) * milk_assim * milk_provisioning * mammary_efficiency *
            (milk_production_rate/D_milk) * t_suckling * feeding_nursing_effort[1] * D_milk *
            (milk_lipids * ED_lipids + milk_protein * ED_protein);
          
          // std::cout<<animal->respnoise << "|"<<milk_assim << "|"<< milk_provisioning<< "|"<<mammary_efficiency
          //   << "|"<< mammary_M<< "|"<< milk_production_rate<< "|"<< t_suckling<< "|"<< 
          // feeding_nursing_effort[1]<< "|"<< D_milk<< "|"<< milk_lipids * ED_lipids + milk_protein * ED_protein <<std::endl;
          // 
          // } # Separation conditional
          
        }
        
        // ------------------------------------------------------------
        // LOCOMOTORY COSTS
        // ------------------------------------------------------------
        
        // Metabolic scalar for RMR based on age class -- (d.u.)
        // 1.4 increase in RMR for calves to account for elevated demands associated with active growth?
        // if(cohortID <= 2) phi_rmr = 1.2; else phi_rmr = 1;
        
        // Resting metabolic rate -- (MJ)
        // Calculated based on lean mass as blubber is more or less metabolically inert (George et al. 2021)
        resting_metab = RMR(animal->lean_mass);
        if(cohortID == 5) resting_metab_calves = RMR(calves[current_animal].lean_mass); 
        
        // Stroke rate during routine swimming -- (d.u.)
        stroke_rate_foraging = rtnorm(0.1635468, 0.004291999, 0, INFINITY);
        stroke_rate = rtnorm(0.1051259, 0.02903964, 0, INFINITY);
        
        // Percentage of time spent gliding -- (d.u.)
        percent_glide_foraging = R::rgamma(17.393772, 0.02006944); // Corresponds to mean of 36, min of 10 and max of 75%
        percent_glide = R::rnorm(0.09, 0.00795); // Corresponds to mean of 9, min of 9 and max of 12%
        
        // Locomotory costs (J converted to MJ)
        // 
        // Based on Noren 2008 (DOI: 10.1111/j.1365-2435.2007.01354.x) and Noren et al. (2011) (DOI: 10.1242/jeb.059121):
        // Apply 17% increase in stroke frequency for lactating females with 0-1 month-old calves swimming in echelon position.
        // Increase stroke rate by 16.2791% (1/0.86) for pregnant females in their last month of pregnancy.

        if(cohortID == 4 & current_month == densitySeq[365]){ 
          // Increase locomotory costs in the last month of pregnancy
          scalar_LC = 1/0.86;
        } else if(cohortID == 5 & current_month == densitySeq[0]){
          // Increase locomotory costs for lactating mothers swimming in echelon position
          scalar_LC = 1.17;
        } else {
          scalar_LC = 1;
        }
        
        LC_tot = locomotor_costs(animal->tot_mass, stroke_rate, stroke_rate_foraging, percent_glide, percent_glide_foraging,
                                 feeding_time, (24*3600) - feeding_time, scalar_LC);
        
        if(cohortID == 5){
          LC_tot_calves = locomotor_costs(calves[current_animal].tot_mass, stroke_rate, stroke_rate_foraging, percent_glide, percent_glide_foraging,
                                          0, 86400, scalar_LC);
        }
        
        // ------------------------------------------------------------
        // FETAL DEVELOPMENT
        // ------------------------------------------------------------

        // Abortion
        if(cohortID == 4){
          
          bc_gest = ((1000*cumgest[current_day])/(ED_lipids*lipid_anacat[1]) + 0.05*animal->tot_mass) / animal->tot_mass;
          if(bc_gest > animal->bc) animal->abort = 1;
          std::cout << animal->bc << " == " << bc_gest << std::endl;
        }
        
        if(cohortID == 4 && animal->abort == 0){ // Pregnant females
          
          // current_day goes from 2 (first day) to 366 (last day)
          
          if(current_day == 365){
            
            fetal_growth_cost = 0.0L;
            placental_cost = 0.0L;
            HIC = 0.0L;
            E_gest = 0.0L;
            
          } else {
          
          fetus_l = fetal_length(current_day - 365, animal->length);
          fetus_m = fetal_mass(current_day - 365, animal->length);
          
          fetus_l_next = fetal_length(current_day - 364, animal->length);
          fetus_m_next = fetal_mass(current_day - 364, animal->length);

          // Daily growth rate of the fetus
          length_increment_fetus = fetus_l_next - fetus_l;
          mass_increment_fetus = fetus_m_next - fetus_m;
          
          // Mass of muscles in fetus -- (kg)
          mass_muscle = fetal_tissue_mass(P_muscle, fetus_l);
          mass_muscle_next = fetal_tissue_mass(P_muscle, fetus_l_next);

          // Mass of viscera in fetus -- (kg)
          mass_viscera = fetal_tissue_mass(P_viscera, fetus_l);
          mass_viscera_next = fetal_tissue_mass(P_viscera, fetus_l_next);

          // Mass of bones in fetus -- (kg)
          mass_bones = fetal_tissue_mass(P_bones, fetus_l);
          mass_bones_next = fetal_tissue_mass(P_bones, fetus_l_next);

          // Mass of blubber in fetus -- (kg)
          mass_blubber = fetal_blubber_mass(fetus_l,
                                            mass_muscle,
                                            mass_viscera,
                                            mass_bones,
                                            blubber_density,
                                            muscle_density,
                                            visceral_density,
                                            bone_density);
          
          mass_blubber_next = fetal_blubber_mass(fetus_l_next,
                                            mass_muscle_next,
                                            mass_viscera_next,
                                            mass_bones_next,
                                            blubber_density,
                                            muscle_density,
                                            visceral_density,
                                            bone_density);

          // Energetic cost of fetal growth during pregnancy -- (MJ) [G]
          FG0 = mass_muscle * (ED_lipids * muscle_lipids[current_animal] + ED_protein * muscle_protein[current_animal]) +
            mass_viscera * (ED_lipids * visceral_lipids[current_animal] + ED_protein * visceral_protein[current_animal]) +
            mass_bones * (ED_lipids * bones_lipids[current_animal] + ED_protein * bones_protein) +
            mass_blubber * (ED_lipids * blubber_lipids[current_animal] +  ED_protein * blubber_protein[current_animal]);
          
          FG1 = mass_muscle_next * (ED_lipids * muscle_lipids[current_animal] + ED_protein * muscle_protein[current_animal]) +
            mass_viscera_next * (ED_lipids * visceral_lipids[current_animal] + ED_protein * visceral_protein[current_animal]) +
            mass_bones_next * (ED_lipids * bones_lipids[current_animal] + ED_protein * bones_protein) +
            mass_blubber_next * (ED_lipids * blubber_lipids[current_animal] +  ED_protein * blubber_protein[current_animal]);
          
          fetal_growth_cost = FG1 - FG0;
          
          // fetal_growth_cost = fetal_growth(fetus_l, 
          //                                  fetus_l_next,
          //                                  P_muscle, 
          //                                  P_viscera, 
          //                                  P_bones,
          //                                  muscle_density, 
          //                                  Dens_v = visceral_density, 
          //                                  Dens_bo = bone_density, 
          //                                  Dens_bl = blubber_density,
          //                                  Lip_m = muscle_lipids, 
          //                                  Lip_v = visceral_lipids, 
          //                                  Lip_bo = bones_lipids, 
          //                                  Lip_bl = blubber_lipids,
          //                                  Pro_m = muscle_protein,
          //                                  Pro_v = visceral_protein,
          //                                  Pro_bo = bones_protein,
          //                                  Pro_bl = blubber_protein,
          //                                  EL = ED_lipids,
          //                                  EP = ED_protein);

          // std::cout << ED_lipids << std::endl;
          // std::cout << ED_protein << std::endl;
          // std::cout << muscle_lipids[current_animal] << std::endl;
          // std::cout << muscle_protein[current_animal] << std::endl;
          // std::cout<< "-------------" << std::endl;
          
          // std::cout << ED_lipids << std::endl;
          // std::cout << ED_protein << std::endl;
          // 
          // std::cout<< "Muscle ==" << std::endl;
          // std::cout << mass_muscle << std::endl;
          // std::cout << muscle_lipids[current_animal] << std::endl;
          // std::cout << muscle_protein[current_animal] << std::endl;
          // std::cout<< "Viscera ==" << std::endl;
          // std::cout << mass_viscera << std::endl;
          // std::cout << visceral_lipids[current_animal] << std::endl;
          // std::cout << visceral_protein[current_animal] << std::endl;
          // std::cout<< "Bones ==" << std::endl;
          // std::cout << mass_bones << std::endl;
          // std::cout << bones_lipids[current_animal] << std::endl;
          // std::cout << bones_protein << std::endl;
          // std::cout<< "blubber ==" << std::endl;
          // std::cout << mass_blubber << std::endl;
          // std::cout << blubber_lipids[current_animal] << std::endl;
          // std::cout << blubber_protein[current_animal] << std::endl;
          // std::cout<< "---------- ==" << std::endl;
          
          // Energetic cost of placental maintenance during pregnancy -- (MJ) [P]
          placental_cost = placental_maintenance(fetal_growth_cost);
          
          // Heat increment of gestation -- (kJ converted to MJ) [Q]
          HIC = heat_gestation(birth_mass[current_animal], mass_increment_fetus)/1000;
          
          // Total cost of gestation -- (kJ)
          E_gest = fetal_growth_cost + placental_cost + HIC;
          
          } // End if current_day == 365
        } // End if cohortID == 4
        
        // ------------------------------------------------------------
        // GROWTH
        // ------------------------------------------------------------
        
        // Cost of growth (adults and juveniles) -- (MJ)
        length_nextday = age2length(animal->age + 1.0L/365.0L, animal->lengthatage);
        lean_mass_nextday = length2mass(length_nextday, animal->massatlength, true);
        mass_increment = (lean_mass_nextday - animal->lean_mass);
        
        E_growth = growth_cost(mass_increment, animal->bc, lean_water, blubber_lipids[current_animal],
                               ED_lipids, ED_protein, deposition_lipids, deposition_protein);
        
        // Cost of growth (calves) -- (MJ)
        
        if(cohortID == 5 && calves[current_animal].alive == 1){
          
          length_nextday_calves = age2length(calves[current_animal].age + 1.0L/365.0L, calves[current_animal].lengthatage);
          lean_mass_nextday_calves = length2mass(length_nextday_calves, calves[current_animal].massatlength, true);
          mass_increment_calves = (lean_mass_nextday_calves - calves[current_animal].lean_mass);
          
          E_growth_calves = growth_cost(mass_increment_calves, calves[current_animal].bc, lean_water, blubber_lipids[current_animal],
                                        ED_lipids, ED_protein, deposition_lipids, deposition_protein);
          
          // Cost of lactation estimated from the energy requirements of the calf
          E_lac = resting_metab_calves + LC_tot_calves + E_growth_calves;
        }
        
        // ------------------------------------------------------------
        // ENERGY EXPENDITURE
        // ------------------------------------------------------------
        
        // Adults and juveniles -- (MJ)
        E_out = (resting_metab + LC_tot + E_gest + E_lac + E_growth);
        if(animal->entangled(0) == 1) E_out = E_out + entanglement_cost;
        
        // Adults and juveniles -- (MJ)
        if(cohortID == 5){
          E_out_calves = (resting_metab_calves + LC_tot_calves + E_growth_calves);
          if(animal->entangled(0) == 1) E_out_calves = E_out_calves + entanglement_cost;
        }
        
        // ------------------------------------------------------------
        // ENERGY ALLOCATION
        // ------------------------------------------------------------
        
        if(growth && animal->alive == 1){
          
          E_balance = E_in < E_out;
          
          // Percent protein breakdown during catabolism
          // protein_breakdown = R::runif(0.40, 0.47);
          
          // Age (+ 1 day) -- Long double precision required here
          animal->age = animal->age + 1.0L/365.0L;
          
          // Body length
          animal->length = age2length(animal->age, animal->lengthatage);
          
          // Lean mass
          animal->lean_mass = lean_mass_nextday;
          
          // Lipid growth or breakdown
          delta_blubber = lipid_anacat[E_balance] * deposition_lipids * (E_in - E_out) / ED_lipids;
          
          animal->fat_mass = animal->fat_mass + delta_blubber;
          animal->tot_mass = animal->fat_mass + animal->lean_mass;
          animal->bc = animal->fat_mass / animal->tot_mass;
          
          // Same for calves
          if(cohortID == 5 && calves[current_animal].alive == 1){
            
            E_balance_calves = E_in_calves < E_out_calves;
            
            // Percent lipid breakdown during catabolism
            lipid_anacat_calves[1] = R::runif(0.53, 0.6);
            
            calves[current_animal].age = calves[current_animal].age + 1.0L/365.0L;
            calves[current_animal].length = age2length(calves[current_animal].age, calves[current_animal].lengthatage);
            calves[current_animal].lean_mass = lean_mass_nextday_calves;
            
            // Lipid growth or breakdown
            delta_blubber_calves = lipid_anacat_calves[E_balance_calves] * deposition_lipids * (E_in_calves - E_out_calves) / ED_lipids;
            
            calves[current_animal].fat_mass = calves[current_animal].fat_mass + delta_blubber_calves;
            calves[current_animal].tot_mass = calves[current_animal].fat_mass + calves[current_animal].lean_mass;
            calves[current_animal].bc = calves[current_animal].fat_mass / calves[current_animal].tot_mass;
            
          }
          
        } // End energy allocation
        
        // Convert back to relevant units
        travel_time = travel_time / 3600;
        feeding_time = feeding_time / 3600;
        resting_time = resting_time / 3600;
        nursing_time = nursing_time / 3600;
        
        // ------------------------------------------------------------
        // STARVATION
        // ------------------------------------------------------------
        
        if(animal->bc < starvation) animal->alive = 0;
        if(calves[current_animal].bc < starvation) calves[current_animal].alive = 0;
        
        } // End if current_day > 0
        
      // }
        
        // ------------------------------------------------------------
        // STORE VALUES in output matrices
        // ------------------------------------------------------------
        
      // // Restart if(alive) conditional otherwise some attributes may still be
      // // recorded despite animals having died of starvation

      } // End if animal->alive
       
      if(animal->alive){
        
        // +++ Animal traits +++
        
        *(attrib_out++) = animal->cohortID;                // Unique cohort identifier
        *(attrib_out++) = animal->alive;                   // Alive or dead
        *(attrib_out++) = animal->age;                     // Age
        *(attrib_out++) = animal->bc;                      // Body condition
        *(attrib_out++) = animal->length;                  // Total body length 
        *(attrib_out++) = animal->lengthatage(0,0);        // Coefficient of the Gompertz growth curve
        *(attrib_out++) = animal->lengthatage(0,1);        // Coefficient of the Gompertz growth curve
        *(attrib_out++) = animal->lengthatage(0,2);        // Coefficient of the Gompertz growth curve
        *(attrib_out++) = animal->tot_mass;                // Total body mass (including blubber)
        *(attrib_out++) = animal->lean_mass;               // Lean body mass (excluding blubber)
        *(attrib_out++) = animal->fat_mass;                // Blubber mass
        *(attrib_out++) = animal->massatlength(0,0);       // Intercept of the logarithmic mass-at-length relationship
        *(attrib_out++) = animal->massatlength(0,1);       // Slope of the logarithmic mass-at-length relationship
        *(attrib_out++) = animal->mouth_ratio;             // Ratio of mouth width to total body length
        *(attrib_out++) = animal->mouth_angle;             // Upper angle of the mouth
        *(attrib_out++) = animal->mouth_width;             // Width of the mouth
        *(attrib_out++) = animal->gsl;                     // Gulf of St Lawrence
        *(attrib_out++) = animal->seus;                    // Southeastern U.S.
        *(attrib_out++) = animal->abort;                   // Abortion (pregnant females)
        
        // +++ Stressors +++
        
        if(stressors){
          
        *(attrib_stressors_out++) = gear_risk;             // Probability of becoming entangled
        *(attrib_stressors_out++) = animal->entangled(0);  // Is the animal entangled?
        *(attrib_stressors_out++) = animal->entangled(1);  // Does the entanglement involve the anterior part of the body?
        *(attrib_stressors_out++) = animal->entangled(2);  // Entanglement severity
        *(attrib_stressors_out++) = animal->entangled(3);  // Entanglement duration
        *(attrib_stressors_out++) = animal->entangled(4);  // Day when entanglement event started
        *(attrib_stressors_out++) = animal->entangled(5);  // Day when entanglement event ended
        *(attrib_stressors_out++) = strike_risk;           // Probability of being struck by vessels
        *(attrib_stressors_out++) = animal->strike;        // Incidence of vessel strike
        *(attrib_stressors_out++) = animal->respnoise;     // Incidence of behavioral response to noise exposure
        *(attrib_stressors_out++) = dB;                    // Noise levels
        *(attrib_stressors_out++) = response_dB;           // Behavioral response threshold
        
        }
        
        // +++ Energy budget +++
        
        *(attrib_E_out++) = E_in - E_out;                  // Net energy (adults and juveniles)
        *(attrib_E_out++) = E_in;                          // Energy gain (adults and juveniles)
        *(attrib_E_out++) = E_out;                         // Energy expenditure (adults and juveniles)
        *(attrib_E_out++) = delta_blubber;                 // Change in blubber mass (adults and juveniles)
        *(attrib_E_out++) = deposition_lipids;             // Deposition efficiency of lipids (adults and juveniles)
        *(attrib_E_out++) = ED_lipids;                     // Energy density of lipids (adults and juveniles)
        *(attrib_E_out++) = lipid_anacat[0];               // Lipid anabolism efficiency (adults and juveniles)
        *(attrib_E_out++) = lipid_anacat[1];               // Lipid catabolism efficiency (adults and juveniles)
        
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
        *(attrib_costs_out++) = scalar_LC;               // Scalar on locomotory costs
        *(attrib_costs_out++) = stroke_rate;             // Stroke rate
        *(attrib_costs_out++) = E_growth;                // Energetic cost associated with growth
        *(attrib_costs_out++) = E_gest;                  // Energetic cost associated with gestation
        *(attrib_costs_out++) = fetal_growth_cost;       // Energetic cost associated with fetal growth
        *(attrib_costs_out++) = placental_cost;          // Energetic cost associated with placental maintenance
        *(attrib_costs_out++) = HIC;                     // Heat increment of gestation
        *(attrib_costs_out++) = E_lac;                   // Energetic cost associated with lactation
        *(attrib_costs_out++) = mass_increment;          // Daily change in mass of adult
        
        // +++ Activity budgets +++
        
        *(attrib_activity_out++) = travel_dist;          // Distance traveled
        *(attrib_activity_out++) = swim_speed;           // Swimming speed
        *(attrib_activity_out++) = travel_time;          // Time spent traveling
        *(attrib_activity_out++) = feeding_time;         // Time spent feeding
        *(attrib_activity_out++) = nursing_time;         // Time spent nursing
        *(attrib_activity_out++) = resting_time;         // Time spent resting
        *(attrib_activity_out++) = n_zeroes;             // Number of behaviors not displayed 
        *(attrib_activity_out++) = t_sum;                // Total time spent in all behaviors
        *(attrib_activity_out++) = t_remain;             // Remaining time
        
        // +++ Fetus development +++
        
        if(cohortID == 4){
          
          *(attrib_fetus_out++) = fetus_l;                               // Fetus length
          *(attrib_fetus_out++) = fetus_m;                               // Fetus mass
          *(attrib_fetus_out++) = mass_increment_fetus;                  // Daily change in mass of fetus
          *(attrib_fetus_out++) = length_increment_fetus;                // Daily change in length of fetus
          *(attrib_fetus_out++) = birth_length[current_animal];          // Expected length at birth
          *(attrib_fetus_out++) = birth_mass[current_animal];            // Expected mass at birth
          *(attrib_fetus_out++) = mass_muscle;                           // Fetal muscle mass
          *(attrib_fetus_out++) = mass_viscera;                          // Fetal visceral mass
          *(attrib_fetus_out++) = mass_bones;                            // Fetal bone mass
          *(attrib_fetus_out++) = mass_blubber;                          // Fetal blubber mass
          *(attrib_fetus_out++) = muscle_lipids[current_animal];         // Lipid content in fetal muscles
          *(attrib_fetus_out++) = muscle_protein[current_animal];        // Protein content in fetal muscles
          *(attrib_fetus_out++) = visceral_lipids[current_animal];       // Lipid content in fetal viscera
          *(attrib_fetus_out++) = visceral_protein[current_animal];      // Protein content in fetal viscera
          *(attrib_fetus_out++) = bones_lipids[current_animal];          // Lipid content in fetal bones
          *(attrib_fetus_out++) = blubber_lipids[current_animal];        // Lipid content in fetal blubber
          *(attrib_fetus_out++) = blubber_protein[current_animal];       // Protein content in fetal blubber
          
        }
        
      } else { 
        
        // If the animal is dead
        // Set all outputs other than location to 0
        // This is required otherwise attributes from the next animal are erroneously recorded
        
        *(attrib_out++) = animal->cohortID;  // Unique cohort identifier
        *(attrib_out++) = animal->alive;     // Alive or dead
        *(attrib_out++) = 0.0L;              // Age
        *(attrib_out++) = 0.0L;              // Body condition
        *(attrib_out++) = 0.0L;              // Total body length 
        *(attrib_out++) = 0.0L;              // Coefficient of the Gompertz growth curve
        *(attrib_out++) = 0.0L;              // Coefficient of the Gompertz growth curve
        *(attrib_out++) = 0.0L;              // Coefficient of the Gompertz growth curve
        *(attrib_out++) = 0.0L;              // Total body mass (including blubber)
        *(attrib_out++) = 0.0L;              // Lean body mass (excluding blubber)
        *(attrib_out++) = 0.0L;              // Blubber mass
        *(attrib_out++) = 0.0L;              // Intercept of the logarithmic mass-at-length relationship
        *(attrib_out++) = 0.0L;              // Slope of the logarithmic mass-at-length relationship
        *(attrib_out++) = 0.0L;              // Ratio of mouth width to total body length
        *(attrib_out++) = 0.0L;              // Upper angle of the mouth
        *(attrib_out++) = 0.0L;              // Width of the mouth
        *(attrib_out++) = 0.0L;              // Gulf of St Lawrence
        *(attrib_out++) = 0.0L;              // Southeastern U.S.
        *(attrib_out++) = animal->abort;              // Abortion (pregnant females)
        
        // +++ Stressors +++
        
        *(attrib_stressors_out++) = 0.0L;                     // Probability of becoming entangled
        *(attrib_stressors_out++) = animal->entangled(0);     // Is the animal entangled?
        *(attrib_stressors_out++) = animal->entangled(1);     // Does the entanglement involve the anterior part of the body?
        *(attrib_stressors_out++) = animal->entangled(2);     // Entanglement severity
        *(attrib_stressors_out++) = animal->entangled(3);     // Entanglement duration
        *(attrib_stressors_out++) = animal->entangled(4);     // Day when entanglement event started
        *(attrib_stressors_out++) = animal->entangled(5);     // Day when entanglement event ended
        *(attrib_stressors_out++) = 0.0L;                     // Probability of being struck by vessels
        *(attrib_stressors_out++) = animal->strike;           // Incidence of vessel strike
        *(attrib_stressors_out++) = animal->respnoise;        // Incidence of behavioral response to noise exposure
        *(attrib_stressors_out++) = 0.0L;                     // Noise levels
        *(attrib_stressors_out++) = 0.0L;                     // Behavioral response threshold
        
        // +++ Energy budget +++
        
        *(attrib_E_out++) = 0.0L;            // Net energy (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Energy gain (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Energy expenditure (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Change in blubber mass (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Deposition efficiency of lipids (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Energy density of lipids (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Lipid anabolism efficiency (adults and juveniles)
        *(attrib_E_out++) = 0.0L;            // Lipid catabolism efficiency (adults and juveniles)         
        
        // +++ Energy intake +++
        
        *(attrib_feeding_out++) = 0.0L;        // Incidence of foraging
        *(attrib_feeding_out++) = 0.0L;        // Prey concentration
        *(attrib_feeding_out++) = 0.0L;        // Minimum prey concentration needed for foraging
        *(attrib_feeding_out++) = 0.0L;        // Energy gain (adults and juveniles)
        *(attrib_feeding_out++) = 0.0L;        // Swimming speed while foraging
        *(attrib_feeding_out++) = 0.0L;        // Capture efficiency
        *(attrib_feeding_out++) = 0.0L;        // % reduction in the gape when entangled
        *(attrib_feeding_out++) = 0.0L;        // Number of daylight hours for given location/season
        *(attrib_feeding_out++) = 0.0L;        // Feeding/nursing effort
        *(attrib_feeding_out++) = 0.0L;        // Target body condition
        *(attrib_feeding_out++) = 0.0L;        // Copepod mass
        *(attrib_feeding_out++) = 0.0L;        // Energy density of Calanus finmarchicus stage V
        *(attrib_feeding_out++) = 0.0L;        // Digestive efficiency
        *(attrib_feeding_out++) = 0.0L;        // Metabolizing efficiency
        *(attrib_feeding_out++) = 0.0L;        // Prey energy content
        
        // +++ Energy costs +++
        
        *(attrib_costs_out++) = 0.0L;           // Resting metabolic rate
        *(attrib_costs_out++) = 0.0L;           // Locomotory costs
        *(attrib_costs_out++) = 0.0L;           // Scalar on locomotory costs
        *(attrib_costs_out++) = 0.0L;           // Stroke rate
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with growth
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with gestation
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with fetal growth
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with placental maintenance
        *(attrib_costs_out++) = 0.0L;           // Heat increment of gestation
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with lactation
        *(attrib_costs_out++) = 0.0L;           // Daily change in mass
        
        // +++ Activity budgets +++
        
        *(attrib_activity_out++) = 0.0L;         // Distance traveled
        *(attrib_activity_out++) = 0.0L;         // Swimming speed
        *(attrib_activity_out++) = 0.0L;         // Time spent traveling
        *(attrib_activity_out++) = 0.0L;         // Time spent feeding
        *(attrib_activity_out++) = 0.0L;         // Time spent nursing
        *(attrib_activity_out++) = 0.0L;         // Time spent resting
        *(attrib_activity_out++) = 0.0L;         // Number of behaviors not displayed
        *(attrib_activity_out++) = 0.0L;         // Total time engaged in any behavior
        *(attrib_activity_out++) = 0.0L;         // Remaining time
        
        // +++ Fetus development +++
        
        if(cohortID == 4){
          
          *(attrib_fetus_out++) = 0.0L;          // Fetus length
          *(attrib_fetus_out++) = 0.0L;          // Fetus mass
          *(attrib_fetus_out++) = 0.0L;          // Fetus daily mass gain
          *(attrib_fetus_out++) = 0.0L;          // Fetus daily length gain
          *(attrib_fetus_out++) = 0.0L;          // Expected length at birth
          *(attrib_fetus_out++) = 0.0L;          // Expected mass at birth
          *(attrib_fetus_out++) = 0.0L;          // Fetal muscle mass
          *(attrib_fetus_out++) = 0.0L;          // Fetal visceral mass
          *(attrib_fetus_out++) = 0.0L;          // Fetal bone mass
          *(attrib_fetus_out++) = 0.0L;          // Fetal blubber mass
          *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal muscles
          *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal muscles
          *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal viscera
          *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal viscera
          *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal bones
          *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal blubber
          *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal blubber
          
        }
        
      } // End if alive
        
        // +++ Calf traits +++
        
        if(cohortID == 5 && calves[current_animal].alive == 1){
          
          *(attrib_calves_out++) = calves[current_animal].cohortID;            // Unique cohort identifier
          *(attrib_calves_out++) = calves[current_animal].alive;               // Alive or dead
          *(attrib_calves_out++) = calves[current_animal].age;                 // Age
          *(attrib_calves_out++) = calves[current_animal].bc;                  // Body condition
          *(attrib_calves_out++) = calves[current_animal].length;              // Body length
          *(attrib_calves_out++) = calves[current_animal].lengthatage(0,0);    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = calves[current_animal].lengthatage(0,1);    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = calves[current_animal].lengthatage(0,2);    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = calves[current_animal].tot_mass;            // Total body mass (including blubber)
          *(attrib_calves_out++) = calves[current_animal].lean_mass;           // Lean body mass (excluding blubber)
          *(attrib_calves_out++) = calves[current_animal].fat_mass;            // Blubber mass
          *(attrib_calves_out++) = calves[current_animal].massatlength(0,0);   // Intercept of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = calves[current_animal].massatlength(0,1);   // Slope of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = calves[current_animal].mouth_ratio;         // Ratio of mouth width to total body length
          *(attrib_calves_out++) = calves[current_animal].mouth_angle;         // Upper angle of the mouth
          *(attrib_calves_out++) = calves[current_animal].mouth_width;         // Width of the mouth
          
          *(attrib_nursing_out++) = milk_assim;                     // Milk assimilation
          *(attrib_nursing_out++) = milk_provisioning;              // Milk provisioning
          *(attrib_nursing_out++) = mammary_M;                      // Mass of mammary glands
          *(attrib_nursing_out++) = milk_production_rate;           // Milk yield
          *(attrib_nursing_out++) = D_milk;                         // Density of milk
          *(attrib_nursing_out++) = t_suckling;                     // Time spent suckling during nursing activities
          *(attrib_nursing_out++) = target_bc[1];                   // Target body condition (mother)
          *(attrib_nursing_out++) = feeding_nursing_effort[1];      // Nursing effort
          *(attrib_nursing_out++) = milk_lipids;                    // Proportion of lipids in milk
          *(attrib_nursing_out++) = milk_protein;                   // Proportion of protein in milk
          *(attrib_nursing_out++) = ED_lipids;                      // Energy density of lipids
          *(attrib_nursing_out++) = ED_protein;                     // Energy density of protein
          
          *(attrib_costs_calves_out++) = resting_metab_calves;           // Resting metabolic rate
          *(attrib_costs_calves_out++) = LC_tot_calves;                  // Locomotory costs
          *(attrib_costs_calves_out++) = scalar_LC;                      // Scalar on locomotory costs
          *(attrib_costs_calves_out++) = E_growth_calves;                // Energetic cost associated with growth
          *(attrib_costs_calves_out++) = mass_increment_calves;          // Daily change in mass
          
          *(attrib_E_calves_out++) = E_in_calves - E_out_calves;    // Net energy (calves)
          *(attrib_E_calves_out++) = E_in_calves;                   // Energy gain (calves)
          *(attrib_E_calves_out++) = E_out_calves;                  // Energy expenditure (calves)
          *(attrib_E_calves_out++) = delta_blubber_calves;          // Change in blubber mass (calves)
          
        } else {
          
          *(attrib_calves_out++) = calves[current_animal].cohortID ; // Unique cohort identifier
          *(attrib_calves_out++) = 0.0L;    // Alive or dead
          *(attrib_calves_out++) = 0.0L;    // Age
          *(attrib_calves_out++) = 0.0L;    // Body condition
          *(attrib_calves_out++) = 0.0L;    // Body length
          *(attrib_calves_out++) = 0.0L;    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;    // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;    // Total body mass (including blubber)
          *(attrib_calves_out++) = 0.0L;    // Lean body mass (excluding blubber)
          *(attrib_calves_out++) = 0.0L;    // Blubber mass
          *(attrib_calves_out++) = 0.0L;    // Intercept of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = 0.0L;    // Slope of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = 0.0L;    // Ratio of mouth width to total body length
          *(attrib_calves_out++) = 0.0L;    // Upper angle of the mouth
          *(attrib_calves_out++) = 0.0L;    // Width of the mouth
          
          *(attrib_nursing_out++) = 0.0L;    // Milk assimilation
          *(attrib_nursing_out++) = 0.0L;    // Milk provisioning
          *(attrib_nursing_out++) = 0.0L;    // Mass of mammary glands
          *(attrib_nursing_out++) = 0.0L;    // Milk yield
          *(attrib_nursing_out++) = 0.0L;    // Density of milk
          *(attrib_nursing_out++) = 0.0L;    // Time spent suckling during nursing activities
          *(attrib_nursing_out++) = 0.0L;    // Target body condition (mother)
          *(attrib_nursing_out++) = 0.0L;    // Nursing effort
          *(attrib_nursing_out++) = 0.0L;    // Proportion of lipids in milk
          *(attrib_nursing_out++) = 0.0L;    // Proportion of protein in milk
          *(attrib_nursing_out++) = 0.0L;    // Energy density of lipids
          *(attrib_nursing_out++) = 0.0L;    // Energy density of protein
          
          *(attrib_costs_calves_out++) = 0.0L;   // Resting metabolic rate
          *(attrib_costs_calves_out++) = 0.0L;   // Locomotory costs
          *(attrib_costs_calves_out++) = 0.0L;   // Scalar on locomotory costs
          *(attrib_costs_calves_out++) = 0.0L;   // Energetic cost associated with growth
          *(attrib_costs_calves_out++) = 0.0L;   // Daily change in mass
          
          *(attrib_E_calves_out++) = 0.0L;   // Net energy (calves)
          *(attrib_E_calves_out++) = 0.0L;   // Energy gain (calves)
          *(attrib_E_calves_out++) = 0.0L;   // Energy expenditure (calves)
          *(attrib_E_calves_out++) = 0.0L;   // Change in blubber mass (calves)
          
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
  results["E_calves"] = attrib_E_calves;
  
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
 //' @name NARW_simulator
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
 Rcpp::List NARW_simulator(
     int cohortID,
     Rcpp::NumericVector seus,
     Rcpp::NumericVector gsl,
     Eigen::MatrixXd support,
     std::vector<Eigen::MatrixXd> densities,
     std::vector<Eigen::MatrixXd> densities_seus,
     std::vector<Eigen::MatrixXd> densities_gsl,
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
     bool stressors,
     bool growth,
     Rcpp::NumericVector cumgest,
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
     latent_environments.emplace_back(densities[*d], densities_seus[*d], densities_gsl[*d], 
                                      prey[*d], fishing[*d], vessels[*d], noise[*d], 
                                      daylight[*d], regions, limits, resolution, *d);
   
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
       latent_animals.emplace_back(Animal(*(xiter++), *(yiter++), cohortID, seus[j], gsl[j]));
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
                  limits, resolution, stepsize, stressors, growth, cumgest, progress);
   
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
                       Eigen::MatrixXd density_seus, 
                       Eigen::MatrixXd density_gsl, 
                       Eigen::MatrixXd prey, 
                       Eigen::MatrixXd fishing, 
                       Eigen::MatrixXd vessels, 
                       Eigen::MatrixXd noise, 
                       Eigen::MatrixXd daylight, 
                       Eigen::MatrixXd regions, 
                       Eigen::VectorXd limits, 
                       Eigen::VectorXd resolution,
                       double x, double y, char layer) {
  Environment e(density, density_seus, density_gsl, prey, fishing, vessels, noise, daylight, regions, limits, resolution, 0);
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