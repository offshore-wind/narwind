// R and 'Eigen' integration using 'Rcpp'. 
// 'Eigen' is a C++ template library 
// for linear algebra: matrices, vectors, 
// numerical solvers and related algorithms.

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include "progressbar.hpp"
#include <Eigen/Eigen>
#include <RcppEigen.h>
#include <Rcpp.h> // For lists
#include <iostream>
#include <cstdlib>
#include "environments.h"
#include "movements.h"
#include "geodesic.h"
#include "bioenergetics.h"
#include <random>
using namespace std;


/** Basic state for an animal **/

// C++ Note: -----------------------------------------------------------------------------------
// Structures (also called structs) are a way to group several related variables into one place. 
// Each variable in the structure is known as a member of the structure.
// Unlike an array, a structure can contain many different data types (int, string, bool, etc.).
// The public keyword is an access specifier. 
// Access specifiers define how the members (attributes and methods) 
// of a class can be accessed. In the example above, the members are public 
// - which means that they can be accessed and modified from outside the code.
// ---------------------------------------------------------------------------------------------
// Using const to declare data members constant, i.e., must have some value before the object is constructed

//-----------------------------------------------------------------------------------
// PARAMETERS 
// -----------------------------------------------------------------------------------
// x: Easting
// y: Northing
// x0: Initial location (easting)
// y0: Initial location (northing)
// cohortID: Unique identifier for the cohort to which an animal belongs.
// seus: Binary variable indicating whether the animal is set to migrate to the Southeast U.S. calving grounds
// gsl: Binary variable indicating whether the animal is set to migrate to the Gulf of St Lawrence feeding grounds
// alive: Binary variable indicating whether the animal is alive (1) or mort (0).
// age: Age of the animal in decimal years
// bc: Body condition, expressed as a relative % of fat
// lengthatage: Coefficients of the Gompertz length-at-age relationship
// length: Total body length in cm
// massatlength: Coefficients of the logarithmic mass-at-length relationship
// lean_mass: Total mass of lean tissues (e.g., muscles, bones) in kg
// tot_mass: Total body mass in kg
// fat_mass: Total blubber mass in kg
// mouth_ratio: Ratio of mouth width to body length (d.u.)
// mouth_angle: Angle between the tip of the baleen plates and the outer edges of the baleen rack (rad)
// mouth_width: Width of the mouth, i.e. body width at 10% of the body length from the snout (m)
// entangled: Entanglement status - vector length 7. See entanglement_event for details.
// strike: Binary variable indicating whether the animal has been struck by a vessel (1) or not (0)
// respnoise: Binary variable indicating whether the animal is displaying a behavioral response to pile-driving noise (1) or not (0)
// abort: Binary variable indicating whether a pregnancy has been aborted (1) or not (0)
// starve: Binary variable indicating whether the animal has died of starvation (1) or not (0)
// mort: Binary variable indicating whether the animal has died of old age
// region: Region in which the animal is located
// calvgrounds: Binary variable indicating whether the animal is in the SEUS calving grounds

struct Animal {
  double x, y, x0, y0;
  int cohortID;
  int seus;
  int gsl;
  int north;
  int alive;
  long double age;
  int sex;
  long double bc;
  Eigen::MatrixXd lengthatage;
  long double length;
  Eigen::MatrixXd massatlength;
  long double lean_mass;
  long double tot_mass;
  long double fat_mass;
  double mouth_ratio;
  double mouth_angle;
  double mouth_width;
  double gape;
  Rcpp::NumericVector entangled;
  int strike;
  int respnoise;
  int abort;
  int starve;
  int mort;
  int region;
  int calvgrounds;
  
  // Member initializer list
  Animal() : 
    x(0), 
    y(0), 
    x0(0),
    y0(0),
    cohortID(0), 
    seus(0),
    gsl(0), 
    north(0),
    alive(cohortID >=1), 
    age(0),
    sex(is_female(cohortID)),
    bc(start_bcondition(cohortID)),
    lengthatage(agL(age)), 
    length(age2length(age, lengthatage)),
    massatlength(mL()), 
    lean_mass(length2mass(length, massatlength)),
    tot_mass(lean_mass/(1-bc)),
    fat_mass(bc*tot_mass),
    mouth_ratio(start_mouth(cohortID, age)), 
    mouth_angle(76.7), 
    mouth_width(length*mouth_ratio), 
    gape(gape_size(length, mouth_width, mouth_angle)),
    entangled(Rcpp::NumericVector::create(0,0,0,0,0,0,0,0)),
    strike(0), 
    respnoise(0),
    abort(0),
    starve(0),
    mort(0),
    region(0),
    calvgrounds(0){ }
  
  Animal(const double & x0, 
         const double & y0,
         int cohort,
         int to_seus,
         int to_gsl) : 
    x(x0), 
    y(y0), 
    cohortID(cohort),
    seus(to_seus),
    gsl(to_gsl),
    north(0),
    alive(cohortID >=1),
    age(start_age(cohort)),
    sex(is_female(cohort)),
    bc(start_bcondition(cohort)),
    lengthatage(agL(age)), 
    length(age2length(age, lengthatage)),
    massatlength(mL()),
    lean_mass(length2mass(length, massatlength)),
    tot_mass(lean_mass/(1-bc)),
    fat_mass(bc*tot_mass),
    mouth_ratio(start_mouth(cohort, age)), 
    mouth_angle(76.7), 
    mouth_width(length*mouth_ratio),
    gape(gape_size(length, mouth_width, mouth_angle)),
    entangled(Rcpp::NumericVector::create(0,0,0,0,0,0,0,0)),
    strike(0), 
    respnoise(0),
    abort(0),
    starve(0),
    mort(0),
    region(0),
    calvgrounds(0){ }
  
};

/**
 //' @title Animal with capacity for attraction to other latent animals
 //'
 //' @description Initialize an animal with latent members, and couple the observed position to one of the latent animals
 //' 
 //' @param init_latent Latent animals to be attracted to
 //' @param init_animal Index of latent animal to use as starting location
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
 * @tparam CoordSize Number of coordinates
 * @param nparam_animal Number of parameters needed to capture individual attributes
 * @param nparam_calves Number of parameters needed to capture calf attributes
 * @param nparam_E Number of parameters needed to capture dynamics of energy gain / expenditure
 * @param nparam_E_calves Number of parameters to capture dynamics of energy gain / expenditure in calves
 * @param nparam_fetus Number of parameters needed to capture fetus growth
 * @param nparam_stressors Number of parameters needed to capture stressor values
 * @param nparam_activity Number of parameters needed to capture activity budgets
 * @param nparam_feeding Number of parameters needed to capture feeding behavior
 * @param nparam_nursing Number of parameters needed to capture nursing behavior
 * @param nparam_costs Number of parameters needed to capture energetic costs
 * @param nparam_costs_calves Number of parameters needed to capture energetic costs in calves
 * @return array of simulated data including time-stamped locations
 */

template<typename AnimalType, 
         typename EnvironmentContainer,
         typename MovementRule, 
         int CoordSize = 7,
         int nparam_animal = 26,
         int nparam_calves = 26,
         int nparam_E = 8,
         int nparam_E_calves = 4,
         int nparam_fetus = 25,
         int nparam_stressors = 24,
         int nparam_activity = 8,
         int nparam_feeding = 18,
         int nparam_nursing = 13,
         int nparam_costs = 15,
         int nparam_costs_calves = 5>

Rcpp::List movesim(
    int cohortID,
    std::vector<AnimalType> & animals,
    std::vector<AnimalType> & calves,
    EnvironmentContainer & environments,
    std::vector<Environment> & layers,
    Rcpp::NumericMatrix doseresp_seed,
    Rcpp::NumericVector doseresp,
    std::vector<std::size_t> densitySeq,
    MovementRule & m,
    Eigen::MatrixXd support,
    Eigen::VectorXd limits, 
    Eigen::VectorXd limits_regions,
    Eigen::VectorXd limits_prey,
    Eigen::VectorXd limits_fishing,
    Eigen::VectorXd limits_vessels,
    Eigen::VectorXd limits_noise,
    Eigen::VectorXd resolution,
    Eigen::VectorXd resolution_regions,
    Eigen::VectorXd resolution_prey,
    Eigen::VectorXd resolution_fishing,
    Eigen::VectorXd resolution_vessels,
    Eigen::VectorXd resolution_noise,
    bool stressors,
    bool growth,
    double prey_scale,
    Rcpp::NumericVector starvation_coefs,
    double starvation_death,
    double starvation_onset,
    Rcpp::NumericVector nursing_cessation,
    double piling_hrs,
    bool progress) {
  
  // Number of animals and time points to simulate
  std::size_t n = animals.size();  // nsim
  std::size_t t = environments.size(); // 458 days (365 days + 3 month burn-in + initial conditions)
  
  // Initialize list in which results will be stored
  Rcpp::List results;
  
  // Initialize and label array in which location information will be stored
  Rcpp::NumericVector out(n*t*CoordSize);
  out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
  std::vector<string> attrib_coords = {"easting", "northing", "region", "inseus", "resid_m", "resid_sd", "pleave"};
  
  // Initialize and label array in which animal attribute data will be stored
  Rcpp::NumericVector attrib(n*t*nparam_animal);
  attrib.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  std::vector<string> attrib_names = {
    "cohort",
    "alive",
    "age", 
    "bc", 
    "length",
    "length_a",
    "length_b",
    "length_c",
    "mass",
    "leanmass",
    "fatmass",
    "mass_a",
    "mass_b",
    "gape",
    "mouth_r",
    "mouth_a",
    "mouth_w",
    "gsl",
    "seus",
    "north",
    "abort",
    "p_starve",
    "starve",
    "died",
    "date_died",
    "p_surv"};
  
  // Create a separate array for dependent calves
  Rcpp::NumericVector attrib_calves(n*t*(nparam_calves));
  attrib_calves.attr("dim") = Rcpp::Dimension(nparam_calves,n,t);
  std::vector<string> attrib_calves_names = {
    "cohort_calf",
    "born",
    "pbirth",
    "dob",
    "date_died_calf",
    "alive_calf",
    "age_calf",
    "bc_calf",
    "length_calf", 
    "La_calf",
    "Lb_calf",
    "Lc_calf",
    "mass_calf",
    "leanmass_calf",
    "fatmass_calf",
    "ma_calf",
    "mb_calf",
    "gape_calf",
    "mouth_r_calf",
    "mouth_a_calf",
    "mouth_w_calf",
    "p_starve_calf",
    "starve_calf",
    "died_calf",
    "date_died_calf",
    "p_surv_calf"};
  
  // Initialize and label array in which energy budget data will be stored
  Rcpp::NumericVector attrib_E(n*t*nparam_E);
  attrib_E.attr("dim") = Rcpp::Dimension(nparam_E,n,t);
  std::vector<string> attrib_E_names = {
    "E_tot",
    "E_in",
    "E_out",
    "delta_fat",
    "lip_anab",
    "lip_catab",
    "EDlip",
    "EDpro"};

  // Initialize and label array in which energy budget data will be stored for calves
  Rcpp::NumericVector attrib_E_calves(n*t*nparam_E_calves);
  attrib_E_calves.attr("dim") = Rcpp::Dimension(nparam_E_calves,n,t);
  std::vector<string> attrib_E_calves_names = {
    "E_tot_calf", 
    "E_in_calf",
    "E_out_calf",
    "delta_fat_calf"
  };

  // Initialize and label array in which fetal growth data will be stored
  Rcpp::NumericVector attrib_fetus(n*t*nparam_fetus);
  attrib_fetus.attr("dim") = Rcpp::Dimension(nparam_fetus,n,t);
  std::vector<string> attrib_fetus_names = {
    "fetus_l",
    "fetus_m",
    "delta_fetus_m",
    "delta_fetus_l",
    "birth_l", 
    "birth_m",
    "muscle_m", 
    "viscera_m", 
    "bones_m",
    "blubber_m", 
    "muscle_lip",
    "muscle_pro",
    "visc_lip", 
    "visc_pro",
    "bone_lip",
    "bone_pro",
    "blubber_lip",
    "blubber_pro",
    "prop_mu",
    "prop_visc",
    "prop_bones",
    "dens_blu",
    "dens_mu",
    "dens_visc",
    "dens_bo"};
  
  // Initialize and label array in which stressor data will be stored
  Rcpp::NumericVector attrib_stressors(n*t*nparam_stressors);
  attrib_stressors.attr("dim") = Rcpp::Dimension(nparam_stressors,n,t);
  std::vector<string> attrib_stressors_names = {
    "gear_risk",
    "entgl_cost",
    "is_entgl",
    "entgl_sev", 
    "entgl_head",
    "gape_perc",
    "entgl_d", 
    "entgl_start",
    "entgl_end", 
    "is_entgl_calf",
    "entgl_head_calf",
    "entgl_sev_calf", 
    "entgl_d_calf", 
    "entgl_start_calf", 
    "entgl_end_calf",
    "strike_risk", 
    "strike",
    "strike_calf", 
    "noise_resp", 
    "noise_lvl",
    "dB_thresh",
    "piling_hrs",
    "nursing_stop",
    "nursing_stop_hrs"};
  
  // Initialize and label array in which activity budget data will be stored
  Rcpp::NumericVector attrib_activity(n*t*nparam_activity);
  attrib_activity.attr("dim") = Rcpp::Dimension(nparam_activity,n,t);
  std::vector<string> attrib_activity_names = {
    "d_travel", 
    "swimspeed",
    "glide",
    "glide_feed",
    "glide_echelon",
    "t_travel",
    "t_feed",
    "t_rest_nurse"};
  
  // Create a separate array for data associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_feeding(n*t*nparam_feeding);
  attrib_feeding.attr("dim") = Rcpp::Dimension(nparam_feeding,n,t);
  std::vector<string> attrib_feeding_names = {
    "feed",
    "preyconc", 
    "minprey",
    "feedspeed",
    "captEff",
    "impedance",
    "feed_effort",
    "eta_lwrBC",
    "eta_upprBC",
    "targetBC", 
    "cop_mass",
    "water_content",
    "cop_kJ",
    "digestEff",
    "metabEff_juv",
    "metabEff_ad",
    "E_cop",
    "calf_feed"};
  
  // Create a separate array for data associated with energy intake in calves
  Rcpp::NumericVector attrib_nursing(n*t*nparam_nursing);
  attrib_nursing.attr("dim") = Rcpp::Dimension(nparam_nursing,n,t);
  std::vector<string> attrib_nursing_names = {
    "t_lac",
    "assim",
    "provision",
    "zeta",
    "milk_drop",
    "eta_milk",
    "mamm_M",
    "mammEff",
    "milk_rate",
    "targetBC_calf",
    "nursing",
    "milk_lip",
    "milk_pro"};

  // Create a separate array for data associated with energy costs in adults/juveniles
  Rcpp::NumericVector attrib_costs(n*t*nparam_costs);
  attrib_costs.attr("dim") = Rcpp::Dimension(nparam_costs,n,t);
  std::vector<string> attrib_costs_names = {
    "rmr",
    "LC",
    "scalar_LC",
    "stroke",
    "stroke_feed",
    "E_growth",
    "E_gest",
    "fgrowth",
    "placenta",
    "hic",
    "E_lac",
    "delta_m",
    "perc_muscle",
    "perc_viscera",
    "perc_bones"};
  
  // Create a separate array for data associated with energy costs in calves
  Rcpp::NumericVector attrib_costs_calves(n*t*nparam_costs_calves);
  attrib_costs_calves.attr("dim") = Rcpp::Dimension(nparam_costs_calves,n,t);
  std::vector<string> attrib_costs_calves_names = {
    "rmr_calf", 
    "LC_calf",
    "scalar_LC",
    "E_growth_calf",
    "delta_m_calf"};
  
  out.attr("dimnames") = Rcpp::List::create(attrib_coords, R_NilValue, R_NilValue);
  attrib.attr("dimnames") = Rcpp::List::create(attrib_names, R_NilValue, R_NilValue);
  attrib_calves.attr("dimnames") = Rcpp::List::create(attrib_calves_names, R_NilValue, R_NilValue);
  attrib_E.attr("dimnames") = Rcpp::List::create(attrib_E_names, R_NilValue, R_NilValue);
  attrib_E_calves.attr("dimnames") = Rcpp::List::create(attrib_E_calves_names,R_NilValue, R_NilValue);
  attrib_fetus.attr("dimnames") = Rcpp::List::create(attrib_fetus_names, R_NilValue, R_NilValue);
  attrib_stressors.attr("dimnames") = Rcpp::List::create(attrib_stressors_names, R_NilValue, R_NilValue);
  attrib_activity.attr("dimnames") = Rcpp::List::create(attrib_activity_names, R_NilValue, R_NilValue);
  attrib_feeding.attr("dimnames") = Rcpp::List::create(attrib_feeding_names, R_NilValue, R_NilValue);
  attrib_nursing.attr("dimnames") = Rcpp::List::create(attrib_nursing_names, R_NilValue, R_NilValue);
  attrib_costs.attr("dimnames") = Rcpp::List::create(attrib_costs_names, R_NilValue, R_NilValue);
  attrib_costs_calves.attr("dimnames") = Rcpp::List::create(attrib_costs_calves_names, R_NilValue, R_NilValue);
  
  // Create iterators
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
  
  // Initialize progress bar
  progressbar bar(t);
  
  // ------------------------------------------------------------
  // DECLARE VARIABLES
  // ------------------------------------------------------------
  
  int current_animal, current_month, current_day;
  int enter_SEUS = 0;
  int depart_SEUS = 0;
  double current_x, current_y;
  Rcpp::NumericVector feeding_nursing_effort (2); // feeding, nursing
  Rcpp::NumericVector target_bc (2); // adjuv, calves
  
  // Stressors +++
  
  double gear_risk = 0.0L;
  double strike_risk = 0.0L;
  double dB = 0.0L;
  double p_surv = 1.0L;
  double p_surv_calf = 1.0L;
  
  // Activity budgets +++
  
  double tmax = 0.0L;
  double feeding_time = 0.0L;
  double travel_time = 0.0L;
  double resting_nursing_time = 0.0L;
  double travel_dist = 0.0L;
  // double daylight_hours;
  
  // Foraging +++
  
  double D_cop = 0.0L;
  bool is_feeding = 0.0L;
  double swim_speed = 0.0L;
  double E_calanus = 0.0L;
  double E_in = 0.0L;
  double p_starve = 0.0L;
  double p_starve_calf = 0.0L;
  int has_starved = 0;
  int has_starved_calf = 0;
  double calf_foraging = 0.0L;
  double min_feed = 0.0L;
  
  // Nursing +++
  
  double prob_birth = 0.0L;
  double milk_assim = 0.0L;
  double milk_provisioning = 0.0L;
  double mammary_M = 0.0L;
  double milk_production_rate = 0.0L;
  double E_in_calves = 0.0L;
  // double t_suckling = 0;
  
  // Energetic costs +++
  
  double resting_metab = 0.0L;
  double resting_metab_calves = 0.0L;
  double stroke_rate_foraging  = 0.0L;
  double stroke_rate = 0.0L;
  double percent_glide_foraging = 0.0L;
  double percent_glide = 0.0L;
  double scalar_LC = 0.0L;
  double scalar_LC_calf = 0.0L;
  double LC_tot = 0.0L;
  double LC_tot_calves = 0.0L;
  double percent_glide_echelon = 0.0L;
  
  // Calf birth & death  +++
  
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
  
  if(cohortID == 4 | cohortID == 5){
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {
      int i = animal - animals.begin();
      birth_length[i] = fetal_length(0, animal->length);
      birth_mass[i] = fetal_mass(0, animal->length);
    }
  }
  
  Rcpp::NumericVector born(n);
  Rcpp::NumericVector dob(n);
  Rcpp::NumericVector date_death(n);
  Rcpp::NumericVector date_death_calf(n);
  double p_leave = 0;
  
  // Pregnancy, fetal growth and lactation +++
  
  double fetus_l = 0.0L;
  double fetus_m = 0.0L;
  double bc_gest = 0.0L;
  double abortion_scalar = 1.0L;
  
  double mass_muscle = 0.0L;
  double mass_muscle_next = 0.0L;
  double mass_viscera = 0.0L;
  double mass_viscera_next = 0.0L;
  double mass_bones = 0.0L;
  double mass_bones_next = 0.0L;
  double mass_blubber = 0.0L;
  double mass_blubber_next = 0.0L;
  double FG0 = 0.0L;
  double FG1 = 0.0L;
  double fetal_growth_cost = 0.0L;
  double placental_cost = 0.0L;
  double mass_increment_fetus = 0.0L;
  double length_increment_fetus = 0.0L;
  double fetus_m_next = 0.0L;
  double fetus_l_next = 0.0L;
  double HIC = 0.0L;
  double E_gest = 0.0L;
  double E_lac = 0.0L;
  
  // Growth +++
  
  double length_nextday = 0.0L;
  double length_nextday_calves = 0.0L;
  double lean_mass_nextday = 0.0L;
  double lean_mass_nextday_calves = 0.0L;
  double lean_mass_increment = 0.0L;
  double lean_mass_increment_calves = 0.0L;
  double E_growth = 0.0L;
  double E_growth_calves = 0.0L;
  double E_out = 0.0L;
  double E_out_calves = 0.0L;
  double E_balance = 0.0L;
  double E_balance_calves = 0.0L;
  double delta_blubber = 0.0L;
  double delta_blubber_calves = 0.0L;
  
  // Behavioral response to pile-driving noise +++
  
  double response_dB = 0.0L;
  
  // ------------------------------------------------------------
  // CONSTANTS
  // ------------------------------------------------------------
  
  //  ****************** Mortality / Stressors ******************
  
  // -------------------------------
  // IMMEDIATE EFFECT OF MINOR ENTANGLEMENTS ON SURVIVAL
  // -------------------------------
  // Definition: Immediate drop survival probability following a minor entanglement event.
  // Notes: Survival is calculated over 3-month periods on the cloglog scale.
  // Source: Pirotta et al. (2023) -- Oikos
  const double entgl_minor = 0.0L;
  
  // -------------------------------
  // IMMEDIATE EFFECT OF MODERATE ENTANGLEMENTS ON SURVIVAL
  // -------------------------------
  // Definition: Immediate drop survival probability following a moderate entanglement event.
  // Notes: Survival is calculated over 3-month periods on the cloglog scale.
  // Source: Pirotta et al. (2023) -- Oikos
  const double entgl_moderate = 0.0L;
  
  // -------------------------------
  // IMMEDIATE EFFECT OF SEVERE ENTANGLEMENTS ON SURVIVAL
  // -------------------------------
  // Definition: Immediate drop survival probability following a severe entanglement event.
  // Notes: Survival is calculated over 3-month periods on the cloglog scale.
  // Source: Pirotta et al. (2023) -- Oikos
  const double entgl_severe = -0.46;
  
  Rcpp::NumericVector entgl_pcoms(3);
  entgl_pcoms(0) = entgl_minor;
  entgl_pcoms(1) = entgl_moderate;
  entgl_pcoms(2) = entgl_severe;
  
  // -------------------------------
  // ENERGETIC COST OF ENTANGLEMENT
  // -------------------------------
  // Definition: Additional energy expenditure associated with carrying fishing gear
  // Notes: Average of values reported in Figure 2a of Van Der Hoop et al. (2017)
  // Units: MJ
  // Source(s): Van Der Hoop et al. (2017)
  const double entanglement_cost = 248.8222913;
  
  //  ****************** Foraging ******************
  
  // -------------------------------
  // ENERGY DENSITY OF PREY
  // -------------------------------
  // Definition: Caloric content of calanoid copepods
  // Notes: Includes entries for Calanus finmarchicus, Calanus finmarchicus (C5),
  // Calanus finmarchicus (C5) / C. hyperboreus (C4), & Calanus finmarchicus / C. hyperboreus
  // Value: meta function applied to spreadsheet of model parameters
  // Units: kJ / g converted to MJ/g
  // Source(s): Comita et al. (1966),	Laurence et al. (1976),
  // Båmstedt et al. (1986),	Michaud et al. (2007),
  // Davies et al. (2012),	McKinstry et al. (2013)
  const double cop_kJ = 23.91922/1000;
  
  // -------------------------------
  // COPEPOD BIOMASS
  // -------------------------------
  // Definition: Mass of prey (wet weight)
  // Value: meta function applied to spreadsheet of model parameters
  // Units: mg converted to g / copepod
  // Source(s): Wishner et al. (1988), Michaud et al. (2007),  McKinstry et al. (2013)
  const double cop_mass = 1.670151/1000;
  
  // -------------------------------
  // PREY WATER CONTENT
  // -------------------------------
  // Definition: Quantity of water contained in prey (% wet weight)
  // Value: meta function applied to spreadsheet of model parameters
  // Units: % (d.u.)
  // Source(s): Ideka et al. (1989, 2006) & Gabrielsen et al. (1991)
  const double water_content = 0.7368484;
  
  // -------------------------------
  // PREY THRESHOLD
  // -------------------------------
  // Definition: Minimum prey concentration below which foraging stops
  // Notes: Excludes Kenney (1986)
  // Value: meta function applied to spreadsheet of model parameters
  // Units: copepods / m3, converted to g / m3
  // Source(s): Murison et al. (1989),	Mayo et al. (1990),	Mayo et al. (1992),
  // Wishner et al. (1995),	Beardsley et al. (1996),	Goodyear et al. (1996),
  // Woodley et al. (1996),	Baumgartner et al. (2003),	DeLorenzo Costa et al. (2006),
  // Pendleton et al. (2009),	Patrician et al. (2010),	McKinstry et al. (2013)
  const double min_calanus = 1066.208 * cop_mass * (1-water_content);

  // -------------------------------
  // CAPTURE EFFICIENCY
  // -------------------------------
  // Definition: Fraction of available prey retained on baleen filaments
  // Notes: Based on controlled lab experiments (flow tank)
  // Value: meta function applied to spreadsheet of model parameters
  // Units: % (d.u.)
  // Source(s): Mayo et al. (2001) & Werth et al. (2012)
  const double capture_efficiency = 0.91575;
  
  // -------------------------------
  // DIGESTIVE EFFICIENCY
  // -------------------------------
  // Definition: Proportion of ingested food absorbed through the gut wall and taken up into the animal's tissues 
  // Notes: Also referred to as assimilation efficiency or apparent digestibility.
  // Value: meta function applied to spreadsheet of model parameters
  // Units: % (d.u.)
  // Source(s): Swaim et al. (2009) & Fortune et al. (2013)
  const double digestive_efficiency = 0.94;
  
  // -------------------------------
  // URINARY EFFICIENCY
  // -------------------------------
  // Definition: Proportion of digestible energy lost as urea and other metabolic end products in the urine.
  // Notes: Assumed for mixed diet of all life stages of C. finmarchicus, Centropages hamatus, Centropages typicus and Pseudocalanus spp.
  // Units: % (d.u.)
  // Source(s): Fortune et al. (2013)
  const double urinary_efficiency = 0.92;
  
  // -------------------------------
  // HEAT INCREMENT OF FEEDING
  // -------------------------------
  // Definition: Proportion of metabolizable energy lost from an increase in heat production associated with food digestion
  // Notes: Value originally borrowed from Steller sea lions
  // Value: meta function applied to spreadsheet of model parameters
  // Units: % (d.u.)
  // Source(s): Fortune et al. (2013)
  Rcpp::NumericVector metabolizing_efficiency(2);
  metabolizing_efficiency(0) = 0.7403977; // Calves and juveniles
  metabolizing_efficiency(1) = 0.875; // Adults
  
  // Percent lipid breakdown during catabolism
  // const double lip_breakdown = 0.57;
  
  // -------------------------------
  // CATABOLISM EFFICIENCY
  // -------------------------------
  // Definition: Efficiency of lipid breakdown during catabolism in a period of energy deficit
  // Units: % (d.u.)
  // Source(s): Barboza et al. (2008)
  Rcpp::NumericVector catabolism_efficiency(2);
  catabolism_efficiency(0) = 1; // During anabolism - see equation
  catabolism_efficiency(1) = 0.80; // Efficiency
  
  // -------------------------------
  // ANABOLISM EFFICIENCY
  // -------------------------------
  // Definition: Efficiency of lipid deposition during anabolism in a period of energy surplus
  // Units: % (d.u.)
  // Source(s): Barboza et al. (2008)
  Rcpp::NumericVector anabolism_efficiency(2);
  anabolism_efficiency(0) = 0.80; // Efficiency
  anabolism_efficiency(1) = 1; // During catabolism -- see equation
  
  // -------------------------------
  // STEEPNESS OF THE FEEDING/NURSING EFFORT CURVE
  // -------------------------------
  // Units: d.u.
  Rcpp::NumericVector eta(2);
  eta(0) = 10; // For lower target BC for all animals
  eta(1) = 30; // For upper target BC for pregnant females and calves
  
  //  ****************** Lactation / Pregnancy ******************
  
  // -------------------------------
  // START OF SIMULATION
  // -------------------------------
  // Definition: Day when the one-year simulation starts, after a 3-month long burn in 
  // Notes: For instance, this is when wind farm stressors are activated
  // or when females become pregnant and gestation begins
  const int start_simulation = 93;
  
  // -------------------------------
  // PROPORTION OF LIPIDS IN MILK
  // -------------------------------
  // Value: Average of values reported for bowhead whales in McHuron et al. (2023)
  // Consistent with data from humpback whales as reported in Oftedal (1997) -- Figure 4
  // Units: % (d.u.)
  // Source(s): See Supplementary Information in McHuron et al. (2023)
  const double milk_lipids = 0.316;
  
  // -------------------------------
  // PROPORTION OF PROTEIN IN MILK
  // -------------------------------
  // Value: Average of values reported for bowhead whales in McHuron et al. (2023)
  // Consistent with data from humpback whales as reported in Oftedal (1997) -- Figure 4
  // Units: % (d.u.)
  // Source(s): See Supplementary Information in McHuron et al. (2023)
  const double milk_protein = 0.082;
  
  // -------------------------------
  // MAMMARY GLAND EFFICIENCY
  // -------------------------------
  // Definition: Efficiency of milk production from nutrients consumed
  // Units: % (d.u.)
  // Source(s): Brody et al. (1968)
  const double mammary_efficiency = 0.95; 
  
  // -------------------------------
  // NONLINEARITY BETWEEN MILK SUPPLY AND BODY CONDITION OF THE MOTHER
  // -------------------------------
  // Units: d.u.
  // const int zeta = -2;
  const double zeta = -7.5;
  
  // -------------------------------
  // AGE AT WHICH MILK COMNSUMPTION STARTS TO DECREASE
  // -------------------------------
  // Notes: Value reflecting the inflection point between the two phases of NARW growth.
  // Unit: days
  // Source(s): Fortune et al. (2020)
  const int milk_decrease = 288;
  
  // -------------------------------
  // NONLINEARITY BETWEEN MILK ASSIMILATION AND CALF AGE
  // -------------------------------
  // Units: d.u.
  // Source(s): Hin et al. (2019) 
  const double eta_milk = 0.9;
  
  // -------------------------------
  // DURATION OF LACTATION
  // -------------------------------
  // Units: days
  const int T_lac = 365;

  // -------------------------------
  // ENERGY DENSITY OF LIPIDS
  // -------------------------------
  // Units: kJ / kg - converted to MJ/kg
  // Source(s): Schmidt-Nielsen et al. (1997)
  const double ED_lipids = 39539.0L/1000;
  
  // -------------------------------
  // ENERGY DENSITY OF PROTEIN
  // -------------------------------
  // Notes: Other papers such as Beltran et al. (2018) have used Unif(17900,23800) 
  // Units: kJ / kg - converted to MJ/kg
  // Source(s): As per Christiansen et al. (2022)
  const double ED_protein = 23640.0L/1000;
  
  // -------------------------------
  // PROPORTION OF BODY VOLUME COMPRISING MUSCLE IN FETUS
  // -------------------------------
  // Source(s): Christiansen et al. (2022)
  const double P_muscle = 0.2820;
  
  // -------------------------------
  // PROPORTION OF BODY VOLUME COMPRISING VISCERA IN FETUS
  // -------------------------------
  // Source(s): Christiansen et al. (2022)
  const double P_viscera = 0.1020;
  
  // -------------------------------
  // PROPORTION OF BODY VOLUME COMPRISING BONES IN FETUS
  // -------------------------------
  // Source(s): Christiansen et al. (2022)
  const double P_bones = 0.1250;
  
  // -------------------------------
  // BLUBBER DENSITY
  // -------------------------------
  // Notes: Value for "Animal fat" in primary reference
  // Units: kg / m3 
  // Source(s): Charrondière et al. (2012)
  const int blubber_density = 700;
  
  // -------------------------------
  // MUSCLE DENSITY
  // -------------------------------
  // Notes: Value for "Cow, lean, no bone, raw" in primary reference
  // Units: kg / m3 
  // Source(s): Charrondière et al. (2012)
  const int muscle_density = 960;
  
  // -------------------------------
  // BONE DENSITY
  // -------------------------------
  // Notes: Value for "Bone and meat meal, loose" in primary reference 
  // Units: kg / m3 
  // Source(s): Charrondière et al. (2012)
  const int bone_density = 720;
  
  // -------------------------------
  // VISCERA DENSITY
  // -------------------------------
  // Notes: Value for "Cow, intestine, raw" in primary reference
  // Units: kg / m3 
  // Source(s): Charrondière et al. (2012)
  const int visceral_density = 930;
  
  // -------------------------------
  // RELATIVE PROPORTIONS OF LEAN TISSUES
  // -------------------------------
  // Notes: Used in calculations of growth costs; standardized to sum to one.
  // Units: % (d.u.)
  // Source(s): Christiansen et al. (2022)
  const double prop_muscle = 0.553797468;
  const double prop_viscera = 0.200949367;
  const double prop_bones = 0.245253165;
  
  // -------------------------------
  // PROPORTIONS OF LIPIDS/PROTEIN IN TISSUES
  // -------------------------------
  // Notes: There are no data on the energy content of different body tissues for NARW
  // Christiansen et al. (2022, MEPS) report values from studies of fin, sei, and minke whales.
  
  const double prop_lipids_in_muscle = 0.114;
  const double prop_protein_in_muscle = 0.221;
  
  const double prop_lipids_in_viscera = 0.758;
  const double prop_protein_in_viscera = 0.037;
  
  // Average of values reported by Montie et al. (2010)
  const double prop_lipids_in_blubber = 0.6116667;
  const double prop_protein_in_blubber = 0.1020;
  
  const double prop_lipids_in_bones = 0.218;
  const double prop_protein_in_bones = 0.248;
  
  //  ****************** Behavior ******************
  
  // -------------------------------
  // TIME SPENT IN ECHELON POSITION
  // -------------------------------
  // Definition: Proportion of time that young calves (< 1 month-old) spend swimming in echelon position
  // Units: % (d.u.)
  // Source(s): Noren et al. (2008)
  const double glide_echelon = 0.36;
  
  Rcpp::NumericMatrix residency_seus(6,2);

  // -------------------------------
  // MEAN RESIDENCY TIME IN SEUS
  // -------------------------------
  // Definition: Average time spent on the calving grounds
  // Units: days
  // Source(s): Krzystan et al. (2018)
  residency_seus(0,0) = 39.6; // Juveniles
  residency_seus(1,0) = 39.6; // Juveniles
  residency_seus(2,0) = 20.4; // Adult males
  residency_seus(3,0) = 41.1; // Pregnant females (non-calving adult)
  residency_seus(4,0) = 87.5; // Lactating females
  residency_seus(5,0) = 41.1; // Resting females (non-calving adult)

  // -------------------------------
  // SD RESIDENCY TIME IN SEUS
  // -------------------------------
  // Definition: Average time spent on the calving grounds
  // Units: days
  // Source(s): Krzystan et al. (2018)
  residency_seus(0,1) = 2.1; // Juveniles
  residency_seus(1,1) = 2.1; // Juveniles
  residency_seus(2,1) = 1.1; // Adult males
  residency_seus(3,1) = 6.3; // Pregnant females (non-calving adult)
  residency_seus(4,1) = 4.9; // Lactating females
  residency_seus(5,1) = 6.3; // Resting females (non-calving adult)
  
  // ------------------------------------------------------------
  // START SIMULATION
  // ------------------------------------------------------------
  
  // RcppThread::ProgressBar pbar(t, 1); // '1' means: print update in increments of 1%
  
  // N days (+1) steps per animal
  for(auto env = environments.begin(); env != env_end; ++env) {
    
    current_day = env - environments.begin();

    if(progress) bar.update();
    // if(progress) pbar++;
    
    // env is a vector iterator, so env - environments.begin() returns the associated index
    current_month = densitySeq[env - environments.begin()];
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {
      
      current_animal = animal - animals.begin();
      
      // Initial location
      if(current_day == 0){
        animal->x0 = animal->x;
        animal->y0 = animal->y;
      }
      
      current_x = animal->x;
      current_y = animal->y;
      animal->region = layers[current_day](current_x, current_y, 'R');
      
      // Mark arrival on calving grounds
      if(animal->region == 9 && animal->calvgrounds == 0){
        animal->calvgrounds = 1;
        enter_SEUS = current_day;
      }
      
      // Mark departure from calving grounds
      if(animal->region != 9 && animal->calvgrounds == 1){
        animal->calvgrounds = 0;
        depart_SEUS = current_day;
      }
      
      if(animal->alive == 0) animal->calvgrounds = 0;
      
      // Reset northward pull after 20 days
      if(animal->north == 1 && current_day == (depart_SEUS + 20)) animal->north = 0;
      
      // Stop second migration to SEUS
      if(current_day > 365 && cohortID == 5){
        animal->seus = 0;
      }
      
      // ------------------------------------------------------------
      // TARGET BC
      // ------------------------------------------------------------
      
      // Note: This is different from the maximum allowable BC (see feeding effort figure)
      
      // Upper target ––– for calves and pregnant females corresponds to 
      // the maximum condition observed in southern right whales (Christiansen, pers. comm + Figure 7D).
      // 
      // Lower target ––– for other animals corresponds to the weighted mean relative blubber mass 
      // reported in bowhead whales and NP right whales in the spreadsheet of model parameters
      // This value aligns with relative fat mass estimated for individual of "average body condition" 
      // in Christiansen et al. (2022) Fig 7
      
      // Pregnant females
      if((animal->cohortID == 4 && current_day >= start_simulation)){
        
        target_bc[0] = 0.6;
        target_bc[1] = 0;
        
      // Lactating mothers and their calves
      } else if(animal->cohortID == 5 && current_day >= start_simulation && born[current_animal] == 1){ 
          
          target_bc[0] = 0.3674076;
          // target_bc[0] = 0.665289;
          
          if(calves[current_animal].alive == 1){
            
            target_bc[1] = 0.6;
            
          } else {
            
            target_bc[1] = 0;  
            
          }
      
      // All other cohorts
      } else { 
        
        target_bc[0] = 0.3674076;
        // target_bc[0] = 0.41718;
        target_bc[1] = 0;
        
      }
      
      // Save positional data
      *(data_out++) = animal->x;
      *(data_out++) = animal->y;
      *(data_out++) = animal->region;
      *(data_out++) = animal->calvgrounds;
      *(data_out++) = residency_seus(cohortID-1,0);
      *(data_out++) = residency_seus(cohortID-1,1);
      *(data_out++) = p_leave;
      
      if(current_day > 0){
        
        // ------------------------------------------------------------
        // STRESSORS
        // ------------------------------------------------------------
        
        if(animal->alive){
          
          if(stressors){
            
            // ------------------------------------------------------------
            // STARVATION MORTALITY
            // ------------------------------------------------------------

            // No mortality during burn-in
            if(current_day >= start_simulation){
            
            // Body condition-dependent mortality, similar to Hin et al. (2019)
            // Probability of the animal starving given its health
            // p_starve = starvation_mortality(animal->bc, starvation_onset, starvation, starvation_scalar);
            p_starve = starvation_mortality(animal->bc, starvation_coefs, starvation_death, starvation_onset);
            p_starve = 1-(std::pow(1-p_starve, 1.0L/365.0L));
            has_starved = R::rbinom(1, p_starve);
            
            // p_starve_calf = starvation_mortality(calves[current_animal].bc, starvation_onset, starvation, starvation_scalar);
            p_starve_calf = starvation_mortality(calves[current_animal].bc, starvation_coefs, starvation_death, starvation_onset);
            p_starve_calf = 1-(std::pow(1-p_starve_calf, 1.0L/365.0L));
            has_starved_calf = R::rbinom(1, p_starve_calf);
            
            if(has_starved){
              
              animal->alive = 0;
              animal->starve = 1;
              date_death[current_animal] = current_day;
              
              if(calves[current_animal].alive == 1){
                calves[current_animal].alive = 0;
                date_death_calf[current_animal] = current_day;
              }
            }
            
            if(calves[current_animal].alive == 1 && born[current_animal] == 1) {
              if(has_starved_calf){
                calves[current_animal].alive = 0;
                calves[current_animal].starve = 1;
                date_death_calf[current_animal] = current_day;
                animal->north = 1;
              }
            }
            
            // ===================================
            // ENTANGLEMENT
            // ===================================
            
            gear_risk = layers[current_day](current_x, current_y, 'F');
            
            // ++++ Juveniles and adults ++++
            
            // An animal can only become entangled if not already carrying gear 
            if(animal->entangled(0) == 0){
              
              // Simulate entanglement event with given probability
              animal->entangled = entanglement_event(gear_risk);
              
              // Start and end dates of entanglement event
              if(animal->entangled(0) == 1){
                animal->entangled(5) = current_day;
                animal->entangled(6) = current_day + animal->entangled(4);
              }
              
            } else {
              
              // Terminate entanglement event
              if(current_day == animal->entangled(6)){
                animal->entangled = entanglement_event(0); // Reset entanglement state
              }
              
            } // End if entangled
            
            
            // ++++ Calves ++++
            
            if(calves[current_animal].alive == 1){
              
              // If the calf is not carrying gear
              if(calves[current_animal].entangled(0) == 0){
                
                calves[current_animal].entangled = entanglement_event(gear_risk);
                
                if(calves[current_animal].entangled(0) == 1){
                  calves[current_animal].entangled(5) = current_day;
                  calves[current_animal].entangled(6) = current_day + calves[current_animal].entangled(4);
                }
                
              } else {
                
                // Terminate entanglement event
                if(current_day == calves[current_animal].entangled(6)){
                  calves[current_animal].entangled = entanglement_event(0); // Reset entanglement state
                }
                
              } 
            } // End if calves entangled
            
              
              // ===================================
              // VESSEL STRIKES
              // ===================================
              
              // Vessel strike risk
              strike_risk = layers[current_day](current_x, current_y, 'V');
              
              // Incidence of a vessel strike -- (d.u.)
              if(animal->strike == 0) animal->strike = R::rbinom(1, strike_risk);
              
              // ++++ Mortality from vessel strikes ++++
              
              if(animal->strike){
                
                animal->alive = 0;
                date_death[current_animal] = current_day;
                
                if(calves[current_animal].alive == 1){
                  calves[current_animal].alive = 0;
                  date_death_calf[current_animal] = current_day;
                  
                }
                
              } else if(animal->strike == 0){
                
                if(calves[current_animal].alive & calves[current_animal].strike == 0) 
                  calves[current_animal].strike = R::rbinom(1, strike_risk);
                
                if(calves[current_animal].alive && calves[current_animal].strike){
                  calves[current_animal].alive = 0;
                  date_death_calf[current_animal] = current_day;
                  animal->north = 1;
                }
                
              }
              
              // ===================================
              // AGE-DEPENDENT MORTALITY
              // ===================================
              
              // Annual survival probability based on Linden (2023) – converted to daily probability
              p_surv = std::pow(survival(animal->age, 0, animal->sex), 1.0L/365);
              p_surv_calf = std::pow(survival(calves[current_animal].age, 0, calves[current_animal].sex), 1.0L/365);
              
              // If the animal is entangled, add the initial effect ("hit") of entanglement on survival
              // on the first day of entanglement.
              // This requires converting to a 3-month survival probability (taken to be 90 days),
              // then adding the gear effect on the complementary log-log scale
              if(animal->entangled(0) == 1 && current_day == animal->entangled(5)){
                p_surv = entanglement_effect(p_surv, animal->entangled(1), 90, entgl_pcoms);
              }
              
              if(calves[current_animal].entangled(0) == 1 && current_day == calves[current_animal].entangled(5)){
                p_surv_calf = entanglement_effect(p_surv_calf, calves[current_animal].entangled(1), 90, entgl_pcoms);
              }
              
              // Run Bernoulli trials
              if(animal->mort == 0) animal->mort = R::rbinom(1, 1-p_surv);
              if(calves[current_animal].mort==0) calves[current_animal].mort = R::rbinom(1, 1-p_surv_calf);
              
              if(animal->alive == 1 && animal->mort == 1){
                
                // Record death and date
                animal->alive = 0;
                date_death[current_animal] = current_day;
                
                if(calves[current_animal].alive == 1){
                  // Kill the calf
                  calves[current_animal].alive = 0;
                  date_death_calf[current_animal] = current_day;
                }
                
              }
              
              if(calves[current_animal].alive == 1 && calves[current_animal].mort == 1){
                
                // Record death and date
                calves[current_animal].alive = 0;
                date_death_calf[current_animal] = current_day;
                
                // Ensure mother leaves calving grounds
                animal->north = 1;
              }
              
            } // End no mortality during burn-in
            
          } // End if stressors
          
        } // End if alive
        
        if(animal->alive){
          
          // Retrieve the seed for the current animal
          int currentseed = doseresp_seed(current_animal, cohortID - 1);
          
          // ------------------------------------------------------------
          // MOVEMENT
          // ------------------------------------------------------------
          
          // Update the animal's state
          m.update(*animal, *env, animal->region, support, 
                   limits, limits_regions, limits_prey,
                   limits_fishing, limits_vessels, limits_noise,
                   resolution, resolution_regions, resolution_prey, 
                   resolution_fishing, resolution_vessels, resolution_noise);
          
          // If pregnant females, resting females, or lactating females in their
          // last trimester move past Cape Hatteras, make them turn around and head north again
          if(cohortID == 6 | cohortID == 4 | (cohortID == 5 && current_day > 365)){
            if(animal->y <= 220) animal->north = 1;
          }
         
          // ------------------------------------------------------------
          // NOISE EXPOSURE & BEHAVIORAL RESPONSES
          // ------------------------------------------------------------
          
          animal->respnoise = 0; // Reset daily
          
          if(current_day >= start_simulation){
          
          // Noise levels +++
          dB = layers[current_day](current_x, current_y, 'N');
          
          // Response threshold (dB)
          response_dB = response_threshold(doseresp, current_day, t, currentseed);
          
          // Behavioral response to pile-driving noise exposure -- (d.u.)
          if(dB >= response_dB) animal->respnoise = 1;
          
          }
          
          // ------------------------------------------------------------
          // BIRTHS
          // ------------------------------------------------------------
          
          // Fixed birth date of Jan 1st
          if(cohortID == 5 && current_day == start_simulation){
            calves[current_animal].alive = 1;
            born[current_animal] = 1;
            dob[current_animal] = current_day;
          }
          
          // if(cohortID == 5 && animal->region == 9 && born[current_animal] == 0){
          //   prob_birth = pbirth(current_day, enter_SEUS, 30);
          //   calves[current_animal].alive = R::rbinom(1, prob_birth);
          // }
          
          // if(calves[current_animal].alive == 1 && born[current_animal] == 0){
          //   born[current_animal] = 1;
          //   dob[current_animal] = current_day;
          // }
          
          // ------------------------------------------------------------
          // Distance traveled
          // ------------------------------------------------------------

          travel_dist = std::sqrt(std::pow(current_x - animal->x, 2) + std::pow(current_y - animal->y, 2));
          
          // ------------------------------------------------------------
          // PREY
          // ------------------------------------------------------------
          
          // Prey concentration in cell at time t -- (g/m3)
          D_cop = prey_scale * layers[current_day](current_x, current_y, 'P');
          
          // Feeding response
          is_feeding = feeding_threshold(min_calanus, D_cop);
          
          // ------------------------------------------------------------
          // ACTIVITY BUDGET
          // ------------------------------------------------------------
          
          // Total time available for foraging activities --(hr converted to s)
          // daylight_hours = layers[current_day](current_x, current_y, 'L');
          
          // Swimming speed (during travel)
          // 4.6 m/s taken as maximum average speed based on satellite tag record described in Mate (1997)
          // Minimum chosen to ensure that travel_time does not exceed 24 hours
          // Mean and SD values derived from meta-analysis of spreadsheet of model parameters
          
          if(is_feeding == 0){
            
            // MIDA and SCOS
            if((animal->region == 7 | animal->region == 8)){
              
              swim_speed = rtnorm(0.9771623, 0.1847643, (travel_dist*1000)/(24*3600), 4.6);
              
            // Elsewhere
            } else {
              
              swim_speed = rtnorm(0.4536798, 0.2290059, (travel_dist*1000)/(24*3600), 4.6);
              
            }
            
          } else if (is_feeding){
            
            swim_speed = rtnorm(0.9066851, 0.1777602, (travel_dist*1000)/(24*3600), 4.6);
            
          }
          
          // Time spent traveling per day -- (hr)
          // Calculated based on distance covered that day and swimming speed
          travel_time = (travel_dist * 1000) / (swim_speed * 3600); // m per hour
          
          // Behavioral modes include: traveling, feeding, and resting/nursing
          
          if(animal->region == 9){ // SEUS
            
            feeding_time = 0;
            
          } else { // Feeding grounds and elsewhere
            
            // tmax = std::min(24 - travel_time, daylight_hours); // Daytime foraging
            tmax = 24 - travel_time; // Round-the-clock foraging
            // Cannot have tmax < lower bound
            if(tmax > 12) min_feed = 12;
            else min_feed = 0.0L;
            
            feeding_time = is_feeding * rtnorm(16.11904, 0.6199398, min_feed, tmax);
 
          }
          
          // ------------------------------------------------------------
          // RESPONSE TO NOISE
          // -----------------------------------------------------------
          
          // Lower the time spent foraging following a behavioral response
          if(animal->respnoise){
            feeding_time -= piling_hrs;
            if(feeding_time < 0) feeding_time = 0;
          }

          resting_nursing_time = 24 - (travel_time + feeding_time);
          
          // ------------------------------------------------------------
          // ENERGY INTAKE - FORAGING
          // ------------------------------------------------------------
          
          // Foraging / nursing response to body condition -- (d.u.)
          // ---------------------
          // Pregnant females
          if((animal->cohortID == 4 && current_day >= start_simulation) | 
             (animal->cohortID == 5 && born[current_animal] == 0)){ 
            
            // Values for 5-parameter logistic curve obtained through optimization using target BC passed to optim_feeding
            // feeding_nursing_effort[0] = scale_effort(animal->bc, 1, 0, logis_slope(1), logis_location(1), logis_asymmetry(1));
            feeding_nursing_effort[0] = feeding_effort(eta(1), target_bc[0], animal->bc);
            
          } else if(animal->cohortID == 5 && born[current_animal] == 1){ // Lactating mothers and their calves
            
            feeding_nursing_effort[0] = feeding_effort(eta(0), target_bc[0], animal->bc);
            feeding_nursing_effort[1] = feeding_effort(eta(1), target_bc[1], calves[current_animal].bc);
            
          } else { // All other cohorts
            
            feeding_nursing_effort[0] = feeding_effort(eta(0), target_bc[0], animal->bc);
            
          }
          
          // Prey energy content -- (MJ / g)
          int meta = animal->cohortID > 2;
          E_calanus = (1-water_content) * cop_kJ * digestive_efficiency * urinary_efficiency * metabolizing_efficiency[meta];
          
          // ENERGY INTAKE (MJ)
          E_in = is_feeding * D_cop * animal->gape * (1 - animal->entangled(3)) *
            swim_speed * capture_efficiency * (feeding_time * 3600) * feeding_nursing_effort[0] *
            E_calanus;
          
          // ------------------------------------------------------------
          // ENERGY INTAKE - SUCKLING
          // ------------------------------------------------------------
          
          if(cohortID == 5){
              
              // Skeleton code --------------------------------------------
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
              // --------------------------------------------
              
              // Milk assimilation -- (d.u.)
              // Varies with calf age
              milk_assim = (calves[current_animal].alive == 1) * 
                milk_assimilation(current_day-start_simulation+1, T_lac, milk_decrease, eta_milk);
              
              // Milk provisioning -- (d.u.)
              // Varies with mother's body condition
              milk_provisioning = (calves[current_animal].alive == 1) * 
                milk_supply(starvation_onset, target_bc[0], animal->tot_mass, animal->fat_mass, zeta);
              
              // Total mass of mammary glands -- (kg)
              mammary_M = mammary_mass(animal->tot_mass);
              
              // Milk production rate by the mother -- (kg per day)
              milk_production_rate = milk_production(mammary_M);
              
              // Time spent nursing/suckling per day -- (hours to sec)
              // Meta:: No filter on time spent nursing, time spent suckling during nursing bouts filtered by species (value for E. australis)
              // t_suckling = rtnorm(0, 3.688672, 0, nursing_time) * 3600;
              // t_suckling = 10*3600;
              
              // ENERGY INTAKE (MJ)
              
              // Nursing cessation coded as reduction in energy intake as time spent nursing
              // is combined with time spent resting
              if(animal->respnoise){
                if(nursing_cessation(0)){
                  feeding_nursing_effort[1] *= ((24-nursing_cessation(1))/24);
                  if(feeding_nursing_effort[1] < 0) feeding_nursing_effort[1] = 0;
                }
              }
              
              E_lac = (calves[current_animal].alive == 1) * 
                milk_provisioning * milk_production_rate * 
                feeding_nursing_effort[1] * milk_assim * 
                (milk_lipids * ED_lipids + milk_protein * ED_protein);
              
              // Calves start independently foraging only once milk consumption begins to decrease
              if(current_day >= start_simulation + milk_decrease){

                calf_foraging = (calves[current_animal].alive == 1) * is_feeding * 
                  D_cop * calves[current_animal].gape * 
                  (1 - calves[current_animal].entangled(3)) * swim_speed * capture_efficiency * 
                  feeding_nursing_effort[1] * (feeding_time * 3600) * (1-milk_assim) *
                  E_calanus;
                
              }
              
              E_in_calves = mammary_efficiency * E_lac + calf_foraging;
              
              // } # Separation conditional
            
          } // End if cohort == 5
          
          // ------------------------------------------------------------
          // LOCOMOTORY COSTS
          // ------------------------------------------------------------
          
          // Metabolic scalar for RMR based on age class -- (d.u.)
          // 1.4 increase in RMR for calves to account for elevated demands associated with active growth?
          // if(cohortID <= 2) phi_rmr = 1.2; else phi_rmr = 1;
          
          // Resting metabolic rate -- (MJ)
          // Calculated based on lean mass as blubber is more or less metabolically inert (George et al. 2021)
          resting_metab = RMR(animal->lean_mass);
          if(cohortID == 5 && calves[current_animal].alive == 1) resting_metab_calves = RMR(calves[current_animal].lean_mass); 
          
          // Stroke rate during routine swimming -- (d.u.)
          stroke_rate_foraging = rtnorm(0.1650784, 0.003930294, 0, INFINITY);
          stroke_rate = rtnorm(0.11375, 0.02678894, 0, INFINITY);
          
          // Percentage of time spent gliding -- (d.u.)
          percent_glide_foraging = R::rgamma(17.393772, 0.02006944); // Corresponds to mean of 36%, min of 10% and max of 75%
          percent_glide = R::rnorm(0.09, 0.00795); // Corresponds to mean of 9, min of 6 and max of 12%
          
          // Locomotory costs (J converted to MJ)
          // ----------------
          // Based on Noren 2008 (DOI: 10.1111/j.1365-2435.2007.01354.x) and 
          // Noren et al. (2011) (DOI: 10.1242/jeb.059121)
          // Apply 17% increase in stroke frequency for lactating females with 0-1 
          // month-old calves swimming in echelon position.
          // Increase stroke rate by 16.2791% (1/0.86) for pregnant females 
          // in their last month of pregnancy.
          
          // Pregnant females in the last month of pregnancy
          if((cohortID == 4 && current_day > start_simulation) |
             (cohortID == 5 && current_day < start_simulation)){
            
            if(current_month == densitySeq.back()){
              scalar_LC = 1/0.86;
            } else {
              scalar_LC = 1.0L;
            }
            
          // Lactating mothers with calf swimming in echelon position
          } else if(cohortID == 5 && current_day >= start_simulation){
            
            if(calves[current_animal].alive == 1 && calves[current_animal].age <= 30/365.0L){
              
              scalar_LC = 1.17;
              scalar_LC_calf = 0.83;
              percent_glide_echelon = glide_echelon;
              
            } else {
              
              scalar_LC = 1.0L;
              scalar_LC_calf = 1.0L;
              percent_glide_echelon = percent_glide;
            }
            
          } else {
            
            scalar_LC = 1.0L;
          }
          
          LC_tot = locomotor_costs(animal->tot_mass, 
                                   travel_dist*1000,
                                   stroke_rate, 
                                   stroke_rate_foraging, 
                                   percent_glide, 
                                   percent_glide_foraging,
                                   feeding_time * 3600, 
                                   travel_time * 3600, 
                                   scalar_LC);
          
          if(cohortID == 5 && calves[current_animal].alive == 1){
            LC_tot_calves = locomotor_costs(calves[current_animal].tot_mass,
                                            travel_dist*1000,
                                            stroke_rate, 
                                            stroke_rate_foraging, 
                                            percent_glide_echelon, 
                                            percent_glide_echelon,
                                            0, 
                                            travel_time * 3600, 
                                            scalar_LC_calf);
          }
          
          // ------------------------------------------------------------
          // FETAL DEVELOPMENT
          // ------------------------------------------------------------
            
            // Nil gestation costs for:
            // (1) Calves, juveniles, adult males
            // (2) Resting females
            // (3) Females that have terminated their pregnancy
            // (4) Pregnant females during burn-in
            // (5) Pregnant females on very last day (birth) 
            // (6) Lactating females after burn-in

            if(cohortID < 4 | cohortID == 6 | animal->abort == 1 |
               (cohortID == 4 && current_day < start_simulation) |
               (cohortID == 4 && current_day == densitySeq.size()) |
               (cohortID == 5 && current_day > start_simulation))
              {
              
              fetal_growth_cost = 0.0L;
              placental_cost = 0.0L;
              HIC = 0.0L;
              E_gest = 0.0L;
              
            } else {
              
              // Define a new day counter to calculate correct gestation costs
              // for lactating females starting the simulation as late pregnant
              int gday = 0;
              
              if(cohortID == 4 && current_day >= start_simulation){
                gday = current_day - start_simulation + 1;
              } else if(cohortID == 5 && current_day < start_simulation){
                gday = current_day + 272;
              }
              
              fetus_l = fetal_length(gday - 365, animal->length);
              fetus_m = fetal_mass(gday - 365, animal->length);
              
              fetus_l_next = fetal_length(gday - 364, animal->length);
              fetus_m_next = fetal_mass(gday - 364, animal->length);
              
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
                                                0,
                                                mass_muscle,
                                                mass_viscera,
                                                mass_bones,
                                                blubber_density,
                                                muscle_density,
                                                visceral_density,
                                                bone_density);
              
              mass_blubber_next = fetal_blubber_mass(fetus_l_next,
                                                     0,
                                                     mass_muscle_next,
                                                     mass_viscera_next,
                                                     mass_bones_next,
                                                     blubber_density,
                                                     muscle_density,
                                                     visceral_density,
                                                     bone_density);
              
              // Energetic cost of fetal growth during pregnancy -- (MJ) [G]
              FG0 = mass_muscle * (ED_lipids * prop_lipids_in_muscle + ED_protein * prop_protein_in_muscle) +
                mass_viscera * (ED_lipids * prop_lipids_in_viscera + ED_protein * prop_protein_in_viscera) +
                mass_bones * (ED_lipids * prop_lipids_in_bones + ED_protein * prop_protein_in_bones) +
                mass_blubber * (ED_lipids * prop_lipids_in_blubber +  ED_protein * prop_protein_in_blubber);
              
              FG1 = mass_muscle_next * (ED_lipids * prop_lipids_in_muscle + ED_protein * prop_protein_in_muscle) +
                mass_viscera_next * (ED_lipids * prop_lipids_in_viscera + ED_protein * prop_protein_in_viscera) +
                mass_bones_next * (ED_lipids * prop_lipids_in_bones + ED_protein * prop_protein_in_bones) +
                mass_blubber_next * (ED_lipids * prop_lipids_in_blubber +  ED_protein * prop_protein_in_blubber);
              
              fetal_growth_cost = FG1 - FG0;
              
              // Energetic cost of placental maintenance during pregnancy -- (MJ) [P]
              placental_cost = placental_maintenance(fetal_growth_cost);
              
              // Heat increment of gestation -- (kJ converted to MJ) [Q]
              HIC = heat_gestation(birth_mass[current_animal], mass_increment_fetus)/1000;
              
              // Total cost of gestation -- (MJ)
              E_gest = fetal_growth_cost + placental_cost + HIC;
              
            } // End if current_day == 365
          
          // ------------------------------------------------------------
          // SOMATIC GROWTH
          // ------------------------------------------------------------
          
          // Cost of growth (adults and juveniles) -- (MJ)
          length_nextday = age2length(animal->age + 1.0L/365.0L, animal->lengthatage);
          lean_mass_nextday = length2mass(length_nextday, animal->massatlength);
          lean_mass_increment = (lean_mass_nextday - animal->lean_mass);
          
          // double sumprop = animal->bc + prop_muscle + prop_viscera + prop_bones;
          // bc_diff = (1 - sumprop)/3;
          
          E_growth = growth_cost(lean_mass_increment, 
                                 ED_lipids, 
                                 ED_protein, 
                                 prop_lipids_in_muscle,
                                 prop_lipids_in_viscera,
                                 prop_lipids_in_bones,
                                 prop_protein_in_muscle,
                                 prop_protein_in_viscera,
                                 prop_protein_in_bones,
                                 prop_muscle,
                                 prop_viscera,
                                 prop_bones);
          
          if(cohortID == 5 && calves[current_animal].alive == 1){
            
            length_nextday_calves = age2length(calves[current_animal].age + 1.0L/365.0L, calves[current_animal].lengthatage);
            lean_mass_nextday_calves = length2mass(length_nextday_calves, calves[current_animal].massatlength);
            lean_mass_increment_calves = (lean_mass_nextday_calves - calves[current_animal].lean_mass);
            
            E_growth_calves = growth_cost(lean_mass_increment_calves,
                                               ED_lipids,
                                               ED_protein,
                                               prop_lipids_in_muscle,
                                               prop_lipids_in_viscera,
                                               prop_lipids_in_bones,
                                               prop_protein_in_muscle,
                                               prop_protein_in_viscera,
                                               prop_protein_in_bones,
                                               prop_muscle,
                                               prop_viscera,
                                               prop_bones);
            
          }
          
          // ------------------------------------------------------------
          // ENERGY EXPENDITURE
          // ------------------------------------------------------------
          
          // Adults and juveniles -- (MJ)
          E_out = (resting_metab + LC_tot + E_gest + E_lac + E_growth);
          if(animal->entangled(0) == 1) E_out = E_out + entanglement_cost;
          
          // Calves -- (MJ)
          if(cohortID == 5 && calves[current_animal].alive == 1){
            E_out_calves = (resting_metab_calves + LC_tot_calves + E_growth_calves);
            if(calves[current_animal].entangled(0) == 1) E_out_calves = E_out_calves + entanglement_cost;
          }
          
          // ------------------------------------------------------------
          // ENERGY ALLOCATION
          // ------------------------------------------------------------
          
          if(animal->alive == 1){
            
            // Age (+ 1 day) -- Long double precision required here
            animal->age = animal->age + 1.0L/365.0L;
            
            // Body length
            animal->length = age2length(animal->age, animal->lengthatage);
            
            // Mouth width
            animal->mouth_width = animal->length * animal->mouth_ratio;
            
            // Mouth gape area -- (m^2)
            animal->gape = gape_size(animal->length, animal->mouth_width, animal->mouth_angle);
              
            // Lean mass
            animal->lean_mass = lean_mass_nextday;
            
            // E_balance = 0 (false) when gains exceed costs (-> anabolism)
            // E_balance = 1 (true) when costs exceed gains (-> catabolism)
            E_balance = E_in < E_out;
            
            // Lipid growth or breakdown
            delta_blubber = (E_in - E_out) * anabolism_efficiency[E_balance] / (ED_lipids * catabolism_efficiency[E_balance]);
            
            // // Costs of fat metabolism
            // double E_fat = fatdeposition_cost(delta_blubber,
            //                                   ED_lipids, 
            //                                   ED_protein,
            //                                   prop_lipids_in_blubber,
            //                                   prop_protein_in_blubber);
            // 
            // E_out = E_out + E_fat;
            // 
            // // E_balance = 0 (false) when gains exceed costs (-> anabolism)
            // // E_balance = 1 (true) when costs exceed gains (-> catabolism)
            // E_balance = E_in < E_out;
            // 
            // // Lipid growth or breakdown
            // delta_blubber = (E_in - E_out) * anabolism_efficiency[E_balance] / (ED_lipids * catabolism_efficiency[E_balance]);
            
            if(growth){
              // Use std::max to ensure that fat_mass never becomes negative
              animal->fat_mass = std::max(0.0L, animal->fat_mass + delta_blubber);
              animal->tot_mass = animal->fat_mass + animal->lean_mass;
              animal->bc = animal->fat_mass / animal->tot_mass;
            }
            
            // Same for calves
            if(cohortID == 5 && calves[current_animal].alive == 1){
              
              E_balance_calves = E_in_calves < E_out_calves;
              
              calves[current_animal].age = calves[current_animal].age + 1.0L/365.0L;
              calves[current_animal].length = age2length(calves[current_animal].age, calves[current_animal].lengthatage);
              calves[current_animal].lean_mass = lean_mass_nextday_calves;
              
              // Update mouth ratio and mouth width for calves > 1 month of age
              if(current_day == start_simulation + 30){
                calves[current_animal].mouth_ratio = start_mouth(0,31.0L/365.0L);
              }
              
              calves[current_animal].mouth_width = calves[current_animal].mouth_ratio*calves[current_animal].length;
              
              calves[current_animal].gape = gape_size(calves[current_animal].length,
                                                      calves[current_animal].mouth_width, 
                                                      calves[current_animal].mouth_angle);
              
              // Lipid growth or breakdown
              delta_blubber_calves = (E_in_calves - E_out_calves) * anabolism_efficiency[E_balance_calves] / (ED_lipids * catabolism_efficiency[E_balance_calves]);
              // delta_blubber_calves = lipid_anacat_calves[E_balance_calves] * (E_in_calves - E_out_calves) / ED_lipids;
              
              if(growth){
                calves[current_animal].fat_mass = calves[current_animal].fat_mass + delta_blubber_calves;
                calves[current_animal].tot_mass = calves[current_animal].fat_mass + calves[current_animal].lean_mass;
                calves[current_animal].bc = calves[current_animal].fat_mass / calves[current_animal].tot_mass;
              }
              
            }
            
          } // End energy allocation
          
          
          // ------------------------------------------------------------
          // ABORTION
          // ------------------------------------------------------------
          
          if((cohortID == 4 && animal->abort == 0) |
             (cohortID == 5 && animal->abort == 0 && current_day <= start_simulation)){
            
            bc_gest = abortion_scalar * (((E_out / (ED_lipids * catabolism_efficiency[1])) + starvation_onset*animal->tot_mass) / animal->tot_mass);
            
            // Abort pregnancy if the animal's reserves are insufficient to sustain its daily energetic costs
            if(bc_gest > animal->bc){
              animal->abort = 1;
              animal->north = 1;
            }
          }
          
          // ------------------------------------------------------------
          // LEAVE CALVING GROUNDS
          // ------------------------------------------------------------

          // Probability of migrating north is conditional on time already spent on calving grounds
          if(animal->north == 0 && animal->region == 9){
            double p_leave = pleave(current_day, enter_SEUS, cohortID, 3, residency_seus);
            animal->north = R::rbinom(1, p_leave);
          }
          
        } // End if alive
      } // End if current_day > 0
      
      // ------------------------------------------------------------
      // STORE VALUES in output matrices
      // ------------------------------------------------------------
      
      // // Restart if(alive) conditional otherwise some attributes may still be
      // // recorded despite animals having died of starvation
      if(current_day == 0){ // Record initial parameter values
        
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
        *(attrib_out++) = animal->gape;                    // Size of the gape
        *(attrib_out++) = animal->mouth_ratio;             // Ratio of mouth width to total body length
        *(attrib_out++) = animal->mouth_angle;             // Upper angle of the mouth
        *(attrib_out++) = animal->mouth_width;             // Width of the mouth
        *(attrib_out++) = animal->gsl;                     // Gulf of St Lawrence
        *(attrib_out++) = animal->seus;                    // Southeastern U.S.
        *(attrib_out++) = animal->north;                   // Migration interruption
        *(attrib_out++) = animal->abort;                   // Abortion (pregnant females)
        *(attrib_out++) = p_starve;                        // Starvation-induced mortality rate
        *(attrib_out++) = animal->starve;                  // Incidence of starvation
        *(attrib_out++) = animal->mort;                    // Age-dependent mortality
        *(attrib_out++) = date_death[current_animal];      // Date of death
        *(attrib_out++) = p_surv;                          // Survival probability
        
        // +++ Stressors +++
        
        *(attrib_stressors_out++) = 0.0L;                  // Probability of becoming entangled
        *(attrib_stressors_out++) = 0.0L;                  // Entanglement cost
        *(attrib_stressors_out++) = animal->entangled(0);  // Is the animal entangled?
        *(attrib_stressors_out++) = animal->entangled(1);  // Entanglement severity
        *(attrib_stressors_out++) = animal->entangled(2);  // Does the entanglement involve the anterior part of the body?
        *(attrib_stressors_out++) = animal->entangled(3);  // Reduction in mouth gape
        *(attrib_stressors_out++) = animal->entangled(4);  // Entanglement duration
        *(attrib_stressors_out++) = animal->entangled(5);  // Day when entanglement event started
        *(attrib_stressors_out++) = animal->entangled(6);  // Day when entanglement event ended
        
        *(attrib_stressors_out++) = calves[current_animal].entangled(0);  // Is the animal entangled? (calves)
        *(attrib_stressors_out++) = calves[current_animal].entangled(1);  // Entanglement severity (calves)
        *(attrib_stressors_out++) = calves[current_animal].entangled(2);  // Does the entanglement involve the anterior part of the body? (calves)
        *(attrib_stressors_out++) = calves[current_animal].entangled(4);  // Entanglement duration (calves)
        *(attrib_stressors_out++) = calves[current_animal].entangled(5);  // Day when entanglement event started (calves)
        *(attrib_stressors_out++) = calves[current_animal].entangled(6);  // Day when entanglement event ended (calves)
        
        *(attrib_stressors_out++) = 0.0L;                                 // Probability of being struck by vessels
        *(attrib_stressors_out++) = animal->strike;                       // Incidence of vessel strike
        *(attrib_stressors_out++) = calves[current_animal].strike;        // Incidence of vessel strike
        
        *(attrib_stressors_out++) = animal->respnoise;     // Incidence of behavioral response to noise exposure
        *(attrib_stressors_out++) = 0.0L;                  // Noise levels
        *(attrib_stressors_out++) = 0.0L;                  // Behavioral response threshold
        *(attrib_stressors_out++) = 0.0L;                  // Duration of pile-driving activities
        *(attrib_stressors_out++) = 0.0L;                  // Nursing cessation
        *(attrib_stressors_out++) = 0.0L;                  // Duration of nursing cessation
        
        // +++ Energy budget +++
        
        *(attrib_E_out++) = 0.0L;                  // Net energy (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Energy gain (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Energy expenditure (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Change in blubber mass (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Lipid anabolism efficiency (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Lipid catabolism efficiency (adults and juveniles)
        *(attrib_E_out++) = 0.0L;                  // Energy density of lipids
        *(attrib_E_out++) = 0.0L;                  // Energy density of protein
        
        // +++ Energy intake +++
        
        *(attrib_feeding_out++) = 0.0L;             // Incidence of foraging
        *(attrib_feeding_out++) = 0.0L;             // Prey concentration
        *(attrib_feeding_out++) = 0.0L;             // Minimum prey concentration needed for foraging
        *(attrib_feeding_out++) = 0.0L;             // Swimming speed while foraging
        *(attrib_feeding_out++) = 0.0L;             // Capture efficiency
        *(attrib_feeding_out++) = 0.0L;             // % reduction in the gape when entangled
        *(attrib_feeding_out++) = 0.0L;             // Feeding/nursing effort
        *(attrib_feeding_out++) = 0.0L;             // Steepness of feeding/nursing effort curve for lower target BC
        *(attrib_feeding_out++) = 0.0L;             // Steepness of feeding/nursing effort curve for higher target BC
        *(attrib_feeding_out++) = 0.0L;             // Target body condition
        *(attrib_feeding_out++) = 0.0L;             // Copepod mass (wet weight)
        *(attrib_feeding_out++) = 0.0L;             // Water content of prey
        *(attrib_feeding_out++) = 0.0L;             // Energy density of Calanus finmarchicus stage V
        *(attrib_feeding_out++) = 0.0L;             // Digestive efficiency
        *(attrib_feeding_out++) = 0.0L;             // Metabolizing efficiency (calves and juveniles)
        *(attrib_feeding_out++) = 0.0L;             // Metabolizing efficiency (adults)
        *(attrib_feeding_out++) = 0.0L;             // Prey energy content
        *(attrib_feeding_out++) = 0.0L;             // Calf foraging
        
        // +++ Energy costs +++
        
        *(attrib_costs_out++) = 0.0L;           // Resting metabolic rate
        *(attrib_costs_out++) = 0.0L;           // Locomotory costs
        *(attrib_costs_out++) = 0.0L;           // Scalar on locomotory costs
        *(attrib_costs_out++) = 0.0L;           // Stroke rate during non-foraging activities
        *(attrib_costs_out++) = 0.0L;           // Stroke rate during foraging activities
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with growth
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with gestation
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with fetal growth
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with placental maintenance
        *(attrib_costs_out++) = 0.0L;           // Heat increment of gestation
        *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with lactation
        *(attrib_costs_out++) = 0.0L;           // Daily change in mass
        *(attrib_costs_out++) = 0.0L;           // Relative % of muscles in lean mass
        *(attrib_costs_out++) = 0.0L;           // Relative % of viscera in lean mass
        *(attrib_costs_out++) = 0.0L;           // Relative % of bones in lean mass
        
        // +++ Activity budgets +++
        
        *(attrib_activity_out++) = 0.0L;         // Distance traveled
        *(attrib_activity_out++) = 0.0L;         // Swimming speed
        *(attrib_activity_out++) = 0.0L;         // Time spent gliding during non-foraging activities
        *(attrib_activity_out++) = 0.0L;         // Time spent gliding during foraging activities
        *(attrib_activity_out++) = 0.0L;         // Time spent gliding during echelon swimming
        *(attrib_activity_out++) = 0.0L;         // Time spent traveling
        *(attrib_activity_out++) = 0.0L;         // Time spent feeding
        *(attrib_activity_out++) = 0.0L;         // Time spent resting + nursing
        
        if(cohortID == 4 | cohortID == 5){
          
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
          *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal bones
          *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal blubber
          *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal blubber
          *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising muscle in fetus
          *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising viscera in fetus
          *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising bones in fetus
          *(attrib_fetus_out++) = 0.0L;          // Density of blubber
          *(attrib_fetus_out++) = 0.0L;          // Density of muscles
          *(attrib_fetus_out++) = 0.0L;          // Density of viscera
          *(attrib_fetus_out++) = 0.0L;          // Density of bones
          
        }
        
        *(attrib_calves_out++) = calves[current_animal].cohortID ; // Unique cohort identifier
        *(attrib_calves_out++) = born[current_animal];             // Birth
        *(attrib_calves_out++) = prob_birth;                       // Probability of Birth
        *(attrib_calves_out++) = dob[current_animal];              // Date of birth
        *(attrib_calves_out++) = date_death_calf[current_animal];  // Date of death
        *(attrib_calves_out++) = 0.0L;                             // Alive or dead
        *(attrib_calves_out++) = 0.0L;                             // Age
        *(attrib_calves_out++) = 0.0L;                             // Body condition
        *(attrib_calves_out++) = 0.0L;                             // Body length
        *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
        *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
        *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
        *(attrib_calves_out++) = 0.0L;                             // Total body mass (including blubber)
        *(attrib_calves_out++) = 0.0L;                             // Lean body mass (excluding blubber)
        *(attrib_calves_out++) = 0.0L;                             // Blubber mass
        *(attrib_calves_out++) = 0.0L;                             // Intercept of the logarithmic mass-at-length relationship
        *(attrib_calves_out++) = 0.0L;                             // Slope of the logarithmic mass-at-length relationship
        *(attrib_calves_out++) = 0.0L;                             // Size of the gape
        *(attrib_calves_out++) = 0.0L;                             // Ratio of mouth width to total body length
        *(attrib_calves_out++) = 0.0L;                             // Upper angle of the mouth
        *(attrib_calves_out++) = 0.0L;                             // Width of the mouth
        *(attrib_calves_out++) = p_starve_calf;                    // Starvation-induced mortality rate
        *(attrib_calves_out++) = calves[current_animal].starve;    // Incidence of starvation
        *(attrib_calves_out++) = calves[current_animal].mort;      // Age-dependent mortality
        *(attrib_calves_out++) = date_death_calf[current_animal];  // Date of death
        *(attrib_calves_out++) = p_surv_calf;                      // Survival probability
        
        *(attrib_nursing_out++) = 0.0L;    // Duration of lactation period
        *(attrib_nursing_out++) = 0.0L;    // Milk assimilation
        *(attrib_nursing_out++) = 0.0L;    // Milk provisioning
        *(attrib_nursing_out++) = 0.0L;    // Non-linearity between milk supply and BC of mother
        *(attrib_nursing_out++) = 0.0L;    // Age (days) at which milk consumption starts to decrease
        *(attrib_nursing_out++) = 0.0L;    // Non-linearity between milk assimilation and calf age
        *(attrib_nursing_out++) = 0.0L;    // Mass of mammary glands
        *(attrib_nursing_out++) = 0.0L;    // Mammary gland efficiency
        *(attrib_nursing_out++) = 0.0L;    // Milk yield
        *(attrib_nursing_out++) = 0.0L;    // Target body condition (mother)
        *(attrib_nursing_out++) = 0.0L;    // Nursing effort
        *(attrib_nursing_out++) = 0.0L;    // Proportion of lipids in milk
        *(attrib_nursing_out++) = 0.0L;    // Proportion of protein in milk
        
        *(attrib_costs_calves_out++) = 0.0L;   // Resting metabolic rate
        *(attrib_costs_calves_out++) = 0.0L;   // Locomotory costs
        *(attrib_costs_calves_out++) = 0.0L;   // Scalar on locomotory costs
        *(attrib_costs_calves_out++) = 0.0L;   // Energetic cost associated with growth
        *(attrib_costs_calves_out++) = 0.0L;   // Daily change in mass
        
        *(attrib_E_calves_out++) = 0.0L;   // Net energy (calves)
        *(attrib_E_calves_out++) = 0.0L;   // Energy gain (calves)
        *(attrib_E_calves_out++) = 0.0L;   // Energy expenditure (calves)
        *(attrib_E_calves_out++) = 0.0L;   // Change in blubber mass (calves)
        
        
      } else {
        
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
          *(attrib_out++) = animal->gape;                    // Size of the gape
          *(attrib_out++) = animal->mouth_ratio;             // Ratio of mouth width to total body length
          *(attrib_out++) = animal->mouth_angle;             // Upper angle of the mouth
          *(attrib_out++) = animal->mouth_width;             // Width of the mouth
          *(attrib_out++) = animal->gsl;                     // Gulf of St Lawrence
          *(attrib_out++) = animal->seus;                    // Southeastern U.S.
          *(attrib_out++) = animal->north;                   // Migration interruption
          *(attrib_out++) = animal->abort;                   // Abortion (pregnant females)
          *(attrib_out++) = p_starve;                        // Starvation-induced mortality rate
          *(attrib_out++) = animal->starve;                  // Incidence of starvation
          *(attrib_out++) = animal->mort;                    // Mortality from other sources
          *(attrib_out++) = date_death[current_animal];      // Date of death
          *(attrib_out++) = p_surv;                          // Survival probability
          
          // +++ Stressors +++
          
          if(stressors){
            
            *(attrib_stressors_out++) = gear_risk;             // Probability of becoming entangled
            *(attrib_stressors_out++) = entanglement_cost;                  // Entanglement cost
            *(attrib_stressors_out++) = animal->entangled(0);  // Is the animal entangled?
            *(attrib_stressors_out++) = animal->entangled(1);  // Entanglement severity
            *(attrib_stressors_out++) = animal->entangled(2);  // Does the entanglement involve the anterior part of the body?
            *(attrib_stressors_out++) = animal->entangled(3);  // Reduction in mouth gape
            *(attrib_stressors_out++) = animal->entangled(4);  // Entanglement duration
            *(attrib_stressors_out++) = animal->entangled(5);  // Day when entanglement event started
            *(attrib_stressors_out++) = animal->entangled(6);  // Day when entanglement event ended
            
            *(attrib_stressors_out++) = calves[current_animal].entangled(0);  // Is the animal entangled? (calves)
            *(attrib_stressors_out++) = calves[current_animal].entangled(1);  // Entanglement severity (calves)
            *(attrib_stressors_out++) = calves[current_animal].entangled(2);  // Does the entanglement involve the anterior part of the body? (calves)
            *(attrib_stressors_out++) = calves[current_animal].entangled(4);  // Entanglement duration (calves)
            *(attrib_stressors_out++) = calves[current_animal].entangled(5);  // Day when entanglement event started (calves)
            *(attrib_stressors_out++) = calves[current_animal].entangled(6);  // Day when entanglement event ended (calves)
            
            *(attrib_stressors_out++) = strike_risk;                           // Probability of being struck by vessels
            *(attrib_stressors_out++) = animal->strike;                        // Incidence of vessel strike
            *(attrib_stressors_out++) = calves[current_animal].strike;         // Incidence of vessel strike (calves)    
            
            *(attrib_stressors_out++) = animal->respnoise;                     // Incidence of behavioral response to noise exposure
            *(attrib_stressors_out++) = dB;                                    // Noise levels
            *(attrib_stressors_out++) = response_dB;                           // Behavioral response threshold
            *(attrib_stressors_out++) = piling_hrs;                            // Duration of pile-driving activities
            *(attrib_stressors_out++) = nursing_cessation(0);                  // Nursing cessation
            *(attrib_stressors_out++) = nursing_cessation(1);                  // Duration of nursing cessation
            
          }
          
          // +++ Energy budget +++
          
          *(attrib_E_out++) = E_in - E_out;                  // Net energy (adults and juveniles)
          *(attrib_E_out++) = E_in;                          // Energy gain (adults and juveniles)
          *(attrib_E_out++) = E_out;                         // Energy expenditure (adults and juveniles)
          *(attrib_E_out++) = delta_blubber;                 // Change in blubber mass (adults and juveniles)
          *(attrib_E_out++) = anabolism_efficiency[0];       // Lipid anabolism efficiency (adults and juveniles)
          *(attrib_E_out++) = catabolism_efficiency[1];      // Lipid catabolism efficiency (adults and juveniles)
          *(attrib_E_out++) = ED_lipids;                     // Energy density of lipids
          *(attrib_E_out++) = ED_protein;                    // Energy density of protein
          
          // +++ Energy intake +++
          
          *(attrib_feeding_out++) = is_feeding;                          // Incidence of foraging
          *(attrib_feeding_out++) = D_cop;                               // Prey concentration
          *(attrib_feeding_out++) = min_calanus;                         // Minimum prey concentration needed for foraging
          *(attrib_feeding_out++) = swim_speed;                          // Swimming speed while foraging
          *(attrib_feeding_out++) = capture_efficiency;                  // Capture efficiency
          *(attrib_feeding_out++) = animal->entangled(3);                // Reduction in mouth gape when entangled
          *(attrib_feeding_out++) = feeding_nursing_effort[0];           // Feeding/nursing effort
          *(attrib_feeding_out++) = eta(0);                              // Steepness of feeding/nursing effort curve for lower target BC
          *(attrib_feeding_out++) = eta(1);                              // Steepness of feeding/nursing effort curve for higher target BC
          *(attrib_feeding_out++) = target_bc[0];                        // Target body condition
          *(attrib_feeding_out++) = cop_mass;                            // Copepod mass (wet weight)
          *(attrib_feeding_out++) = water_content;                       // Water content of prey
          *(attrib_feeding_out++) = cop_kJ;                              // Energy density of Calanus finmarchicus stage V
          *(attrib_feeding_out++) = digestive_efficiency;                // Digestive efficiency
          *(attrib_feeding_out++) = metabolizing_efficiency(0);          // Metabolizing efficiency (calves and juveniles)
          *(attrib_feeding_out++) = metabolizing_efficiency(1);          // Metabolizing efficiency (adults)
          *(attrib_feeding_out++) = E_calanus;                           // Prey energy content
          *(attrib_feeding_out++) = calf_foraging;                       // Calf foraging
          
          // +++ Energy costs +++
          
          *(attrib_costs_out++) = resting_metab;           // Resting metabolic rate
          *(attrib_costs_out++) = LC_tot;                  // Locomotory costs
          *(attrib_costs_out++) = scalar_LC;               // Scalar on locomotory costs
          *(attrib_costs_out++) = stroke_rate;             // Stroke rate during non-foraging activities
          *(attrib_costs_out++) = stroke_rate_foraging;    // Stroke rate during foraging activities
          *(attrib_costs_out++) = E_growth;                // Energetic cost associated with growth
          *(attrib_costs_out++) = E_gest;                  // Energetic cost associated with gestation
          *(attrib_costs_out++) = fetal_growth_cost;       // Energetic cost associated with fetal growth
          *(attrib_costs_out++) = placental_cost;          // Energetic cost associated with placental maintenance
          *(attrib_costs_out++) = HIC;                     // Heat increment of gestation
          *(attrib_costs_out++) = E_lac;                   // Energetic cost associated with lactation
          *(attrib_costs_out++) = lean_mass_increment;     // Daily change in mass of adult
          *(attrib_costs_out++) = prop_muscle;             // Relative % of muscles in lean mass
          *(attrib_costs_out++) = prop_viscera;            // Relative % of viscera in lean mass
          *(attrib_costs_out++) = prop_bones;              // Relative % of bones in lean mass
          
          // +++ Activity budgets +++
          
          *(attrib_activity_out++) = travel_dist;               // Distance traveled
          *(attrib_activity_out++) = swim_speed;                // Swimming speed
          *(attrib_activity_out++) = percent_glide;            // % Time spent gliding during non-foraging activities
          *(attrib_activity_out++) = percent_glide_foraging;    // % Time spent gliding during foraging activities
          *(attrib_activity_out++) = percent_glide_echelon;     // % Time spent gliding during echelon swimming
          *(attrib_activity_out++) = travel_time;               // Time spent traveling
          *(attrib_activity_out++) = feeding_time;              // Time spent feeding
          *(attrib_activity_out++) = resting_nursing_time;      // Time spent resting + nursing
          
          // +++ Fetus development +++
          
          if(cohortID == 4 | cohortID == 5){
            
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
            *(attrib_fetus_out++) = prop_lipids_in_muscle;                 // Lipid content in fetal muscles
            *(attrib_fetus_out++) = prop_protein_in_muscle;                // Protein content in fetal muscles
            *(attrib_fetus_out++) = prop_lipids_in_viscera;                // Lipid content in fetal viscera
            *(attrib_fetus_out++) = prop_protein_in_viscera;               // Protein content in fetal viscera
            *(attrib_fetus_out++) = prop_lipids_in_bones;                  // Lipid content in fetal bones
            *(attrib_fetus_out++) = prop_protein_in_bones;                 // Protein content in fetal bones
            *(attrib_fetus_out++) = prop_lipids_in_blubber;                // Lipid content in fetal blubber
            *(attrib_fetus_out++) = prop_protein_in_blubber;               // Protein content in fetal blubber
            *(attrib_fetus_out++) = P_muscle;          // % of body volume comprising muscle in fetus
            *(attrib_fetus_out++) = P_viscera;          // % of body volume comprising viscera in fetus
            *(attrib_fetus_out++) = P_bones;          // % of body volume comprising bones in fetus
            *(attrib_fetus_out++) = blubber_density;          // Density of blubber
            *(attrib_fetus_out++) = muscle_density;          // Density of muscles
            *(attrib_fetus_out++) = visceral_density;          // Density of viscera
            *(attrib_fetus_out++) = bone_density;          // Density of bones
            
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
          *(attrib_out++) = 0.0L;              // Size of the gape
          *(attrib_out++) = 0.0L;              // Ratio of mouth width to total body length
          *(attrib_out++) = 0.0L;              // Upper angle of the mouth
          *(attrib_out++) = 0.0L;              // Width of the mouth
          *(attrib_out++) = animal->gsl;       // Gulf of St Lawrence
          *(attrib_out++) = animal->seus;      // Southeastern U.S.
          *(attrib_out++) = animal->north;     // Migration interruption
          *(attrib_out++) = animal->abort;     // Abortion (pregnant females)
          *(attrib_out++) = p_starve;          // Starvation-induced mortality rate
          *(attrib_out++) = animal->starve;    // Incide of starvation
          *(attrib_out++) = animal->mort;      // Mortality from other sources
          *(attrib_out++) = date_death[current_animal];       // Date of death
          *(attrib_out++) = p_surv;                          // Survival probability
          
          // +++ Stressors +++
          
          *(attrib_stressors_out++) = 0.0L;                     // Probability of becoming entangled
          *(attrib_stressors_out++) = 0.0L;                  // Entanglement cost
          *(attrib_stressors_out++) = animal->entangled(0);  // Is the animal entangled?
          *(attrib_stressors_out++) = animal->entangled(1);  // Entanglement severity
          *(attrib_stressors_out++) = animal->entangled(2);  // Does the entanglement involve the anterior part of the body?
          *(attrib_stressors_out++) = animal->entangled(3);  // Reduction in mouth gape
          *(attrib_stressors_out++) = animal->entangled(4);  // Entanglement duration
          *(attrib_stressors_out++) = animal->entangled(5);  // Day when entanglement event started
          *(attrib_stressors_out++) = animal->entangled(6);  // Day when entanglement event ended
          
          *(attrib_stressors_out++) = calves[current_animal].entangled(0);  // Is the animal entangled? (calves)
          *(attrib_stressors_out++) = calves[current_animal].entangled(1);  // Entanglement severity (calves)
          *(attrib_stressors_out++) = calves[current_animal].entangled(2);  // Does the entanglement involve the anterior part of the body? (calves)
          *(attrib_stressors_out++) = calves[current_animal].entangled(4);  // Entanglement duration (calves)
          *(attrib_stressors_out++) = calves[current_animal].entangled(5);  // Day when entanglement event started (calves)
          *(attrib_stressors_out++) = calves[current_animal].entangled(6);  // Day when entanglement event ended (calves)
          
          *(attrib_stressors_out++) = 0.0L;                                  // Probability of being struck by vessels
          *(attrib_stressors_out++) = animal->strike;                        // Incidence of vessel strike
          *(attrib_stressors_out++) = calves[current_animal].strike;         // Incidence of vessel strike (calves)    
          
          *(attrib_stressors_out++) = animal->respnoise;        // Incidence of behavioral response to noise exposure
          *(attrib_stressors_out++) = 0.0L;                     // Noise levels
          *(attrib_stressors_out++) = 0.0L;                     // Behavioral response threshold
          *(attrib_stressors_out++) = 0.0L;                     // Duration of pile-driving activities
          *(attrib_stressors_out++) = 0.0L;                     // Nursing cessation
          *(attrib_stressors_out++) = 0.0L;                     // Duration of nursing cessation
          
          // +++ Energy budget +++
          
          *(attrib_E_out++) = 0.0L;            // Net energy (adults and juveniles)
          *(attrib_E_out++) = 0.0L;            // Energy gain (adults and juveniles)
          *(attrib_E_out++) = 0.0L;            // Energy expenditure (adults and juveniles)
          *(attrib_E_out++) = 0.0L;            // Change in blubber mass (adults and juveniles)
          *(attrib_E_out++) = 0.0L;            // Lipid anabolism efficiency (adults and juveniles)
          *(attrib_E_out++) = 0.0L;            // Lipid catabolism efficiency (adults and juveniles)      
          *(attrib_E_out++) = 0.0L;            // Energy density of lipids
          *(attrib_E_out++) = 0.0L;            // Energy density of protein
          
          // +++ Energy intake +++
          
          *(attrib_feeding_out++) = 0.0L;        // Incidence of foraging
          *(attrib_feeding_out++) = 0.0L;        // Prey concentration
          *(attrib_feeding_out++) = 0.0L;        // Minimum prey concentration needed for foraging
          *(attrib_feeding_out++) = 0.0L;        // Swimming speed while foraging
          *(attrib_feeding_out++) = 0.0L;        // Capture efficiency
          *(attrib_feeding_out++) = 0.0L;        // % reduction in the gape when entangled
          *(attrib_feeding_out++) = 0.0L;        // Feeding/nursing effort
          *(attrib_feeding_out++) = 0.0L;        // Steepness of feeding/nursing effort curve for lower target BC
          *(attrib_feeding_out++) = 0.0L;        // Steepness of feeding/nursing effort curve for higher target BC
          *(attrib_feeding_out++) = 0.0L;        // Target body condition
          *(attrib_feeding_out++) = 0.0L;        // Copepod mass (wet weight)
          *(attrib_feeding_out++) = 0.0L;        // Water content of prey
          *(attrib_feeding_out++) = 0.0L;        // Energy density of Calanus finmarchicus stage V
          *(attrib_feeding_out++) = 0.0L;        // Digestive efficiency
          *(attrib_feeding_out++) = 0.0L;        // Metabolizing efficiency (calves and juveniles)
          *(attrib_feeding_out++) = 0.0L;        // Metabolizing efficiency (adults)
          *(attrib_feeding_out++) = 0.0L;        // Prey energy content
          *(attrib_feeding_out++) = 0.0L;        // Calf foraging
          
          // +++ Energy costs +++
          
          *(attrib_costs_out++) = 0.0L;           // Resting metabolic rate
          *(attrib_costs_out++) = 0.0L;           // Locomotory costs
          *(attrib_costs_out++) = 0.0L;           // Scalar on locomotory costs
          *(attrib_costs_out++) = 0.0L;           // Stroke rate during non-foraging activities
          *(attrib_costs_out++) = 0.0L;           // Stroke rate during foraging activities
          *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with growth
          *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with gestation
          *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with fetal growth
          *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with placental maintenance
          *(attrib_costs_out++) = 0.0L;           // Heat increment of gestation
          *(attrib_costs_out++) = 0.0L;           // Energetic cost associated with lactation
          *(attrib_costs_out++) = 0.0L;           // Daily change in mass
          *(attrib_costs_out++) = 0.0L;           // Relative % of muscles in lean mass
          *(attrib_costs_out++) = 0.0L;           // Relative % of viscera in lean mass
          *(attrib_costs_out++) = 0.0L;           // Relative % of bones in lean mass
          
          // +++ Activity budgets +++
          
          *(attrib_activity_out++) = 0.0L;         // Distance traveled
          *(attrib_activity_out++) = 0.0L;         // Swimming speed
          *(attrib_activity_out++) = 0.0L;        // Time spent gliding during non-foraging activities
          *(attrib_activity_out++) = 0.0L;         // Time spent gliding during foraging activities
          *(attrib_activity_out++) = 0.0L;         // Time spent gliding during echelon swimming
          *(attrib_activity_out++) = 0.0L;         // Time spent traveling
          *(attrib_activity_out++) = 0.0L;         // Time spent feeding
          *(attrib_activity_out++) = 0.0L;         // Time spent resting + nursing
          
          // +++ Fetus development +++
          
          if(cohortID == 4 | cohortID == 5){
            
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
            *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal bones
            *(attrib_fetus_out++) = 0.0L;          // Lipid content in fetal blubber
            *(attrib_fetus_out++) = 0.0L;          // Protein content in fetal blubber
            *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising muscle in fetus
            *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising viscera in fetus
            *(attrib_fetus_out++) = 0.0L;          // % of body volume comprising bones in fetus
            *(attrib_fetus_out++) = 0.0L;          // Density of blubber
            *(attrib_fetus_out++) = 0.0L;          // Density of muscles
            *(attrib_fetus_out++) = 0.0L;          // Density of viscera
            *(attrib_fetus_out++) = 0.0L;          // Density of bones
            
          }
          
        } // End if alive
        
        // +++ Calf traits +++
        
        if(cohortID == 5 && calves[current_animal].alive == 1){
          
          *(attrib_calves_out++) = calves[current_animal].cohortID;            // Unique cohort identifier
          *(attrib_calves_out++) = born[current_animal];                       // Birth
          *(attrib_calves_out++) = prob_birth;                                 // Probability of Birth
          *(attrib_calves_out++) = dob[current_animal];                        // Date of birth
          *(attrib_calves_out++) = date_death_calf[current_animal];            // Date of death
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
          *(attrib_calves_out++) = calves[current_animal].gape;                // Size of the gape
          *(attrib_calves_out++) = calves[current_animal].mouth_ratio;         // Ratio of mouth width to total body length
          *(attrib_calves_out++) = calves[current_animal].mouth_angle;         // Upper angle of the mouth
          *(attrib_calves_out++) = calves[current_animal].mouth_width;         // Width of the mouth
          *(attrib_calves_out++) = p_starve_calf;                              // Starvation-induced mortality rate
          *(attrib_calves_out++) = calves[current_animal].starve;              // Incidence of starvation
          *(attrib_calves_out++) = calves[current_animal].mort;                // Mortality from other sources
          *(attrib_calves_out++) = date_death_calf[current_animal];            // Mortality date
          *(attrib_calves_out++) = p_surv_calf;                                // Survival probability
          
          *(attrib_nursing_out++) = T_lac;                                     // Duration of lactation period
          *(attrib_nursing_out++) = milk_assim;                                // Milk assimilation
          *(attrib_nursing_out++) = milk_provisioning;                         // Milk provisioning
          *(attrib_nursing_out++) = zeta;                                      // Non-linearity between milk supply and BC of mother
          *(attrib_nursing_out++) = milk_decrease;                             // Age (days) at which milk consumption starts to decrease
          *(attrib_nursing_out++) = eta_milk;                                  // Non-linearity between milk assimilation and calf age
          *(attrib_nursing_out++) = mammary_M;                                 // Mass of mammary glands
          *(attrib_nursing_out++) = mammary_efficiency;                        // Mammary gland efficiency
          *(attrib_nursing_out++) = milk_production_rate;                      // Milk yield
          *(attrib_nursing_out++) = target_bc[1];                              // Target body condition
          *(attrib_nursing_out++) = feeding_nursing_effort[1];                 // Nursing effort
          *(attrib_nursing_out++) = milk_lipids;                               // Proportion of lipids in milk
          *(attrib_nursing_out++) = milk_protein;                              // Proportion of protein in milk
          
          *(attrib_costs_calves_out++) = resting_metab_calves;           // Resting metabolic rate
          *(attrib_costs_calves_out++) = LC_tot_calves;                  // Locomotory costs
          *(attrib_costs_calves_out++) = scalar_LC_calf;                 // Scalar on locomotory costs
          *(attrib_costs_calves_out++) = E_growth_calves;                // Energetic cost associated with growth
          *(attrib_costs_calves_out++) = lean_mass_increment_calves;     // Daily change in mass
          
          *(attrib_E_calves_out++) = E_in_calves - E_out_calves;    // Net energy (calves)
          *(attrib_E_calves_out++) = E_in_calves;                   // Energy gain (calves)
          *(attrib_E_calves_out++) = E_out_calves;                  // Energy expenditure (calves)
          *(attrib_E_calves_out++) = delta_blubber_calves;          // Change in blubber mass (calves)
          
        } else {
          
          *(attrib_calves_out++) = calves[current_animal].cohortID ; // Unique cohort identifier
          *(attrib_calves_out++) = born[current_animal];             // Birth
          *(attrib_calves_out++) = prob_birth;                       // Probability of Birth
          *(attrib_calves_out++) = dob[current_animal];              // Date of birth
          *(attrib_calves_out++) = date_death_calf[current_animal];  // Date of death
          *(attrib_calves_out++) = 0.0L;                             // Alive or dead
          *(attrib_calves_out++) = 0.0L;                             // Age
          *(attrib_calves_out++) = 0.0L;                             // Body condition
          *(attrib_calves_out++) = 0.0L;                             // Body length
          *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;                             // Coefficient of the Gompertz growth curve
          *(attrib_calves_out++) = 0.0L;                             // Total body mass (including blubber)
          *(attrib_calves_out++) = 0.0L;                             // Lean body mass (excluding blubber)
          *(attrib_calves_out++) = 0.0L;                             // Blubber mass
          *(attrib_calves_out++) = 0.0L;                             // Intercept of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = 0.0L;                             // Slope of the logarithmic mass-at-length relationship
          *(attrib_calves_out++) = 0.0L;                             // Size of the gape
          *(attrib_calves_out++) = 0.0L;                             // Ratio of mouth width to total body length
          *(attrib_calves_out++) = 0.0L;                             // Upper angle of the mouth
          *(attrib_calves_out++) = 0.0L;                             // Width of the mouth
          *(attrib_calves_out++) = p_starve_calf;                    // Starvation-induced mortality rate
          *(attrib_calves_out++) = calves[current_animal].starve;    // Incidence of starvation
          *(attrib_calves_out++) = calves[current_animal].mort;      // Age-dependent mortality
          *(attrib_calves_out++) = date_death_calf[current_animal];  // Date of death
          *(attrib_calves_out++) = p_surv_calf;                      // Survival probability
          
          *(attrib_nursing_out++) = 0.0L;   // Duration of lactation period
          *(attrib_nursing_out++) = 0.0L;   // Milk assimilation
          *(attrib_nursing_out++) = 0.0L;   // Milk provisioning
          *(attrib_nursing_out++) = 0.0L;   // Non-linearity between milk supply and BC of mother
          *(attrib_nursing_out++) = 0.0L;   // Age (days) at which milk consumption starts to decrease
          *(attrib_nursing_out++) = 0.0L;   // Non-linearity between milk assimilation and calf age
          *(attrib_nursing_out++) = 0.0L;   // Mass of mammary glands
          *(attrib_nursing_out++) = 0.0L;   // Mammary gland efficiency
          *(attrib_nursing_out++) = 0.0L;   // Milk yield
          *(attrib_nursing_out++) = 0.0L;   // Target body condition (mother)
          *(attrib_nursing_out++) = 0.0L;   // Nursing effort
          *(attrib_nursing_out++) = 0.0L;   // Proportion of lipids in milk
          *(attrib_nursing_out++) = 0.0L;   // Proportion of protein in milk
          
          *(attrib_costs_calves_out++) = 0.0L;  // Resting metabolic rate
          *(attrib_costs_calves_out++) = 0.0L;  // Locomotory costs
          *(attrib_costs_calves_out++) = 0.0L;  // Scalar on locomotory costs
          *(attrib_costs_calves_out++) = 0.0L;  // Energetic cost associated with growth
          *(attrib_costs_calves_out++) = 0.0L;  // Daily change in mass
          
          *(attrib_E_calves_out++) = 0.0L;  // Net energy (calves)
          *(attrib_E_calves_out++) = 0.0L;  // Energy gain (calves)
          *(attrib_E_calves_out++) = 0.0L;  // Energy expenditure (calves)
          *(attrib_E_calves_out++) = 0.0L;  // Change in blubber mass (calves)
          
        } 
        
      } // End if current_day == 0
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

//      std::vector<std::size_t> latentDensitySeq,
      
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
//' @param stepsize Radius of the proposal circle for movement (Michelot, 2019)
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
     std::vector<Eigen::MatrixXd> densities_gom,
     std::vector<std::size_t> densitySeq,
     std::vector<Eigen::MatrixXd> prey,
     std::vector<Eigen::MatrixXd> fishing,
     std::vector<Eigen::MatrixXd> vessels,
     std::vector<Eigen::MatrixXd> noise,
     Rcpp::NumericMatrix doseresp_seed,
     Rcpp::NumericVector doseresp,
     Eigen::MatrixXd regions,
     Eigen::VectorXd limits, 
     Eigen::VectorXd limits_regions,
     Eigen::VectorXd limits_prey,
     Eigen::VectorXd limits_fishing,
     Eigen::VectorXd limits_vessels,
     Eigen::VectorXd limits_noise,
     Eigen::VectorXd resolution,
     Eigen::VectorXd resolution_regions,
     Eigen::VectorXd resolution_prey,
     Eigen::VectorXd resolution_fishing,
     Eigen::VectorXd resolution_vessels,
     Eigen::VectorXd resolution_noise,
     std::size_t M,
     Eigen::MatrixXd xinit, 
     Eigen::MatrixXd yinit,
     bool stressors,
     bool growth,
     double prey_scale,
     Rcpp::NumericVector starvation_coefs,
     double starvation_death,
     double starvation_onset,
     Rcpp::NumericVector nursing_cessation,
     double piling_hrs,
     bool progress
 ) {
   
   // ------------------------------------------------------------------------
   // 0-index variables
   // ------------------------------------------------------------------------
   
   // Coming from R, we assume latentDensitySeq is 1-indexed
   // .begin = Returns an iterator pointing to the first element in the sequence
   // .end = Returns an iterator pointing to the last element in the sequence
   // for(auto d = latentDensitySeq.begin(); d != latentDensitySeq.end(); ++d)
   //   --(*d);
   
   // Ditto for densitySeq
   for(auto d = densitySeq.begin(); d != densitySeq.end(); ++d)
     --(*d);
   
   // ------------------------------------------------------------------------
   // Initialize environments for latent animals
   // ------------------------------------------------------------------------
   // emplace_back: inserts a new element at the end of a vector
   
   // std::vector<Environment> latent_environments;
   // auto dend  = latentDensitySeq.end();
   // 
   // for(auto d = latentDensitySeq.begin(); d != dend; ++d)
   //   latent_environments.emplace_back(densities[*d], densities_seus[*d], densities_gsl[*d],
   //                                    densities_gom[*d],
   //                                    prey[*d], fishing[*d], vessels[*d], noise[*d],
   //                                    daylight[*d], regions,
   //                                    limits, limits_daylight, limits_regions,
   //                                    limits_prey, limits_fishing,
   //                                    limits_vessels, limits_noise,
   //                                    resolution, resolution_daylight, resolution_regions,
   //                                    resolution_prey, resolution_fishing,
   //                                    resolution_vessels, resolution_noise,
   //                                    *d);
   
   std::vector<Environment> latent_environments;
   auto dend = densitySeq.end();
   
   for(auto d = densitySeq.begin(); d != dend; ++d){
     auto day_of_year = d - densitySeq.begin();
     // auto month_of_year = densitySeq[day_of_year];

     latent_environments.emplace_back(densities[*d],
                                      densities_seus[*d],
                                      densities_gsl[*d],
                                      densities_gom[*d],
                                      prey[*d],
                                      fishing[*d],
                                      vessels[*d],
                                      noise[day_of_year],
                                      regions,
                                      limits,
                                      limits_regions,
                                      limits_prey,
                                      limits_fishing,
                                      limits_vessels,
                                      limits_noise,
                                      resolution,
                                      resolution_regions,
                                      resolution_prey,
                                      resolution_fishing,
                                      resolution_vessels,
                                      resolution_noise,
                                      *d);
   }
   
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
   LatentMvmt rm(M);
   // LatentMvmt rm(M, stepsize); // Availability radius model
   
   // Initialize observed movement rule
   CoupledRandomWalk<
     CouplingAnimal, LatentAttractor, LatentMvmt, std::vector<Environment>
   > crm(rm, latent_environments, M);
   // > crm(rm, latent_environments, M, stepsize); // Availability radius model
   
   // ------------------------------------------------------------------------
   // Run simulations
   // ------------------------------------------------------------------------
   
   return movesim(cohortID, animals, calves, environments, 
                  latent_environments, doseresp_seed, doseresp, densitySeq, crm, support, 
                  limits, limits_regions, 
                  limits_prey, limits_fishing,
                  limits_vessels, limits_noise,
                  resolution, resolution_regions,
                  resolution_prey, resolution_fishing,
                  resolution_vessels, resolution_noise, 
                  // stepsize,
                  stressors, growth, prey_scale, 
                  starvation_coefs, starvation_death, starvation_onset,
                  nursing_cessation, piling_hrs, progress);
   
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
                       Eigen::MatrixXd density_gom, 
                       Eigen::MatrixXd prey, 
                       Eigen::MatrixXd fishing, 
                       Eigen::MatrixXd vessels, 
                       Eigen::MatrixXd noise, 
                       Eigen::MatrixXd regions, 
                       Eigen::VectorXd limits, 
                       Eigen::VectorXd limits_regions,
                       Eigen::VectorXd limits_prey,
                       Eigen::VectorXd limits_fishing,
                       Eigen::VectorXd limits_vessels,
                       Eigen::VectorXd limits_noise,
                       Eigen::VectorXd resolution,
                       Eigen::VectorXd resolution_regions,
                       Eigen::VectorXd resolution_prey,
                       Eigen::VectorXd resolution_fishing,
                       Eigen::VectorXd resolution_vessels,
                       Eigen::VectorXd resolution_noise,
                       double x, double y, char layer) {
  Environment e(density, density_seus, density_gsl, density_gom, prey, fishing, vessels, noise, regions, 
                limits, limits_regions, limits_prey, limits_fishing, 
                limits_vessels, limits_noise, resolution, resolution_regions, 
                resolution_prey, resolution_fishing,
                resolution_vessels, resolution_noise, 0);
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