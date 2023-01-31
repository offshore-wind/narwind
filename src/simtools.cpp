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
  float length;
  float mouth_ratio;
  float mouth_angle;
  float mouth_width;
  float mass;
  float fatmass;
  int gear;
  
  // Member initializer list - age initialized by default at 1 day old
  Animal() : x(0), y(0), cohortID(0), alive(1), age(0), bc(R::runif(0.3,0.5)),
  length(age2length(age)), mouth_ratio(start_mouth(cohortID, age)), mouth_angle(72), mouth_width(length*mouth_ratio),
  mass(length2mass(length,age)), fatmass(bc*mass), gear(0) { }
  
  // x: Easting
  // y: Northing
  // cohortID: Unique identifier for the cohort to which an animal belongs.
  // age: Age of the animal in years
  // bc: Body condition, expressed as a relative % of fat
  // length: Total body length in cm
  // mass: Total body mass in kg
  // fatmass: Total blubber mass in kg
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
    bc(start_percfat()),
    length(age2length(age)),
    mouth_ratio(start_mouth(cohort, age)), 
    mouth_angle(72), 
    mouth_width(length*mouth_ratio),
    mass(length2mass(length, age)),
    fatmass(bc*mass),
    gear(start_entangled()){ }
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
         int CoordSize = 2,
         int nparam_animal = 13,
         int nparam_feeding = 17,
         int nparam_nursing = 15,
         int nparam_costs = 3>

Rcpp::List movesim(
    int cohortID,
    std::vector<AnimalType> & animals,
    std::vector<AnimalType> & calves,
    EnvironmentContainer & environments,
    std::vector<Environment> & layers,
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
    Rcpp::CharacterVector::create("easting", "northing"), R_NilValue, R_NilValue
  );
  
  // Initialize and label array in which bioenergetic parameters will be stored
  Rcpp::NumericVector attrib(n*t*nparam_animal);
  attrib.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  attrib.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort", "alive", "age", "bc", "length", "mouth_r", "mouth_a", "mouth_w", "mass", "fatmass", "gear", "strike", "noise"),
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for dependent calves
  Rcpp::NumericVector attrib_calves(n*t*nparam_animal);
  attrib_calves.attr("dim") = Rcpp::Dimension(nparam_animal,n,t);
  attrib_calves.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("cohort", "alive", "age", "bc", "length", "mouth_r", "mouth_a", "mouth_w", "mass", "fatmass", "gear", "strike", "noise"), 
    R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_feeding(n*t*nparam_feeding);
  attrib_feeding.attr("dim") = Rcpp::Dimension(nparam_feeding,n,t);
  attrib_feeding.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_in", "feed", "preyconc", "minprey", "gape", "speed", "captEff", "impedance", "daylight",
                                  "tfeed", "feed_effort", "targetBC", "cop_mass", "cop_kJ", "digestEff", "metabEff", "E_cop"), 
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_nursing(n*t*nparam_nursing);
  attrib_nursing.attr("dim") = Rcpp::Dimension(nparam_nursing,n,t);
  attrib_nursing.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_in", "sep", "t_sep", "assim", "provision", "mamm_M", "milk_rate", "Dmilk", "t_suckling",
                                  "targetBC", "nursing", "milk_lip", "milk_pro", "ED_lip", "ED_pro"),
                                  R_NilValue, R_NilValue
  );
  
  // Create a separate array for variables associated with energy intake in adults/juveniles
  Rcpp::NumericVector attrib_costs(n*t*nparam_costs);
  attrib_costs.attr("dim") = Rcpp::Dimension(nparam_costs,n,t);
  attrib_costs.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("E_out", "p1", "p2"),
                                  R_NilValue, R_NilValue
  );
  
  // Loop over time points and environments: One environment per time point
  Rcpp::NumericVector::iterator data_out = out.begin();
  Rcpp::NumericVector::iterator attrib_out = attrib.begin();
  Rcpp::NumericVector::iterator attrib_calves_out = attrib_calves.begin();
  Rcpp::NumericVector::iterator attrib_feeding_out = attrib_feeding.begin();
  Rcpp::NumericVector::iterator attrib_nursing_out = attrib_nursing.begin();
  Rcpp::NumericVector::iterator attrib_costs_out = attrib_costs.begin();
  
  // std::advance(attrib_out, nparam); // Skip update of parameters at first time point as already initialized
  
  auto env_end = environments.end();
  auto animals_end = animals.end();
  
  // Inititalize progress bar
  progressbar bar(t);
  
  // std::cout << cohortID << std::endl;
  
  // double age, length, mtot, mfat, bc;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> binary_distribution(0,1);
  
  // 365 steps per animal
  for(auto env = environments.begin(); env != env_end; ++env) {
    
    int current_day = env - environments.begin();
    
    if(progress) bar.update();
    
    // env is a vector iterator, so env - environments.begin() returns the associated index
    int current_month = densitySeq[env - environments.begin()];
    
    // std::cout << current_month << " " << environments[env - environments.begin()].id << std::endl;
    
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {
      
      // Save current coordinates (fills by column)
      
      *(data_out++) = animal->x;
      *(data_out++) = animal->y;
      
      // std::cout << current_day << " (" << current_month << ")" << std::endl;
      
      int current_animal = animal - animals.begin();
      // std::cout << animals[current_animal].cohortID << std::endl;
      
      // ------------------------------------------------------------
      // Retrieve environment values (prey + stressors)
      // ------------------------------------------------------------
      
      double D_cop = layers[current_month](animal->x, animal->y, 'P');
      // double dB = layers[current_month](animal->x, animal->y, 'N');
      // double gear_risk = layers[current_month](animal->x, animal->y, 'F');
      // double strike_risk = layers[current_month](animal->x, animal->y, 'V');

      // ------------------------------------------------------------
      // ENERGY INTAKE - FORAGING
      // ------------------------------------------------------------
      
      // Behavioral response to pile-driving noise exposure -- (d.u.)
      // TODO: Update this line with the dose-response function obtained during the expert elicitation.
      int resp_noise = binary_distribution(generator);
      
      // Incidence of a vessel strike -- (d.u.)
      // TODO: Update this line with correct strike process -> leads to death
      int vessel_strike = binary_distribution(generator);
      
      // Minimum prey concentration (Calanus finmarchicus stage V) -- (copepods/cubic m)
      // Meta: From year 1987, by species
      double min_calanus = rtnorm(1785, 718.77, 0, INFINITY);
      
      // Determine whether the prey concentration is sufficient to support foraging -- (d.u.)
      bool is_feeding;
      if(resp_noise){
        is_feeding = 0;
      } else {
        is_feeding = feeding_threshold(min_calanus, D_cop);
      }

      // Mouth gape area -- (m2)
      double mouth_gape = gape_size(animal->length, animal->mouth_width, animal->mouth_angle);
      
      // Swim speed during foraging -- (m/s)
      // Meta: No temporal filter, all records
      double swim_speed_foraging = R::rnorm(0.89940, 0.18310);
      
      // Capture efficiency -- (%)
      // Meta: No temporal filter, all records
      // draw("beta", 0.91575, 0.0392)
      double capture_efficiency = R::rbeta(45.06241, 4.145791);
      
      // Percent reduction in the gape during an entanglement event
      // Meta: Parameter unknown
      double gape_reduction = 0.0;
      if(animal->gear == 1) gape_reduction = R::runif(0,0.5);
      
      // Total time available for foraging activities --(hr converted to s)
      double daylight_hours = layers[current_month](animal->x, animal->y, 'L') * 3600;
      
      // Time spent feeding per day -- (hr converted to s)
      // Meta: No temporal filter, by age class, fill NAs with "adult"
      double feeding_time = rtnorm(12.768, 2.9166, 0, daylight_hours) * 3600;
      
      // Foraging/nursing response to body condition -- (d.u.)
      // Parameters are estimated by:
      // (1) running the meta-analysis function on target relative blubber mass (see separate R script)
      // (2) drawing the resulting curve and recording the maximum value from the associated distribution
      // (3) using the optim_feeding function with bounds of c(0, max) to find the parameters of the logistic function
      Rcpp::NumericVector feeding_nursing_effort (2); // feeding, nursing
      Rcpp::NumericVector target_bc (2); // adjuv, calves
      
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
      
      // Average mass of the prey -- (mg to g)
      // Meta: From year 2000, all records
      double cop_mass = rtnorm(0.38224, 0.095553, 0, INFINITY)/1000;
      
      // Energy density of Calanus stage V -- (kJ / g)
      // Meta: No temporal filter, by species, combine all records for C. finmarchicus
      double cop_kJ = R::rnorm(24.054, 3.3191);
      
      // Digestive efficiency -- (%)
      // Meta: No temporal filter, all records
      double digestive_efficiency = R::rnorm(0.88504, 0.0045112);
      
      // Metabolizing efficiency (1 - heat increment of feeding) -- (%)
      // Meta: No temporal filter, by age_class, combining juveniles and calves
      double metabolizing_efficiency;
      
      if(animal->cohortID <= 2){ // Juveniles
        metabolizing_efficiency = rtnorm(0.740, 0.046, 0, 1);
      } else { // Adults
        metabolizing_efficiency = rtnorm(0.875, 0.036, 0, 1);
      }
      
      // Prey energy content
      double E_calanus = Econtent_cop(cop_mass, cop_kJ, digestive_efficiency, metabolizing_efficiency);
      
      // ENERGY INTAKE (kj to MJ)
      double E_in = (1 - resp_noise) * is_feeding * D_cop * mouth_gape * (1 - gape_reduction) *
        swim_speed_foraging * capture_efficiency * feeding_time * feeding_nursing_effort[0] *
        E_calanus / 1000;
      
      // std::cout << animal->x << ";" << animal->y << " - " << current_month << " - " << D_cop << " - " << min_calanus << " " << is_feeding << " " << E_in<< std::endl;
      
      // ------------------------------------------------------------
      // ENERGY INTAKE - SUCKLING
      // ------------------------------------------------------------
      
      // Incidence of mother-calf separation -- (d.u.)
      // Reported to be between 10 and 40% of sightings as per Hamilton et al. (2022) - take mean = 25%
      bool mumcalf_separation = R::rbinom(1, 0.25);
      double separation_duration = 0;

      // Separation duration -- (days)
      // Using half-normal with mean of 5.9 and max of 23, as per Hamilton et al. (2022)
      // Coded as truncated Normal centered on zero
      if(mumcalf_separation){
        separation_duration = rtnorm(0, 7.44375, 0, 23);
      }

      // Age at which milk consumption starts to decrease -- (days)
      // Fortune et al. (2020)
      double milk_decrease = 288;

      // Non-linearity between milk assimilation and calf age -- (d.u.)
      // Hin et al. (2019)
      double eta_milk = 0.9;

      // Duration of the lactation -- (days)
      int T_lac = 365;

      // Milk assimilation -- (d.u.)
      // Varies with calf age
      double milk_assim = milk_assimilation(current_day, T_lac, milk_decrease, eta_milk);

      // Starvation threshold expressed as relative blubber mass -- (d.u.)
      // As per Pirotta et al. (2018) - to test extreme conditions of leanness
      float starvation = 0.05;

      // Non-linearity between milk supply and the body condition of the mother
      int zeta = -2;

      // Milk provisioning -- (d.u.)
      // Varies with mother's body condition
      double milk_provisioning = milk_supply(starvation, target_bc[0], animal->mass, animal->fatmass, zeta);

      // Mammary gland efficiency -- (%)
      // As per Brody 1968
      float mammary_efficiency = 0.9;

      // Total mass of mammary glands -- (kg)
      double mammary_M = mammary_mass(animal->mass);

      // Milk production rate by the mother -- (cubic m per kg per sec)
      double milk_production_rate = milk_production(mammary_M);

      // Density of milk -- (kg/m3)
      // Meta: No temporal filter, all records
      double D_milk = R::rnorm(1.02, 0.01) * 1000; // L to cubic m

      // Time spent nursing/suckling per day -- (hours to sec)
      // Meta:: No filter on time spent nursing, time spent suckling during nursing bouts filtered by species (value for E. australis)
      double t_suckling = rtnorm(0, 3.688672, 0, 24) *3600;

      // Proportion of lipids in milk -- (%)
      // Meta: All records – as only value for closest relative (Eubalaena mysticetus) deemed low
      // Set maximum to 60% as reported in the literature (White 1953)
      double milk_lipids = rtnorm(0.365680, 0.11932, 0, 0.51);

      // Proportion of protein in milk -- (%)
      double milk_protein = rtnorm(0.13650, 0.031562, 0, 1);

      if((milk_lipids + milk_protein) > 1){
        std::cout << "Inconsistent milk composition" << std::endl;
        break;
      }

      // Energy density of lipids -- (kJ / kg)
      double ED_lipids = 39300;

      // Energy density of protein -- (kJ / kg)
      // mean(c(23600,18000))
      double ED_protein = 20800;

      // ENERGY INTAKE (kj to MJ)
      double E_in_calves = (1 - resp_noise) * (1 - mumcalf_separation) * milk_assim * milk_provisioning * mammary_efficiency *
        mammary_M * milk_production_rate * t_suckling * feeding_nursing_effort[1] * D_milk *
        (milk_lipids * ED_lipids + milk_protein * ED_protein) / 1000;
      
      // ------------------------------------------------------------
      // STORE VALUES in output matrices
      // ------------------------------------------------------------
      
      // +++ Animal traits +++

      *(attrib_out++) = animal->cohortID;
      *(attrib_out++) = animal->alive;
      *(attrib_out++) = animal->age;
      *(attrib_out++) = animal->bc;
      *(attrib_out++) = animal->length;
      *(attrib_out++) = animal->mouth_ratio;
      *(attrib_out++) = animal->mouth_angle;
      *(attrib_out++) = animal->mouth_width;
      *(attrib_out++) = animal->mass;
      *(attrib_out++) = animal->fatmass;
      *(attrib_out++) = animal->gear;
      *(attrib_out++) = vessel_strike;
      *(attrib_out++) = resp_noise;
      
      // +++ Energy intake +++
      
      *(attrib_feeding_out++) = E_in;
      *(attrib_feeding_out++) = is_feeding;
      *(attrib_feeding_out++) = D_cop;
      *(attrib_feeding_out++) = min_calanus;
      *(attrib_feeding_out++) = mouth_gape;
      *(attrib_feeding_out++) = swim_speed_foraging;
      *(attrib_feeding_out++) = capture_efficiency;
      *(attrib_feeding_out++) = gape_reduction;
      *(attrib_feeding_out++) = daylight_hours;
      *(attrib_feeding_out++) = feeding_time;
      *(attrib_feeding_out++) = feeding_nursing_effort[0];
      *(attrib_feeding_out++) = target_bc[0];
      *(attrib_feeding_out++) = cop_mass;
      *(attrib_feeding_out++) = cop_kJ;
      *(attrib_feeding_out++) = digestive_efficiency;
      *(attrib_feeding_out++) = metabolizing_efficiency;
      *(attrib_feeding_out++) = E_calanus;
      
      // +++ Calf traits +++
      
      if(cohortID == 5){
        
        *(attrib_calves_out++) = calves[current_animal].cohortID;
        *(attrib_calves_out++) = calves[current_animal].alive;
        *(attrib_calves_out++) = calves[current_animal].age;
        *(attrib_calves_out++) = calves[current_animal].bc;
        *(attrib_calves_out++) = calves[current_animal].length;
        *(attrib_calves_out++) = calves[current_animal].mouth_ratio;
        *(attrib_calves_out++) = calves[current_animal].mouth_angle;
        *(attrib_calves_out++) = calves[current_animal].mouth_width;
        *(attrib_calves_out++) = calves[current_animal].mass;
        *(attrib_calves_out++) = calves[current_animal].fatmass;
        *(attrib_calves_out++) = calves[current_animal].gear;
        *(attrib_calves_out++) = vessel_strike;
        *(attrib_calves_out++) = resp_noise;
        
        *(attrib_nursing_out++) = E_in_calves;
        *(attrib_nursing_out++) = mumcalf_separation;
        *(attrib_nursing_out++) = separation_duration;
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
      }
      
      // +++ Energy costs +++
      
      *(attrib_costs_out++) = 0;
      *(attrib_costs_out++) = 0;
      *(attrib_costs_out++) = 0;
      
      // // Save current coordinates (fills by column)
      // 
      // *(data_out++) = animal->x;
      // *(data_out++) = animal->y;
      
      double x_0 = animal->x;
      double y_0 = animal->y;
      
      bool calc_geo = false;
      double x_1, y_1;
      int gd;
      int *gdptr;
      gdptr = &gd;
      
      // std::cout << layers[current_month](x_0, y_0, 'P') << std::endl;
      
      // Only activate geodesic calculations in some areas
      if((x_0 >= 685 & x_0 <= 1420 & y_0 >= 1070 & y_0 <= 1865) ||
         (x_0 >= 520 & x_0 <= 730 & y_0 >= 825 & y_0 <= 960)) calc_geo = true;
      
      if(!calc_geo){
        
        // Update the animal's state
        m.update(*animal, *env);
        
      } else {
        
        while(calc_geo == true){
          
          // Update the animal's state
          m.update(*animal, *env);
          
          x_1 = animal->x;
          y_1 = animal->y;
          
          // Calculate geodesic distance
          *gdptr = geoD(support, x_0, y_0, x_1, y_1, limits, resolution);
          
          if(*gdptr >= 0 & *gdptr < stepsize){
            
            calc_geo = false;
            
          } else {
            animal->x = x_0;
            animal->y = y_0;
          }
          
        }
        
      }
      
    }
  }
  
  results["locs"] = out;
  results["attrib"] = attrib;
  results["attrib_calves"] = attrib_calves;
  results["in.kj"] = attrib_feeding;
  results["in.kj_calves"] = attrib_nursing;
  results["out.kj"] = attrib_costs;
  results["out.kj_calves"] = attrib_costs;
  
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
     std::vector<Eigen::MatrixXd> daylight,
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
     latent_environments.emplace_back(densities[*d], prey[*d], fishing[*d], vessels[*d], noise[*d], daylight[*d],
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
                  latent_environments, densitySeq, crm, support, 
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
                       Eigen::VectorXd limits, 
                       Eigen::VectorXd resolution,
                       double x, double y, char layer) {
  Environment e(density, prey, fishing, vessels, noise, daylight, limits, resolution, 0);
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