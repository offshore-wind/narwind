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
 * Behaviors include: traveling (0), resting (1), nursing (2), feeding (3).
 */
struct Animal {
    double x, y;
    int behavior;
    // Member initializer list
    Animal() : x(0), y(0), behavior(0) { }
    Animal(const double & x0, const double & y0) : x(x0), y(y0) { }
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
     CouplingAnimal() { }
     CouplingAnimal(
         const std::vector<Animal> & init_latent,
         std::size_t init_animal
     ) : Animal(init_latent[init_animal].x, init_latent[init_animal].y),
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
template<typename AnimalType, 
         typename EnvironmentContainer,
         typename MovementRule, 
         int CoordSize = 2>
Rcpp::NumericVector movesim(
    std::vector<AnimalType> & animals,
    EnvironmentContainer & environments,
    MovementRule & m,
    Eigen::MatrixXd support,
    Eigen::VectorXd limits, 
    Eigen::VectorXd resolution,
    double stepsize,
    bool progress) {

    // Number of animals and time points to simulate
    std::size_t n = animals.size();
    std::size_t t = environments.size();
    
    // Initialize and label output
    Rcpp::NumericVector out(n*t*CoordSize);
    out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
    out.attr("dimnames") = Rcpp::List::create(
        Rcpp::CharacterVector::create("easting", "northing"), R_NilValue,
        R_NilValue
    );

    // Loop over time points and environments: One environment per time point
    Rcpp::NumericVector::iterator data_out = out.begin();
    
    auto env_end = environments.end();
    auto animals_end = animals.end();
    
    // Initialize progress bar
    progressbar bar(t);
    
    // 365 steps per animal
    for(auto env = environments.begin(); env != env_end; ++env) {
      
      if(progress) bar.update();
      
        for(auto animal = animals.begin(); animal != animals_end; ++animal) {
          
            // Save current coordinates
            *(data_out++) = animal->x;
            *(data_out++) = animal->y;
            
            // Update coordinates
            m.update(*animal, *env);
        }
    }

    return out;
}

template<typename AnimalType, 
         typename EnvironmentContainer,
         typename MovementRule, 
         int CoordSize = 2>
Rcpp::NumericVector movesim_geo(
    std::vector<AnimalType> & animals,
    EnvironmentContainer & environments,
    std::vector<PreyEnv> & prey,
    std::vector<std::size_t> densitySeq,
    MovementRule & m,
    Eigen::MatrixXd support,
    Eigen::VectorXd limits, 
    Eigen::VectorXd resolution,
    double stepsize,
    bool progress) {
  
  // Number of animals and time points to simulate
  std::size_t n = animals.size();
  std::size_t t = environments.size();
  
  // Initialize and label output
  Rcpp::NumericVector out(n*t*CoordSize);
  out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
  out.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create("easting", "northing"), R_NilValue,
    R_NilValue
  );
  
  // Loop over time points and environments: One environment per time point
  Rcpp::NumericVector::iterator data_out = out.begin();
  
  auto env_end = environments.end();
  auto animals_end = animals.end();

  progressbar bar(t);
  int current_month;
  
  // 365 steps per animal
  for(auto env = environments.begin(); env != env_end; ++env) {
    
    if(progress) bar.update();
    // env is a vector iterator, so env - environments.begin() returns the associated index
    current_month = densitySeq[env - environments.begin()];
    
    for(auto animal = animals.begin(); animal != animals_end; ++animal) {

      // Save current coordinates
      *(data_out++) = animal->x;
      *(data_out++) = animal->y;
      // *(data_out++) = 99;
      
      double x_0 = animal->x;
      double y_0 = animal->y;
      
      bool calc_geo = false;
      double x_1, y_1;
      int gd;
      int *gdptr;
      gdptr = &gd;

      // std::cout << current_month << std::endl;
      // std::cout << prey[current_month](x_0, y_0) << std::endl;
      
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
         
         // D = abs(x_1-x_0) + abs(y_1-y_0);
         
         // Calculate geodesic distance
         *gdptr = geoD(support, x_0, y_0, x_1, y_1, limits, resolution);
         
         // std::cout << gd  << " points(" << x_0 << "," << y_0 << ") " << "points(" << x_1 << "," << y_1 << ")\n";
         // if(abs(gd - D) < 2){
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
  
  return out;
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

Rcpp::NumericVector simAnnualCoupled(
    bool geodesic,
    Eigen::MatrixXd support,
    std::vector<Eigen::MatrixXd> densities,
    std::vector<std::size_t> densitySeq,
    std::vector<std::size_t> latentDensitySeq,
    std::vector<Eigen::MatrixXd> prey,
    Eigen::VectorXd limits, 
    Eigen::VectorXd resolution,
    std::size_t M,
    double stepsize, 
    Eigen::MatrixXd xinit, 
    Eigen::MatrixXd yinit,
    bool progress
) {

    // 0-index latentDensitySeq since, coming from R, we assume it is 1-indexed
    // .begin = Returns an iterator pointing to the first element in the sequence
    // .end = Returns an iterator pointing to the last element in the sequence
    for(auto d = latentDensitySeq.begin(); d != latentDensitySeq.end(); ++d)
        --(*d);

    // 0-index densitySeq since, coming from R, we assume it is 1-indexed
    for(auto d = densitySeq.begin(); d != densitySeq.end(); ++d)
        --(*d);

    // Initialize environments for latent animals
    std::vector<Environment> latent_environments;
    auto dend  = latentDensitySeq.end();
    // emplace_back: inserts a new element at the end of a vector
    for(auto d = latentDensitySeq.begin(); d != dend; ++d)
      latent_environments.emplace_back(densities[*d], 
                                       limits, 
                                       resolution, 
                                       *d);
    
    // Initialize prey base environments
    std::vector<PreyEnv> prey_environments;
    auto pend  = latentDensitySeq.end();
    for(auto p = latentDensitySeq.begin(); p != pend; ++p)
      prey_environments.emplace_back(prey[*p], limits, resolution, *p);
    
    // Initialize environments with attraction to latent animals
    std::vector<LatentAttractor> environments;
    auto dattrend = densitySeq.end();
    for(auto d = densitySeq.begin(); d != dattrend; ++d)
        environments.emplace_back(*d);

    // Initialize population
    std::vector<CouplingAnimal> animals;
    double *xiter = xinit.data();
    double *yiter = yinit.data();
    for(std::size_t j = 0; j < xinit.cols(); ++j) {
        std::vector<Animal> latent_animals;
        for(std::size_t i = 0; i < xinit.rows(); ++i) {
            latent_animals.emplace_back(Animal(*(xiter++), *(yiter++)));
        }
        animals.push_back(CouplingAnimal(latent_animals, densitySeq[0]));
    }

    // Initialize latent movement rule
    // typedef gives a type a new name
    typedef ReweightedRandomWalk<Animal, Environment> LatentMvmt;
    LatentMvmt rm(M, stepsize);

    // Initialize observed movement rule
    CoupledRandomWalk<
        CouplingAnimal, LatentAttractor, LatentMvmt, std::vector<Environment>
    > crm(rm, latent_environments, M, stepsize);
    
    if (geodesic) {
      return movesim_geo(animals, environments, prey_environments, densitySeq, crm, support, limits, resolution, stepsize, progress);
    } else {
      return movesim(animals, environments, crm, support, limits, resolution, stepsize, progress);
    }

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
                       Eigen::VectorXd limits, 
                       Eigen::VectorXd resolution,
    double x, double y) {
  Environment e(density, limits, resolution, 0);
    return e(x, y);
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
// 

