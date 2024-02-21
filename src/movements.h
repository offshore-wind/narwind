#ifndef MOVEMENTS_H
#define MOVEMENTS_H

#include <RcppEigen.h>
#include "geodesic.h"
#include "bioenergetics.h"
#include <math.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

// Standard deviations of the Normal half-step function
double sigma_move(int r){
  
  // Standard deviation of Normal half-step density given by 
  // sigma = mean step length / sqrt(pi)
  
  // Scale parameter of corresponding Rayleigh distribution given by
  // lambda = sigma * sqrt(2)
  
  double sigma_travel = 45/std::sqrt(M_PI);
  double sigma_feed = 6.75/std::sqrt(M_PI); 
  double sigma_nurse = 18.5/std::sqrt(M_PI);
  
  // Travel
  if(r == 7 | r == 8){ // MIDA and SCOS
    return sigma_travel;
  } else if(r == 9){ // SEUS
    return sigma_nurse;
  } else { // BOF, CABOT, CCB, GOM, GSL, SNE
    return sigma_feed;
  }
}

/**
 * Locally Gibbs movement, following Michelot et al. (2019)
 *
 * The implementation is abstracted so that the proposal weights come from
 * an overloaded operator() function within EnvironmentType
 *
 * @tparam AnimalType
 * @tparam EnvironmentType
 */
template<typename AnimalType, typename EnvironmentType>
class ReweightedRandomWalk {
  
private:
  
  // double stepsize_sq;
  std::size_t m;
  Eigen::MatrixXd proposals;
  Eigen::VectorXd weights; // Column vector
  // Eigen::VectorXd distances; // Column vector
  
  // CDF sampling from weights, assuming they are already normalized
  // It is a type able to represent the size of any object in bytes: 
  // size_t is the type returned by the sizeof operator and is widely 
  // used in the standard library to represent sizes and counts.
  std::size_t sampleInd() {
    double u = R::runif(0,1);
    double q = 0;
    std::size_t k = 0;
    double *w = weights.data();
    while(q<u) {
      q += *(w++);
      ++k;
    }
    return k - 1;
  }
  
public:
  
  // ReweightedRandomWalk(std::size_t n_proposals, double r) : stepsize_sq(r*r),
  ReweightedRandomWalk(std::size_t n_proposals) :  m(n_proposals)
  {
    // Resize the containers
    proposals.resize(2, m);
    weights.resize(m);
  }
  
  void update(AnimalType & animal, 
              EnvironmentType & environment, 
              int region,
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
              Eigen::VectorXd resolution_noise){
    
    // double d, d_0, a;
    double x0, y0, xn, yn;
    bool sampling = true;
    
    // int n_interm = 100;
    
    while(sampling) {
        
        x0 = animal.x + R::rnorm(0, sigma_move(region));
        y0 = animal.y + R::rnorm(0, sigma_move(region));
        
      //   // Availability radius model of Michelot et al. (2019)
      //   // runif run here in scalar mode - generates a single value
      //   // The square root is needed to sample points uniformly within a circle using inverse transform sampling
      //   // This is because the desired probability density function for point generation grows linearly with the
      //   // circle's radius, i.e., PDF(x) = 2x for a unit circle.
      //   // The corresponding CDF is the integral of this function, i.e., x2
      //   // The inverse CDF is therefore sqrt()
      //   // https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
      //   
      //   d_0 = std::sqrt(R::runif(0,stepsize_sq));
      //   a = R::runif(-M_PI,M_PI);
      // 
      //   // Disc center
      //   x0 = animal.x + d_0 * std::cos(a);
      //   y0 = animal.y + d_0 * std::sin(a);
      
      while(environment(x0,y0,'D')==0){
          
          x0 = animal.x + R::rnorm(0, sigma_move(region));
          y0 = animal.y + R::rnorm(0, sigma_move(region));
          
        //   // Availability radius model of Michelot et al. (2019)
        //   // runif run here in scalar mode - generates a single value
        //   // The square root is needed to sample points uniformly within a circle using inverse transform sampling
        //   // This is because the desired probability density function for point generation grows linearly with the
        //   // circle's radius, i.e., PDF(x) = 2x for a unit circle.
        //   // The corresponding CDF is the integral of this function, i.e., x2
        //   // The inverse CDF is therefore sqrt()
        //   // https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
        //   
        //   d_0 = std::sqrt(R::runif(0,stepsize_sq));
        //   a = R::runif(-M_PI,M_PI);
        //   
        //   // Disc center
        //   x0 = animal.x + d_0 * std::cos(a);
        //   y0 = animal.y + d_0 * std::sin(a);
        
      }
      
      // Sample m points in the circle, and compute their weights
      // Use asterisks to create pointers
      double *prop = proposals.data();
      double *w = weights.data();
      
      Rcpp::NumericVector distances(m);
      Rcpp::NumericVector indices(m);
      
      for(std::size_t i=0; i<m; ++i) {
        
          xn = x0 + R::rnorm(0, sigma_move(region));
          yn = y0 + R::rnorm(0, sigma_move(region));
          
          // // Availability radius model of Michelot et al. (2019)
          // 
          // d = std::sqrt(R::runif(0,stepsize_sq));
          // a = R::runif(-M_PI,M_PI);
          // 
          // // Calculate new x,y
          // xn = x0 + d * std::cos(a);
          // yn = y0 + d * std::sin(a);

        
        // Store x,y proposal
        *(prop++) = xn;
        *(prop++) = yn;
        
        // Retrieve density value (weight) at new x,y
        char lyr;
        if(animal.seus == 1 & animal.gsl == 0){
          lyr = 'S';
        } else if(animal.seus == 0 & animal.gsl == 1){
          lyr = 'G';
        } else {
          lyr = 'D';
        }
        

        *(w++) = environment(xn, yn, lyr);
        distances(i) = std::pow(xn - animal.x0, 2) + std::pow(yn - animal.y0, 2);
        
      }
      
      double d_min = Rcpp::which_min(distances);

      // Select one point
      double totalWeight = weights.sum();
      Rcpp::checkUserInterrupt();
      
      if(totalWeight > 0) {
        
        // weights /= totalWeight is equivalent to weights = weights / totalWeight
        weights /= totalWeight;
        std::size_t k;
        
        if(animal.north == 1){
          k = d_min;
        } else {
          k = sampleInd();
        }

        animal.x = proposals(0,k);
        animal.y = proposals(1,k);
        
        sampling = false;
        
      } // End if(totalWeight)
    } // End while(sampling)
    
  } // End update
};

/**
 * Animal movement with individual movement coupled to latent animals.
 *
 * Each latent animal moves according to a single latent environment and
 * movement type.
 *
 * @tparam AnimalType
 * @tparam EnvironmentType
 * @tparam LatentMovementType
 */
template<typename AnimalType, 
         typename EnvironmentType,
         typename LatentMovementType, 
         typename LatentEnvironmentContainer>

class CoupledRandomWalk {
  
private:
  
  // double stepsize_sq;
  std::size_t m;
  
  LatentMovementType & latent_mvmt;
  LatentEnvironmentContainer & latent_envs;
  
public:
  
  /**
   * @param base_mvmt Movement rule for updating latent animals
   * @param base_environments Environments for updating latent animals
   */
  CoupledRandomWalk(
    LatentMovementType & base_mvmt,
    LatentEnvironmentContainer & base_environments,
    std::size_t n_proposals
  ) : m(n_proposals), latent_mvmt(base_mvmt),
  latent_envs(base_environments) { }
  
  // ) : stepsize_sq(r*r), m(n_proposals), latent_mvmt(base_mvmt),
  
  void update(AnimalType & animal, EnvironmentType & environment, int region, 
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
              Eigen::VectorXd resolution_noise) {
    
    // Rcpp::NumericVector out (3);
    // Rcpp::NumericVector latent_coords (3);
    
    // Update latent animals
    auto latent_env = latent_envs.begin();
    auto latent_animal = animal.latent.begin();
    auto latent_end = animal.latent.end();
    
    while(latent_animal != latent_end) {
      latent_mvmt.update(*(latent_animal++), *(latent_env++), region, support, limits, 
                         limits_regions, limits_prey,
                         limits_fishing, limits_vessels, limits_noise,
                         resolution, resolution_regions, resolution_prey,
                         resolution_fishing, resolution_vessels, resolution_noise);
      // latent_coords = latent_mvmt.update(*(latent_animal++), *(latent_env++), TRUE, region);
      // animal.latent[latent_animal - animal.latent.begin()].x = latent_coords[0];
      // animal.latent[latent_animal - animal.latent.begin()].y = latent_coords[1];
    }
    
    // Environment defines latent animal toward which movement is attracted
    auto active_latent = animal.latent[environment.id];
    
    int from_region = latent_envs[environment.id](animal.x, animal.y, 'R');
    
    // If close enough, animal is coupled to the active latent animal
    double dist_to_latent;
    
    if(animal.north == 1){
      
      dist_to_latent = std::pow(animal.x - animal.x0, 2) + std::pow(animal.y - animal.y0, 2);
      
    } else {
      
      dist_to_latent = std::pow(animal.x - active_latent.x, 2) + std::pow(animal.y - active_latent.y, 2);
      
      if(dist_to_latent < 100) {
        
        animal.x = active_latent.x;
        animal.y = active_latent.y;
        
      } 
    }

    
    // Otherwise, select point within step size closest to latent animal
    // double d_0, d, a;
    double x0, y0, xn, yn, xnew, ynew;
    double dmin = std::numeric_limits<double>::infinity();
    bool sampling = true;
    
    while(sampling) {
      
        x0 = animal.x + R::rnorm(0, sigma_move(region));
        y0 = animal.y + R::rnorm(0, sigma_move(region)); 
        
        // Availability radius model of Michelot
      
        // d_0 = std::sqrt(R::runif(0,stepsize_sq));
        // a = R::runif(-M_PI,M_PI);
        // 
        // // Disc center
        // x0 = animal.x + d_0 * std::cos(a);
        // y0 = animal.y + d_0 * std::sin(a);

      while(latent_envs[environment.id](x0, y0, 'D')==0){

          x0 = animal.x + R::rnorm(0, sigma_move(region));
          y0 = animal.y + R::rnorm(0, sigma_move(region)); 

          // d_0 = std::sqrt(R::runif(0,stepsize_sq));
          // a = R::runif(-M_PI,M_PI);
          // 
          // // Disc center
          // x0 = animal.x + d_0 * std::cos(a);
          // y0 = animal.y + d_0 * std::sin(a);
          // 

      }
      
      // Sample m points in the circle, keep closest valid point
      for(std::size_t i=0; i<m; ++i) {
        
          xn = x0 + R::rnorm(0, sigma_move(region));
          yn = y0 + R::rnorm(0, sigma_move(region));
          
          // d = std::sqrt(R::runif(0,stepsize_sq));
          // a = R::runif(-M_PI,M_PI);
          // 
          // xn = x0 + d * std::cos(a);
          // yn = y0 + d * std::sin(a);
        
        
        char lyr;
        if(animal.seus == 1 & animal.gsl == 0){
          lyr = 'S';
        } else if(animal.seus == 0 & animal.gsl == 1){
          lyr = 'G';
        } else {
          lyr = 'D';
        }
        
        // Only consider proposals with non-zero mass
        if(latent_envs[environment.id](xn, yn, lyr) > 0) {
          
          bool calc_geo = false;
          
          // If in Canadian waters between U.S. border and entrance to GSL or
          // in Cape Cod Bay region
          if((x0 >= 685 & x0 <= 1420 & y0 >= 1070 & y0 <= 1865) ||
             (x0 >= 520 & x0 <= 730 & y0 >= 825 & y0 <= 960)) calc_geo = true;
          
          if(animal.north == 1){
            
            dist_to_latent = std::pow(xn - animal.x0, 2) + std::pow(yn - animal.y0, 2);
            
          } else {
            
            if(calc_geo){
              
              dist_to_latent = geoD(support, xn, yn, active_latent.x, active_latent.y, limits, resolution);
              
            } else {
              
              dist_to_latent = std::pow(xn - active_latent.x, 2) + std::pow(yn - active_latent.y, 2);
            }
          }
          
          bool across_land = false;
          
          // Keep track of the point that gets closest to latent
          if(dist_to_latent < dmin) {
            
            int to_region = latent_envs[environment.id](xn, yn, 'R');
          
          // Stop movements between GSL <--> SCOS <--> BOF
          if((from_region == 8 && to_region == 6) |
             (from_region == 6 && to_region == 8) |
             (from_region == 1 && to_region == 6) |
             (from_region == 6 && to_region == 1) |
             (from_region == 1 && to_region == 8) |
             (from_region == 8 && to_region == 1)) across_land = true;

          if(across_land == 0){
            // if(animal.day == 101) std::cout << latent_envs[environment.id](xn, yn, 'R') << animal.x << " " << animal.y << " " <<  dmin << " " << xn <<  " " << yn << std::endl;
            xnew = xn;
            ynew = yn;
            dmin = dist_to_latent;
          }
          }
          
        } // End conditional non-zero mass
      } // End loop over candidate locations
      
      // CheckUserInterrupt() checks if the ‘ctrl + c’ button was pressed
      // and stops execution if so
      Rcpp::checkUserInterrupt();
      sampling = !std::isfinite(dmin);
    }
    
    animal.x = xnew;
    animal.y = ynew;
    
  } // End update
  
};

#endif