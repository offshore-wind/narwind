#ifndef MOVEMENTS_H
#define MOVEMENTS_H

#include <RcppEigen.h>
#include "geodesic.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

/**
 * Fixed-rate uniform random walk in R^2
 * @tparam AnimalType
 * @tparam EnvironmentType
 */
template<typename AnimalType, typename EnvironmentType>
struct RandomMovement {
  void update(AnimalType & a, EnvironmentType & e) {
    a.x += R::runif(-M_PI,M_PI);
    a.y += R::runif(-M_PI,M_PI);
  }
};

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
  
  double stepsize_sq;
  std::size_t m;
  Eigen::MatrixXd proposals;
  Eigen::VectorXd weights; // Column vector
  
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
  
  ReweightedRandomWalk(std::size_t n_proposals, double r) : stepsize_sq(r*r),
  m(n_proposals)
  {
    // Resize the containers
    proposals.resize(2, m);
    weights.resize(m);
  }
  
  void update(AnimalType & animal, EnvironmentType & environment) {
    double d, a, x0, y0, xn, yn;
    bool sampling = true;
    
    // Constraint on movements between regions
    // bool anticosti_1_from, anticosti_1_to = false;
    // bool anticosti_2_from, anticosti_2_to = false;
    // bool fundy_1_from, fundy_1_to = false;
    // bool fundy_2_from, fundy_2_to = false;
    // bool fundy_3_from, fundy_3_to = false;
    // bool fundy_4_from, fundy_4_to = false;
    // 
    // if(animal.x >= 1042 & animal.x <= 1198 & animal.y >= 1742 & animal.y <= 1782) anticosti_1_from = true;
    // if(animal.x >= 1030 & animal.x <= 1082 & animal.y >= 1742 & animal.y <= 1809) anticosti_2_from = true;
    // if(animal.x >= 1048 & animal.x <= 1162 & animal.y >= 1285 & animal.y <= 1345) fundy_1_from = true;
    // if(animal.x >= 1048 & animal.x <= 1162 & animal.y >= 1285 & animal.y <= 1345) fundy_2_from = true;
    // if(animal.x >= 1035 & animal.x <= 1085 & animal.y >= 1380 & animal.y <= 1422) fundy_3_from = true;
    // if(animal.x >= 1002 & animal.x <= 1175 & animal.y >= 1278 & animal.y <= 1380) fundy_4_from = true;
    
    
    while(sampling) {
      
      // runif run here in scalar mode - generates a single value
      d = std::sqrt(R::runif(0,stepsize_sq));
      a = R::runif(-M_PI,M_PI);
      
      // Disc center
      x0 = animal.x + d * std::cos(a);
      y0 = animal.y + d * std::sin(a);
      
      // Sample m points in the circle, and compute their weights
      // Use asterisks to create pointers
      double *prop = proposals.data();
      double *w = weights.data();
      
      for(std::size_t i=0; i<m; ++i) {
        
        // anticosti_1_to = false;
        // anticosti_2_to = false;
        // fundy_1_to = false;
        // fundy_2_to = false;
        // fundy_3_to = false;
        // fundy_4_to = false;
        
        // Randomly sample a distance and angle
        d = std::sqrt(R::runif(0,stepsize_sq));
        a = R::runif(-M_PI,M_PI);
        
        // Calculate new x,y
        xn = x0 + d * std::cos(a);
        yn = y0 + d * std::sin(a);
        
        // Store x,y proposal
        *(prop++) = xn;
        *(prop++) = yn;

        // if(xn >= 1042 & xn <= 1198 & yn >= 1782 & yn <= 1862) anticosti_1_to = true;
        // if(xn >= 1042 & xn <= 1160 & yn >= 1809 & yn <= 1862) anticosti_2_to = true;
        // if(xn >= 1002 & xn <= 1082 & yn >= 1342 & yn <= 1410) fundy_1_to = true;
        // if(xn >= 1082 & xn <= 1202 & yn >= 1182 & yn <= 1290) fundy_2_to = true;
        // if(xn >= 1042 & xn <= 1222 & yn >= 1380 & yn <= 1482) fundy_3_to = true;
        // if(xn >= 1042 & xn <= 1222 & yn >= 1380 & yn <= 1482) fundy_4_to = true;
        
        // if((anticosti_1_from && anticosti_1_to) || (anticosti_2_from && anticosti_2_to) ||
        //    (fundy_1_from && fundy_1_to) || (fundy_2_from && fundy_2_to) || (fundy_3_from && fundy_3_to) ||
        //    (fundy_4_from && fundy_4_to)){
        //   
        //   *(w++) = 0;
        //   
        // } else {
        //   
        //   *(w++) = environment(xn, yn);
        //   
        // }
        
        *(w++) = environment(xn, yn);
        
        // Constraint on movements between regions
        // Regions are:
        // --------------------------
        // [1] Lower Bay of Fundy
        // [2] Upper Bay of Fundy
        // [3] Cabot Strait
        // [4] Cape Cod Bay
        // [5] Gulf of Maine
        // [6] Gulf of St Lawrence
        // [7] Mid Atlantic
        // [8] Scotian Shelf
        // [9] Southeast U.S.
        // [10] Southern New England
        // int r_to = environment(xn, yn)[1];
        // 
        // if((r_from == 2 & r_to == 6) || (r_from == 6 & r_to == 2) ||
        //    (r_from == 8 & r_to == 6) || (r_from == 6 & r_to == 8) ||
        //    (r_from == 8 & r_to == 2) || (r_from == 2 & r_to == 8)){
        //   *(w++) = 0;
        // } else {
        //   *(w++) = environment(xn, yn)[0];
        // }
        
        // Retrieve density value (weight) at new x,y
        // *(w++) = environment(xn, yn);
        
      }
      
      // Select one point
      double totalWeight = weights.sum();
      Rcpp::checkUserInterrupt();
      
      if(totalWeight > 0) {
        
        // weights /= totalWeight is equivalent to weights = weights / totalWeight
        weights /= totalWeight;
        std::size_t k = sampleInd();
        
        animal.x = proposals(0,k);
        animal.y = proposals(1,k);
        
        sampling = false;
        
      }
    }
  }
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
template<typename AnimalType, typename EnvironmentType,
         typename LatentMovementType, typename LatentEnvironmentContainer>
class CoupledRandomWalk {
  
private:
  
  double stepsize_sq;
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
    std::size_t n_proposals,
    double r
  ) : stepsize_sq(r*r), m(n_proposals), latent_mvmt(base_mvmt),
  latent_envs(base_environments) { }
  
  void update(AnimalType & animal, EnvironmentType & environment) {
    
    // Update latent animals
    auto latent_env = latent_envs.begin();
    auto latent_animal = animal.latent.begin();
    auto latent_end = animal.latent.end();
    while(latent_animal != latent_end) {
      latent_mvmt.update(*(latent_animal++), *(latent_env++));
    }
    
    // Environment defines latent animal toward which movement is attracted
    auto active_latent = animal.latent[environment.id];
    
    // If close enough, animal is coupled to the active latent animal
    double dist_to_latent = std::pow(animal.x - active_latent.x, 2) +
      std::pow(animal.y - active_latent.y, 2);
    
    if(dist_to_latent < stepsize_sq) {
      animal.x = active_latent.x;
      animal.y = active_latent.y;
      return;
    }
    
    // Constraint on movements between regions
    // bool anticosti_1_from, anticosti_1_to = false;
    // bool anticosti_2_from, anticosti_2_to = false;
    // bool fundy_1_from, fundy_1_to = false;
    // bool fundy_2_from, fundy_2_to = false;
    // bool fundy_3_from, fundy_3_to = false;
    // bool fundy_4_from, fundy_4_to = false;
    // 
    // if(animal.x >= 1042 & animal.x <= 1198 & animal.y >= 1742 & animal.y <= 1782) anticosti_1_from = true;
    // if(animal.x >= 1030 & animal.x <= 1082 & animal.y >= 1742 & animal.y <= 1809) anticosti_2_from = true;
    // if(animal.x >= 1048 & animal.x <= 1162 & animal.y >= 1285 & animal.y <= 1345) fundy_1_from = true;
    // if(animal.x >= 1048 & animal.x <= 1162 & animal.y >= 1285 & animal.y <= 1345) fundy_2_from = true;
    // if(animal.x >= 1035 & animal.x <= 1085 & animal.y >= 1380 & animal.y <= 1422) fundy_3_from = true;
    // if(animal.x >= 1002 & animal.x <= 1175 & animal.y >= 1278 & animal.y <= 1380) fundy_4_from = true;
    
    // bool area_from = false;
    // bool area_to = false;
    
    // bool bof_from = false;
    // bool bof_to = false;
    // 
    // bool gsl_from = false;
    // bool gsl_to = false;
    // 
    // bool cod_from = false;
    // bool cod_to = false;
    // 
    // bool codsouth_from = false;
    // bool codsouth_to = false;
    
    // // Bay of Fundy
    // if(animal.x >= 1012.593 & animal.x <= 1176.901 & animal.y >= 1278.953 & animal.y <= 1389.127) bof_from = true;
    // 
    // // Gulf of St Lawrence (west and east)
    // if(animal.x >= 979.2029 & animal.x <= 1176.901 & animal.y >= 1389.127 & animal.y <= 1544.332) gsl_from = true;
    // if(animal.x >= 1176.901 & animal.x <= 1313.094 & animal.y >= 1378.905 & animal.y <= 1544.332) gsl_from = true;
    // 
    // // Cape Cod Bay (inside and south)
    // if(animal.x >= 602.6252 & animal.x <= 664.1906 & animal.y >= 883.3870 & animal.y <= 930.7010) cod_from = true;
    // if(animal.x >= 530 & animal.x <= 676.6508 & animal.y >= 800 & animal.y <= 883.3873) codsouth_from = true;
    
    
    
    // Otherwise, select point within stepsize closest to latent animal
    double d, a, x0, y0, xn, yn, xnew, ynew;
    double dmin = std::numeric_limits<double>::infinity();
    bool sampling = true;
    while(sampling) {
      
      d = std::sqrt(R::runif(0,stepsize_sq));
      a = R::runif(-M_PI,M_PI);
      
      // Disc center
      x0 = animal.x + d * std::cos(a);
      y0 = animal.y + d * std::sin(a);
      
      // Sample m points in the circle, keep closest valid point
      for(std::size_t i=0; i<m; ++i) {
        
        // anticosti_1_to = false;
        // anticosti_2_to = false;
        // fundy_1_to = false;
        // fundy_2_to = false;
        // fundy_3_to = false;
        // fundy_4_to = false;
        
        // bof_to = false;
        // gsl_to = false;
        // cod_to = false;
        // codsouth_to = false;
        
        d = std::sqrt(R::runif(0,stepsize_sq));
        a = R::runif(-M_PI,M_PI);
        
        xn = x0 + d * std::cos(a);
        yn = y0 + d * std::sin(a);
        
        // if(xn >= 1042 & xn <= 1198 & yn >= 1782 & yn <= 1862) anticosti_1_to = true;
        // if(xn >= 1042 & xn <= 1160 & yn >= 1809 & yn <= 1862) anticosti_2_to = true;
        // if(xn >= 1002 & xn <= 1082 & yn >= 1342 & yn <= 1410) fundy_1_to = true;
        // if(xn >= 1082 & xn <= 1202 & yn >= 1182 & yn <= 1290) fundy_2_to = true;
        // if(xn >= 1042 & xn <= 1222 & yn >= 1380 & yn <= 1482) fundy_3_to = true;
        // if(xn >= 1042 & xn <= 1222 & yn >= 1380 & yn <= 1482) fundy_4_to = true;
        // 
        // if((anticosti_1_from && anticosti_1_to) || (anticosti_2_from && anticosti_2_to) ||
        //    (fundy_1_from && fundy_1_to) || (fundy_2_from && fundy_2_to) || (fundy_3_from && fundy_3_to) ||
        //    (fundy_4_from && fundy_4_to)){
        //   
        //   // Do nothing
        //   
        // } else {
        //   
        //   // std::cout << animal.x << "," << animal.y << " " << xn << "," << yn << " " <<
        //   //   bof_from << "," << gsl_from << " " << bof_to << "," << gsl_to << " " << latent_envs[environment.id](xn, yn) << "\n";
        //   
        //   // Only consider proposals with non-zero mass
        //   if(latent_envs[environment.id](xn, yn) > 0) {
        //     
        //     dist_to_latent = std::pow(xn - active_latent.x, 2) + std::pow(yn - active_latent.y, 2);
        //     
        //     // Keep track of the point that gets closest to latent
        //     if(dist_to_latent < dmin) {
        //       xnew = xn;
        //       ynew = yn;
        //       dmin = dist_to_latent;
        //     }
        //   }
        //   
        // }
        
        // Constraint on movements between regions
        // int r_to = latent_envs[environment.id](xn, yn)[1];
        // 
        // if((r_from == 2 & r_to == 6) || (r_from == 6 & r_to == 2) ||
        //    (r_from == 8 & r_to == 6) || (r_from == 6 & r_to == 8) ||
        //    (r_from == 8 & r_to == 2) || (r_from == 2 & r_to == 8)){
        //  
        //  // Do nothing
        //  
        // } else {
        //   
        //   // Only consider proposals with non-zero mass
        //   if(latent_envs[environment.id](xn, yn)[0] > 0) {
        //     
        //     dist_to_latent = std::pow(xn - active_latent.x, 2) +
        //       std::pow(yn - active_latent.y, 2);
        //     
        //     // Keep track of the point that gets closest to latent
        //     if(dist_to_latent < dmin) {
        //       xnew = xn;
        //       ynew = yn;
        //       dmin = dist_to_latent;
        //     }
        //   }  
        //   
        // }
        
        // Only consider proposals with non-zero mass
        if(latent_envs[environment.id](xn, yn) > 0) {

          dist_to_latent = std::pow(xn - active_latent.x, 2) +
            std::pow(yn - active_latent.y, 2);

          // Keep track of the point that gets closest to latent
          if(dist_to_latent < dmin) {
            xnew = xn;
            ynew = yn;
            dmin = dist_to_latent;
          }
        }
        
      }
      
      Rcpp::checkUserInterrupt();
      sampling = !std::isfinite(dmin);
    }
    animal.x = xnew;
    animal.y = ynew;
  }
  
};

#endif
