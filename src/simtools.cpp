// R and 'Eigen' integration using 'Rcpp'. 
// 'Eigen' is a C++ template library 
// for linear algebra: matrices, vectors, 
// numerical solvers and related algorithms.

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppEigen.h>
#include "environments.h"
#include "movements.h"

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
 * basic state for an animal
 */
struct Animal {
    double x, y;
    Animal() : x(0), y(0) { }
    Animal(const double & x0, const double & y0) : x(x0), y(y0) { }
};

//' @name Animal
//' @title Animal with capacity for attraction to other latent animals
//'
//' @description Initialize an animal with latent members, and couple the observed position to one of the latent animals
//' 
//' @param init_latent latent animals to be attracted to
//' @param init_animal index of latent animal to use as starting location

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
template<typename AnimalType, typename EnvironmentContainer,
         typename MovementRule, int CoordSize = 2>
Rcpp::NumericVector movesim(
    std::vector<AnimalType> & animals,
    EnvironmentContainer & environments,
    MovementRule & m
) {

    // number of animals and timepoints to simulate
    std::size_t n = animals.size();
    std::size_t t = environments.size();

    // initialize and label output
    Rcpp::NumericVector out(n*t*CoordSize);
    out.attr("dim") = Rcpp::Dimension(CoordSize,n,t);
    out.attr("dimnames") = Rcpp::List::create(
        Rcpp::CharacterVector::create("easting", "northing"), R_NilValue,
        R_NilValue
    );

    // loop over timepoints and environments: one environment per timepoint
    Rcpp::NumericVector::iterator data_out = out.begin();
    auto env_end = environments.end();
    auto animals_end = animals.end();
    for(auto env = environments.begin(); env != env_end; ++env) {
        for(auto animal = animals.begin(); animal != animals_end; ++animal) {
            // save current coords
            *(data_out++) = animal->x;
            *(data_out++) = animal->y;
            // update the animal's state
            m.update(*animal, *env);
        }
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector simBasic(
    std::size_t t, std::size_t n, Eigen::MatrixXd cov1, Eigen::MatrixXd cov2,
    Eigen::VectorXd beta, Eigen::VectorXd lim, Eigen::VectorXd res, double allr
) {

    // initialize population
    std::vector<Animal> animals(n);

    // initialize environment container
    RSFEnvironment env(cov1, cov2, beta, lim, res);
    ConstEnvironment<RSFEnvironment> environments(env, t);

    // initialize movement rule
    // RandomMovement<Animal, RSFEnvironment> rm;
    ReweightedRandomWalk<Animal, RSFEnvironment> rm(100, allr);

    return movesim(animals, environments, rm);
}

/**
 * Simulate
 *
 * @param densities array of density maps, stored in column-major order with
 *   indices [xind, yind, map] for projected spaces (xind is index for a
 *   projected lon coordinate, yind is index for a projected lat)
 * @param densitySeq 1-indexed vector specifying which density maps should be
 *   used to simulate movement at each timepoint; number of timepoints to
 *   simulate comes from length of vector densitySeq
 * @param limits vector (xmin, xmax, ymin, ymax) of spatial coordinate extents
 *   for spatial densities
 * @param resolution vector (xres, yres) of spatial step sizes for spatial
 *   densities
 * @param M specify the number of proposals used in the importance sampler for
 *   movement (Michelot, 2019)
 * @param stepsize the radius of the proposal circle for movement
 *   (Michelot, 2019)
 * @param xinit initial x coordinates for animals to simulate
 * @param yinit initial y coordinates for animals to simulate
 * @return
 */
// [[Rcpp::export]]
Rcpp::NumericVector simAnnual(
    std::vector<Eigen::MatrixXd> densities, std::vector<std::size_t> densitySeq,
    Eigen::VectorXd limits, Eigen::VectorXd resolution, std::size_t M,
    double stepsize, std::vector<double> xinit, std::vector<double> yinit
) {

    // 0-index densitySeq since, coming from R, we assume it is 1-indexed
    for(auto d = densitySeq.begin(); d != densitySeq.end(); ++d)
        --(*d);

    // initialize environments
    std::vector<Environment> environments;
    auto dend  = densitySeq.end();
    for(auto d = densitySeq.begin(); d != dend; ++d)
        environments.emplace_back(densities[*d], limits, resolution, *d);

    // initialize population
    std::vector<Animal> animals;
    auto xend = xinit.end();
    auto xiter = xinit.begin();
    auto yiter = yinit.begin();
    while(xiter != xend) {
        animals.push_back(Animal(*(xiter++), *(yiter++)));
    }

    // initialize movement rule
    ReweightedRandomWalk<Animal, Environment> rm(M, stepsize);

    return movesim(animals, environments, rm);
}


//' Simulate right whale movements
//' 
//' @param densities array of density maps, stored in column-major order with 
//' indices [xind, yind, map] for projected spaces (xind is index for a projected lon coordinate, 
//' yind is index for a projected lat)
//' @param densitySeq 1-indexed vector specifying which density maps should be
//' used to simulate movement at each timepoint; number of timepoints to
//' simulate comes from length of vector densitySeq
//' @param latentDensitySeq 1-indexed vector specifying which density maps should
//' be used to simulate movement for each latent animal
//' @param limits vector (xmin, xmax, ymin, ymax) of spatial coordinate extents
//' for spatial densities
//' @param resolution vector (xres, yres) of spatial step sizes for spatial
//' densities
//' @param M specify the number of proposals used in the importance sampler for
//' movement (Michelot, 2019)
//' @param stepsize the radius of the proposal circle for movement
//' (Michelot, 2019)
//' @param xinit matrix of initial x coordinates for latent animals (nlatent x n)
//' @param yinit matrix of initial y coordinates for latent animals (nlatent x n)
//' @return List of coordinates
// [[Rcpp::export]]

Rcpp::NumericVector simAnnualCoupled(
    std::vector<Eigen::MatrixXd> densities, std::vector<std::size_t> densitySeq,
    std::vector<std::size_t> latentDensitySeq,
    Eigen::VectorXd limits, Eigen::VectorXd resolution, std::size_t M,
    double stepsize, Eigen::MatrixXd xinit, Eigen::MatrixXd yinit
) {

    // 0-index latentDensitySeq since, coming from R, we assume it is 1-indexed
    for(auto d = latentDensitySeq.begin(); d != latentDensitySeq.end(); ++d)
        --(*d);

    // 0-index densitySeq since, coming from R, we assume it is 1-indexed
    for(auto d = densitySeq.begin(); d != densitySeq.end(); ++d)
        --(*d);

    // initialize environments for latent animals
    std::vector<Environment> latent_environments;
    auto dend  = latentDensitySeq.end();
    for(auto d = latentDensitySeq.begin(); d != dend; ++d)
        latent_environments.emplace_back(densities[*d], limits, resolution, *d);

    // initialize environments with attraction to latent animals
    std::vector<LatentAttractor> environments;
    auto dattrend = densitySeq.end();
    for(auto d = densitySeq.begin(); d != dattrend; ++d)
        environments.emplace_back(*d);

    // initialize population
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

    // initialize latent movement rule
    typedef ReweightedRandomWalk<Animal, Environment> LatentMvmt;
    LatentMvmt rm(M, stepsize);

    // initialize observed movement rule
    CoupledRandomWalk<
        CouplingAnimal, LatentAttractor, LatentMvmt, std::vector<Environment>
    > crm(rm, latent_environments, M, stepsize);

    return movesim(animals, environments, crm);
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
double evalEnvironment(
    Eigen::MatrixXd density, Eigen::VectorXd limits, Eigen::VectorXd resolution,
    double x, double y
) {
    Environment e(density, limits, resolution, 0);
    return e(x, y);
}
