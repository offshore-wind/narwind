#ifndef ENVIRONMENTS_H
#define ENVIRONMENTS_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

class LatentAttractor {
  
public:
  
  // Index of latent animal toward which movement is attracted
  const std::size_t &id;
  LatentAttractor(const std::size_t & animal_id) : id(animal_id) { }
};

// Environment that captures features for NARW simulation

class Environment {

    private:

        using MatrixType = Eigen::MatrixXd;

        // Population density raster (in projected coordinates)
        const MatrixType &m_density;
        
    public:
      
      // Spatial extents and resolution of density raster
      const Eigen::VectorXd &m_limits, &m_resolution;
      
      // Numeric id for map
      const std::size_t id;
      
      // Define a constructor for Environment - gets called every time en environment is constructed
      // This is used to initialize the class members
      // The code after the : is a member initializer list
      // m_density(density) is the same as m_density = density; and ensures that the value
      // passed to the <density> parameter during initialization is assigned to class member m_density
      Environment(
        const MatrixType & density,
        const Eigen::VectorXd & limits,
        const Eigen::VectorXd & resolution, 
        const std::size_t & mapid
      ) : 
        m_density(density),
        m_limits(limits), 
        m_resolution(resolution),
        id(mapid) { }

        /**
         * Return the density value at the specified coordinates
         * @param x
         * @param y
         * @return Density value
         */
        // Eigen::VectorXd operator()(const double & x, const double & y) {
        double operator()(const double & x, const double & y) {
          
          // Return 0 if outside the study area boundaries
            if(x < m_limits[0] || x > m_limits[1] || y < m_limits[2] || 
               y > m_limits[3]){
              return 0;
              // return Eigen::Matrix2d::Zero().row(0);
            }


            // Phil (Nov 16, 2022) - changed floor to round
            // Retrieve row and col of cell corresponding to x,y
            std::size_t i = std::round((x-m_limits[0])/m_resolution[0]);
            std::size_t j = std::round((m_limits[3]-y)/m_resolution[1]);

            double d = m_density(i,j);
            
            if(!std::isfinite(d)) {
              return 0;
              // return Eigen::Matrix2d::Zero().row(0);
            }

            return d;
            // Eigen::VectorXd out(2);
            // out[0] = d;
            // out[1] = r;
            // return out;
        }

        void Print(){
          std::cout << "I am happy" << std::endl;
        }
      
};

// Prey environment

class PreyEnv {
  
private:
  
  using MatrixType = Eigen::MatrixXd;
  
  // Population density raster (in projected coordinates)
  const MatrixType &m_prey;
  
public:
  
  // Spatial extents and resolution of density raster
  const Eigen::VectorXd &m_limits, &m_resolution;
  
  // Numeric id for map
  const std::size_t id;
  
  // Initialize class members
  PreyEnv(
    const MatrixType & prey,
    const Eigen::VectorXd & limits,
    const Eigen::VectorXd & resolution, 
    const std::size_t & mapid
  ) : 
    m_prey(prey),
    m_limits(limits), 
    m_resolution(resolution),
    id(mapid) { }
  
  /**
   * Return the prey concentration at the specified coordinates
   * @param x
   * @param y
   * @return Prey concentration
   */
  double operator()(const double & x, const double & y) {
    
    // Return 0 if outside the study area boundaries
    if(x < m_limits[0] || x > m_limits[1] || y < m_limits[2] || y > m_limits[3]) return 0;
    
    // Retrieve row and col of cell corresponding to x,y
    std::size_t i = std::round((x-m_limits[0])/m_resolution[0]);
    std::size_t j = std::round((m_limits[3]-y)/m_resolution[1]);
    
    double p = m_prey(i,j);
    if(!std::isfinite(p)) return 0;
    return p;
  }
  
};

/**
 * Wrap an environment object inside a container that will appear to iterate
 * over n_elem copies of the wrapped environment
 *
 * @tparam EnvironmentType Type for wrapped environment object
 */
template<typename EnvironmentType>
struct ConstEnvironment {

    using MatrixType = Eigen::MatrixXd;
  
    // Iterator keeps count of iterations, but always points to the same object
    struct Iterator {

        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = EnvironmentType;
        using pointer           = value_type*;
        using reference         = value_type&;

        Iterator(pointer ptr, std::size_t i) : m_ptr(ptr), ind(i) { }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment: constant return until end
        Iterator& operator++() { ++ind; return *this; }

        // Postfix increment: constant return
        Iterator operator++(int) { ++ind; return this; }

        friend bool operator== (const Iterator& a, const Iterator& b) {
            return a.ind == b.ind;
        };
        friend bool operator!= (const Iterator& a, const Iterator& b) {
            return a.ind != b.ind;
        };

        private:
            pointer m_ptr;
            std::size_t ind;
    };

    EnvironmentType & m_env;
    std::size_t n_elem;
    const Eigen::VectorXd &m_limits, &m_resolution;

    /**
     * @param environment Environment object iterator should point to
     * @param n Number of virtual copies to "create"
     */
    ConstEnvironment(EnvironmentType & environment, 
                     std::size_t n, 
                     const Eigen::VectorXd &lim,
                     const Eigen::VectorXd &res) :
      m_env(environment), n_elem(n), m_limits(lim), m_resolution(res) { }
    
    std::size_t size() { return n_elem; }
    Iterator begin() { return Iterator(&m_env, 0); }
    Iterator end()   { return Iterator(&m_env, n_elem); }
};

#endif
