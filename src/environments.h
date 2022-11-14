#ifndef ENVIRONMENTS_H
#define ENVIRONMENTS_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

class LatentAttractor {
    public:
        // index of latent animal toward which movement is attracted
        const std::size_t &id;
        LatentAttractor(const std::size_t & animal_id) : id(animal_id) { }
};

/**
 * Environment that captures features for NARW simulation
 */
class Environment {

    private:

        using MatrixType = Eigen::MatrixXd;

        // population density raster (in projected coordinates)
        const MatrixType &m_density;
        // spatial extents and resolution of density raster
        const Eigen::VectorXd &m_limits, &m_resolution;

    public:

        // numeric id for map
        const std::size_t id;

        Environment(
            const MatrixType & density, const Eigen::VectorXd & limits,
            const Eigen::VectorXd & resolution, const std::size_t & mapid
        ) : m_density(density), m_limits(limits), m_resolution(resolution),
            id(mapid) { }

        /**
         * Return the density value at the specified coordinates
         * @param x
         * @param y
         * @return
         */
        double operator()(const double & x, const double & y) {
            if(x < m_limits[0] || x > m_limits[1] || y < m_limits[2] ||
               y > m_limits[3])
                return 0;

            std::size_t i = std::floor((x-m_limits[0])/m_resolution[0]);
            std::size_t j = std::floor((m_limits[3]-y)/m_resolution[1]);

            double d = m_density(i,j);

            if(!std::isfinite(d))
                return 0;

            return d;
        }

};

/**
 * dummy environment for use when simulation does not require environment state
 */
struct EmptyEnvironment { };

/**
 * Environment that provides a resource selection function, as in Theo's example
 */
class RSFEnvironment {

    private:

        using MatrixType = Eigen::MatrixXd;

        const MatrixType &m_cov1, &m_cov2;
        const Eigen::VectorXd &m_beta, &m_lim, &m_res;

    public:

        const std::size_t id = 0;

        /**
         * @param cov1 Matrix of values for the first resource
         * @param cov2 Matrix of values for the second resource
         * @param beta Vector of resource selection coefficients
         * @param lim Four values: x min, x max, y min, y max
         * @param res Two values: resolution in x, and resolution in y
         */
        RSFEnvironment(
            const MatrixType &cov1, const MatrixType &cov2,
            const Eigen::VectorXd &beta, const Eigen::VectorXd &lim,
            const Eigen::VectorXd &res
        ) : m_cov1(cov1), m_cov2(cov2), m_beta(beta), m_lim(lim), m_res(res) { }

        /**
         * Evaluate the resource selection function at the specified coordinates
         *
         * @param x
         * @param y
         * @return
         */
        double operator()(const double &x, const double &y) {
            if(x < m_lim[0] || x > m_lim[1] || y < m_lim[2] || y > m_lim[3])
                return 0;

            std::size_t i = std::floor((m_lim[3]-y)/m_res[1]);
            std::size_t j = std::floor((x-m_lim[0])/m_res[0]);

            double res1 = m_cov1(i,j);
            double res2 = m_cov2(i,j);

            if(!std::isfinite(res1) || !std::isfinite(res2))
                return 0;

            double reg = m_beta[0]*res1 + m_beta[1]*res2;

            return std::exp(reg);
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

    // iterator keeps count of iterations, but always points to the same object
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

    /**
     * @param environment Environment object iterator should point to
     * @param n Number of virtual copies to "create"
     */
    ConstEnvironment(EnvironmentType & environment, std::size_t n) :
        m_env(environment), n_elem(n) { }

    std::size_t size() { return n_elem; }

    Iterator begin() { return Iterator(&m_env, 0); }
    Iterator end()   { return Iterator(&m_env, n_elem); }
};

#endif
