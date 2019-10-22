#include <Eigen/Core>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <random>
#include "walker.h"

namespace walkers {
    //Initializers
    //walker::walker() {};

    walker::walker(const uint& npoints, const uint& l0) : npoints_(npoints), l0_(l0)  {
        x_.resize(npoints_,3);
        x_.row(0) = vector3_t::Zero().transpose();
        x_.row(1) = rand_sphere().transpose();
        l_ = vector_t::Ones(npoints_-1)*187;
    };

    //Setters
    void walker::set_monomers(const uint& npoints) {npoints_ = npoints;};
    void walker::set_bond_length(const uint& l0) {l0_ = l0;};

    //Getters
    const matrix3_t& walker::get_coord(void) const {return x_;};
    const vector_t& walker::get_lengths(void) const {return l_;};

    //Utility Functions
    vector3_t walker::get_com(void) const {return x_.colwise().mean();};

    vector3_t walker::rand_sphere(void) {
        //Random number generation
        std::random_device rd{};
        std::mt19937 generator{rd()};
        std::uniform_real_distribution<f_type> distribution(0.0,1.0);

        vector3_t vec = vector3_t::Zero();
        f_type theta = 2.0 * M_PI * distribution(generator);
        f_type phi = acos(1.0 - 2.0 * distribution(generator));
        vec << sin(phi)*cos(theta),
               sin(phi)*sin(theta),
               cos(phi);
        return vec*l0_;
    };

    f_type walker::dist(const int& a, const int& b){
       vector3_t vec = vector3_t::Zero(); 
       vec = x_.row(a) - x_.row(b);
       return vec.norm();
    };

    f_type walker::get_rg(void) const{
        f_type rg = 0.0;
        vector3_t com = get_com();
        rg = (x_.rowwise()-com.transpose()).rowwise().squaredNorm().sum()/float(x_.rows());
        return rg;
    };

    f_type walker::get_gene_length(void) const{return l_.sum();};

    int walker::nearest_neighbor(const vector3_t& v, f_type tol) {
        Eigen::ArrayXXf tmp(x_.rows(),3);
        tmp = (x_.rowwise()-v.transpose()).rowwise().norm();
        return (tmp < tol).count();
    };

    //Chain Growth
    void walker::chain_growth(const f_type tol, const uint max_trial) {

        if (tol < 0.0) 
            throw std::invalid_argument("Tolerance for nearest neighbor must be greater than zero");
        if (max_trial < 0)
            throw std::invalid_argument("Maximum number of trials must be greater than zero");

        bool chain_bool = false;
        while (!chain_bool)
            chain_bool = chain_test(tol,max_trial);

        return;
    };

    bool walker::chain_test(const f_type tol, const uint max_trial) {
        uint jdx = 1;
        for (size_t i=2; i<npoints_; i++) {
            for (size_t j=0; j<max_trial; j++) {
                vector3_t x = x_.row(i-1).transpose()+rand_sphere();
                int count = nearest_neighbor(x,tol);
                if (count == 0) {
                    x_.row(i) = x.transpose();
                    break;
                }
                jdx++;
            }
            if (jdx == max_trial) {
                return false;
            }
        }
        return true;
    };
}

