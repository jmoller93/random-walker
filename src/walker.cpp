#include <Eigen/Core>
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
    vector3_t walker::get_com(void) {return x_.colwise().mean();};

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

    f_type walker::dist(const int& a, const int& b) {
       vector3_t vec = vector3_t::Zero(); 
       vec = x_.row(a) - x_.row(b);
       return vec.norm();
    };

    f_type walker::get_rg(void) {
        f_type rg = 0.0;
        vector3_t com = get_com();
        rg = (x_.rowwise()-com.transpose()).rowwise().squaredNorm().sum()/float(x_.rows());
        return rg;
    };

    f_type walker::get_gene_length(void) {return l_.sum();};

    int walker::nearest_neighbor(const vector3_t& v) {
        f_type tol = 5.0;
        Eigen::ArrayXXf tmp(x_.rows(),3);
        tmp = (x_.rowwise()-v.transpose()).rowwise().norm();
        return (tmp < tol).count();
    };

    //Chain Growth
    bool walker::chain_growth(void) {
        int jdx = 1;
        const int max_trial=100; 
        for (size_t i=2; i<npoints_; i++) {
            for (size_t j=0; j<max_trial; j++) {
                vector3_t x = x_.row(i-1).transpose()+rand_sphere();
                int count = nearest_neighbor(x);
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

//using namespace walkers;
//int main(int argc, char* argv) {
//    walker sarw = walker(100000,10);
//    bool testBool = false;
//    while(!testBool) {
//        testBool = sarw.chain_growth();
//    }
//    std::cout<<sarw.get_coord()<<std::endl;
//    std::cout<<sarw.get_rg()<<std::endl;
//    //std::cout<<sarw.rand_sphere()<<std::endl;
//    //std::cout<<sarw.get_lengths()<<std::endl;
//    //std::cout<<sarw.get_glength()<<std::endl;
//    //std::cout<<sarw.get_com().transpose()<<std::endl;
//    return 0;
//}

