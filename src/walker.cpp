#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include "walker.h"
#include "utils.h"

namespace walkers {
    //Initializers
    //walker::walker() {};

    walker::walker(const uint& npoints, const uint& l0) : npoints_(npoints), l0_(l0)  {
        x_.resize(npoints_,3);
        x_.row(0) = vector3_t::Zero().transpose();
        x_.row(1) = rand_sphere().transpose();
        l_ = vector_t::Ones(npoints_-1)*187;
    };

    //Initialize by dists
    walker::walker(void) {};
    // Based on https://stackoverflow.com/questions/25389480    
    walker::walker(const std::string& filename)
    {
        #ifdef _OPENMP
        if(num_threads) omp_set_num_threads(num_threads);
        #endif

        std::ifstream in(filename, std::ios::in | std::ios::binary);
        
        // Read in distance matrix. 
        {
            matrix_t::Index rows = 0, cols = 0;
            in.read((char*) (&rows), sizeof(matrix_t::Index));
            in.read((char*) (&cols), sizeof(matrix_t::Index));
            x_.resize(rows, cols);
            in.read((char*) x_.data(), x_.size()*sizeof(matrix_t::Scalar));
        }  
        in.close();
    }

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
        vector3_t vec = vector3_t::Zero();
        f_type theta = 2.0 * M_PI * random_uniform_float(1.0);
        f_type phi = acos(1.0 - 2.0 * random_uniform_float(1.0));
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

    bool walker::nearest_neighbor(const vector3_t& v, uint idx, f_type tol) {
        Eigen::ArrayXXf tmp(idx,3);
        tmp = (x_.block(0,0,idx,3).rowwise()-v.transpose()).rowwise().norm();
        return (tmp < tol).any();
    };

    //Chain Growth
    //Grow according to a Rosenbluth-weighted random walk
    void walker::chain_growth(const f_type tol, const uint max_trial) {
        if (tol < 0.0) 
            throw std::invalid_argument("Tolerance for nearest neighbor must be greater than zero");
        if (max_trial < 0)
            throw std::invalid_argument("Maximum number of trials must be greater than zero");

        uint idx = 2;
        while (idx<npoints_) {
            matrix3_t x_test;
            if (chain_test(idx,tol,max_trial,x_test)) {
                x_.row(idx) = x_test.row(random_uniform_int(x_test.rows()-1));
                idx++;
            }
            else {
                idx = 2;
                x_.block(idx,0,npoints_-2,3) = matrix3_t::Zero(npoints_-2,3);
            }
        }
        return;
    };

    //Grow according to a distance array
    void walker::chain_growth(const vector_t dists) {

        vector3_t tmp = vector3_t::Zero();
        x_.resize(dists.rows()+1,3);
        x_.row(0) = vector3_t(0,0,0).transpose();
        l_.resize(dists.rows());

        for (size_t i=0; i<dists.rows(); i++) {
            tmp += vector3_t(dists(i),0,0);
            x_.row(i+1) = tmp;
        }
    };

    bool walker::chain_test(uint idx, const f_type tol, const uint max_trial, matrix3_t& x_test) {
        uint count = 0;
        for (size_t j=0; j<max_trial; j++) {
            vector3_t x = x_.row(idx-1).transpose()+rand_sphere();
            if (!nearest_neighbor(x,idx,tol)) {
				if (count == 0) {
					x_test.resize(1,3);
					x_test.row(0) = x;
				}
				else {
					x_test.conservativeResize(x_test.rows()+1,x_test.cols());
					x_test.row(x_test.rows()-1) = x;
				}
                count++;
            }
            if (count == 0) 
                return false;
        }
        return true;
    };

    // Based on https://stackoverflow.com/questions/25389480
    void walker::save(const std::string& filename) const
    {
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        
        // Write out distance matrix. 
        {
            matrix_t::Index rows = x_.rows(), cols = x_.cols();
            out.write((char*) (&rows), sizeof(matrix_t::Index));
            out.write((char*) (&cols), sizeof(matrix_t::Index));
            out.write((char*) x_.data(), x_.size()*sizeof(matrix_t::Scalar));
        } 
        out.close();
    };
}

