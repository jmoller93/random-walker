#pragma once

#include <string>
#include "types.h"

namespace walkers {
    //Walker class initialization
    class walker {

        private:
            uint npoints_; // Number of points in the walker
            uint l0_;      // Bond distance
            matrix3_t x_;  // Positions
            vector_t l_;   // Linker lengths

            //Utility Functions
            vector3_t rand_sphere(void);  // Returns a random point on the surface of a sphere
            int nearest_neighbor(const vector3_t&); // Calculates the distance to the nearest neighbor
        public:
            //Initializers
            //walker();
            walker(const uint&, const uint&);
            //walker(const std::string&); // Initialize by file to be implemented

            //Setters
            void set_monomers(const uint&);
            void set_bond_length(const uint&);

            //Getters
            const matrix3_t& get_coord(void) const;
            //const matrix_t& get_dists(void) const;
            const vector_t& get_lengths(void) const;

            //Growth of chain
            bool chain_growth(void);      // Grow chain

            //Utility Functions
            vector3_t get_com(void);      // Gets center of mass
            f_type get_rg(void);          // Get the radius of gyration
            f_type get_gene_length(void);     // Gets the genomic length 
            f_type dist(const int&, const int&); // Calculate the euclidean distance of two points
    };
}
