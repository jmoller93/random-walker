#pragma once

#include <Eigen/Core>
#include <random>
#include "types.h"

namespace walkers
{
    extern int random_uniform_int(uint max_int) { 
        //Random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::uniform_int_distribution<int> dist(0,max_int);
        return dist(gen);
    };

    extern f_type random_uniform_float(f_type max_float) {
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::uniform_real_distribution<f_type> dist(0.0,1.0);
        return dist(gen);
    };


	//Adapted from https://stackoverflow.com/questions/13290395/ 
	extern void remove_row(matrix3_t& matrix, unsigned int rowToRemove) {
		uint numRows = matrix.rows()-1;
		uint numCols = matrix.cols();

		if( rowToRemove < numRows )
			matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

		matrix.conservativeResize(numRows,numCols);
	};
}
