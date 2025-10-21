#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "basinsCUDA.cuh"
#include "cudaMacros.cuh"

#include <iomanip>
#include <string>

namespace Basins {

    void basinsOfAttraction_2(
        const double tMax,                              
        const int nPts,                                 
        const double h,                                 
        const int amountOfInitialConditions,            
        const double* initialConditions,                
        const double* ranges,                           
        const int* indicesOfMutVars,                   
        const int writableVar,                         
        const double maxValue,                          
        const double transientTime,                    
        const double* values,                           
        const int amountOfValues,                       
        const int preScaller,                          
        const double eps,                              
        std::string OUT_FILE_PATH);                    

}