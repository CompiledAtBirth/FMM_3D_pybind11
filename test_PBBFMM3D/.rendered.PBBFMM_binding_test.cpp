


//path to pybind11/include needed in cfg['include_dirs'] if it is not contained in the CPATH environment variable

#define _USE_MATH_DEFINES
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include "../PBBFMM3D/include/bbfmm3d.hpp"
#include "stdint.h"
#include <cmath>

namespace py = pybind11;

//The request method of a pybind11::array returns:
//struct buffer_info {
//    void *ptr;
//    size_t itemsize;
//    std::string format;
//    int ndim;
//    std::vector<size_t> shape;
//    std::vector<size_t> strides;
//};


//Converts input file into vector of custom vector3 type ; the "c_style" protection ensures row-major packing
void loc2vector(py::array_t<double, py::array::c_style | py::array::forcecast> locations, int N, std::vector<vector3>& source, std::vector<vector3>& target){

    py::buffer_info loc = locations.request();
    double *ptrLoc = (double *) loc.ptr;
    
    
    //Input array is 3*N, so second row starts at pointer N for instance
    for(int k = 0; k < N ; k++){
        
        source[k].x = ptrLoc[k];
        source[k].y = ptrLoc[k+N];   
        source[k].z = ptrLoc[k+2*N];    
    }
    
    for(int n = 0 ; n < N ; n++){
    
        target[n].x = source[n].x;
        target[n].y = source[n].y;
        target[n].z = source[n].z;
    
    }
}

void chargesSet(py::array_t<double, py::array::c_style | py::array::forcecast> weights, int N, std::vector<double>& charges){
    
    py::buffer_info ch = weights.request();
    double *ptrCh = (double *) ch.ptr;
    
    for(int j = 0 ; j < N; j++){
        charges[j] = ptrCh[j];
    }

}

//MAIN FUNCTION -> MODULE
py::array_t<double, py::array::c_style | py::array::forcecast> pbbfmm_3D(py::array_t<double, py::array::c_style | py::array::forcecast> locations, py::array_t<double, py::array::c_style | py::array::forcecast> weights, int interpolation_order, int tree_level, double eps){

    int nCols = 1; //always one set of charges
    int use_chebyshev = 1; //use Chebyshev interpolation scheme
    double L = 1.0/M_PI; //kspace samples live in [-1/(2pi), 1/(2pi)]^3
    std::cout << "L = " << L << std::endl;
    
    //Pybind11 buffers
    py::buffer_info loc = locations.request();
    int N = loc.shape[1];
    
    std::vector<vector3> source(N); // Position array for the source points
    
    std::vector<vector3> target(N);  // Position array for the target points

    std::vector<double> charges(N*nCols); // charges array
    
    chargesSet(weights, N, charges);
    
    //std::cout<<"Sum of all charges is :"<<std::endl;
    //std::cout<< std::accumulate(charges.begin(),charges.end(),0) << std::endl;
    
    loc2vector(locations, N, source, target);
    
    ///////////DEBUG///////////////
    
    //std::cout << "target0 x = " << target[0].x << std::endl;
    //std::cout << "target0 y = " << target[0].y << std::endl;
    //std::cout << "target0 z = " << target[0].z << std::endl;
    
    //std::cout << "source0 x = " << source[0].x << std::endl;
    //std::cout << "source0 y = " << source[0].y << std::endl;
    //std::cout << "source0 z = " << source[0].z << std::endl;
    
    //const char* outname;
    //outname = "../data/castedTarget_radialIO3434.txt";
    
    //ofstream outfile;
    //outfile.open(outname);
    
    //for(int n = 0; n < N*nCols ; n++){
    
    //    outfile << target[n].x;
    //    outfile << "," << target[n].y << "," << target[n].z ;
    //    if(n != N-1){
    //        outfile << '\n';
    //    }
    //}
    //outfile.close();
    //////////////////////////////////////
    
    std::vector<double> output(N*nCols);   // output array (PBBFMM calculation)    
    
    for (int i = 0; i < N*nCols; i++)
        output[i] = 0;                          //do we need to initialize?
        
     /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    kernel_Laplacian Atree(L,tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();
    
    /*****      FMM Computation     *******/  
    H2_3D_Compute<kernel_Laplacian> compute(Atree, target, source, charges, nCols, output);
    
    std::cout << "output[0] = " << output[0] << " (binding)"<< std::endl;
    
    //output tunrned into py::array_t
    
    py::array_t<double, py::array::c_style | py::array::forcecast> Qh({N, nCols});
    py::buffer_info out = Qh.request();
    double *ptrOut = (double *) out.ptr;
    
    for(int j = 0 ; j < N*nCols ; j++){
        ptrOut[j] = output[j];
    }
    
    //Clean memory
    
    target.clear();
    source.clear();
    charges.clear();
    
    target.shrink_to_fit();
    source.shrink_to_fit();
    charges.shrink_to_fit();
    
    return Qh;
}

PYBIND11_MODULE(PBBFMM_binding_test, m) {
    m.doc() = "Calculate the potentials with PBBFMM"; // optional

    m.def("pbbfmm_3D", &pbbfmm_3D, py::return_value_policy::reference);
}
