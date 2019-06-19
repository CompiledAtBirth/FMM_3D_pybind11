#define _USE_MATH_DEFINES
#include <iostream>
#include "../PBBFMM3D/include/bbfmm3d.hpp"
#include "stdint.h"
#include <cmath>


void split(const string& s, char* del, std::vector<vector3>& source, std::vector<vector3>& target){
    
    size_t i = 0;
    size_t j = s.find(del);
    vector3 v;
    
    v.x = atof(s.substr(0, j).c_str()); //c_str() turns string into char*
    i = ++j;
    j = s.find(del, j);
    v.y = atof(s.substr(i, j-i).c_str());
    v.z = atof(s.substr(j+1, s.length()-j-1).c_str());
    
    source.push_back(v);
    target.push_back(v);
}

//Transform input into a vector of struct vector3, and returns the size by the same occasion   
int txt2Vec(const char* filename, char* del, std::vector<vector3>& source, std::vector<vector3>& target){

    string line;
    ifstream inFile;
    inFile.open(filename);
    int N = 0;
    
    while (getline(inFile,line))
    {
     split(line, del, source, target);
     N++;
    }
    
    inFile.close();
    return N;
}

//Perfomes the FMM on the input distribution
int main(){
    
    const char* filename;
    filename = "../data/radialIO3D_3434.txt";
    char del[] = ",";
    
    int use_chebyshev = 1; //use Chebyshev interpolation scheme
    double L = 1.0/M_PI; //kspace samples live in [-1/(2pi), 1/(2pi)]^3
    std::cout << "L = " << L << std::endl;
    int nCols = 1;
    int interpolation_order = 5;
    int tree_level = 6;
    double eps = 1e-5;
    
    std::vector<vector3> source;
    std::vector<vector3> target;
    
    int N = txt2Vec(filename, del, source, target);
    
    ////////////////DEBUG////////////////
    
    const char* castname;
    castname = "../data/castedSourceDirect_radialIO3434.txt";
    
    ofstream castfile;
    castfile.open(castname);
    
    for(int n = 0; n < N*nCols ; n++){
    
        castfile << source[n].x;
        castfile << "," << source[n].y << "," << source[n].z ;
        if(n != N-1){
            castfile << '\n';
        }
    }
    castfile.close();
    
    //Unitary set of charges for testing
    
    std::vector<double> charges(N*nCols);
    
    for(int k = 0 ; k < N*nCols ; k++){    
        charges[k] = 1.0;
    }
    
    std::vector<double> output(N*nCols);
    
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
    
    std::cout << "output[0] = " << output[0] <<std::endl;
    
    //write potentials into text file
    
    const char* outname;
    outname = "../data/QH_radialIO3434_pbbfmm.txt";
    
    ofstream outfile;
    outfile.open(outname);
    
    for(int n = 0; n < N-1; n++){
    
        outfile << output[n];
        outfile << '\n'; 
    }
    
    outfile << output[N-1];
    outfile.close();

    return 0;
}
