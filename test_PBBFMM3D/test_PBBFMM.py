import numpy as np
import cppimport

if __name__ == "__main__":

    pbbfmm = cppimport.imp("PBBFMM_binding_test")
    N = 10000
    kshot = (np.random.rand(3,N) - 0.5)*1/np.pi
    weights = np.random.rand(N,1)
    
    interpolation_order = 5
    tree_depth = 4
    eps = 1e-5
    
    
    Qh = pbbfmm.pbbfmm_3D(kshot, weights, interpolation_order, tree_depth, eps)
    print("shape of output is " + str(Qh.shape))
