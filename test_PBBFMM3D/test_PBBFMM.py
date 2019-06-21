import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import cppimport

#%%  QH calculated with pybind11 & PBBFMM3D
shot = np.loadtxt("../data/distribution_radial3D.txt", delimiter = ',', dtype = np.float64)

pfmm = cppimport.imp("PBBFMM_binding_test")

interpolation_order = 5 #Number of Chebyshev nodes for interpolations (between 3 and 10)
tree_depth = 6 #in PBBFMM3D, the user controls the tree depth
eps = 1e-5 #target precision for compressed M2L operator
H = np.ones(shot.shape[0]) #set of unitary charges for testing
Qh = pfmm.pbbfmm_3D(shot.T, H, interpolation_order, tree_depth, eps) ##pbbfmm_3D takes 3*N arrays
np.savetxt("../data/QH_radial3D_binding.txt", Qh, delimiter = ",")

fig2 = plt.figure()
plt.plot(Qh[:,0], linewidth = 0.2)
plt.title("QH product with pybind11 and PBBFMM3D")
plt.show()


#%% QH calculated "directly" in PBBFMM3D
Qh_bb = np.loadtxt("../data/QH_radial3D_pbbfmm.txt")

fig3 = plt.figure()
plt.plot(Qh_bb, linewidth = 0.2)
plt.xlabel("Particles indices", fontsize = 18)
plt.ylabel("Potential", fontsize = 18)
plt.title("QH product, PBBFMM3D directly")
plt.show()


#%% Some figures to diplay calculation differences
diffQh = Qh[:,0] - Qh_bb
ratioQh = 100*diffQh/Qh_bb

fig4 = plt.figure()
plt.plot(diffQh, linewidth = 0.4)
plt.title("Difference between PBBFMM and PBBFMM + pybind11")
plt.show()
#plt.savefig("../Results/QH_radial3D_difference.png")

fig5 = plt.figure()
plt.plot(ratioQh, linewidth = 0.4)
plt.title("Relative error betwwen PBBFMM and pybind11 + PBBFMM3D")
plt.ylabel("Relative error in %", fontsize = 18)
plt.show()

#%% Checking if anything happens in the casting of the vector in the C++ script
castedSourceDirect = np.loadtxt("../data/castedSourceDirect_radial3D.txt", delimiter = ",", dtype = np.float64)

#%% Display shot (computationally heavy) for visual check

#fig = plt.figure()
#ax = fig.gca(projection = "3d")
#
#unit = 1/4 ; tick = np.arange(-0.5, 0.5 + unit, unit)
#label_pi = [r"$-\frac{1}{2\pi}$", r"$-\frac{1}{4\pi}$", r"$0$", r"$+\frac{1}{4\pi}$", r"$+\frac{1}{2\pi}$"]
#
#ax.scatter(castedSourceDirect[:,0], castedSourceDirect[:,1], castedSourceDirect[:,2], marker = ".", s = 0.1)
#
#ax.set_xticks(tick/np.pi)
#ax.set_xlabel(r"$k_x$", fontsize = 18)
#ax.set_xticklabels(label_pi)
#ax.set_yticks(tick/np.pi)
#ax.set_ylabel(r"$k_y$", fontsize = 18)
#ax.set_yticklabels(label_pi)
#ax.set_zticks(tick/np.pi)
#ax.set_zlabel(r"$k_z$", fontsize = 18)
#ax.set_zticklabels(label_pi)
#plt.title("Distribution of particles")
#plt.show()

#%%
#if __name__ == "__main__":
#
#    pbbfmm = cppimport.imp("PBBFMM_binding_test")
#    N = 10000
#    kshot = (np.random.rand(3,N) - 0.5)*1/np.pi
#    weights = np.random.rand(N,1)
#    
#    interpolation_order = 5
#    tree_depth = 4
#    eps = 1e-5
#    
#    
#    Qh = pbbfmm.pbbfmm_3D(kshot, weights, interpolation_order, tree_depth, eps)
#    print("shape of output is " + str(Qh.shape))
