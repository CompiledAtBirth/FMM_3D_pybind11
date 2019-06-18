import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import cppimport

#%%  QH calculated with pybind11 & PBBFMM3D
shot = np.loadtxt("../data/radialIO3D_3434.txt", delimiter = ',')

pfmm = cppimport.imp("PBBFMM_binding_test")

interpolation_order = 5 #Number of Chebyshev nodes for interpolations (between 3 and 10)
tree_depth = 6
eps = 1e-5
H = np.ones(shot.shape[0]) #set of unitary charges
Qh = pfmm.pbbfmm_3D(shot.T, H, interpolation_order, tree_depth, eps) ##bbfmm_2D takes 2*N arrays
np.savetxt("Qh_radialIO3434_binding.txt", Qh, delimiter = ",")

fig2 = plt.figure()
plt.plot(Qh[:,0], linewidth = 0.2)
plt.title("Qh product with pybind11 and PBBFMM3D")
plt.show()

#%% Qh calculated "directly" in PBBFMM3D
Qh_bb = np.loadtxt("../data/QH_radialIO3434_pbbfmm.txt")

fig3 = plt.figure()
plt.plot(Qh_bb, linewidth = 0.2)
plt.xlabel("Particles indices", fontsize = 18)
plt.ylabel("Potential", fontsize = 18)
plt.title("QH product, PBBFMM3D directly")
plt.title("First QH product, bbfmm no binding")
plt.show()

diffQh = Qh[:,0] - Qh_bb
ratioQh = 100*diffQh/Qh_bb

fig4 = plt.figure()
plt.plot(diffQh, linewidth = 0.4)
plt.title("Difference between PBBFMM calculation and calcutation with pybind11")
plt.show()

fig5 = plt.figure()
plt.plot(ratioQh, linewidth = 0.4)
plt.title("Relative error betwwen direct pbbfmm calculation and pybind11")
plt.ylabel("Relative error in %", fontsize = 18)
plt.show()
#%% Checking if anything happens in the casting
castedSourceDirect = np.loadtxt("../data/castedSourceDirect_radialIO3434.txt", delimiter = ",")


#%% Display shot (computationally heavy)

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