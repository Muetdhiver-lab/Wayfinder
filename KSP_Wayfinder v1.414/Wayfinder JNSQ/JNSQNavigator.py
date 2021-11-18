from pykep import planet, DEG2RAD, epoch, AU
import pykep
import sys
sys.path.append("../WayfinderCore")
print(pykep.__version__)
from _JNSQ_System import Moho, Eve, Kerbin, Duna

KAU = 37525647898 #m

Planets = [Moho,Eve,Kerbin,Duna]



def calc_SOI(body):
    SOI = body.orbital_elements[0]*(body.mu_self/body.mu_central_body)**0.4
    print(SOI)
    
    


def plot_innerKerbol(epoch = epoch(0)):
    
    """
    Plots the Galilean Moons of Jupiter at epoch

    USAGE: plot_moons(epoch = epoch(34654, epoch.epoch_type.MJD)):
    * epoch: the epoch one wants the galilean moons to be plotted at
	
    """
    from pykep.orbit_plots import plot_planet
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    print(epoch)
    fig = plt.figure()
    fig.set_size_inches(10,6)
    #ax1 = fig.gca(projection='3d')

	 #plot_planet(ax,Moho,color = 'y', units = 1.0, t0 = epoch, legend=True)
	 #plot_planet(ax,Eve,color = 'm', units = 1.0, t0 = epoch, legend=True)
	 #lot_planet(ax=ax1,Kerbin,color = 'c', units = 1.0, t0, legend=True)
    #plot_planet(ax = ax1,Duna,color = 'r', units = 1.0, t0 = epoch, legend=True)
    plot_planet(Moho, t0=epoch, color = 'y', legend=True, units=KAU)
    plot_planet(Eve, t0=epoch, color = 'm',   legend=True, units=KAU)
    plot_planet(Kerbin, t0=epoch, color = 'c', legend=True, units=KAU)
    plot_planet(Duna, t0=epoch, color = 'r',   legend=True, units=KAU) 
    ax1.set_xlim3d(-1,1)
    ax1.set_ylim3d(-1,1)
    ax1.set_zlim3d(-1,1)
    plt.show()
    
    
plot_innerKerbol(epoch(100))

for Planet in Planets :
    calc_SOI(Planet)