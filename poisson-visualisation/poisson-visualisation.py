#----------------------------------------
#   Visualisation:  2D Poisson
#   Authors:        Azin and Alexander
#----------------------------------------


#--------------------
#   Packages
#--------------------
import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib import gridspec
from matplotlib import ticker
from matplotlib import cm




#------------------------------
#   Open files   >   2D array
#------------------------------
approx_file = open("/Users/azinshahiri/Desktop/mpi-projects/mpi-poisson/mpi-poisson/approximation.txt",'r')
approx = [[float(number) for number in line.split()] for line in approx_file]
approx = np.array(approx)
approx_file.close()

diff_file = open("/Users/azinshahiri/Desktop/mpi-projects/mpi-poisson/mpi-poisson/difference.txt",'r')
diff = [[float(number) for number in line.split()] for line in diff_file]
diff = np.array(diff)
diff_file.close()




#--------------------------------------------------
#   Figure:		Poisson equation:	Jacobi method
#--------------------------------------------------
poissonplots = pyplot.figure("Poisson equation: Jacobi method", figsize=(11,7))
spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[1, 1], height_ratios=[3,2.2], wspace=0.5, hspace=0.2)


#------------------------------
#   Surfaceplots:   create mesh
#------------------------------
N = 20
X = np.arange(0, N+2, 1)
Y = np.arange(0, N+2, 1)
X, Y = np.meshgrid(X, Y)


#--------------------
#   Surfaceplots:
#--------------------
surfaceplot = poissonplots.add_subplot(spec[0], projection='3d')
surfaceplot.set_title("Approximation: u_approx")
surfaceplot.dist=8
surfaceplot.elev=20
# surfaceplot.zaxis.set_major_locator(ticker.LinearLocator(3))
surfaceplot.plot_surface(X,Y,approx, cmap=cm.coolwarm, linewidth=0)


surfaceplot = poissonplots.add_subplot(spec[1], projection='3d')
surfaceplot.set_title("Difference: u_diff")
surfaceplot.dist=8
surfaceplot.elev=20
surfaceplot.plot_surface(X,Y,diff, cmap=cm.coolwarm, linewidth=0)


#--------------------
#   Colour plot
#--------------------
colourplot = poissonplots.add_subplot(spec[2])
colourplot.imshow(approx, cmap=cm.coolwarm)
colourplot.invert_yaxis()


colourplot = poissonplots.add_subplot(spec[3])
colourplot.imshow(diff, cmap=cm.coolwarm)
colourplot.invert_yaxis()


#--------------------------------------------------
#	Show:	Poisson equation:	Jacobi method
#--------------------------------------------------
pyplot.show()



#--------------------
#   END-OF-FILE
#--------------------