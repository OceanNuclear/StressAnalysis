#!/home/ocean/anaconda3/bin/python3
from numpy import cos, arccos, sin, arctan, tan, pi, sqrt; from numpy import array as ary; import numpy as np; tau = 2*pi
import matplotlib.pyplot as plt; import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from quat import *
def RotateSigma(stressTensor,theta, THETA, PHI):
	q = ThreeVar_Q([theta, THETA, PHI])
	R = QuatToR(q)
	#(R)(sigma)(R^T)
	return np.matmul(R,np.matmul(stressTensor,R.T))
#sphereFillingCurve
def sphereCoveringCurve(n, numSteps):	#Generate a curve that fills all of a unit sphere. Used in animating the pole figure for an arbitrary input of quaternion.
	steps = np.linspace(0, 1, numSteps)
	theta= arccos(steps*2 -1)
	phi = tau*steps*n	#n is the winding number, i.e. number of longitudinal cycles to go through.
	return ary([theta, phi], dtype=float)

sigmaxx = 3
sigmayy = 2	#the program rotates about the y axis by phi. Therefore if sigmaxx==sigmazz then we have isotropic stress in the zx plane; rotating about phi does not vary it.
sigmazz = 3	#
shear_xy= 3
shear_yz= 2
shear_zx= 1
#Number of steps to be shown:
numPoints=30	#for rotation around the chosen asix.
small_theta_lim = tau
phiLim = pi	#rotates the axis around half of the greater circle.
numLines = 500
windingNum=4
selectedLines=range(numLines)
#selectedLines=[399,398]
#Make the original stress state.
sigma = np.zeros([3,3])
sigma[0][0] = sigmaxx
sigma[1][1] = sigmayy
sigma[2][2] = sigmazz
sigma[0][1] = shear_xy
sigma[1][0] = shear_xy
sigma[0][2] = shear_zx
sigma[2][0] = shear_zx
sigma[1][2] = shear_yz
sigma[2][1] = shear_yz

fig = plt.figure()
ax = fig.gca(projection='3d')
ColorMixer1 = np.transpose(	#Make the color vary
	ary([np.linspace(0,0,numLines), #Starts off as red
	np.linspace(1,0,numLines),
	np.linspace(0,1,numLines),	#ends up as blue
	np.linspace(0.8,0.8,numLines)])
	)
ColorMixer2 = np.transpose(	#Make the color vary
	ary([np.linspace(1,0,numLines), #start off as orange
	np.linspace(1,1,numLines),	#ends up as green
	np.linspace(0,1,numLines),	
	np.linspace(0.8,0.8,numLines)])	#transparency =.8
	)

THETA,PHI = sphereCoveringCurve(windingNum,numLines)
if True:
	theta=np.linspace(0,small_theta_lim, numPoints)

	for axis in range(numLines):
		xStressList,yStressList,zStressList = [],[],[]
		for n in range(numPoints):
			R = RotateSigma(sigma, theta[n], THETA[axis], PHI[axis])
			xStressList.append(R[0,0])
			yStressList.append(R[1,1])
			zStressList.append(R[2,2])
		if axis in selectedLines: ax.plot(xStressList, yStressList,zStressList,	c=ColorMixer1[axis])

	for axis in range(numLines):
		shearxyList,shearyzList,shearzxList = [],[],[]
		for n in range(numPoints):
			[[xStress,shear1,shear3],[shear2,yStress,shear5],[shear4,shear6,zStress]] = RotateSigma(sigma, theta[n], THETA[axis], PHI[axis])
			if (shear1-shear2)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			if (shear3-shear4)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			if (shear5-shear6)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			shearxyList.append(shear1)
			shearyzList.append(shear5)
			shearzxList.append(shear3)
		if axis in selectedLines: ax.plot(shearxyList,shearyzList,shearzxList,	c=ColorMixer2[axis])
if False:
	for step in range(numLines):
		x,y,z = 3*ary(spherical_cartesian(THETA[step],PHI[step]))
		ax.scatter(x,y,z, c = ColorMixer2[step])
ax.set_xlim([-3,3])
ax.set_xlabel("x-stress")
ax.set_ylim([-3,3])
ax.set_ylabel("y-stress")
ax.set_zlim([-3,3])
ax.set_zlabel("z-stress")
# ax.set_aspect(1)
plt.show()

'''	#discarded code
if False:#loop to calculate the principle stresses:
	for phi in np.linspace(0, phiLim, numPoints):
		xStressList,yStressList,zStressList = [],[],[]
		for theta in np.linspace(0,tau,numPoints):
			[[xStress,shear1,shear3],[shear2,yStress,shear5],[shear4,shear6,zStress]] = RotateSigma(sigma, theta, phi)
			xStressList.append(xStress)
			yStressList.append(yStress)
			zStressList.append(zStress)
		xStressList.append(xStressList[0])
		yStressList.append(yStressList[0])
		zStressList.append(zStressList[0])
		#ax.plot(xStressList,yStressList,zStressList,	c=ColorMixer1[int(numPoints*phi/phiLim -1)])

if False:#loop to calculate teh shear stresses		
	for phi in np.linspace(0, phiLim, numPoints):
		shearxyList,shearyzList,shearzxList = [],[],[]
		for theta in np.linspace(0,tau,numPoints):
			[[xStress,shear1,shear3],[shear2,yStress,shear5],[shear4,shear6,zStress]] = RotateSigma(sigma, theta, phi)
			if (shear1-shear2)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			if (shear3-shear4)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			if (shear5-shear6)>1e-14: print( FloatingPointError("sigmaxy = ",shear1," is not equal to shearyx = ",shear2))
			shearxyList.append(shear1)
			shearyzList.append(shear5)
			shearzxList.append(shear3)
		shearxyList.append(shearxyList[0])
		shearyzList.append(shearyzList[0])
		shearzxList.append(shearzxList[0])
		ax.plot(shearxyList,shearyzList,shearzxList,	c=ColorMixer1[int(numPoints*phi/phiLim -1)])
'''