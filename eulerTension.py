#Calculates tension in CNT using Euler Beam analtyic expression with
#constant T and recursively approximating T
# Note TPrime_0 is residual, built in CNT tension (due to fabrication), and
# is given such that T_0 = TPrime_0 * E*momentInertia/L^2

# To Do: 
#  - Fix K_elec component for high dangerZ (TODO: This was a MATLAB error, see if python shares it)
#   ie, for xi=7.95e9, L=1e-6,
#   (sinh(xi*L)/(cosh(xi*L)-1))*(cosh(xi*myz)-1)-sinh(xi*myz)+xi*myz-xi*myz^2/L
#   throws NaN

import numpy as np
import math

# Load physical, CNT, and SMM constants
import myConstants as mc

def eulerTension(L=None,diameter=None,Q=None,Vg=None,Vbias=None,h=None,cRatio=None,dBdz=None,Sz=None,TPrime_0=None,zPrime_0=None,stepCount=None,*args,**kwargs):
	# Numeric Settings
	spatialElements=100000
	epsilon=1e-10
	xiL_UpperLimit=25
	xiL_LowerLimit=0.03

	# Read in keyword inputs
	effectiveCapLengthRatio=1
	for key, value in kwargs.items():
		if "effCapLenRatio" in key:
			effectiveCapLengthRatio = value
		else:
			raise ValueError('An unrecognized key was given %s'.format(key))

	# Make sure inputs are floats where needed (note Python3 gets rid of this need...)
	L = float(L)
	diameter = float(diameter)
	Q = float(Q)
	Vg = float(Vg)
	Vbias = float(Vbias)
	h = float(h)
	cRatio = float(cRatio)
	dBdz = float(dBdz)
	TPrime_0 = float(TPrime_0)

	# Calculated parameters
	rOut=diameter/2
	rIn=rOut - mc.wallThickness
	mCNT=mc.rhoA * 0.735 * math.pi * diameter * L
	A=math.pi * (rOut ** 2)
	momentInertia=(math.pi / 4) * (rOut ** 4)
	lengthDensity=math.pi * diameter * mc.rhoA
	Cg=4 * math.pi * mc.epsilon0 * L * effectiveCapLengthRatio / (2 * np.log(2 * h / rOut))
	C=cRatio * Cg
	CL=(C - Cg) / 2
	a=zPrime_0 * L # The position of the point force
	K_electric=(1.0 / (4 * math.pi * mc.epsilon0 * L ** 2 * h))*((1 / cRatio) ** 2.0) * (C * Vg - CL * Vbias) ** 2
	F_mag=mc.g * mc.muB * dBdz * Sz
	T_0=TPrime_0 * mc.E * momentInertia / L ** 2

	# If initial T estimate would give xi < 1, use a low T limit starting guess
	TInit=T_0 + (mc.E * A / 24 * (K_electric ** 2 * L ** 2 + 3 * K_electric * L * F_mag + 3 * F_mag ** 2)) ** (1./3)
	if abs(np.sqrt(TInit / (mc.E * momentInertia)) * L) < 1:
		TInit=T_0 + K_electric ** 2 * L ** 6 * A / (60480 * mc.E * momentInertia ** 2) + F_mag ** 2 * L ** 4 * A / (30720 * mc.E * momentInertia ** 2)

	# Determine the end condition on T iteration
	threshold=abs(epsilon * TInit)

	# Generate posSpace and determine special location indices
	z=np.linspace(0,L,spatialElements)
	dz=z[2] - z[1]
	dxdz=np.zeros((spatialElements,1))
	aIndex=int(np.floor(spatialElements * zPrime_0))
	halfIndex=int(np.floor(spatialElements / 2))
	quarterIndex=int(np.floor(spatialElements / 4))
	threeQuarterIndex=spatialElements - quarterIndex

	# Initialize outputs
	T = np.zeros((stepCount,1))

	# Loop until error threshold or max iteration reached
	counter=1
	escapeIn=99999999
	for j in range(0,stepCount):
		print('Run {}'.format(j))
		xi = np.sqrt(TInit/(mc.E*momentInertia))
		dangerZ = xiL_UpperLimit / abs(xi) # position at which xiL issues could take effect

		# Calculate modeshape
		x=(K_electric * L / (2 * TInit * xi)) * ((np.sinh(xi * L) / (np.cosh(xi * L) - 1)) * (np.cosh(xi * z) - 1) - np.sinh(xi * z) + xi * z - xi * z ** 2 / L)
		sinhka=np.sinh(xi * a)
		sinhkL=np.sinh(xi * L)
		sinhkLa=np.sinh(xi * (L - a))
		coshka=np.cosh(xi * a)
		coshkL=np.cosh(xi * L)
		coshkLa=np.cosh(xi * (L - a))
		FPrime=F_mag / (mc.E * momentInertia)
		k=xi
		sigma1=FPrime * (sinhka - sinhkL + sinhkLa + a * k + k * (L - a) * coshkL - L * k * coshkLa)
		sigma2=L * k * sinhkL - 2 * coshkL + 2
		sigma3=FPrime * (coshkL - coshka + coshkLa - k * (L - a) * sinhkL - 1)
		c1=sigma3 / (k * sigma2)
		c2=sigma1 / (k * sigma2)
		c3=- sigma3 / (k ** 2 * sigma2)
		c4=- sigma1 / (k ** 3 * sigma2)
		z1=z[0:aIndex-1]
		z2=z[aIndex:]
		x[0:aIndex-1]=x[0:aIndex-1] + c1 * np.sinh(k * z1) / k ** 2 + c2 * np.cosh(k * z1) / k ** 2 + c3 * z1 + c4
		x[aIndex:]=x[aIndex:] + FPrime * np.sinh(k * (z2 - a)) / k ** 3 - FPrime * (z2 - a) / k ** 2 + c1 * np.sinh(k * z2) / k ** 2 + c2 * np.cosh(k * z2) / k ** 2 + c3 * z2 + c4

		# Determine slope
		dxdz[0]=(x[1] - x[0]) / dz
		for count in range(1,(spatialElements - 1)):
			dxdz[count]=(x[count + 1] - x[count - 1]) / (2 * dz)
		dxdz[spatialElements-1]=(x[spatialElements-1] - x[spatialElements - 2]) / dz
		
		# Calculate tension
		T[j]=T_0 + mc.E * A / (2) * sum(dxdz ** 2) / spatialElements

		# Convergence speed is increased by only partially updating tension guess
		TInit=(0.6 + 0.37 * j / stepCount) * TInit + (0.4 - 0.37 * j / stepCount) * T[j]
		T[j]=TInit

		# If we have reached error threshold, do a couple more loops then exit
		if j > 2 and abs(T[j] - T[j - 1]) < threshold and escapeIn > 10:
			escapeIn=3
		if escapeIn < 1:
			break
		# If the tension is going not moving, average the last few runs and restart
		if counter > 25 and abs(T[j] - T[j - 24]) < (10000.0 * threshold):
			counter=1
			TInit=np.mean(T[j - 10:j])
			T[j]=TInit
		counter=counter + 1
		escapeIn=escapeIn - 1

	# Trim T to only include updated values
	T = T[0:counter]
	maxx=max(abs(x))
#	return T,maxx,x,z,dxdz,K_electric,F_mag
	return T, K_electric, F_mag

