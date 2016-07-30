# Class for solving euler beam equations designed to determine resonant frequencies of CNT resonators with both a point and uniform applied force

import numpy as np
import math

# Load physical, CNT, and SMM constants
import myConstants as mc

class eulerBeamSolver:
	def __init__(self,L,diameter,Q,Vg,VgAC,Vbias,h,cRatio,dBdz,Sz,T_0,pos):
		# The inputs can be singletons or tuples of values
		self.L = L
		self.diameter = diameter
		self.Q = Q
		self.Vg = Vg
		self.VgAC = VgAC
		self.Vbias = Vbias
		self.h = h
		self.cRatio = cRatio
		self.dBdz = dBdz
		self.Sz = Sz
		self.T_0 = T_0
		self.pos = pos
		
		# Default run parameters
		self.fStepCount = 10000
		self.TStepCount = 500
		
	def calcParameters(self):
		self.rOut = self.diameter/2
		self.rIn=rOut - mc.wallThickness
		self.area=math.pi * (rOut ** 2)
		self.mCNT=mc.rhoA * 0.735 * math.pi * self.diameter * self.L
		self.momentInertia=(math.pi / 4) * (rOut ** 4)
		self.lengthDensity=math.pi * diameter * mc.rhoA
		self.Cg=4 * math.pi * mc.epsilon0 * self.L * effectiveCapLengthRatio / (2 * np.log(2 * self.h / self.rOut))
		self.C = self.cRatio * self.Cg
		self.CL=(C - Cg) / 2
		a=zPrime_0 * L # The position of the point force
		K_electric=(1.0 / (4 * math.pi * mc.epsilon0 * L ** 2 * h))*((1 / cRatio) ** 2.0) * (C * Vg - CL * Vbias) ** 2
		F_mag=mc.g * mc.muB * dBdz * Sz
		
	def eulerFreq(self):		
		# Run parameters
		epsilon=0
		escapeIn = -9999
		spatialElements=1000000
		z=np.linspace(0,1,spatialElements)
		
		# Begin calculation
		freqStart=self.quickFreq()
	#    threshold=epsilon * freqStart
		threshold=1
		slowFactor=1 # A scaling factor that will slow down corrections when the solution starts to oscillate
		errorSlow=1
		slowFactorCountdown=2000 # Number of loops before slowFactor can be used, this is reset in the slowFactor region to slowing "period"
		slowFactorPeriod=2000
		
		# Scale the variables (as done in thesis)
		nu=(1 / L ** 2) * np.sqrt(mc.E * momentInertia / (mCNT / L))
		omegaPrimeInit=freqStart / nu
		TPrime=T * L ** 2 / (mc.E * momentInertia)
		kPrime=np.sqrt(TPrime)
		lDcPrime=K_elec * L ** 4 / (mc.E * momentInertia * rOut)
		lMagPrime=F_mag * L ** 3 / (mc.E * momentInertia * rOut)
		
		# Store some repeatedly used calculations
		a=pos
		k=kPrime
		sinhka=np.sinh(k * a)
		sinhk=np.sinh(k)
		sinhkLa=np.sinh(k * (1 - a))
		coshka=np.cosh(k * a)
		coshk=np.cosh(k)
		coshkLa=np.cosh(k * (1 - a))
		sigma1=(coshk - coshka + coshkLa - k * (1 - a) * sinhk - 1)
		sigma2=(sinhka - sinhk + sinhkLa + a * k + k * (1 - a) * coshk - k * coshkLa)
		sigma3=k * sinhk - 2 * coshk + 2
		
		# Compute the DC function derivatives needed later
		d2udz0=((lDcPrime / (2 * kPrime ** 2)) * ((kPrime * np.sinh(kPrime)) / (np.cosh(kPrime) - 1) - 2) + lMagPrime * sigma2 / (kPrime * sigma3))
		d2udz1=((lDcPrime / (2 * kPrime ** 2)) * ((kPrime * np.sinh(kPrime)) / (np.cosh(kPrime) - 1) - 2) + lMagPrime / (kPrime * sigma3) * (sigma1 * sinhk + sigma2 * coshk) + lMagPrime / kPrime * sinhkLa)
		d3udz0=- lDcPrime / 2 + lMagPrime * sigma1 / sigma3
		d3udz1=lDcPrime / 2 + lMagPrime / sigma3 * (sigma1 * coshk + sigma2 * sinhk) + lMagPrime * coshkLa
		d3udzAdiscont=lMagPrime
		d5udzAdiscont=kPrime ** 2 * lMagPrime

		# Initialize storage and counter
		resFreq = np.zeros((stepCount,1))
		counter = 0
		
		for j in range(0,stepCount):
			print('Run {}'.format(j))
	#		if kPrime > 19:
	#			resFreq[0]=0
	#			break
			kPlus=1 / 2 * np.sqrt(-2 * TPrime + 2 * np.sqrt(TPrime ** 2 + 4 * omegaPrimeInit ** 2))
			kMinus=1 / 2 * np.sqrt(2 * TPrime + 2 * np.sqrt(TPrime ** 2 + 4 * omegaPrimeInit ** 2))
			cospa=np.cos(kPlus * a)
			sinpa=np.sin(kPlus * a)
			cospLa=np.cos(kPlus * (1 - a))
			sinpLa=np.sin(kPlus * (1 - a))
			coshma=np.cosh(kMinus * a)
			sinhma=np.sinh(kMinus * a)
			coshmLa=np.cosh(kMinus * (1 - a))
			sinhmLa=np.sinh(kMinus * (1 - a))
			TacC1=TacC1_calc(kPrime,kPlus,lDcPrime,lMagPrime,a)
			TacC2=TacC2_calc(kPrime,kPlus,lDcPrime,lMagPrime,a)
			TacC3=TacC3_calc(kPrime,kMinus,lDcPrime,lMagPrime,a)
			TacC4=TacC4_calc(kPrime,kMinus,lDcPrime,lMagPrime,a)
			TacB1=TacB1_calc(kPrime,kPlus,lDcPrime,lMagPrime,a)
			TacB2=TacB2_calc(kPrime,kPlus,lDcPrime,lMagPrime,a)
			TacB3=TacB3_calc(kPrime,kMinus,lDcPrime,lMagPrime,a)
			TacB4=TacB4_calc(kPrime,kMinus,lDcPrime,lMagPrime,a)
			TacB5=TacE_calc(kPrime,lDcPrime,lMagPrime,a)
			line1=np.array([1,0,1,0,0,0,0,0,d2udz0])
			line2=np.array([0,kPlus,0,kMinus,0,0,0,0,d3udz0])
			line3=np.array([0,0,0,0,1,0,1,0,d2udz1])
			line4=np.array([0,0,0,0,0,- kPlus,0,- kMinus,d3udz1])
			line5=np.array([cospa,sinpa,coshma,sinhma,- cospLa,- sinpLa,- coshmLa,- sinhmLa,0])
			line6=np.array([- kPlus * sinpa,kPlus * cospa,kMinus * sinhma,kMinus * coshma,- kPlus * sinpLa,kPlus * cospLa,kMinus * sinhmLa,kMinus * coshmLa,- d3udzAdiscont])
			line7=np.array([- kPlus ** 2 * cospa,- kPlus ** 2 * sinpa,kMinus ** 2 * coshma,kMinus ** 2 * sinhma,kPlus ** 2 * cospLa,kPlus ** 2 * sinpLa,- kMinus ** 2 * coshmLa,- kMinus ** 2 * sinhmLa,0])
			line8=np.array([kPlus ** 3 * sinpa,- kPlus ** 3 * cospa,kMinus ** 3 * sinhma,kMinus ** 3 * coshma,kPlus ** 3 * sinpLa,- kPlus ** 3 * cospLa,kMinus ** 3 * sinhmLa,kMinus ** 3 * coshmLa,-d5udzAdiscont])
			line9=np.array([TacC1,TacC2,TacC3,TacC4,TacB1,TacB2,TacB3,TacB4,TacB5 + omegaPrimeInit ** 2])
			solverMatrix=np.vstack((line1,line2,line3,line4,line5,line6,line7,line8,line9))
			soln=np.linalg.det(solverMatrix)
			if j == 0:
				initSoln=(soln)
			if np.abs(soln) < epsilon:
				resFreq[j]=omegaPrimeInit * nu
				break
			if slowFactorCountdown < 1 and 0.95 * np.mean(np.abs(resFreq[j - 50:j - 45] - resFreq[j - 1])) < np.mean(np.abs(resFreq[j - 25:j - 20] - resFreq[j - 1])):
				slowFactor=slowFactor * 0.1
				slowFactorCountdown=slowFactorPeriod
			
			if np.abs(TPrime) > 12:
				resFreq[j]=omegaPrimeInit * nu * (1 + (0.05 - 0.03 * j / stepCount) * 5.0 * slowFactor * errorSlow * soln / (np.abs(initSoln) ** (12 / 11)))
			elif np.abs(TPrime) > 5:
				resFreq[j]=omegaPrimeInit * nu * (1 + (0.05 - 0.03 * j / stepCount) * 0.5 * slowFactor * errorSlow * soln / (np.abs(initSoln) ** (12 / 11)))
			elif np.abs(TPrime) > 1:
				resFreq[j]=omegaPrimeInit * nu * (1 + (0.05 - 0.03 * j / stepCount) * 0.05 * slowFactor * errorSlow * soln / (np.abs(initSoln) ** (12 / 11)))
			elif np.abs(TPrime) > 1e-05:
				resFreq[j]=omegaPrimeInit * nu * (1 + (0.05 - 0.03 * j / stepCount) * slowFactor * errorSlow * 0.001 * soln / max_(np.abs(initSoln),1))
			else:
				resFreq[j]=omegaPrimeInit * nu * (1 + (0.05 - 0.03 * j / stepCount) * slowFactor * errorSlow * 0.001 * soln / max_(np.abs(initSoln),1))
				
			if np.abs(np.log10((resFreq[j] / freqStart))) > 1:
				if j < 4:
					resFreq[j]=freqStart
				else:
					resFreq[j]=resFreq[j - 3]
				errorSlow=errorSlow * 0.1
			omegaPrimeInit=resFreq[j][0] / nu
			if j > 5 and np.mean(np.abs(resFreq[j - 4:j - 1] - resFreq[j])) < threshold and escapeIn < - 1:
				escapeIn=3
			if escapeIn < 1 and escapeIn > - 1:
				break

			counter = counter + 1
			escapeIn=escapeIn - 1
			slowFactorCountdown=slowFactorCountdown - 1

		# Trim the resFreq to just the updated results
		resFreq = resFreq[0:counter]

		return resFreq

	def quickFreq(self)
		tL = self.L # L is referenced a lot, so store it
		rOut=self.diameter / 2
		rIn=self.rOut - mc.wallThickness
		mCNT=mc.rhoA * 0.735 * math.pi * self.diameter * tL
		momentInertia=(math.pi / 4) * (self.rOut ** 4)
		epsilon=1e-12
		xi=np.sqrt(self.T / (mc.E * momentInertia))
		if xi * tL < 1:
			freqStart=np.sqrt(mc.E * momentInertia / (mCNT / tL)) * (22.38 / tL ** 2 + 0.28 * xi ** 2)
		else:
			freqStart=np.sqrt(mc.E * momentInertia / (mCNT / tL)) * (2 * math.pi / tL ** 2 + math.pi * xi / tL)
		threshold=epsilon * freqStart
		escapeIn=9999999
		counter = 0
		resFreq = np.zeros((stepCount,1))

		for k in range(0,stepCount):
			_lambda=np.sqrt(mCNT / (tL * mc.E * momentInertia)) * freqStart
			yPlus=(L / np.sqrt(2)) * np.sqrt(np.sqrt(xi ** 4 + 4 * _lambda ** 2) + xi ** 2)
			yMinus=(L / np.sqrt(2)) * np.sqrt(np.sqrt(xi ** 4 + 4 * _lambda ** 2) - xi ** 2)
			soln=np.cosh(yPlus) * np.cos(yMinus) - (yPlus ** 2 - yMinus ** 2) / (2 * yPlus * yMinus) * np.sinh(yPlus) * np.sin(yMinus)
			resFreq[k]=freqStart * (1 - (0.05 - 0.04 * k / stepCount) * (soln - 1) / np.exp(yPlus))
			freqStart=resFreq[k]
			if k > 2 and np.abs(resFreq[k] - resFreq[k - 1]) < threshold and escapeIn > 10:
				escapeIn=3
			if escapeIn < 1:
				break

			counter = counter + 1
			escapeIn=escapeIn - 1

		# Trim the result
		resFreq = resFreq[0:counter]	

		return resFreq

	def TacB4_calc(k=None,m=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		sigma10=k * np.sinh(k) - 2 * np.cosh(k) + 2
		sigma9=k ** 2 - m ** 2
		sigma8=np.cosh(m * (a - 1))
		sigma7=np.sinh(m * (a - 1))
		sigma6=k ** 2 * m * sigma10
		sigma5=k ** 2 * sigma9 * sigma10
		sigma3=k * (np.sinh(k) - np.sinh(a * k) * sigma8) + m * sigma7 * np.cosh(a * k)
		sigma4=k * (np.cosh(k) - np.cosh(a * k) * sigma8) + m * sigma7 * np.sinh(a * k)
		sigma2=np.sinh(k * (a - 1))
		sigma1=np.cosh(k * (a - 1))
		DCPart=sigma3 / (2 * k ** 2 * sigma9) - (a * m * sigma7 - sigma8 + 1) / (k ** 2 * m ** 2) + sigma7 / (2 * k ** 2 * m) - (np.sinh(k) * sigma4) / (2 * k ** 2 * (np.cosh(k) - 1) * sigma9)
		magLine1=(k * sigma2 - m * sigma7) / (k ** 2 * sigma9) - sigma7 / (k ** 2 * m) + sigma3 / sigma5 + sigma7 / sigma6 + (sigma1 * sigma4) / (k * sigma9 * sigma10) + (np.sinh(k) * sigma3) / (k * sigma9 * sigma10) - sigma1 * sigma3 / sigma5 + sigma2 * sigma4 / sigma5 + sigma7 * (np.cosh(a * k) - np.cosh(k)) / sigma6 + sigma7 * np.sinh(k) / (k * m * sigma10) - a * sigma4 / (k * sigma9 * sigma10) - sigma7 * sigma1 / sigma6 + np.cosh(a * k) * sigma3 / sigma5 - np.sinh(a * k) * sigma4 / sigma5 - np.cosh(k) * sigma4 / (k * sigma9 * sigma10) - np.cosh(k) * sigma3 / sigma5 + np.sinh(k) * sigma4 / sigma5
		magLine2=(a * np.cosh(k) * sigma4) / (k * sigma9 * sigma10) - (a * np.sinh(k) * sigma3) / (k * sigma9 * sigma10) - a * sigma7 * np.sinh(k) / (k * m * sigma10)
		TacB4=fDc * 4 * m * DCPart + fMag * 4 * m * (magLine1 + magLine2)
		return TacB4

	def TacB3_calc(k=None,m=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		s10=k * np.sinh(k) - 2 * np.cosh(k) + 2
		s9=k ** 2 - m ** 2
		s8=np.sinh(m * (a - 1))
		s7=np.cosh(m * (a - 1))
		s6=k ** 2 * m * s10
		s5=k ** 2 * s9 * s10
		s4=m * (np.cosh(k) - np.cosh(a * k) * s7) + k * s8 * np.sinh(a * k)
		s3=m * (np.sinh(k) - np.sinh(a * k) * s7) + k * s8 * np.cosh(a * k)
		s2=np.cosh(k * (a - 1))
		s1=(np.sinh(m * (a - 1) / 2)) ** 2
		DcPart=(s8 - m * (a * s7 - 1)) / (k ** 2 * m ** 2) + s1 / (k ** 2 * m) - s4 / (2 * k ** 2 * s9) + np.sinh(k) * s3 / (2 * k ** 2 * (np.cosh(k) - 1) * s9)
		magLine1=2 * s1 / (k ** 2 * m) + (s4 - s2 * s4 + np.sinh(k) * s3 + np.sinh(k * (a - 1)) * s3 + np.cosh(a * k) * s4 - np.sinh(a * k) * s3 - np.cosh(k) * s4) / s5 - m * (s2 - s7) / (k ** 2 * s9) + (- 2 * s1 - s1 * 2 * np.cosh(a * k) + s1 * np.cosh(k) * 2 + 2 * s1 * s2) / s6
		magLine2=(s2 * s3 - a * s3 - np.cosh(k) * s3 + np.sinh(k) * s4 + a * np.cosh(k) * s3 - a * np.sinh(k) * s4) / (k * s9 * s10) + (- 2 * s1 * np.sinh(k) + a * 2 * s1 * np.sinh(k)) / (k * m * s10)
		TacB3=- fDc * 4 * m * DcPart + fMag * 4 * m * (magLine1 + magLine2)
		return TacB3

	def TacB2_calc(k=None,p=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		C=fDc
		M=fMag
		s24=np.cos(p * (1 - a))
		s23=np.sin(p * (1 - a))
		s22=8 * C * k ** 2 * s24
		s21=8 * C * p ** 2 * s24
		s20=8 * C * a * k ** 2 * p * s23
		s19=4 * M * a * k * p ** 3 * s23
		s18=8 * C * a * p ** 3 * s23
		s17=2 * C * k * p ** 3 * s23
		s16=4 * C * k ** 2 * p * s23
		s15=4 * C * k * p ** 2 * s24
		s14=4 * M * k * p ** 2 * s24
		s13=4 * C * p ** 3 * s23
		s12=4 * M * k ** 2 * p * s23
		s11=8 * C * k ** 2
		s10=8 * C * p ** 2
		s9=2 * k ** 4 * p
		s8=k - 2 * a * k
		s7=4 * M * k * p ** 2
		s6=2 * C * k ** 2 * p ** 2
		s5=4 * M * k ** 2 * p ** 2
		s4=2 * k ** 2 * p ** 3
		s3=4 * M * a * k ** 2 * p ** 2
		s2=np.cosh(k * (1 - a))
		s1=np.sinh(k * (1 - a))
		denom=(- s9 - s4) * np.cosh(k) + (k ** 5 * p + k ** 3 * p ** 3) * np.sinh(k) + s9 + s4
		Ccoshk=s22 - s10 - s6 - s11 + s21 + s13 - s18 + s16 - s12 - s3 - 2 * M * k ** 2 * p ** 2 * s24 - s20
		Csinhk=4 * C * k ** 3 * (1 - s24) + 8 * C * k * p ** 2 + s7 - s15 + s14 - s17 - 2 * C * k ** 3 * p * s23 - 2 * M * k * p ** 3 * s23 + 4 * C * a * k * p ** 3 * s23 + 4 * C * a * k ** 3 * p * s23 + s19 + 4 * M * a * k ** 3 * p * s23
		Ccoshak=s5 + s13 - s12 + 2 * C * k ** 2 * p ** 2 * s24 + 4 * M * a * k ** 2 * p ** 2 * s24
		Csinhak=- s7 - s15 - s14 - s17 - s19
		line01=s11 + s10 - s6 - s5 - s22 - s21 - s13 - 4 * C * p ** 3 * s2 * s23 + s18 - s16 - 4 * M * k * p ** 2 * s1 + s12 + s3 + 4 * M * k ** 2 * p ** 2 * s2 * s24 - 2 * M * k ** 2 * p ** 2 * np.cosh(s8) * s24 - 4 * C * k * p ** 2 * s24 * s1 + s20 - 4 * M * k * p ** 2 * s24 * s1
		line2=4 * M * k ** 2 * p * s2 * s23 + 2 * C * k * p ** 3 * s1 * s23 + 4 * M * k * p ** 3 * s1 * s23 - 2 * M * k * p ** 3 * np.sinh(s8) * s23 + 2 * C * k ** 2 * p ** 2 * s2 * s24 - 4 * M * a * k * p ** 3 * s1 * s23 - 4 * M * a * k ** 2 * p ** 2 * s2 * s24
		TacB2=(Ccoshk * np.cosh(k) + Csinhk * np.sinh(k) + Ccoshak * np.cosh(a * k) + Csinhak * np.sinh(a * k) + line01 + line2) / denom
		return TacB2

	def TacB1_calc(k=None,p=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		C=fDc
		M=fMag
		s20=np.sin(p * (1 - a))
		s19=np.cos(p * (1 - a))
		s18=8 * C * k ** 2 * s20
		s17=8 * C * p ** 2 * s20
		s16=8 * C * a * p ** 3 * s19
		s15=2 * C * k * p ** 3 * s19
		s14=4 * C * k ** 2 * p * s19
		s13=4 * C * k * p ** 2 * s20
		s12=4 * M * k * p ** 2 * s20
		s11=8 * C * a * k ** 2 * p * s19
		s10=4 * M * a * k * p ** 3 * s19
		s9=4 * C * p ** 3 * s19
		s8=4 * M * k ** 2 * p * s19
		s7=2 * k ** 4 * p
		s6=k - 2 * a * k
		s5=4 * C * k ** 2 * p
		s4=2 * k ** 2 * p ** 3
		s3=4 * M * k ** 2 * p
		s2=np.sinh(k * (1 - a))
		s1=np.cosh(k * (1 - a))
		denom=(- s7 - s4) * np.cosh(k) + (k ** 5 * p + k ** 3 * p ** 3) * np.sinh(k) + s7 + s4
		Ccoshk=s9 + s5 - s18 - s17 + s3 - s16 + s14 - s8 + 2 * M * k ** 2 * p ** 2 * s20 - s11
		Csinhk=4 * C * k ** 3 * s20 - 2 * C * k ** 3 * p - s15 - 2 * C * k ** 3 * p * s19 - 4 * M * a * k ** 3 * p - 2 * M * k * p ** 3 * s19 + s13 - s12 + 4 * C * a * k * p ** 3 * s19 + 4 * C * a * k ** 3 * p * s19 + s10 + 4 * M * a * k ** 3 * p * s19
		Ccoshak=s9 + s3 - s8 - 2 * C * k ** 2 * p ** 2 * s20 - 4 * M * a * k ** 2 * p ** 2 * s20
		Csinhak=s13 - s15 + s12 - s10
		line0=s18 - s5 - s9 + s17 - s3
		line1=- 4 * C * p ** 3 * s1 * s19 + s16 - s14 - 4 * M * k ** 2 * p * s1 + s8 - 2 * C * k ** 2 * p ** 2 * s1 * s20 - 4 * M * k ** 2 * p ** 2 * s1 * s20 + 2 * M * k ** 2 * p ** 2 * np.cosh(s6) * s20 + s11 + 4 * M * k ** 2 * p * s1 * s19 + 2 * C * k * p ** 3 * s19 * s2 + 4 * M * k * p ** 3 * s19 * s2
		line2=- 2 * M * k * p ** 3 * s19 * np.sinh(s6) + 4 * C * k * p ** 2 * s2 * s20 + 4 * M * k * p ** 2 * s2 * s20 - 4 * M * a * k * p ** 3 * s19 * s2 + 4 * M * a * k ** 2 * p ** 2 * s1 * s20
		TacB1=(Ccoshk * np.cosh(k) + Csinhk * np.sinh(k) + Ccoshak * np.cosh(a * k) + Csinhak * np.sinh(a * k) + line0 + line1 + line2) / denom
		return TacB1

	def TacC4_calc(k=None,m=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		s7=k * np.sinh(k) - 2 * np.cosh(k) + 2
		s6=k ** 2 - m ** 2
		s5=k ** 2 * m * s7
		s4=k ** 2 * s6 * s7
		s3=np.cosh(k * (a - 1))
		s2=k * np.cosh(a * m) * np.sinh(a * k) - m * np.cosh(a * k) * np.sinh(a * m)
		s1=k * (np.cosh(a * k) * np.cosh(a * m) - 1) - m * np.sinh(a * k) * np.sinh(a * m)
		DcPart=(a * m * np.sinh(a * m) - np.cosh(a * m) + 1) / (k ** 2 * m ** 2) + s2 / (2 * k ** 2 * s6) - np.sinh(a * m) / (2 * k ** 2 * m) - np.sinh(k) * s1 / (2 * k ** 2 * (np.cosh(k) - 1) * s6)
		magLine1=(- np.sinh(a * m) + np.sinh(a * m) * s3 - np.cosh(a * k) * np.sinh(a * m) + np.sinh(a * m) * np.cosh(k)) / s5 + (s2 - np.sinh(a * k) * s1 + np.cosh(a * k) * s2 + np.sinh(k) * s1 + np.sinh(k * (a - 1)) * s1 - np.cosh(k) * s2 - s2 * s3) / s4
		magLine2=(- np.sinh(a * m) * np.sinh(k) + a * np.sinh(a * m) * np.sinh(k)) / (k * m * s7) + (- np.cosh(k) * s1 + s3 * s1 + np.sinh(k) * s2 - a * s1 - a * np.sinh(k) * s2 + a * np.cosh(k) * s1) / (k * s6 * s7)
		TacC4=fDc * (- 4 * m) * DcPart + fMag * (- 4 * m) * (magLine1 + magLine2)
		return TacC4

	def TacC3_calc(k=None,m=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		s8=k * np.sinh(k) - 2 * np.cosh(k) + 2
		s7=k ** 2 - m ** 2
		s6=k ** 2 * m * s8
		s5=k ** 2 * s7 * s8
		s4=np.cosh(k * (a - 1))
		s3=(np.sinh(a * m / 2)) ** 2
		s2=k * np.cosh(a * k) * np.sinh(a * m) - m * np.cosh(a * m) * np.sinh(a * k)
		s1=m * (np.cosh(a * k) * np.cosh(a * m) - 1) - k * np.sinh(a * k) * np.sinh(a * m)
		DcPart=(np.sinh(a * m) - a * m * np.cosh(a * m)) / (k ** 2 * m ** 2) + s1 / (2 * k ** 2 * s7) + s3 / (k ** 2 * m) + np.sinh(k) * s2 / (2 * k ** 2 * (np.cosh(k) - 1) * s7)
		magLine1=(- s1 - np.cosh(a * k) * s1 + np.cosh(k) * s1 - np.sinh(a * k) * s2 + s4 * s1 + np.sinh(k) * s2 + np.sinh(k * (a - 1)) * s2) / s5 + (- s3 * 2 - np.cosh(a * k) * 2 * s3 + s3 * np.cosh(k) * 2 + s3 * s4 * 2) / s6
		magLine2=(- a * s2 - np.sinh(k) * s1 - np.cosh(k) * s2 + s4 * s2 + a * np.sinh(k) * s1 + a * np.cosh(k) * s2) / (k * s7 * s8) + (- s3 * 2 * np.sinh(k) + a * s3 * 2 * np.sinh(k)) / (k * m * s8)
		TacC3=fDc * 4 * m * DcPart + fMag * (- 4 * m) * (magLine1 + magLine2)
		return TacC3

	def TacC2_calc(k=None,p=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		C=fDc
		M=fMag
		s21=np.sinh(k * (1 - a))
		s20=np.cosh(k * (1 - a))
		s19=4 * M * k * p ** 2 * s21
		s18=4 * M * k ** 2 * p ** 2 * s20
		s17=4 * C * k ** 3
		s16=2 * k ** 4 * p
		s15=k - 2 * a * k
		s14=8 * C * a * p ** 3
		s13=2 * C * k * p ** 3
		s12=4 * C * k * p ** 2
		s11=4 * C * k ** 2 * p
		s10=8 * C * a * k ** 2 * p
		s9=4 * M * a * k * p ** 3
		s8=2 * k ** 2 * p ** 3
		s7=4 * C * p ** 3
		s6=4 * M * k ** 2 * p
		s5=8 * C * k ** 2
		s4=8 * C * p ** 2
		s3=4 * M * k * p ** 2
		s2=2 * C * k ** 2 * p ** 2
		s1=4 * M * a * k ** 2 * p ** 2
		line1=(s1 + s2) * np.cos(a * p) * np.cosh(a * k) + (s3 - s12 - s17) * np.cos(a * p) * np.sinh(k) + (- s12 - s3) * np.cos(a * p) * np.sinh(a * k) + (- 2 * M * k ** 2 * p ** 2 + s5 + s4) * np.cos(a * p) * np.cosh(k)
		line2=(2 * C * k ** 2 * p ** 2 * s20 - s4 - 4 * C * k * p ** 2 * s21 - s19 - s5 + s18 - 2 * M * k ** 2 * p ** 2 * np.cosh(s15) - 4 * M * a * k ** 2 * p ** 2 * s20) * np.cos(a * p) + (s6 - s7) * np.sin(a * p) * np.cosh(a * k)
		line3=(s13 + 2 * C * k ** 3 * p + 2 * M * k * p ** 3 + 4 * M * k ** 3 * p - 4 * C * a * k * p ** 3 - 4 * C * a * k ** 3 * p - s9 - 4 * M * a * k ** 3 * p) * np.sin(a * p) * np.sinh(k) + (s13 + s9) * np.sin(a * p) * np.sinh(a * k)
		line4=(s14 - s7 - s11 - s6 + s10) * np.sin(a * p) * np.cosh(k) + (s7 + 4 * C * p ** 3 * s20 - s14 + s11 + s6 - s10 - 4 * M * k ** 2 * p * s20 - 2 * C * k * p ** 3 * s21 - 4 * M * k * p ** 3 * s21 + 2 * M * k * p ** 3 * np.sinh(s15) + 4 * M * a * k * p ** 3 * s21) * np.sin(a * p)
		line5=(s17 + 8 * C * k * p ** 2 + s3) * np.sinh(k) - s3 * np.sinh(a * k) + (s1 - s4 - s2 - 4 * M * k ** 2 * p ** 2 - s5) * np.cosh(k) + s5 + s4 - s2 - s19 + s18 - s1
		denom=(k ** 5 * p + k ** 3 * p ** 3) * np.sinh(k) + (- s16 - s8) * np.cosh(k) + s16 + s8
		TacC2=(line1 + line2 + line3 + line4 + line5) / denom
		return TacC2

	def TacC1_calc(k=None,p=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		s14=np.sinh(k * (a - 1))
		s13=2 * k * p ** 2 * np.sin(a * p) * s14
		s12=k * (2 * a - 1)
		s11=2 * k ** 2 * p
		s10=2 * k ** 2 * p * np.cosh(k)
		s9=2 * k ** 2 * p * np.cos(a * p)
		s8=k * p ** 3 * np.cos(a * p) * np.sinh(k)
		s7=2 * k ** 2 * p * np.cos(a * p) * np.cosh(k)
		s6=2 * k * p ** 2 * np.sin(a * p) * np.sinh(k)
		s5=2 * a * k * p ** 3 * np.cos(a * p) * np.sinh(k)
		s4=2 * a * k ** 3 * p * np.cos(a * p) * np.sinh(k)
		s3=2 * k * p ** 2 * np.sinh(a * k) * np.sin(a * p)
		s2=k ** 2 * p * (k ** 2 + p ** 2) * (k * np.sinh(k) - 2 * np.cosh(k) + 2)
		s1=np.cosh(k * (a - 1))
		dcLine1=2 * p ** 3 * np.cos(a * p) - s11 + 4 * k ** 2 * np.sin(a * p) + 4 * p ** 2 * np.sin(a * p) - 4 * a * p ** 3 * np.cos(a * p) + s9 - 2 * p ** 3 * np.cosh(a * k) * np.cos(a * p) + s10 - k ** 3 * p * np.sinh(k) - 2 * p ** 3 * np.cos(a * p) * np.cosh(k) - 4 * k ** 2 * np.sin(a * p) * np.cosh(k) - 4 * p ** 2 * np.sin(a * p) * np.cosh(k) + 2 * k ** 3 * np.sin(a * p) * np.sinh(k) + 2 * p ** 3 * np.cos(a * p) * s1
		dcLine2=s8 + k ** 3 * p * np.cos(a * p) * np.sinh(k) + s6 + k * p ** 3 * np.cos(a * p) * s14 - k ** 2 * p ** 2 * np.cosh(a * k) * np.sin(a * p) - s13 - k ** 2 * p ** 2 * np.sin(a * p) * s1 - 4 * a * k ** 2 * p * np.cos(a * p) + k * p ** 3 * np.cos(a * p) * np.sinh(a * k) + s3 + 4 * a * p ** 3 * np.cos(a * p) * np.cosh(k) - s7 + 4 * a * k ** 2 * p * np.cos(a * p) * np.cosh(k) - s5 - s4
		magLine1=s11 + 2 * k ** 2 * p * np.cosh(a * k) - s9 - s10 + 2 * k ** 3 * p * np.sinh(k) - 2 * k ** 2 * p * s1 - s8 - 2 * k ** 3 * p * np.cos(a * p) * np.sinh(k) + 2 * k ** 2 * p * np.cos(a * p) * s1 + s6 - 2 * k * p ** 3 * np.cos(a * p) * s14 + s13 - k ** 2 * p ** 2 * np.sin(a * p) * np.cosh(k) + k * p ** 3 * np.sinh(s12) * np.cos(a * p)
		magLine2=2 * k ** 2 * p ** 2 * np.sin(a * p) * s1 - 2 * k ** 2 * p * np.cosh(a * k) * np.cos(a * p) - 2 * a * k ** 3 * p * np.sinh(k) - s3 - k ** 2 * p ** 2 * np.cosh(s12) * np.sin(a * p) + s7 - 2 * a * k * p ** 3 * np.cos(a * p) * np.sinh(a * k) + s5 + s4 + 2 * a * k * p ** 3 * np.cos(a * p) * s14 + 2 * a * k ** 2 * p ** 2 * np.cosh(a * k) * np.sin(a * p) - 2 * a * k ** 2 * p ** 2 * np.sin(a * p) * s1
		TacC1=fDc * 2 * (dcLine1 + dcLine2) / s2 - fMag * 2 * (magLine1 + magLine2) / s2
		return TacC1

	def TacE_calc(k=None,fDc=None,fMag=None,a=None,*args,**kwargs):
		C=fDc
		M=fMag
		TacE=(- (32 * np.cosh(k) - 8 * np.cosh(2 * k) + 10 * k * np.sinh(2 * k) - 3 * k ** 3 * np.sinh(k) - 4 * k ** 2 * np.cosh(2 * k) + (k ** 3 * np.sinh(2 * k)) / 2 - 20 * k * np.sinh(k) + k ** 4 + 4 * k ** 2 * np.cosh(k) + k ** 4 * np.cosh(k) - 24) / (k ** 4 * (np.cosh(k) - 1) * (4 * np.cosh(k) - 4 * k * np.sinh(k) + k ** 2 + k ** 2 * np.cosh(k) - 4))) * C ** 2 + (- (4 * k * np.sinh(2 * k) - 6 * k ** 3 * np.sinh(k) - 8 * k ** 2 * np.cosh(k * (a - 1)) + 4 * k ** 2 * np.cosh(k * (a - 2)) - 6 * k ** 3 * np.sinh(k * (a - 1)) + 2 * k ** 3 * np.sinh(k * (a + 1)) + 2 * k ** 3 * np.sinh(k * (a - 2)) - (k ** 3 * np.sinh(2 * k * (a - 1))) / 2 - 4 * k ** 2 * np.cosh(2 * k) + (3 * k ** 3 * np.sinh(2 * k)) / 2 + 12 * k * np.sinh(a * k) + 2 * k ** 3 * np.sinh(k * (2 * a - 1)) - 4 * a * k ** 4 - 8 * k * np.sinh(k) + 4 * k ** 2 * np.cosh(a * k) - 4 * k ** 4 * np.cosh(a * k) + k ** 4 * np.cosh(2 * a * k) - 12 * k * np.sinh(k * (a - 1)) - 4 * k * np.sinh(k * (a + 1)) + 4 * k * np.sinh(k * (a - 2)) + 2 * k ** 3 * np.sinh(a * k) - (3 * k ** 3 * np.sinh(2 * a * k)) / 2 - 4 * k ** 2 + 3 * k ** 4 + 4 * a ** 2 * k ** 4 + 8 * k ** 2 * np.cosh(k) - 12 * a * k ** 2 * np.cosh(a * k) + 6 * a * k ** 4 * np.cosh(a * k) - a * k ** 4 * np.cosh(2 * a * k) + 2 * a ** 2 * k ** 3 * np.sinh(2 * k) - 8 * a * k ** 3 * np.sinh(a * k) + 2 * a * k ** 3 * np.sinh(2 * a * k) + 4 * a * k ** 4 * np.cosh(k) + 4 * a * k ** 3 * np.sinh(k) + 12 * a * k ** 2 * np.cosh(k * (a - 1)) + 4 * a * k ** 2 * np.cosh(k * (a + 1)) - 4 * a * k ** 2 * np.cosh(k * (a - 2)) - 2 * a * k ** 4 * np.cosh(k * (a - 1)) - 2 * a * k ** 4 * np.cosh(k * (a + 1)) - 2 * a * k ** 4 * np.cosh(k * (a - 2)) + a * k ** 4 * np.cosh(2 * k * (a - 1)) - 2 * a ** 2 * k ** 4 * np.cosh(a * k) + 16 * a * k ** 3 * np.sinh(k * (a - 1)) - 8 * a * k ** 3 * np.sinh(k * (a - 2)) + 2 * a * k ** 3 * np.sinh(2 * k * (a - 1)) + 12 * a ** 2 * k ** 3 * np.sinh(a * k) - 4 * a ** 2 * k ** 4 * np.cosh(k) - 2 * a * k ** 3 * np.sinh(2 * k) - 4 * a ** 2 * k ** 3 * np.sinh(k) - 2 * a ** 2 * k ** 4 * np.cosh(k * (a - 1)) + 2 * a ** 2 * k ** 4 * np.cosh(k * (a + 1)) + 2 * a ** 2 * k ** 4 * np.cosh(k * (a - 2)) - 4 * a * k ** 3 * np.sinh(k * (2 * a - 1)) - 12 * a ** 2 * k ** 3 * np.sinh(k * (a - 1)) - 4 * a ** 2 * k ** 3 * np.sinh(k * (a + 1)) + 4 * a ** 2 * k ** 3 * np.sinh(k * (a - 2))) / (k ** 4 * (np.cosh(k) - 1) * (4 * np.cosh(k) - 4 * k * np.sinh(k) + k ** 2 + k ** 2 * np.cosh(k) - 4))) * M ** 2 + (4 * C * M * k ** 2 - 2 * C * M * k ** 4 + 8 * C * M * k * np.sinh(k) - 4 * C * M * k ** 2 * np.cosh(a * k) + 2 * C * M * k ** 4 * np.cosh(a * k) + 12 * C * M * k * np.sinh(k * (a - 1)) + 4 * C * M * k * np.sinh(k * (a + 1)) - 4 * C * M * k * np.sinh(k * (a - 2)) + C * M * k ** 3 * np.sinh(a * k) - 8 * C * M * k ** 2 * np.cosh(k) - 2 * C * M * k ** 4 * np.cosh(k) - 4 * C * M * k * np.sinh(2 * k) + 6 * C * M * k ** 3 * np.sinh(k) + 8 * C * M * k ** 2 * np.cosh(k * (a - 1)) - 4 * C * M * k ** 2 * np.cosh(k * (a - 2)) + C * M * k ** 4 * np.cosh(k * (a - 1)) + C * M * k ** 4 * np.cosh(k * (a + 1)) + 3 * C * M * k ** 3 * np.sinh(k * (a - 1)) - 3 * C * M * k ** 3 * np.sinh(k * (a + 1)) - C * M * k ** 3 * np.sinh(k * (a - 2)) + 4 * C * M * k ** 2 * np.cosh(2 * k) - C * M * k ** 3 * np.sinh(2 * k) - 12 * C * M * k * np.sinh(a * k) - 12 * C * M * a * k ** 2 * np.cosh(k * (a - 1)) - 4 * C * M * a * k ** 2 * np.cosh(k * (a + 1)) + 4 * C * M * a * k ** 2 * np.cosh(k * (a - 2)) + C * M * a * k ** 4 * np.cosh(k * (a - 1)) - C * M * a * k ** 4 * np.cosh(k * (a + 1)) + C * M * a * k ** 4 * np.cosh(k * (a - 2)) - 4 * C * M * a * k ** 3 * np.sinh(k * (a - 1)) + 4 * C * M * a * k ** 3 * np.sinh(k * (a + 1)) + 4 * C * M * a * k ** 3 * np.sinh(k * (a - 2)) + 12 * C * M * a * k ** 2 * np.cosh(a * k) - C * M * a * k ** 4 * np.cosh(a * k) - 4 * C * M * a * k ** 3 * np.sinh(a * k)) / (k ** 4 * (np.cosh(k) - 1) * (4 * np.cosh(k) - 4 * k * np.sinh(k) + k ** 2 + k ** 2 * np.cosh(k) - 4))
		return TacE

	def TacEmag_calc(a=None,k=None,lMagPrime=None,lDcPrime=None,*args,**kwargs):
		sigma1=np.sinh(k * (a - 1))
		sigma2=np.sinh(k * (a - 2))
		sigma3=np.sinh(k * (a + 1))
		sigma4=np.cosh(k * (a - 2))
		sigma5=np.cosh(k * (a - 1))
		sigma6=np.cosh(k * (a + 1))
		sigma7=np.sinh(k * (2 * a - 1))
		sigma8=np.sinh(2 * a * k)
		sigma9=np.cosh(2 * a * k)
		sigma10=2 * k * (a - 1)
		sinhk=np.sinh(k)
		sinhak=np.sinh(a * k)
		sinh2k=np.sinh(2 * k)
		coshk=np.cosh(k)
		coshak=np.cosh(a * k)
		cosh2k=np.cosh(2 * k)
		line1=8 * k - 8 * sinh2k - 24 * sinhak + 16 * sinhk + 24 * sigma1 + 8 * sigma3 - 8 * sigma2 + 12 * k ** 2 * sinhk + 12 * k ** 2 * sigma1 - 4 * k ** 2 * sigma3 - 4 * k ** 2 * sigma2 + k ** 2 * np.sinh(sigma10) - 8 * k * coshak - 3 * k ** 2 * sinh2k - 4 * k ** 2 * sigma7 + 8 * a * k ** 3 - 16 * k * np.cosh(k) + 16 * k * sigma5 - 8 * k * sigma4
		line2=8 * k ** 3 * coshak - 2 * k ** 3 * sigma9 - 4 * k ** 2 * sinhak + 3 * k ** 2 * sigma8 - 6 * k ** 3 - 8 * a ** 2 * k ** 3 + 8 * k * cosh2k - 24 * a * k * sigma5 - 8 * a * k * sigma6 + 8 * a * k * sigma4 - 12 * a * k ** 3 * coshak + 2 * a * k ** 3 * sigma9 - 4 * a ** 2 * k ** 2 * sinh2k + 16 * a * k ** 2 * sinhak - 4 * a * k ** 2 * sigma8 - 8 * a * k ** 3 * coshk
		line3=- 8 * a * k ** 2 * sinhk + 4 * a * k ** 3 * sigma5 + 4 * a * k ** 3 * sigma6 + 4 * a * k ** 3 * sigma4 - 2 * a * k ** 3 * np.cosh(sigma10) + 4 * a ** 2 * k ** 3 * coshak - 32 * a * k ** 2 * sigma1 + 16 * a * k ** 2 * sigma2 - 4 * a * k ** 2 * np.sinh(sigma10) - 24 * a ** 2 * k ** 2 * sinhak + 8 * a ** 2 * k ** 3 * coshk + 24 * a * k * coshak + 4 * a * k ** 2 * sinh2k
		line4=8 * a ** 2 * k ** 2 * sinhk + 4 * a ** 2 * k ** 3 * (sigma5 - sigma6 - sigma4) + 8 * a * k ** 2 * sigma7 + 24 * a ** 2 * k ** 2 * sigma1 + 8 * a ** 2 * k ** 2 * (sigma3 - sigma2)
		if k > 0.05:
			denom=2 * k ** 3 * (coshk - 1) * (4 * coshk - 4 * k * sinhk + k ** 2 + k ** 2 * coshk - 4)
		else:
			denom=k ** 6 / 72 + k ** 8 / 1440 + k ** 10 / 67200
		part1=lMagPrime ** 2 * (line1 + line2 + line3 + line4) / denom
		if k > 0.05:
			line5=2 * k - 2 * sinh2k - 6 * sinhak + 4 * sinhk + 6 * sigma1 + 2 * sigma3 - 2 * sigma2 + 2 * k ** 2 * sinhk + k ** 2 * (sigma1 - sigma3) - 3 * k * coshak - 4 * k * coshk + 3 * k * sigma5 + k * sigma6 - k * sigma4 + 2 * k * coshk ** 2 - 6 * a * k * sigma5 - 2 * a * k * sigma6 + 2 * a * k * sigma4 - a * k ** 2 * sinhak - a * k ** 2 * sigma1 + a * k ** 2 * sigma3 + a * k ** 2 * sigma2 + 6 * a * k * coshak
		else:
			line5=a ** 2 * (a - 1) ** 2 * (k ** 9 / 72 + k ** 11 * (4 * a ** 2 - 4 * a + 7) / 4320)
		denom2=k ** 3 * (coshk - 1) * (k * sinhk - 2 * coshk + 2)
		part2=- lMagPrime * lDcPrime * line5 / denom2
		TacEmag=part1 + part2
		return TacEmag