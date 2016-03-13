module Gauss where

import Mol
import Data.Number.Erf --error function
import Data.Vector as V

-- Here all importand integrals involving gaussian functions will be evaluated.

--function for computing the distance between to vecotrs: |r_A - r-B|^2
distance :: Floating a => Vector a -> Vector a -> a
distance rA rB = V.sum $ V.map (**2) $ V.zipWith (-) rA rB


frac :: Fractional a => a -> a -> a
frac alpha beta = alpha*beta/(alpha + beta)

--calculates weighted radius of two gaussians
--input: alpha, beta rA, rB
weightedR :: Fractional b => b -> Vector b -> b -> Vector b -> Vector b
weightedR a rA b rB =  V.map(/(a+b)) $ V.zipWith (+) (V.map (b*)  rB) (V.map (a*)  rA)


--Overlap integral for s orbitals
--solves <1s,alpha, A | 1s, beta, B>
overlaps :: Floating a => a -> a -> Vector a -> Vector a -> a
overlaps alpha beta rA rB = prefactor * exp(exponent)
	where 	prefactor = (pi/(alpha + beta))**(3/2)
		exponent = - ((frac alpha beta) * distance rA rB)

--Kinetic integral
--solves <1s,alpha, A | - Laplace | 1s, beta, B>
kinetic :: Floating a => a -> a -> Vector a -> Vector a -> a
kinetic alpha beta rA rB = prefactor * exp(exponent)
	where 	prefactor = (frac alpha beta)
			 *(6 -4 *(frac alpha beta) * distance rA rB) 
			 *(pi/(alpha + beta))**(3/2)
		exponent = - ((frac alpha beta) * distance rA rB)


--Nuclear interaction integral
--solves <1s,alpha, A | - Z/r_C | 1s, beta, B>
nuclear :: Erf a => a -> a -> Vector a -> Vector a -> Vector a -> a -> a
nuclear alpha beta rA rB rC z = pref1 * pref2 * exp(exponent)
	where 	pref1 = -2*pi*z/(alpha + beta)
		pref2 = 1/arg * erf arg
		arg   = sqrt((distance (weightedR alpha rA beta rB) rC )   * (alpha*beta))
		exponent = - ((frac alpha beta) * (distance rA rB))

--two-electron integral
--important!
--solves <1s,alpha, A ; 1s, beta, B | 1s,gamma, C ; 1s, delta, D >
