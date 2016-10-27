module Gauss where


import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Data hiding (linspace)


-- Here all importand integrals involving gaussian functions will be evaluated.

--function for computing the distance between to vectors: |r_A - r-B|^2
distance :: Vector Double -> Vector Double -> Double
distance rA rB = norm2(rA - rB)


frac ::  Double -> Double -> Double
frac alpha beta = alpha*beta/(alpha + beta)

--calculates weighted radius of two gaussians
--input: alpha, beta rA, rB
center :: Double -> Vector Double -> Double -> Vector Double -> Vector Double
center a rA b rB = ((scalar b * rB) * (scalar a * rA))/scalar (a+b)

erf :: Double -> Double
erf x = 2/sqrt pi * sum [x/(2*n+1) * product [-x^2/k |  k <- [1..n] ] | n <- [0..100]]

--error function calculation
f_0 x = sqrt pi/2 * 1/sqrt x * (sqrt x) -- erf function missing

--Overlap integral for s orbitals
--calculates <1s,alpha, A | 1s, beta, B>
overlaps :: Double -> Double -> Vector Double -> Vector Double -> Double
overlaps alpha beta rA rB = prefactor * exp exponent
	where 	prefactor = (pi/(alpha + beta))**(3/2)
		exponent = - (frac alpha beta * distance rA rB)

--Kinetic integral
--calculates <1s,alpha, A | - Laplace | 1s, beta, B>
kinetic :: Double -> Double -> Vector Double -> Vector Double -> Double
kinetic alpha beta rA rB = prefactor * exp(exponent)
	where 	prefactor = (frac alpha beta)
			 *(6 -4 *(frac alpha beta) * distance rA rB)
			 *(pi/(alpha + beta))**(3/2)
		exponent = - ((frac alpha beta) * distance rA rB)



--Nuclear interaction integral
--calculates <1s,alpha, A | - Z/r_C | 1s, beta, B>
nuclear :: Double -> Double -> Vector Double -> Vector Double -> Vector Double -> Double -> Double
nuclear alpha beta rA rB rC z = pref1 * pref2 * exp(exponent)
	where 	pref1 = -2*pi*z/(alpha + beta)
		pref2 = f_0 arg
		arg   = (distance (center alpha rA beta rB) rC )   * (alpha*beta)
		exponent = - ((frac alpha beta) * (distance rA rB))

--two-electron integral
--important!
--calculates <1s,alpha, A ; 1s, beta, B | 1s,gamma, C ; 1s, delta, D >
twoelectron :: Double -> Double -> Double -> Double -> Vector Double -> Vector Double -> Vector Double -> Vector Double -> Double
twoelectron a b g d rA rB rC rD = pref * exp(exponent) * (f_0 arg)
	where	pref = 			(2*pi**(5/2))/
				( (a+g)*(b+d)*(a+b+g+d)**(1/2)   )
		exponent = 	-a*g/(a+g) * (distance rA rC)
				-b*d/(b+d) * (distance rB rD)
		arg = (a+g)*(b+d)/(a+b+g+d) * ( distance (center a rA g rC) (center b rB d rD) )




--Overlap Matrix S_pq = <p|q>
--overlap p q  = t
--where
--	t = buildMatrix (n*n) (n*n) (\(i,j) -> overlaps alpha beta rA rB )










