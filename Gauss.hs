module Gauss where

import Mol
import Data.Number.Erf --error function

-- Here all importand integrals involving gaussian functions will be evaluated.

--recursive function for computing the distance between to vecotrs: |r_A - r-B|^2
distance :: Floating a => Int -> [a] -> [a] -> a
distance n a b
	| n == 0 = 0
	| otherwise  = (a !! n - b !! n)**2 + distance (n-1) a b


--Overlap integral for s orbitals
--solves <1s,alpha, A | 1s, beta, B>
overlaps :: Floating a => a -> a -> [a] -> [a] -> a
overlaps alpha beta rA rB = prefactor * exp(exponent)
	where 	prefactor = (pi/(alpha + beta))**(3/2)
		exponent = - (alpha*beta/(alpha + beta) * distance 2 rA rB)

--Kinetic integral
--solves <1s,alpha, A | - Laplace | 1s, beta, B>


--Nuclear interaction integral
--solves <1s,alpha, A | - Z/r_C | 1s, beta, B>


--two-electron integral
--important!
--solves <1s,alpha, A ; 1s, beta, B | 1s,gamma, C ; 1s, delta, D >
