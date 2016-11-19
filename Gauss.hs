module Gauss where

import Data
import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Data hiding (linspace)
import Data.Maybe

{-- | Notice: Formulas of Gaussian Integrals are taken from:
	Fundamentals of Molecular Integral Evaluation
	by Justin T. Fermann and Edward F. Valeev
--}



---------------------
---------------------
--Helper Functions---
---------------------
---------------------


-- |Calculates the Doublefactorial f(x) = x!!
factorial2 :: (Eq a, Num a) => a -> a
factorial2 0 = 1
factorial2 1 = 1
factorial2 n = n * factorial2 (n-2)

-- |Calculates the binomialcoefficient n over k
--binom :: Integer -> Integer -> Double
binom n 0 = 1
binom 0 k = 0
binom n k = binom (n-1) (k-1) * n `div` k 

-- |Computes the distance between to vectors: |r_A - r-B|^2
distance :: Vector Double -> Vector Double -> Double
distance rA rB = norm2(rA - rB)

-- |Calculates
frac ::  Double -> Double -> Double
frac alpha beta = alpha*beta/(alpha + beta)

--calculates weighted radius of two gaussians
--input: alpha, beta rA, rB
center :: Double -> Vector Double -> Double -> Vector Double -> Vector Double
center a rA b rB = ((scalar b * rB) * (scalar a * rA))/scalar (a+b)

-- |errof function calculation
erf :: Double -> Double
erf x = 2/sqrt pi * sum [x/(2*n+1) * product [-x^2/k |  k <- [1..n] ] | n <- [0..100]]

-- |error function calculation
f_0 :: Double -> Double
f_0 x = sqrt pi/2 * 1/sqrt x * (sqrt x) -- erf function missing







--------------------------
---Important Functions----
--------------------------


-- |Calculates normalization factor of contracted Gaussian with arbitrary angular momentum
normfactor :: Ctr -> Double
normfactor contr = (prefactor * summation)**(-1.0/2.0)
	--See also Eq. 2.11
	where
		(l, m, n) = lmncontract contr
		--(l,m,n) = (fromIntegral l_, fromIntegral m_, fromIntegral n_) --This looks quite dirty
		n_sum = lengthcontract contr -1
		a = toList $ coefflist contr
		alp = [alpha prim | prim <- (gaussians contr)]
		mom = fromIntegral $ momentumctr contr
		prefactor = 1.0/(2.0**mom)*pi**(3.0/2.0)* factorial2 (2*l-1) * factorial2 (2*m-1)* factorial2 (2*n-1)
		summation = sum $ concat $ [[(a !! i)*(a !! j)/((alp !! i +alp !! j)**(mom + 3.0/2.0)) | i <-[0..n_sum]]| j <- [0..n_sum]]



-- |Calculates f_k (l1 l2 pa_x pb_x) used in Gaussian Product Theorem
f :: Int -> Int -> Int -> Double -> Double -> Double
f k l1_ l2_ pa_x pb_x = sum $ concat $ [[pa_x**(fromIntegral (l1_-i))  *  pb_x**(fromIntegral(l2_-j)) * fromIntegral ((binom l1_ i) * (binom l2_ j))| i <- [0..l1_], i+j <= k]| j <-[0..l2_]]
--See also Eq. 2.45


-- |Evaluates a given function for two contractions
zipContractionWith :: (PG -> PG -> Double) -> Ctr -> Ctr -> Double
zipContractionWith zipfunction (Ctr pglist1 coeffs1) (Ctr pglist2 coeffs2) = pr1 * pr2 * value
	where
		coefflist1 = toList coeffs1
		coefflist2 = toList coeffs2
		n1 = length coefflist1 -1
		n2 = length coefflist2 -1
		pr1 = normfactor (Ctr pglist1 coeffs1) --I am not so sure about this
		pr2 = normfactor (Ctr pglist2 coeffs2) --I am not so sure about this
		value = sum $ [(coefflist1 !! i)* (coefflist2 !! j)* (zipfunction (pglist1 !! i) (pglist2 !! j)) | i <- [0..n1], j <- [0..n2]]



-- |Construct a matrix out of a list of Contractions and a function for two primitive gaussians
constr_matrix :: [Ctr] -> (PG -> PG -> Double) -> Matrix Double
constr_matrix contractionlist function = buildMatrix (length contractionlist) (length contractionlist) (\(i,j) -> zipContractionWith function (contractionlist !! i) (contractionlist !! j) )



-- |Calculates overlap of two primitive gaussians
s_12 :: PG -> PG -> Double
s_12 (PG lmn1 alpha1 pos1) (PG lmn2 alpha2 pos2) = prefactor * (i 0) * (i 1) * (i 2) --pref * I_x * I_y * I_z
	where
		(l1, m1, n1) = (lmn1 !! 0, lmn1 !! 1, lmn1 !! 2)
		(l2, m2, n2) = (lmn2 !! 0, lmn2 !! 1, lmn2 !! 2)
		g = alpha1 + alpha2
		prefactor = exp (-alpha1*alpha2 * (distance pos1 pos2) /g)
		p = center alpha1 pos1 alpha2 pos2
		pa = toList $ p  - pos1
		pb = toList $ p  - pos2
		range = [0..(fromIntegral (l1+l2) /2)] :: (Enum a, Fractional a) => [a]
		-- See also Eq. 3.15
		i k = sum $ [f (round (2*i)) l1 l2 (pa !! k) (pb !! k) * factorial2 (2.0* i -1.0) * ((2*g)** (- i)) * (pi/g)**(1/2) | i <- range ]

-- |Calculates kinetic energy of two primitive gaussians
t_12 :: PG -> PG -> Double
t_12 (PG lmn1 alpha1 pos1) (PG lmn2 alpha2 pos2) =  (i 0) + (i 1) + (i 2) -- I_x + I_y + I_z
	where
		--lifts or lowers l, m or n, indicated by index k for a lmn-list, respectively
		pp1 list k = [if i == k then list !! i + 1 else list !! i | i<-[0..2]]	--lifts
		mm1 list k = [if i == k then list !! i - 1 else list !! i | i<-[0..2]]	--lowers
		--overlaps where parts of l, m or n ar lifted +1 or lowered -1
		m1m1 k = s_12 (PG (mm1 lmn1 k) alpha1 pos1) (PG (mm1 lmn2 k) alpha2 pos2) -- <-1|-1>_k
		p1p1 k = s_12 (PG (pp1 lmn1 k) alpha1 pos1) (PG (pp1 lmn2 k) alpha2 pos2) -- <+1|+1>_k
		p1m1 k = s_12 (PG (pp1 lmn1 k) alpha1 pos1) (PG (mm1 lmn2 k) alpha2 pos2) -- <+1|-1>_k
		m1p1 k = s_12 (PG (mm1 lmn1 k) alpha1 pos1) (PG (pp1 lmn2 k) alpha2 pos2) -- <-1|+1>_k
		--See also Eq.4.13
		i k = 1/2 * (fromIntegral (lmn1 !! k)) * (fromIntegral (lmn2 !! k)) * (m1m1 k) 
			+ 2 * alpha1 * alpha2 * (p1p1 k)
			- alpha1 * (fromIntegral (lmn2 !! k)) * (p1m1 k)  
			- alpha2 * (fromIntegral (lmn1 !! k)) * (m1p1 k)




	 

-- Here all importand integrals involving gaussian functions will be evaluated.

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












