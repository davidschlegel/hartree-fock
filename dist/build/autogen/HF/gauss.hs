
{-|

 Module      :  HF.Gauss
 Copyright   :  Copyright (c) David Schlegel
 License     :  BSD
 Maintainer  :  David Schlegel
 Stability   :  experimental
 Portability :  Haskell

Gaussian Integral evaluation plays an important role in quantum chemistry. Here functions will be provided to compute the most important integrals involving gaussian-type orbitals (GTOs).

An unnormalized primitive cartesian Gaussian Orbital has the form

<<centered_gaussian.svg centered_gaussian>> 

where /A/ denotes the center of the orbital.

An unnormalized contracted Gaussian Orbital is a linear combination of primitive Gaussian:

<<centered_contracted_gaussian.svg centered_contracted_gaussian>>

"Gauss" uses the datatype contstructors 'PG' (Primitive Gaussian) and 'Ctr' (Contraction), defined in "HF.Data".

Notice: Formulas of Gaussian Integrals are taken from:
	Fundamentals of Molecular Integral Evaluation
	by Justin T. Fermann and Edward F. Valeev

-}



module HF.Gauss (
-- * Basic functions
distance, frac, center, itersum, erf, f, factorial2, general_erf,
-- * Normalization
normCtr, normPG,
-- * Integral Evaluation
s_12, t_12,
-- * Contraction operations
zipContractionWith, constr_matrix


) where

import HF.Data
import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Data hiding (linspace)
import Data.Maybe
import Control.DeepSeq
import Debug.Trace



---------------------
---------------------
--Helper Functions---
---------------------
---------------------


itersum function range
	| range == [] = 0
	| otherwise = (function $ head range) + itersum function (tail range)



-- |Calculates the factorial f(x) = x!
factorial n
    | n < 0     = 0
    | otherwise = product [1..n]

-- |Calculates the Doublefactorial f(x) = x!!
factorial2 :: (Show a, Eq a, Num a) => a -> a
factorial2 0 = 1
factorial2 1 = 1
factorial2 (-1) = 1
factorial2 n = n * factorial2 (n-2)

-- |Calculates the binomialcoefficient n over k
--binom :: Integer -> Integer -> Double
binom n 0 = 1
binom 0 k = 0
binom n k = binom (n-1) (k-1) * n `div` k 

-- |Computes the distance between to vectors:
--
-- <<norm.svg norm>>
distance :: Vector Double -> Vector Double -> Double
distance rA rB = norm2(rA - rB)

-- | Calculates effective coefficient:
--
-- >>> frac alpha beta = alpha*beta/(alpha + beta)
frac ::  Double -> Double -> Double
frac alpha beta = alpha*beta/(alpha + beta)

-- | Calculates weighted radius of two gaussians
center :: Double -- ^ alpha
			-> Vector Double -- ^ rA
			-> Double -- ^ beta
			-> Vector Double -- ^ rB
			-> Vector Double
center a rA b rB = ((scalar b * rB) * (scalar a * rA))/scalar (a+b)

-- |error function calculation
erf :: Double -> Double
erf x = 2/sqrt pi * sum [x/(2*n+1) * product [-x^2/k |  k <- [1..n] ] | n <- [0..100]]

--error function calculation
f_0 :: Double -> Double
f_0 = erf 

-- | General error function. See Eq. 4.118 Thijssen.
general_erf nu u
	| nu == 0 = sqrt pi * erf (sqrt u) * 1/(2*sqrt u)
	| otherwise = 1/(2*u) * (2*nu - 1) * (general_erf (nu-1) u ) - exp (-u)






--------------------------
---Important Functions----
--------------------------


-- |Calculates normalization factor of contracted Gaussian with arbitrary angular momentum
normCtr :: Ctr -> Double
normCtr contr = (prefactor * summation)**(-1.0/2.0)
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


-- |Calculates normalization factor of primitive Gaussian with arbitrary angular momentum
normPG :: PG -> Double
normPG (PG lmn alpha pos) = (2*alpha/pi)**(3/4) * (a/b)**(1/2)
	where
		(l, m, n) = (fromIntegral $ lmn !! 0, fromIntegral $ lmn !! 1, fromIntegral $ lmn !! 2)
		a = (8*alpha)**(l+m+n) * factorial l * factorial m * factorial n
		b = factorial (2*l) * factorial (2*m) * factorial (2*n)


-- |Calculates f_k (l1 l2 pa_x pb_x) used in Gaussian Product Theorem
f :: Int -> Int -> Int -> Double -> Double -> Double
--trace ("f" ++ " " ++  show k ++ " " ++ show l1_ ++ " " ++ show l2_ ++ " " ++ show pa_x ++ " " ++ show pb_x) $ 
f k l1_ l2_ pa_x pb_x = sum $ concat $ [[pa_x**(fromIntegral (l1_-i))  *  pb_x**(fromIntegral(l2_-j)) * fromIntegral ((binom l1_ i) * (binom l2_ j))| i <- [0..l1_], i+j <= k]| j <-[0..l2_]]
--See also Eq. 2.45


-- |Evaluates a given function for two contractions. The operation is symmetric in contraction arguments.
zipContractionWith :: (PG -> PG -> Double) -- ^ Function of two primitive Gaussians 
							-> Ctr -- ^ Contraction
							-> Ctr -- ^ Contraction 
							-> Double
zipContractionWith zipfunction (Ctr pglist1 coeffs1) (Ctr pglist2 coeffs2) = pr1 * pr2 * value
	where
		coefflist1 = toList coeffs1
		coefflist2 = toList coeffs2
		n1 = length coefflist1 -1
		n2 = length coefflist2 -1
		pr1 = normCtr (Ctr pglist1 coeffs1) --I am not so sure about this
		pr2 = normCtr (Ctr pglist2 coeffs2) --I am not so sure about this
		value = sum $ [(coefflist1 !! i)* (coefflist2 !! j)* (zipfunction (pglist1 !! i) (pglist2 !! j)) | i <- [0..n1], j <- [0..n2]]



-- |Construct a matrix out of a list of Contractions and a function for two primitive gaussians
constr_matrix :: [Ctr] -- ^ List of Contraction
						-> (PG -> PG -> Double) -- ^ Function of two primitive Gaussians
						-> Matrix Double -- ^ Matrix of dimensions (length [Ctr]) x (length [Ctr])
constr_matrix contractionlist function = buildMatrix (length contractionlist) (length contractionlist) (\(i,j) -> zipContractionWith function (contractionlist !! i) (contractionlist !! j) )



-- |Calculates overlap integral of two primitive gaussians

-- | <<s_12.svg s_12>>
s_12 :: PG -> PG -> Double
s_12 (PG lmn1 alpha1 pos1) (PG lmn2 alpha2 pos2) = trace ("s_12:   range " ++ show range  ++ "   f 0 =" ++ show (factor 0) ++ "   f 1 = " ++ show (factor 1)++ "   f 2 =" ++ show (factor 2)) $ prefactor * (factor 0) * (factor 1) * (factor 2) --pref * I_x * I_y * I_z
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
		function k i = force $ f (round (2*i)) l1 l2 (pa !! k) (pb !! k) *  (factorial2 (2.0* i -1.0)) * ((2*g)** (- i)) * (pi/g)**(1/2)
		factor k = itersum (function k) range
		--i k = sum $ [f (round (2*i)) l1 l2 (pa !! k) (pb !! k) * factorial2 (2.0* i -1.0) * ((2*g)** (- i)) * (pi/g)**(1/2) | i <- range ]



-- |Calculates kinetic energy integral of two primitive gaussians

-- | <<t_12.svg t_12>>
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

{-
-- |Calculates nuclear attraction integral of two primitive gaussians
v_12 :: Double -> Vecotr Double -> PG -> PG -> Double
v_12 z r_c (PG lmn1 alpha1 pos1) (PG lmn2 alpha2 pos2) = z*pi*pre/g * 2* summation
	where
		ax, ay, az = (fromList pos1) !! 0, (fromList pos1) !! 1, (fromList pos1) !! 2
		bx, by, bz = (fromList pos2) !! 0, (fromList pos2) !! 1, (fromList pos2) !! 2
		r_cx, r_cy, r_cz = (fromList r_c) !! 0, (fromList r_c) !! 1, (fromList r_c) !! 2
		pref = exp (-alpha1*alpha2 * (distance pos1 pos2) /g)
		g = alpha1 + alpha2
		p = 1/g * (alpha1*pos1 + alpha2*pos2)
		px, py, pz = (fromList p) !! 0, (fromList p) !! 1, (fromList p) !! 2
		eta = (alpha1 + alpha2)/g
		mu_x = l1 + l2 - 2*(ijk1 + ijk2) - (opq1 + opq2)
		mu_y = l1 + l2 - 2*(j1 + j2) - (p1 + p2)
		mu_z = l1 + l2 - 2*(k1 + k2) - (q1 + q2)
		nu = mu_x + mu_y + mu_z - (u + v + w)
		summation = general_erf nu (g * distance p r_c) * a_x * a_y * a_z

		--ranges:
		ijk1_range = [0..floor $ (fromIntegral l1 /2)] :: (Integral a, Enum a) => [a]
		ijk2_range = [0..floor $ (fromIntegral l2 /2)] :: (Integral a, Enum a) => [a]
		opq1_range ijk1 = [0..floor $ (fromIntegral l1 - 2*ijk1)] :: (Integral a, Enum a) => [a]
		opq2_range ijk2 = [0..floor $ (fromIntegral l2 - 2*ijk2)] :: (Integral a, Enum a) => [a]
		rst_range opq1 opq2 = [0..floor $ fromIntegral (opq1+l2) /2] :: (Integral a, Enum a) => [a]
		u_range mu_x = [0..floor $ (fromIntegral mu_x /2)] :: (Integral a, Enum a) => [a]



		a_x = (-1)**(l1+l2) * (factorial l1) * (factorial l2) * 
		(-1)**(opq2 + r) * factorial (opq1 + opq2) * 1/(4**(ijk1+ijk2+r) * factorial ijk1 * factorial ijk2 * factorial opq1 * factorial opq2 * factorial r) * 
		alpha1**(opq2-ijk1-r) * alpha2**(opq1-ijk2-r) * (ax-bx)**(opq1+opq2-2*r) * 1/(factorial (l1-2*ijk1-opq1) * factorial (l2-2*ijk2-opq2) * factorial (opq1-2*opq2-2*r)) * 
		* [ (-1)**u * factorial mu_x * (px -r_cx)**(mu_x - 2*u) * 1/ ((4)**u * factorial u * factorial (mu_x-2*u) * g**(opq1 + opq2 - r + u)) | u <- [0..1/2 * mu_x]]
-}

--v_12_test ::  Double -> Vector Double -> PG -> PG -> [Double]
{- v_12_test z r_c (PG lmn1 alpha1 pos1) (PG lmn2 alpha2 pos2) = summation 0
	where
		--(ax, ay, az) = ((toList pos1) !! 0, (toList pos1) !! 1, (toList pos1) !! 2)
		--(bx, by, bz) = ((toList pos2) !! 0, (toList pos2) !! 1, (toList pos2) !! 2)
		--(r_cx, r_cy, r_cz) = ((toList r_c) !! 0, (toList r_c) !! 1, (toList r_c) !! 2)
		--lmn1_ = fromList lmn1
		--lmn2_ = fromList lmn2
		
		pref = exp (-alpha1*alpha2 * (distance pos1 pos2) /g)
		g = alpha1 + alpha2
		p = scale (1/g) (add (scale alpha1 pos1) (scale alpha2 pos2))
		--(px, py, pz) = ((toList p) !! 0, (toList p) !! 1, (toList p) !! 2)

		--f_nu ijk1 ijk2 opq1 opq2 = general_erf (nu ijk1 ijk2 opq1 opq2) (g * distance p r_c)
		
		eta = (alpha1 + alpha2)/g
		mu ijk1 ijk2 opq1 opq2 s = fromIntegral (atIndex lmn1 s) + fromIntegral (atIndex lmn2 s) - 2*(ijk1 + ijk2) - (opq1 + opq2)
		--nu ijk1 ijk2 opq1 opq2 = (mu ijk1 ijk2 opq1 opq2 s) + (mu ijk1 ijk2 opq1 opq2 s) + (mu ijk1 ijk2 opq1 opq2 s) - (u v q)

		--ranges alle ranges (X,Y,Z) mit einer klappe schlagen, dafür index s:
		ijk1_range s = [0..fromIntegral $ floor $ fromIntegral (lmn1 !! s) /2] :: [Double]
		ijk2_range s = [0..fromIntegral $ floor $ fromIntegral (lmn2 !! s) /2] :: [Double]
		opq1_range ijk1 s = [0..fromIntegral $ floor $ (fromIntegral (lmn1 !! s) - 2*ijk1)] :: [Double]
		opq2_range ijk2 s = [0..fromIntegral $ floor $ (fromIntegral (lmn2 !! s) - 2*ijk2)] :: [Double]
		rst_range opq1 opq2 s = [0..fromIntegral $ floor $ (opq1+opq2) /2] :: [Double]
		uvw_range ijk1 ijk2 opq1 opq2 s = [0..fromIntegral $ floor $ ((mu ijk1 ijk2 opq1 opq2 s) /2)] :: [Double]

	
		--a_x :: Num a => a -> a -> a -> a -> a -> Double

		--(-1)**(l1+l2) * (factorial fromIntegral (lmn1 !! s)) * (factorial fromIntegral (lmn2 !! s)) *
		function ijk1 ijk2 opq1 opq2 rst s = (-1)**(opq2 + rst) * factorial (opq1 + opq2) * 1/(4**(ijk1+ijk2+rst) * factorial ijk1 * factorial ijk2 * factorial opq1 * factorial opq2 * factorial rst) * alpha1**(opq2-ijk1-rst) * alpha2**(opq1-ijk2-rst) * (atIndex pos1 s - atIndex pos2 s)**(opq1+opq2-2*rst) * 1/(factorial ((lmn1 !! s)-2*ijk1-opq1) * factorial ((lmn2 !! s)-2*ijk2-opq2) * factorial (opq1-2*opq2-2*rst)) * sum $ [ (-1)**uvw * factorial (mu ijk1 ijk2 opq1 opq2 s) * ((atIndex p s) - (atIndex r_c s))**((mu ijk1 ijk2 opq1 opq2 s) - 2*uvw) * 1/ ((4)**uvw * factorial uvw * factorial ((mu ijk1 ijk2 opq1 opq2 s)-2*uvw) * g**(opq1 + opq2 - rst + uvw)) | uvw <- (uvw_range ijk1 ijk2 opq1 opq2 s)]  

		summation s = (-1)**(lmn1 !! s + lmn2 !! s) * (factorial fromIntegral (lmn1 !! s)) * (factorial fromIntegral (lmn2 !! s)) * sum $ concat $ concat $ [[[function ijk1 ijk2 opq1 opq2 rst s | rst <- (rst_range opq1 opq2 s)] | opq1 <- (opq1_range ijk1 s), opq2 <- (opq2_range ijk2 s)]| ijk1 <- ijk1_range s, ijk2 <- ijk2_range s]
-}












-----------------
--Deprecated-----
-----------------

	 

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












