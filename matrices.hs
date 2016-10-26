--module Matrices where

import Gauss
import Mol
import Read_basis
import Numeric.Container
import Data.Vector hiding (sum, length, fromList)


--Get orbitals from i-indexed Atom for a given molecule
-- Input: 	molecule
--				index of atom
--Output: List of contracted Orbitals
getorbitals :: Mol -> Int -> [Orbital]
getorbitals mol i = orbitals $ (fst ((config mol) !! i))


--Get list of coefficients out of a given set (list of Orbitals)
--Input:		List of Orbitals
--Output:	List of Coefficient Vectors 
getcoeffs :: [Orbital] -> [Data.Vector.Vector Double]
getcoeffs set = [coeffs i | i <- set]


--Get list of exponents out of a given set (list of Orbitals)
--Input:		List of Orbitals
--Output:	List of Exponent Vectors 
getexps :: [Orbital] -> [Data.Vector.Vector Double]
getexps set = [exponents i | i <- set]


--Get position vector for i-th Atom
getposofatom :: Mol -> Int -> Numeric.Container.Vector Double
getposofatom mol i = geometry mol !! i


--element wise zipped list (first coeff, second exp)
zippedVector :: Data.Vector.Vector Double -> Data.Vector.Vector Double -> [(Double, Double)]
zippedVector a b = Data.Vector.toList $ Data.Vector.zipWith (,) a b


--Converts the indices used for matrix calculations for a given molecule to a multiindex (m,n)
-- Where 	m = atom index
--				n = set index for given atom
indextomultiindex :: Mol -> Int -> (Int, Int)
indextomultiindex mol i = rekursion i 0
	where rekursion i n
			| length (getorbitals mol n) > i 	= (n,i)
 			| otherwise 															= rekursion (i-length (getorbitals mol n)) (n+1)



-- Gives amount of contracted gaussian orbitals for a given molecule
nbasisfunctions :: Mol -> Int
nbasisfunctions mol = sum [length ( orbitals ( fst ( (config mol) !!(i)))) | i <- [0.. (length (config mol)) -1]]



--Calculates the S-Matrix that contains the overlap of the orbitals
--Input: Molecule
--Output: S_Matrix, property: symmetric
s_pq :: Mol -> Matrix Double
--Building the Matrix: Size = nbasisfuctions of molecule
s_pq mol  = buildMatrix (nbasisfunctions mol) (nbasisfunctions mol) (\(i,j) -> calculateoverlap mol i j)
	where

		-- Actual function that calculates the overlap of the orbitals
		--Input: 	Molecule
		--				Indices i j
		--Output:	Calculated value of overlap (Gaussian Integrals)
		calculateoverlap :: Mol -> Int -> Int -> Double
		calculateoverlap mol i j = sum overlaplist
			where	
				--Getting atomindex and contracted orbital index a,k and b,l, respectively
				(a,k) = indextomultiindex mol i
				(b,l) = indextomultiindex mol j
				--Getting coefficients
				cs_i = (getcoeffs $ getorbitals mol a) !! k
				cs_j = (getcoeffs $ getorbitals mol b) !! l
				--Getting coefficients
				es_i = (getexps $ getorbitals mol a) !! k
				es_j = (getexps $ getorbitals mol b) !! l
				--Getting positions of corresponding atoms
				posi = getposofatom mol a
				posj = getposofatom mol b
				--Zip coeffs and exps into List of Tuples: (Could possibly be improved)
				zippedList_i = zippedVector cs_i es_i
				zippedList_j = zippedVector cs_j es_j
				--actually calculate the overlaps
				--Output: 	List of overlap integrals, whereby each coeffs and exps of index i are
				--				combined in every possible way with coeffs and exps of index j
				--				Length of List should be orderofcontraction_i x orderofcontractin_j 
				overlaplist = [(fst i)*(fst j)* ( overlaps (snd i) (snd j) posi posj) | i <- zippedList_i , j <- zippedList_j]





