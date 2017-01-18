--IMPORTANT: This moduel is deprecated and no longer in use!!!


module HF.Matrices where

import HF.Data
import Hf.Gauss
import Numeric.Container
import Data.Maybe

------------------------------------------------------------------------------------
-----------IMPORTANT: This moduel is deprecated and no longer in use!!!-------------
------------------------------------------------------------------------------------


--element wise zipped list (first coeff, second exp)
zippedVector :: Vector Double -> Vector Double -> [(Double, Double)]
zippedVector a b = zipWith (,) (toList a) (toList b)


--Converts the indices used for matrix calculations for a given molecule to a multiindex (m,n)
-- Where 	m = atom index
--				n = set index for given atom
indextomultiindex :: Mol -> Int -> (Int, Int)
indextomultiindex mol i = rekursion i 0
	where rekursion i n
			| length (getorbitals mol n) > i 	= (n,i)
 			| otherwise 															= rekursion (i-length (getorbitals mol n)) (n+1)




--Calculates the S-Matrix that contains the overlap of the orbitals
--Input: Molecule
--Output: S_Matrix, property: symmetric
s_mat :: Mol -> Matrix Double
--Building the Matrix: Size = nbasisfuctions of molecule
s_mat mol  = buildMatrix (nbasisfunctions mol) (nbasisfunctions mol) (\(i,j) -> calculateoverlap mol i j)
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



--Calculates the kinetic-Matrix that contains the kinetic interactions of the orbitals
--Input: Molecule
--Output: 	kinetic_Matrix, property: symmetric
--				<alpha, A|-Nabla^2|beta,B>
--Reference: Eq. 4.101 p. 74 Thijssen
kin_mat :: Mol -> Matrix Double
kin_mat mol  = buildMatrix (nbasisfunctions mol) (nbasisfunctions mol) (\(i,j) -> calculatekinetic mol i j)
	where

		-- Actual function that calculates the kinetic interactions of the orbitals
		--Input: 	Molecule
		--				Indices i j
		--Output:	Calculated value of kinetic interactions (Gaussian Integrals)
		calculatekinetic :: Mol -> Int -> Int -> Double
		calculatekinetic mol i j = sum kineticlist
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
				--Output: 	List of kinetic integrals, whereby each coeffs and exps of index i are
				--				combined in every possible way with coeffs and exps of index j
				--				Length of List should be orderofcontraction_i x orderofcontractin_j 
				kineticlist = [(fst i)*(fst j)* ( kinetic (snd i) (snd j) posi posj) | i <- zippedList_i , j <- zippedList_j]



--Calculates the nuclear-Matrix that contains the nuclear attractions of the orbitals
--Input: Molecule
--Output: 	nuclear_Matrix, property: symmetric
--				<alpha, A|-Z/R_C^2|beta,B>
--Reference: Eq. 4.101 p. 74 Thijssen
nuc_mat_sng :: Mol -> Double -> Vector Double -> Matrix Double
nuc_mat_sng mol z r = buildMatrix (nbasisfunctions mol) (nbasisfunctions mol) (\(i,j) -> calculatenuclear mol i j z r)
	where

		-- Actual function that calculates the nuclear attractions of the orbitals
		--Input: 	Molecule
		--				Indices i j
		--Output:	Calculated value of nuclear attractions (Gaussian Integrals)
		calculatenuclear :: Mol -> Int -> Int -> Double -> Vector Double -> Double
		calculatenuclear mol i j z r = sum nuclearlist
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
				--Output: 	List of nuclear integrals, whereby each coeffs and exps of index i are
				--				combined in every possible way with coeffs and exps of index j
				--				Length of List should be orderofcontraction_i x orderofcontractin_j 
				nuclearlist = [(fst i)*(fst j)* ( nuclear (snd i) (snd j) posi posj r z) | i <- zippedList_i , j <- zippedList_j]



--Calculates the sum over all nuclear-Matrices for every atomic number and position
--Input: Molecule
--Output: 	total nuclear_Matrix, property: symmetric
--				Sum_{n} <alpha, A|-Z_n/R_n^2|beta,B>
--Reference: Eq. 4.101 p. 74 Thijssen
nuc_mat :: Mol -> Matrix Double
nuc_mat mol = addmat mol
	where
		--Number of Atoms in mol
		n = length $ config mol
		--get i-th Atomstring
		atomstring i = atomname $ fst ((config mol) !! i)
		--look it up in dictionary
		atomicnumber string = fromJust ( lookup string (zipWith (,) atomstrings numbers))
		pos i = geometry mol !! i 
		--sum over all atoms
		addmat mol = sum [nuc_mat_sng mol (atomicnumber $ atomstring i) (pos i) | i <- [0..n-1]]


--Construct h-Matrix
--Input:	Molecule
--Output: h-Matrix
h_mat :: Mol -> Matrix Double
-- (-)sign is already included in calculation of nuc_mat
h_mat mol =  0.5*kin_mat mol + nuc_mat mol






















