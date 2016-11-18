module Data where

-----------------------------------------------------------
---- ALL data type declarations used in this software -----
---Important helper functions for the used data types------
---Important Conversion functions between data types-------
-----------------------------------------------------------

import Numeric.Container hiding (linspace)
import Data.Maybe





----------------
--Dictionaries--
----------------

dict3 = [("S", 0), ("P", 1), ("D", 2), ("F", 3)]






------------------------
----- DATA TYPES -------
------------------------

--data type Orbital: Contains all necessary information about a certain orbital
data Orbital = Orbital { description :: (String, Int) 	-- eg ("S", 3)
			, numbering :: Vector Double 	-- array of the numbering of gaussians
		 	, exponents :: Vector Double 	-- exponents of gaussians
			, coeffs :: Vector Double    	-- coefficients of gaussians
			} deriving (Show)

--data type Atom: Collection of several (contracted) Orbitals
data Atom = Atom {atomname :: String	-- name of the Atom e.g. "CARBON"
		, orbitals :: [Orbital] -- array of the orbitals corresponding to an atom
		 } deriving (Show)


--data type Orbital: Contains all necessary information about a certain orbital
data Mol = Mol {  molname :: String		--molecule name , eg "H20" etc.
		, config :: [(Atom , Vector Double)] --configuration
		} deriving (Show)


--data type primitive Gaussian or short: Primitives
data PG = PG { lmn :: [Int]	--[l,m,n] for angular momentum
										, alpha :: Double --exponent 
										, position :: Vector Double --position vector
		} deriving (Eq, Show)

--data type Contraction
data Ctr = Ctr{ gaussians :: [PG]	--list of primitives
											, coefflist :: Vector Double --coefficient vector
											} deriving (Eq, Show)

------------------------------
------------------------------
----- DATA TYPE Functions ----
------------------------------
------------------------------



--------------------------
-- Orbital functions -----
--------------------------

--Get list of coefficients out of a given set (list of Orbitals)
--Input:		List of Orbitals
--Output:	List of Coefficient Vectors 
getcoeffs :: [Orbital] -> [Vector Double]
getcoeffs set = [coeffs i | i <- set]


--Get list of exponents out of a given set (list of Orbitals)
--Input:		List of Orbitals
--Output:	List of Exponent Vectors 
getexps :: [Orbital] -> [Vector Double]
getexps set = [exponents i | i <- set]


--see also below: getorbitals





--------------------------
---- Mol functions -------
--------------------------

geometry :: Mol -> [Vector Double]
geometry mol = geom k
	where 	k = length $ config mol
		geom n 
			| n == 0 = []
			| otherwise  = snd ((config mol) !! (k-n)) : geom (n-1)

natoms :: Mol -> Int
natoms mol = length $ config mol


--Get orbitals from i-indexed Atom for a given molecule
-- Input: 	molecule
--				index of atom
--Output: List of contracted Orbitals
getorbitals :: Mol -> Int -> [Orbital]
getorbitals mol i = orbitals $ (fst ((config mol) !! i))



--Get position vector for i-th Atom
getposofatom :: Mol -> Int -> Vector Double
getposofatom mol i = geometry mol !! i


-- Gives amount of contracted gaussian orbitals for a given molecule
nbasisfunctions :: Mol -> Int
nbasisfunctions mol = sum [length ( orbitals ( fst ( (config mol) !!(i)))) | i <- [0.. (length (config mol)) -1]]



--------------------------
------- PG functions -----
--------------------------


--gives angular momentum quantum number
momentum :: PG -> Int
momentum gauss = sum $ lmn gauss








--------------------------
------CONVERSION----------
--------------------------
---- Mol -> [Ctr] --------
--------------------------


--Converts molecule datatype to a list of Ctr Datatype
--Used for calculations using gaussian functions
mol_to_gaussians :: Mol -> [Ctr]
mol_to_gaussians mol = concat $ toCtr mol
	where
		--i indicates i-th Atom in Mol
		--j indicates j-th Orbital in Atom

		--These functions could also be made public
		pos mol i = snd $ config mol !! i 
		atm mol i = fst $ config mol !! i	
		getorb mol i j =  orbitals (atm mol i) !! j
		getLnumber mol i j = fromJust $ lookup (fst $ description $ orbitals (atm mol i) !! j) dict3
		norbs mol i = length $ orbitals $ atm mol i
		natms mol = length $ config mol


		llist l
			| l == 0 = [[0,0,0]]
			| l == 1 = [[0,0,1],[0,1,0],[1,0,0]]
			| l == 2 = [[1,1,0],[1,0,1],[0,1,1]] --contstructing the D-type needs still to be done
			| otherwise = []

		--gives Ctr for given molecule for i-th atom and j-th orbital and lmn config
		toprimitives mol i j lmn = Ctr ([PG lmn k (pos mol i) | k <- toList $ exponents $ getorb mol i j] ) (coeffs $ getorb mol i j)


		toCtr mol = concat ( [[[toprimitives mol i j lmn| lmn <- llist $ getLnumber mol i j] | j <-[0..(norbs mol i -1)]] | i <- [0..(natms mol -1)]  ] )

