{-|

 Module      :  Gauss
 Copyright   :  Copyright (c) David Schlegel
 License     :  BSD
 Maintainer  :  David Schlegel
 Stability   :  experimental
 Portability :  Haskell

 Module "Data" provides all important datastructures and corresponding functions as well as dictionary lists.
-}


module Data (
-- * Definitions
-- | To provide a handy overview about angular momentum, element names and corresponding cardinal numbers, we give these list definitions which are beeing used with lookup tables.
momentumdict, numbers, atomsymbs, atomstrings,
-- * Constructors
-- |
--		* 'Orbital', defining a certain orbital; dependencies: None
--		* 'Atom', defining a collection of several (contracted) Orbitals; dependencies: 'Orbital'
--		* 'Mol', defining a Strucure containing multiple Atoms; dependencies: 'Atom', 'Orbital'
--		* 'PG', defining a primitive Gaussian function; dependencies: None
--		* 'Ctr', defining a set of contractions of primitive Gaussians; dependencies: 'PG'

Orbital(..), Atom(..), Mol(..), PG(..), Ctr(..),
-- | 

-- * Functions
-- ** Orbital functions
getcoeffs, getexps,
-- ** Mol functions
geometry, natoms, getorbitals, getposofatom, nbasisfunctions,
-- ** PG functions
momentumpg, 
-- ** Ctr functions
momentumctr, lengthcontract, lmncontract,
-- ** Conversion functions
mol_to_gaussians	
			) where

import Numeric.Container hiding (linspace)
import Data.Maybe


-- | Dictionary for angular momentum: [( \"S \", 0), ( \"P \", 1), ( \"D \", 2), ( \"F \", 3)]
momentumdict = [("S", 0), ("P", 1), ("D", 2), ("F", 3)]

-- | Ordered list of cardinal numbers: [1..118]
numbers = [1..118] :: [Double]

-- | Ordered list of element symbols: [\"H \", \"He \", \"Li \", \"Be \", \"B \", \"C \", \"N \", \"O \", \"F \", ... ]
atomsymbs = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"]

-- | Ordered list of element names: [ \"HYDROGEN \", \"HELIUM \", \"LITHIUM \", \"BERYLLIUM \", \"BORON \", \"CARBON \", \"NITROGEN \", \"OXYGEN \", \"FLUORINE \", ... ]
atomstrings = ["HYDROGEN","HELIUM","LITHIUM","BERYLLIUM","BORON","CARBON","NITROGEN","OXYGEN","FLUORINE","NEON","SODIUM","MAGNESIUM","ALUMINUM","SILICON","PHOSPHORUS","SULFUR","CHLORINE","ARGON","POTASSIUM","CALCIUM","SCANDIUM","TITANIUM","VANADIUM","CHROMIUM","MANGANESE","IRON","COBALT","NICKEL","COPPER","ZINC","GALLIUM","GERMANIUM","ARSENIC","SELENIUM","BROMINE","KRYPTON","RUBIDIUM","STRONTIUM","YTTRIUM","ZIRCONIUM","NIOBIUM","MOLYBDENUM","TECHNETIUM","RUTHENIUM","RHODIUM","PALLADIUM","SILVER","CADMIUM","INDIUM","TIN","ANTIMONY","TELLURIUM","IODINE","XENON","CESIUM","BARIUM","LANTHANUM","CERIUM","PRASEODYMIUM","NEODYMIUM","PROMETHIUM","SAMARIUM","EUROPIUM","GADOLINIUM","TERBIUM","DYSPROSIUM","HOLMIUM","ERBIUM","THULIUM","YTTERBIUM","LUTETIUM","HAFNIUM","TANTALUM","TUNGSTEN","RHENIUM","OSMIUM","IRIDIUM","PLATINUM","GOLD","MERCURY","THALLIUM","LEAD","BISMUTH","POLONIUM","ASTATINE","RADON","FRANCIUM","RADIUM","ACTINIUM","THORIUM","PROTACTINIUM","URANIUM","NEPTUNIUM","PLUTONIUM","AMERICIUM","CURIUM","BERKELIUM","CALIFORNIUM","EINSTEINIUM","FERMIUM","MENDELEVIUM","NOBELIUM","LAWRENCIUM","RUTHERFORDIUM","DUBNIUM","SEABORGIUM","BOHRIUM","HASSIUM","MEITNERIUM","DARMSTADTIUM","ROENTGENIUM","COPERNICIUM","UNUNTRIUM","FLEROVIUM","UNUNPENTIUM","LIVERMORIUM","UNUNSEPTIUM","UNUNOCTIUM"]







-- |Defines all information about a certain Gaussian Orbital
data Orbital = Orbital { 
			-- | orbital-type (angular momentum) and contraction length  
			description :: (String, Int)
			-- | numbering of gaussians
			, numbering :: Vector Double
			-- | exponents of gaussians  
		 	, exponents :: Vector Double
			-- | coefficients of gaussians  
			, coeffs :: Vector Double  
			} deriving (Eq, Show)

-- |Defines a collection of several (contracted) Orbitals
data Atom = Atom {atomname :: String	-- ^ name of the Atom e.g. \"CARBON\"  
		, orbitals :: [Orbital] -- ^ list of the orbitals corresponding to an atom   
		 } deriving (Eq, Show)


-- |Defines all necessary information about a Strucure containing multiple Atoms
data Mol = Mol {  molname :: String		-- ^ molecule name , eg \"H20\" etc.  
		, config :: [(Atom , Vector Double)] -- ^ configuration  
		} deriving (Eq, Show)


-- |Defines a primitive Gaussian, which has the form 

-- | <<primitive_gaussian.svg primitive_gaussian>>

-- | where l+m+n=L defines the angular momentum.
data PG = PG { lmn :: [Int]	-- ^ \[l,m,n\] angular momentum information  
										, alpha :: Double -- ^ exponent  
										, position :: Vector Double -- ^ position vector  
		} deriving (Eq, Show)

-- |Defines a Contraction of Gaussians, which has the form

-- | <<contraction.svg contraction>>
data Ctr = Ctr{ gaussians :: [PG]	-- ^ list of primitives  
											, coefflist :: Vector Double -- ^ coefficient vector  
											} deriving (Eq, Show)




-- | Gives list of coefficients out of a given Orbital set (list of Orbitals)
getcoeffs :: [Orbital] -> [Vector Double]
getcoeffs set = [coeffs i | i <- set]


-- |Get list of exponents out of a given set (list of Orbitals)
getexps :: [Orbital] -> [Vector Double]
getexps set = [exponents i | i <- set]



-- | Gives list of the positions for a given Mol. Length of list ist the number of Atoms in the given Mol.
geometry :: Mol -> [Vector Double]
geometry mol = geom k
	where 	k = length $ config mol
		geom n 
			| n == 0 = []
			| otherwise  = snd ((config mol) !! (k-n)) : geom (n-1)

-- | Gives the number of Atoms in a given Mol
natoms :: Mol -> Int
natoms mol = length $ config mol


-- |Gives list of orbitals from i-indexed Atom for a given Mol
getorbitals :: Mol -> Int -> [Orbital]
getorbitals mol i = orbitals $ (fst ((config mol) !! i))



-- |Gives position vector for i-th Atom
getposofatom :: Mol -> Int -> Vector Double
getposofatom mol i = geometry mol !! i


-- |Gives amount of contracted gaussian orbitals for a given Mol
nbasisfunctions :: Mol -> Int
nbasisfunctions mol = sum [length ( orbitals ( fst ( (config mol) !!(i)))) | i <- [0.. (length (config mol)) -1]]



-- |Gives angular momentum quantum number
momentumpg :: PG -> Int
momentumpg gauss = sum $ lmn gauss

-- |Gives angular momentum quantum number
momentumctr :: Ctr -> Int
momentumctr ctr = momentumpg (gaussians ctr !! 0)

-- |Gives length of contracttion
lengthcontract :: Ctr -> Int
lengthcontract ctr = length $ gaussians ctr

-- |Gives l m n out of a contraction
lmncontract :: (Num t, Num t1, Num t2) => Ctr -> (t, t1, t2)
lmncontract ctr = (fromIntegral $ a !! 0, fromIntegral $ a !! 1, fromIntegral $ a !! 2)
	where a = lmn (gaussians ctr !! 0)










-- |Converts molecule datatype to a list of Ctr Datatype
-- Used for calculations using gaussian functions
mol_to_gaussians :: Mol -> [Ctr]
mol_to_gaussians mol = concat $ toCtr mol
	where
		--i indicates i-th Atom in Mol
		--j indicates j-th Orbital in Atom

		--These functions could also be made public
		pos mol i = snd $ config mol !! i 
		atm mol i = fst $ config mol !! i	
		getorb mol i j =  orbitals (atm mol i) !! j
		getLnumber mol i j = fromJust $ lookup (fst $ description $ orbitals (atm mol i) !! j) momentumdict
		norbs mol i = length $ orbitals $ atm mol i
		natms mol = length $ config mol


		llist l
			| l == 0 = [[0,0,0]] --S
			| l == 1 = [[0,0,1],[0,1,0],[1,0,0]] --P
			| l == 2 = [[1,1,0],[1,0,1],[0,1,1],[2,0,0],[0,2,0],[0,0,2]] --D
			| l == 3 = [[1,1,1],[2,1,0],[1,2,0],[0,1,2],[0,2,1],[1,0,2],[2,0,1],[3,0,0],[0,3,0],[0,0,3]] --F
			| otherwise = error "Angular momentum not in range"

		--gives Ctr for given molecule for i-th atom and j-th orbital and lmn config
		toprimitives mol i j lmn = Ctr ([PG lmn k (pos mol i) | k <- toList $ exponents $ getorb mol i j] ) (coeffs $ getorb mol i j)


		toCtr mol = concat ( [[[toprimitives mol i j lmn| lmn <- llist $ getLnumber mol i j] | j <-[0..(norbs mol i -1)]] | i <- [0..(natms mol -1)]  ] )

