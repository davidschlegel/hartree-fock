{-|

 Module      :  HF.Gauss
 Copyright   :  Copyright (c) David Schlegel
 License     :  BSD
 Maintainer  :  David Schlegel
 Stability   :  experimental
 Portability :  Haskell

 In "ReadWrite" important functions are provided to

	* Read input (Geometry information, Basis set information) and convert to internal datastructures, provided in "Data"
	* Write output (Geometry, structure and basis set information)
-}
module HF.ReadWrite (
-- * Read
-- ** Basis set data

{-| Basis set data can be read using the <https://bse.pnl.gov/bse EMSL> database.  This is currently done by using the <https://github.com/TApplencourt/EMSL_Basis_Set_Exchange_Local EMSL_Basis_Set_Exchange_Local Python API>.

__Create a 'Mol' from basis set files:__

Assuming we have a basis set file \"C_STO_3G.dat\" (created with 'generateFile'), we can construct a 'Mol' datastructure by

>>> let carbon = Mol "SINGLE-Carbon" [(getAtomData "C_STO_3G.dat", fromList [1,0,0])]
 -}
generateFile, getAtomData,
-- ** Geometry
{-| To provide geometry information, it might be useful to read geometry information from an external file, which was created by another programm (e.g. Molden, for which it has been tested). 
This makes it easy to design a geometric structure with preguessed bond lengths and dihedral angles, etc.

__Examplary procedure__ to integrate external geometry into a 'Mol' datastructure:

Suppose, we have a file "CO2.xyz" containing the cartesian geometric information.

	* Read the geometry information, using 'getgeom':

			>>> let list = getgeom "CO2.xyz"

	* With 'constr_set_from_file', construct basis set data for the different elements which will be saved to files:

			>>> constr_set_from_file list "STO-6G"

	* With 'get_mol_from_files', read the basis set data to construct the complete molecule:

			>>> let mol = get_mol_from_files list "STO-6G" "CO2"
-}
getgeom, constr_set_from_file, get_mol_from_files,
-- * Write
molInfoPrint, molGeomPrint
) where

import HF.Data
import System.Process
import System.IO
import System.IO.Unsafe
import Data.List
import Data.List.Split
import Numeric.Container
import Text.Printf
import Data.Maybe


--helper function
--compares two elements
areTheySame x y | x == y = []
                | otherwise = [y]

--helper function
--removes a certain caracter from a list
removeItem _ [] = []
removeItem x (y:ys) = areTheySame x y ++ removeItem x ys






{-| Reads basis set data from a file into 'Atom' datastructure.

For creating a basis set file, see 'generateFile'.

NOTE: The basis set file has to exist before executing this function.

Example:

>>> getAtomData "C_STO_3G.dat"
-}
getAtomData :: [Char] -- ^ filename of the file where the basis set for a certain element is stored
					-> Atom -- ^ Atom datastructure
getAtomData filename = Atom atomname orbitals
	where 
	      	--descriptions: string array of orbital description 
	      descriptions =  [splitRow !! i | i <-[0..((length splitGeneral)-1)], (length (splitRow !! i) == 2)]
		--orbitallist: string array of remaining orbital data (numbering, exponents, coeffs)
              orbitallist = [splitRow !! i | i <-[0..((length splitGeneral)-1)],  (length (splitRow !! i) >= 3)]
	      	-- splitGeneral splits the string into array elements seperated by "\n"
	      splitGeneral = removeItem "" $ split (dropDelims $ oneOf "\n") (unsafePerformIO (readFile filepath))
	      	--e.g. filepath = "EMSL_Basis_Set_Exchange_Local/data/data2.txt"
	      filepath = "EMSL_Basis_Set_Exchange_Local/data/" ++ filename
	      	--splitRow removes empty spaces
	      splitRow = [removeItem "" $ split (dropDelims $ oneOf " ") (splitGeneral !! i) | i <- [0..((length splitGeneral)-1)]]
	      
	      atomname = (splitRow !! 0) !! 0 --first entry of splitRow
	      orbitals = setAtomDatab 0 --start of recursive function


              --function that gives the description of a certain orbital (given by index) from the descriptions array
	      --output: ((Orbital Letter, e.g. "S", "P", "D"), (Number of contraction coeffs for given basis))
	      get_description :: Int -> ([Char], Int)
	      get_description index = (x, y)
		where 	x  = ((descriptions !! index) !! 0) --S/P/D ..etc
			y  = (read ((descriptions !! index) !! 1) :: Int) --number of coefficients

	      --index1: numbering (0), exponents(1), coeffs (2)
	      --index2: number of orbital (0...norb-1)
	      --output: List of corresponding Doubles for given atom	
	      --get_orbtial :: Int -> Int -> [Double]
	      get_orbital index1 index2 = [read ((orbitallist !! j) !! index1) :: Double | j <- [l..l+(k-1)]]
		where l =  sum [snd (get_description i) | i <- [0..index2-1] ]
	      	      k =  snd ( get_description index2)


		--helper function to set the orbitals
		--output: List of all Orbitals for given atom
	      setAtomDatab :: Int -> [Orbital]
	      setAtomDatab n
		|  n >= length descriptions    = [] --base case
		| otherwise  = (Orbital description numbering exponents coeffs) : setAtomDatab (n+1) --recursion step
		where description = get_description n
	 	      numbering   = fromList $ get_orbital 0 n
	  	      exponents	  = fromList $ get_orbital 1 n
	 	      coeffs	  = fromList $ get_orbital 2 n 


{-| Generates a file with orbital data from the <https://github.com/TApplencourt/EMSL_Basis_Set_Exchange_Local EMSL_Basis_Set_Exchange_Local Python API>. 

Therefore the EMSL_Basis_Set_Exchange_Local folder has to be in the same working directory as from which this program is beeing executed.

NOTE: Not for every basis type and element a basis set exists. For a list of possible Atoms for a given basis type, run 

@ ./EMSL_api.py list_atoms --basis=\<basistype\>@

Example:

>>> generateFile "STO-3G" "C" "C_STO_3G.dat" 
-}
generateFile:: [Char] -- ^ basis, eg. \"STO-6G"\"
				-> [Char] -- ^ elementsymbol, eg. \"C\"
				-> [Char] -- ^ output filename
     			-> IO (Maybe Handle, Maybe Handle, Maybe Handle, ProcessHandle) -- ^ output file is saved to @ EMSL_Basis_Set_Exchange_Local\/data\/ @
generateFile basis atmname filename = do
	createProcess (shell string) {cwd = Just "EMSL_Basis_Set_Exchange_Local/", std_out = CreatePipe}
	where string = "python2.7 EMSL_api.py get_basis_data --basis '" ++ basis ++"' --atom " ++ atmname ++ " --treat_l > data/" ++ filename
	      


{-| Show all Important Information about a given Mol datatype.

__Output:__ 

		Molecule Name

		Number of Atoms

		For each Atom:

		* Atomname
		* Number of Orbitals
		* contracted set descriptions
		* geometry

__Usage Examples:__

>>> let mol = Mol "SINGLE-Carbon" [(getAtomData "C_STO_3G.dat", fromList [1,0,0])]

* Normal output: 

>>> putStrLn $ molInfoprint mol
Molecule name:   SINGLE-Carbon
Number of Atoms:  1
	Atom:    CARBON
	 	 Number of Orbitals:	3
	 	 contracted sets:	[("S",3),("S",3),("P",3)]
	 	 geometry:	 	[1.0,0.0,0.0]


* Save output to file:

>>> writeFile filename $ molInfoPrint mol
-}
molInfoPrint :: Mol -> String
molInfoPrint mol = "\n" ++ "Molecule name:   " ++ molname mol ++ "\n" ++ "Number of Atoms:  " ++ show k ++ "\n" ++ printAtoms mol k
		where
			k = length $ config mol
			printAtoms mol 0 = "\n"
		 	printAtoms mol n = "\tAtom:    " ++ atomname ( fst ( (config mol) !! (n-1) )) ++ "\n"
									++ "\t \t Number of Orbitals:" ++  "\t" ++ show (length ( orbitals ( fst ( (config mol) !!(k-n)))) ) ++ "\n"
									++ "\t \t contracted sets:" ++ "\t" ++ show( [description k | k <- orbitals (fst ( (config mol) !! (k-n)))]) ++ "\n"
									++ "\t \t geometry:" ++ "\t \t" ++ show( toList $ snd ( (config mol) !!(k-n)) )  ++ "\n"
									++ printAtoms mol (n-1)



{-| Show atoms with corresponding positions
Especially useful for viewing the molecule in a different graphical program, s.a. molden

__Usage:__

* Normal output: 

>>> putStrLn $ molGeomprint mol
    1
 C	1.000000	0.000000	0.000000

* Save output to file:

>>> writeFile filename $ molGeomPrint mol
-}
molGeomPrint :: Mol -> String
molGeomPrint mol = "    " ++ show (natoms mol) ++ "\n\n" ++ printAtoms mol k
	where
		--Little helper for printing lists tab-separated
		printList []	= ""
		printList list = "\t" ++ (printf "%.6f" (list !! 0) :: String) ++ printList (snd ( Prelude.splitAt 1 list))

		k = length $ config mol
		printAtoms mol 0 = []
	 	printAtoms mol n = " " ++ fromJust ( lookup (atomname ( fst ( (config mol) !! (n-1)))) (zipWith (,) atomstrings atomsymbs))  ++ printList (toList ( (geometry mol) !! (k-n))) ++ "\n" ++ printAtoms mol (n-1)



{-| Reads a file containing geometry information and outputs a list of element symbols, e.g. \"C\" with corresponding positions.

Especially nice, when creating the geometry with a different programm, e.g. Molden.

NOTE: The file should have the cartesian XYZ-file format that is used in Molden.-}
getgeom :: [Char] -- ^ directory
			-> [([Char], Vector Double)] -- ^ List with element symbols and corresponding position vectors
getgeom string = splitrows (n-1)
	where
		list = drop 2 ( removeItem "" $ split (dropDelims $ oneOf "\n") (unsafePerformIO (readFile string)) )
		splitRow = [removeItem "" $ split (dropDelims $ oneOf " ") (list !! i) | i <- [0..((length list)-1)]]

		n = length list

		tolistelement i = (a,b)
			where 
				a = (splitRow !! i) !! 0
				b = fromList [read ((splitRow !! i) !! j) :: Double | j <- [1..3]]

		splitrows i
			| i == -1 = []
			| otherwise = (tolistelement i):(splitrows (i-1))



{-| For a list of element symbols with corresponding positions (@\[(\[Char\], Vector Double)\]@) and a given basis set, e.g. \"STO-3G\", save basis set data to multiple files.
 
		IMPORTANT: The files are saved in the format @\<basis\>_\<element_name\>.dat@.

To construct a list of element symbols with corresponding positions (@\[(\[Char\], Vector Double)\]@) from a file, see 'getgeom'.
-}
constr_set_from_file :: [([Char], Vector Double)] -> [Char] -> IO ()
constr_set_from_file list basis = rekursion n
	where
		n = length list -1
		rekursion 0 = do return ()
		rekursion i = do 
			generateFile basis (fst (list !! i)) (basis ++ "_" ++ fst (list !! i) ++ ".dat")
			rekursion (i-1)


{-| For a list of element symbols with corresponding positions a given basis set and a molecule name, __read the existing basis set files__ (which can be created by using 'constr_set_from_file') and construct a 'Mol' datastructure. NOTE the difference between this function and 'constr_set_from_file'!
-}
get_mol_from_files :: [([Char], Vector Double)] -- ^ list of element symbols with corresponding positions
						-> [Char] -- ^ basis set name, eg. \"STO-3G\"
						-> [Char] -- ^ molname, eg. \"CO2\"
						-> Mol
get_mol_from_files list basis name = Mol name [(getAtomData (basis ++ "_" ++ fst (list !! i) ++ ".dat"), snd (list !! i))| i <- [0..(length list -1)] ]


