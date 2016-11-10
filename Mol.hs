module Mol where

import Read_basis
import Data.Maybe
import Numeric.Container hiding (linspace)
import System.IO
import System.IO.Unsafe
import Data.List
import Data.List.Split
import Text.Printf



--data type Orbital: Contains all necessary information about a certain orbital
data Mol = Mol {  molname :: String		--molecule name , eg "H20" etc.
		, config :: [(Atom , Vector Double)] --configuration
		} deriving (Show)

geometry :: Mol -> [Vector Double]
geometry mol = geom k
	where 	k = length $ config mol
		geom n 
			| n == 0 = []
			| otherwise  = snd ((config mol) !! (k-n)) : geom (n-1)

natoms :: Mol -> Int
natoms mol = length $ config mol



--Show all Important Information about the given Molecule
--Input: Mol
--Output: 
--	Molecule Name
--	Number of Atoms
--	For each Atom:
--		Atomname
--		Number of Orbitals
--		contracted set descriptions
--		geometry
molinfo :: Mol -> IO ()
molinfo mol = do
	putStrLn ""
	putStrLn ("Molecule name:   " ++ molname mol)
	putStrLn ("Number of Atoms:  " ++ show k)
	printAtoms mol k
		where 	k = length $ config mol
			printAtoms mol 0 = putStrLn " "
		 	printAtoms mol n = do 
				putStrLn ("\tAtom:    " ++ atomname ( fst ( (config mol) !! (n-1))))
				putStrLn  ("\t \t Number of Orbitals:" ++  "\t" ++ show (length ( orbitals ( fst ( (config mol) !!(k-n)))) ))
				putStrLn ("\t \t contracted sets:" ++ "\t" ++ show( [description k | k <- orbitals (fst ( (config mol) !! (k-n)))]))
				putStrLn ("\t \t geometry:" ++ "\t \t" ++ show( toList $ snd ( (config mol) !!(k-n)) ) )
				--putStrLn (str ++ orb)
				printAtoms mol (n-1)


molsafe :: Mol -> String -> IO ()
molsafe mol filename = writeFile filename $ molprint mol

--Show atoms with corresponding positions
--Especially for viewing the molecule in a different graphical program, s.a. molden
molprint :: Mol -> String
molprint mol = "    " ++ show (natoms mol) ++ "\n\n" ++ printAtoms mol k
	where
		--Little helper for printing lists tab-separated
		printList []	= ""
		printList list = "\t" ++ (printf "%.6f" (list !! 0) :: String) ++ printList (snd ( Prelude.splitAt 1 list))
		--Dictionary
		dict = [("CARBON", "C"), ("NITROGEN", "N"), ("OXYGEN", "O"), ("HYDROGEN", "H")] 	
		k = length $ config mol
		printAtoms mol 0 = []
	 	printAtoms mol n = " " ++ fromJust ( lookup (atomname ( fst ( (config mol) !! (n-1)))) dict)  ++ printList (toList ( (geometry mol) !! (k-n))) ++ "\n" ++ printAtoms mol (n-1)

showgeom mol = putStrLn $ molprint mol


--Function to read Geometry from a file
--Especially nice, when creating the geometry with a different programm, e.g. molden
--Input: directory
--Output: List with Atomstring and corresponding position vector
getgeom :: [Char] -> [([Char], Vector Double)]
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



constr_set_from_file list basis = rekursion n
	where
		n = length list -1
		rekursion 0 = do return ()
		rekursion i = do 
			generateFile basis (fst (list !! i)) (basis ++ "_" ++ fst (list !! i) ++ ".dat")
			rekursion (i-1)

get_mol_from_files :: [([Char], Vector Double)] -> [Char] -> [Char] -> Mol
get_mol_from_files list basis name = Mol name [(getAtomData (basis ++ "_" ++ fst (list !! i) ++ ".dat"), snd (list !! i))| i <- [0..(length list -1)] ]



