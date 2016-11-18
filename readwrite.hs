module ReadWrite where

import Data
import System.Process
import System.IO
import System.IO.Unsafe
import Data.List
import Data.List.Split
import Numeric.Container
import Text.Printf
import Data.Maybe

----------------
--Dictionaries--
----------------

dict1 = [("CARBON", "C"), ("NITROGEN", "N"), ("OXYGEN", "O"), ("HYDROGEN", "H")]


--helper function
--compares two elements
areTheySame x y | x == y = []
                | otherwise = [y]

--helper function
--removes a certain caracter from a list
removeItem _ [] = []
removeItem x (y:ys) = areTheySame x y ++ removeItem x ys






--main function: setAtomData
--input: filatomname of the textfile where the basis set of the given atom is stored
--output: Atom; data type described above
getAtomData :: [Char] -> Atom
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


--generates a file with orbital data from EMSL_api.py
--input: basis, atmname, filename
generateFile:: [Char] -> [Char] -> [Char]
     -> IO (Maybe Handle, Maybe Handle, Maybe Handle, ProcessHandle)
generateFile basis atmname filename = do
	createProcess (shell string) {cwd = Just "EMSL_Basis_Set_Exchange_Local/", std_out = CreatePipe}
	where string = "python2.7 EMSL_api.py get_basis_data --basis " ++ basis ++" --atom " ++ atmname ++ " --treat_l > data/" ++ filename
	      







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

		k = length $ config mol
		printAtoms mol 0 = []
	 	printAtoms mol n = " " ++ fromJust ( lookup (atomname ( fst ( (config mol) !! (n-1)))) dict1)  ++ printList (toList ( (geometry mol) !! (k-n))) ++ "\n" ++ printAtoms mol (n-1)

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

--This creates basis set files for a given list of atomstrings with corresponding position vector
--see also getgeom
--constr_set_from_file :: [([Char], b)] -> [Char] -> IO ()
constr_set_from_file list basis = rekursion n
	where
		n = length list -1
		rekursion 0 = do return ()
		rekursion i = do 
			generateFile basis (fst (list !! i)) (basis ++ "_" ++ fst (list !! i) ++ ".dat")
			rekursion (i-1)

get_mol_from_files :: [([Char], Vector Double)] -> [Char] -> [Char] -> Mol
get_mol_from_files list basis name = Mol name [(getAtomData (basis ++ "_" ++ fst (list !! i) ++ ".dat"), snd (list !! i))| i <- [0..(length list -1)] ]


