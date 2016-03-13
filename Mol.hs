module Mol where

import Read_basis
import Data.Vector hiding ((++), sum, length)

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
				putStrLn ("\t \t geometry:" ++ "\t \t" ++ show( snd ( (config mol) !!(k-n)) ) )
				--putStrLn (str ++ orb)
				printAtoms mol (n-1)
				
