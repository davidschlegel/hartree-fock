import System.Process
import System.IO
import System.IO.Unsafe
--import Text.PrettyPrint.Boxes
import Data.List
import Data.List.Split



--helper function
--compares two elements
areTheySame x y | x == y = []
                | otherwise = [y]

--helper function
--removes a certain caracter from a list
removeItem _ [] = []
removeItem x (y:ys) = areTheySame x y ++ removeItem x ys

--data type Orbital: Contains all necessary information about a certain orbital
data Orbital = Orbital { description :: (String, Int) 	-- eg ("S", 3)
			, numbering :: [Double] 	-- array of the numbering of gaussians
		 	, exponents :: [Double] 	-- exponents of gaussians
			, coeffs :: [Double]    	-- coefficients of gaussians
			} deriving (Show)

--data type Atom: Collection of several (contracted) Orbitals
data Atom = Atom {name :: String	-- name of the Atom e.g. "CARBON"
		, orbitals :: [Orbital] -- array of the orbitals corresponding to an atom
		 } deriving (Show)



--main function: setAtomData
--input: filname of the textfile where the basis set of the given atom is stored
--output: Atom; data type described above
setAtomData :: [Char] -> Atom
setAtomData filename = Atom name orbitals
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
	      
	      name = (splitRow !! 0) !! 0 --first entry of splitRow
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
	 	      numbering   = get_orbital 0 n
	  	      exponents	  = get_orbital 1 n
	 	      coeffs	  = get_orbital 2 n 



generateFile:: [Char] -> [Char] -> [Char]
     -> IO (Maybe Handle, Maybe Handle, Maybe Handle, ProcessHandle)
generateFile basis atomname filename = do
	createProcess (shell string) {cwd = Just "EMSL_Basis_Set_Exchange_Local/", std_out = CreatePipe}
	where string = "python EMSL_api.py get_basis_data --basis " ++ basis ++" --atom " ++ atomname ++ " --treat_l > data/" ++ filename
	      

--main :: IO ()
--main = do

--createProcess (shell "python EMSL_api.py get_basis_data --basis STO-3G --atom C --treat_l > data/data2.txt") {cwd = Just "/media/david/Data3/Uni-Goettingen/SoSe 15/Physik/Haskell/QM_Project/EMSL_Basis_Set_Exchange_Local/", std_out = CreatePipe}

--print "Start"

