import Control.Parallel

import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Util
import Numeric.LinearAlgebra.Data hiding (linspace)
import Numeric.LinearAlgebra.Algorithms 
import Numeric.LinearAlgebra.Devel
import Data.List hiding (find)
import DSP.Basic


--gridsize
n = 75

--index function
index :: (Int, Int) -> Int -> Double
index (i,j) n
	| i == j       =  6
	| i == j + n*n = -1
	| i == j - n*n = -1
	| i == j + n   = -1
	| i == j - n   = -1
	| i == j + 1   = -1
	| i == j - 1   = -1
	| otherwise    =  0

-- Laplace Matrix in 3 dimensions
t = buildMatrix (n*n*n) (n*n*n) (\(i,j) -> index (i,j) n)
eigen = eigenvalues t		

--Calculate corresponding Eigenfunction to n-th energy-eigenvalue
psi :: Int -> Vector Double
psi x = nullVector (a x)
	where a x = t - (mapMatrix (* (realPart (eigen @> x))) (ident (n*n*n)) )

--For further image processing
v = asColumn (psi 56)

main :: IO ()
main = do
-- write out matrix to visualize with python
saveMatrix "vmat.txt" "%f" v
