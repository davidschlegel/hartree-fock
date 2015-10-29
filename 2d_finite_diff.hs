import Control.Parallel

import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Util
import Numeric.LinearAlgebra.Data hiding (linspace)
import Numeric.LinearAlgebra.Algorithms 
import Numeric.LinearAlgebra.Devel
import Data.List hiding (find)
import DSP.Basic


--For Charts
--import Graphics.Rendering.Chart hiding (Matrix, Vector)
--import Data.Colour
--import Data.Colour.Names
--import Data.Default.Class
--import Graphics.Rendering.Chart.Backend.Cairo
--import Control.Lens 


--gridsize
n = 60

--index function
index :: (Int, Int) -> Int -> Double
index (i,j) n
	| i == j       =  4
	| i == j + n   = -1
	| i == j - n   = -1
	| i == j + 1   = -1
	| i == j - 1   = -1
	| otherwise    =  0

-- Laplace Matrix in 2 dimensions
t :: Matrix Double
t = buildMatrix (n*n) (n*n) (\(i,j) -> index (i,j) n)
	

-- Compute eigenvalues and corresponding eigenvectors 
e = sort(toList(mapVector realPart (fst(eig(t)))))
psi = toColumns (mapMatrix realPart (snd(eig(t))))


--For further image processing
grid :: Vector Double
grid = fromList (linspace 0.0 1.0 n)

ygrid = fromColumns $ replicate n grid
xgrid = fromRows    $ replicate n grid
-- eigenfunction matrix of n=12 for plotting
z1grid = reshape n $ mapVector (\x -> abs(x)^2) (psi !! 12)

main :: IO ()
main = do
-- write out matrices to visualize with python
saveMatrix "xmat.txt" "%f" xgrid
saveMatrix "ymat.txt" "%f" ygrid
saveMatrix "zmat.txt" "%f" z1grid

