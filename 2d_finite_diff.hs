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

--first blockmatrix
block1 = buildMatrix n n (\(i,j) -> if (i,j) == (i,i) then 4 else if (i,j) == (i, i-1) then -1 else if (i,j) == (i,i+1) then -1 else 0) 
block1_Sparse = mkSparse [((i,j), (block1 @@> (i,j))) | (i,j) <- (find (/= 0) block1)]

--second blockmatrix
block2 = -1 * ident n
block2_Sparse = mkSparse [((i,j), (block2 @@> (i,j))) | (i,j) <- (find (/= 0) block2)]

--function mapping three matrices to a tridiagonal matrix of desired dimension
tridiag ::  Int -> Matrix Double -> Matrix Double -> Matrix Double -> Matrix Double
tridiag n l m r = a + b + c
	where 
		sizeb = fst (size l)
		sizec = fst (size r)
		a = diagBlock (replicate n m)
		b = (zeros sizeb (sizeb^2)) === ((diagBlock (replicate (n-1) l)) |||  (zeros (sizeb^2 - sizeb) sizeb ))
		c = ((zeros (sizec^2 -n) sizec ) ||| (diagBlock (replicate (n-1) r)))  === (zeros sizec (sizec^2))
		
		
t = tridiag n block2 block1 block2		
e = sort(toList(mapVector realPart (fst(eig(t)))))
psi = toColumns (mapMatrix realPart (snd(eig(t))))

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

