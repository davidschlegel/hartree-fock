import Numeric.Container hiding (linspace)
import Numeric.LinearAlgebra.Util
import Numeric.LinearAlgebra.Algorithms 
import Numeric.LinearAlgebra.Devel
import Data.List
import DSP.Basic


--For Charts
import Graphics.Rendering.Chart hiding (Matrix, Vector)
import Data.Colour
import Data.Colour.Names
import Data.Default.Class
import Graphics.Rendering.Chart.Backend.Cairo
import Control.Lens


{- 1.1.1 Method of Finite Differences
# Eigenvalue problem in one dimension. -}


--Create the tridiagonal differential matrix
--Possibly neccessary to convert into a sparse matrix 
d :: Matrix Double
d = buildMatrix 100 100 (\(i,j) -> if (i,j) == (i,i) then 2 else if (i,j) == (i, i-1) then -1 else if (i,j) == (i,i+1) then -1 else 0) 

--Create the potential matrix
--Possibly neccessary to convert into a sparse matrix 
v :: Matrix Double
v = buildMatrix 100 100 (\(i,j) -> if (i,j) == (0,0) then 1e100 else if (i,j) == (99,99) then 1e100 else 0) 

--The energy eigenvalues in a sorted List:
e = sort(toList(mapVector realPart (fst(eig(d)))))

--The corresponding eigenfunctions
psi = toColumns (mapMatrix realPart (snd(eig(d))))


--Plot
chart = toRenderable layout
  where
	
	--Grid for Visualization
	x = linspace 0.0 1.0 100 
	--List of Tuples [(Grid, psi) ] for Plotting
	points1 = zip ( x) (map (\x -> abs(x)^2) (toList (psi !! 1)))
	points2 = zip ( x) (map (\x -> abs(x)^2) (toList (psi !! 2)))
	points3 = zip ( x) (map (\x -> abs(x)^2) (toList (psi !! 3)))
	points4 = zip ( x) (map (\x -> abs(x)^2) (toList (psi !! 0)))  
	psi1 = plot_lines_values .~ [points1]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ "psi1"
              $ def
	psi2 = plot_lines_values .~ [points2]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "psi2"
              $ def
			  
	psi3 = plot_lines_values .~ [points3]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ "psi3"
              $ def		  
			  
	psi4 = plot_lines_values .~ [points4]
              $ plot_lines_style  . line_color .~ opaque black
              $ plot_lines_title .~ "psi4"
              $ def		  

	layout = layout_title .~ "Eigenfunctions"
           $ layout_plots .~ [toPlot psi1,
                              toPlot psi2,
							  toPlot psi3,
							  toPlot psi4]
           $ def
 
 
main = renderableToFile def{_fo_format=PDF} "eigenfunctions1.pdf" chart