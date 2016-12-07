{-|

 Module      :  HF
 Copyright   :  Copyright (c) David Schlegel
 License     :  BSD
 Maintainer  :  David Schlegel
 Stability   :  experimental
 Portability :  Haskell

 Module "HF".
-}

module HF (
{-| Module "HF.Data" provides all important datastructures and corresponding functions as well as dictionary lists.	-} 
module HF.Data, 
{-|In "ReadWrite" important functions are provided to

	* Read input (Geometry information, Basis set information) and convert to internal datastructures, provided in "Data"
	* Write output (Geometry, structure and basis set information)
-}
module HF.ReadWrite, 
{-|Gaussian Integral evaluation plays an important role in quantum chemistry. In "HF.Gauss", functions will be provided to compute the most important integrals involving gaussian-type orbitals (GTOs). -}
module HF.Gauss)where
{-
--Total Packages included
import System.Process
import System.IO
import System.IO.Unsafe
import Data.List
import Data.List.Split
import Numeric.Container hiding (linspace)
import Data.Maybe
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data hiding (linspace)
import Text.Printf
-}

--Some Dictionaries
--dict1 = [("CARBON", "C"), ("NITROGEN", "N"), ("OXYGEN", "O"), ("HYDROGEN", "H")]
--dict2 = [("CARBON", 6.0), ("NITROGEN", 7.0), ("OXYGEN", 8.0), ("HYDROGEN", 1.0)]
--dict3 = [("S", 0), ("P", 1), ("D", 2), ("F", 3)]

import HF.Data

import HF.ReadWrite


import HF.Gauss


--import Matrices
import Numeric.Container hiding (linspace)



