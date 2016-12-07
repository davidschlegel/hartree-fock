module Paths_hartree_fock (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName, getSysconfDir
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
import Prelude

catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
catchIO = Exception.catch

version :: Version
version = Version [0,1,0,0] []
bindir, libdir, datadir, libexecdir, sysconfdir :: FilePath

bindir     = "/home/david/.cabal/bin"
libdir     = "/home/david/.cabal/lib/x86_64-linux-ghc-7.10.3/hartree-fock-0.1.0.0-EbwprckSmo4Ieza9qLBpL8"
datadir    = "/home/david/.cabal/share/x86_64-linux-ghc-7.10.3/hartree-fock-0.1.0.0"
libexecdir = "/home/david/.cabal/libexec"
sysconfdir = "/home/david/.cabal/etc"

getBinDir, getLibDir, getDataDir, getLibexecDir, getSysconfDir :: IO FilePath
getBinDir = catchIO (getEnv "hartree_fock_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "hartree_fock_libdir") (\_ -> return libdir)
getDataDir = catchIO (getEnv "hartree_fock_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "hartree_fock_libexecdir") (\_ -> return libexecdir)
getSysconfDir = catchIO (getEnv "hartree_fock_sysconfdir") (\_ -> return sysconfdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
