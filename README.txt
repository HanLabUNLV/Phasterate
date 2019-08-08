---------------------------------------------------------------------------
PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS
---------------------------------------------------------------------------

QUICK START - INSTALLING PHAST

-Compiling from Source
    Part 1 - Installing Clapack - (If you already have Clapack installed, skip to Part 2)
    Go into the directory (i.e. 'cd CLAPACK-3.2.1')
	and type 
    'cp make.inc.example make.inc && make f2clib && make blaslib && make lib'
    Note: Building Clapack can take several minutes depending on your system.
    There are several online resources for installing Clapack.

    Part 2 - Installing Phast
    1. Change directory to 'phast-1.3/src/' and run 'make CLAPACKPATH=../../CLAPACK-3.2.1'
    2. The Phast binaries should be created in the '../bin/' directory

    Part 3 - Installing Phasterate
    1. Change directory to Phasterate and run make


The Phast package should compile cleanly in most standard linux or
linux-like environments (including MacOS).  If you encounter problems
compiling, please report them to ADDEMAILHERE.  We'll do our
best to help you work around them and to avoid similar problems in the
future.

NOTES

    - If CLAPACK is used, PHAST also depends on the "F2C" (Fortran to C)
      package and on an implementation of the "BLAS" (Basic Linear Algebra
      Subroutines).  By default it uses the versions of these that come
      with CLAPACK.  The default BLAS implementation seems to be fine for
      normal usage.

    - The software requires GNU Make, some standard UNIX tools (e.g.,
      sed, ar, and ln), and a getopt implementation that supports long
      options (e.g., GNU getopt or BSD getopt).  These should be
      available on most UNIX systems, on Mac OS X, and via the Cygwin
      toolkit for Windows.

    - It's possible to compile the software without LAPACK by commenting
      out both the VECLIB and CLAPACKPATH lines in src/make-include.mk.  In
      this case, some programs will be usable, but programs that require
      matrix diagonalization will abort at the critical point of calling a
      LAPACK routine.




ACKNOWLEDGEMENTS

PHAST makes use of the CLAPACK linear algebra library
(http://www.netlib.org/clapack/) and the PCRE regular expression library
(http://www.pcre.org).  We thank the authors of these packages for making
them freely available to the community.
