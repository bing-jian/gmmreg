This directory hosts Python implementation of the robust point set registration
algorithm described in the paper
Robust Point Set Registration Using Gaussian Mixture Models,
Bing Jian and Baba C. Vemuri,
IEEE Transactions on Pattern Analysis and Machine Intelligence, 2011, 33(8),
pp. 1633-1645.
(An earlier conference version of this work appeared in the proceedings of ICCV'05
A Robust Algorithm for Point Set Registration Using Mixture of Gaussians,
Bing Jian and Baba C. Vemuri.)

It is part of software package which can be freely downloaded from
https://github.com/bing-jian/gmmreg


------------------------------------------
| Requirements for Using the Python Code |
------------------------------------------

Platforms:
The Python code should work on all platforms which are supported by Python, NumPy and SciPy.

Requirements:
(1) Python (version >=2.4 )
    Most versions of Linux and Mac OS already have python preinstalled. Windows users can download
    Python from http://www.python.org for free. Please make sure that Python is put in your path,
    so that typing 'python' in the system command line starts the Python interpreter.
(2) NumPy (version >=1.0) and SciPy (version >=0.6.0)
    In addition to the standard python installations, both NumPy and Scipy have to be installed so
    that neither typing 'import numpy' nor 'import scipy' in the python interpreter raise an error.
    For how to obtain and install NumPy/SciPy, please refer to http://www.scipy.org/
(3) Matplotlib [optional] (version >= 0.91.4)
    Matplotlib is only required for plotting the point sets and visualizing the registration results.
    Note that the axes3d support has been removed since Matplotlib v0.98.1. So Matplotlib v0.91.4 has
    to be used in order to plot 3D data.
(4) ConfigObj [optional]

-------------------------------------
| Install the gmmreg Python package |
-------------------------------------

In the system command line, go to the gmmreg/Python directory where you can find setup.py,
type the following command:
   python setup.py build

On Windows, you may receive an error message like below:
************************************************************************************
error: Python was built with Visual Studio 2003;
  extensions must be built with a compiler than can generate compatible binaries.
  Visual Studio 2003 was not found on this system. If you have Cygwin installed,
you can try compiling with MingW32, by passing "-c mingw32" to setup.py.
************************************************************************************
If you do have Cygwin installed, please first try
  python setup.py build -c mingw32
If this succeeds, then proceed with
  python setup.py install --skip-build
You may also install the package to a prefix of your choice by specifying --prefix
or --install-lib with the install command. For detailed help information, type
  python setup.py install --help

Upon the successful installation of gmmreg package, you should be able to "import gmmreg"
in the python interpreter.

http://stackoverflow.com/questions/2817869/error-unable-to-find-vcvarsall-bat
http://stackoverflow.com/questions/12418735/pip-fails-to-install-pil-or-pillow-with-mt-exe-error

---------------------------------
| How to use the gmmreg package |
---------------------------------

Data preparation:

To use this program, you need to have the input point set data ready, and prepare a configuration file (.ini)
to specify path and various parameters. The details of the input data and configuration file will be covered
in the rest of this README.

Example usage with the provided sample data and configuration file:
(1) Assume the program resides in c:\gmmreg.
(2) Assume you have successfully installed gmmreg package
(3) Test the Python code in the Python interpreter or an interactive Python shell. For example, if
IPython (http://ipython.scipy.org/) is installed, which is strongly recommended, then you can test
the following sequence of commands in IPython:
  cd c:\gmmreg\data
  import gmmreg
  gmmreg.test('fish_partial.ini')


---------------------------------------------------------
| About Configuration File and Input/Output Data Format |
---------------------------------------------------------

Configuration file:
We use the standard ini file to store the information about where to read input data, user-specified model
parameters and output files. Please refer to the example ini files in the data subdirectory.

Input data:
The input and output point sets are all stored in plain text files formatted as matrices where each row
represents the coordinates of one point.


-----------------------------------------------------------------
| Last modified: 10/30/2016                                     |
| If you have any questions, please contact bing.jian@gmail.com |
-----------------------------------------------------------------
