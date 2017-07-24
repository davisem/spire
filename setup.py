__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__modname__ = "setup.py"

from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext

setup(name = "c_hasher",
	ext_modules = cythonize("hashers/chash.pyx"),
)