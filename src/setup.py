from distutils.core import Extension, setup
from Cython.Build import cythonize


ext = Extension(name="edge_generation", sources=["edge_generation.pyx"])
setup(ext_modules=cythonize(ext))