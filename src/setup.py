from distutils.core import Extension, setup
from Cython.Build import cythonize


ext = Extension(name="edge_generation", sources=["edge_generation.pyx"],extra_compile_args=["-O3","-fopenmp"],extra_link_args=['-fopenmp'],language="c++")
setup(ext_modules=cythonize(ext))