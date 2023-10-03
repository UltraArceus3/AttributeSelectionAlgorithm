
from cython.parallel cimport prange
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.unordered_set cimport unordered_set
from libcpp.string cimport string
from libcpp.map cimport map
from libc.stdio cimport printf
from cython.operator cimport dereference

cdef extern from "edge_gen.cpp":
    cdef vector[int] EdgeGeneration(vector[vector[string]] de_duplicate,vector[int] samp,int blocking_attribute, map[string,vector[int]] dict_blocks, int threshold)


def edge_generation_c(vector[vector[string]] de_duplicate,vector[int] samp,int blocking_attribute,map[string,vector[int]] dict_blocks,int threshold):
    
    printf("Works till Cython call\n")

    cdef vector[int] output = EdgeGeneration(de_duplicate,samp,blocking_attribute,dict_blocks,threshold)
    return output