"""
computes distance matrix and nuclear-nuclear repulsion
"""
import math
import numpy

def distance( X, Z ):
    natoms = len(Z)
    R = numpy.zeros( ( natoms, natoms) )
    Vnn = 0.0
    for a in range( natoms ):
        for b in range( natoms ):
            if a == b:
                R[ a, b ] = 0.0
            elif b > a:
                r = 0.0
                for k in range(3):
                    r = r + ( X[b,k] - X[a,k] ) * ( X[b,k] - X[a,k] )
                r = math.sqrt(r)
                R[ a, b ] = r
            else:
                R[ a, b ] = R[ b, a ]
                Vnn = Vnn + Z[a] * Z[b] / R[ a, b ]
    return R, Vnn
