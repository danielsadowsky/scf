"""
transforms one-electron hamiltonians (transform.hamiltonian) and two-electron integrals (transform.integrals) from GTO basis to MO basis 
"""

import math
import numpy
from numpy import linalg

def hamiltonian( nbasis, T, Vne, C ):

    # TRANSFORM HAMILTONIAN FROM GTO BASIS, H, TO MO BASIS, H_MO:
    H = T + Vne
    H_MO = numpy.zeros( ( nbasis, nbasis ) )
    for r in range(nbasis):
        for s in range(r,nbasis):
            H_MO[r,s] = numpy.dot( numpy.asarray( [C[:,r]] ), numpy.dot( H, numpy.asarray( [C[:,s]] ).transpose() ) )
            H_MO[s,r] = H_MO[r,s]

    return H_MO

def integrals( nbasis, Q, C_a, C_b ):

    # TRANSFORM 4-CENTER INTEGRALS FROM GTO BASIS, Q, TO MO BASIS, Q_MO:
    temp1 = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )
    temp2 = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )
    temp3 = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )
    Q_MO = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )

    # USE r,s,t,u FOR MOs AND i, j, k, l FOR GTOs.
    for r in range(nbasis):
        for i in range(nbasis):
            temp1[r,:,:,:] = temp1[r,:,:,:] + C_a[i,r] * Q[i,:,:,:]    
        for s in range(nbasis):
            for j in range(nbasis):
                temp2[r,s,:,:] = temp2[r,s,:,:] + C_a[j,s] * temp1[r,j,:,:]
            for t in range(nbasis):
                for k in range(nbasis):
                    temp3[r,s,t,:] = temp3[r,s,t,:] + C_b[k,t] * temp2[r,s,k,:]
                for u in range(nbasis):
                    for l in range(nbasis):
                        Q_MO[r,s,t,u] = Q_MO[r,s,t,u] + C_b[l,u] * temp3[r,s,t,l]

    # BRUTE FORCE METHOD FOR TESTING (SCALES O^8)
    #Q_MO = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )
    #for r in range(nbasis):
    #    for s in range(nbasis):
    #        for t in range(nbasis):
    #            for u in range(nbasis):
    #                for i in range(nbasis):
    #                    for j in range(nbasis):
    #                        for k in range(nbasis):
    #                            for l in range(nbasis):
    #                                Q_MO[r,s,t,u] = Q_MO[r,s,t,u] + C[i,r] * C[j,s] * C[k,t] * C[l,u] * Q[i,j,k,l]

    return Q_MO
