"""
runs either a restricted (scf.rhf) or unrestricted (scf.uhf) Hartree-Fock calculation and returns orbital energies, coefficients, and ground state energy
"""

import math
import numpy
from numpy import linalg 

def rhf( S, T, Vne, Q, C, run_settings ):
    # run_settings = [ nalpha, nbeta, orth, maxscf, econv, variables.conv ]
    nalpha = run_settings[0]
    orth = run_settings[2]
    maxscf = run_settings[3]
    econv = run_settings[4]
    print_conv = run_settings[5]
    nbasis = len(S)

    # COMPUTE TRANSFORMATION MATRIX, X:
    if orth == 'symmetric':
        S_eigenvalues, S_eigenvectors = linalg.eigh(S)
        X = numpy.dot( S_eigenvectors, numpy.dot( S_eigenvalues**(-0.5) * numpy.identity(nbasis), S_eigenvectors.transpose() ) ) 
    elif orth == "canonical":
        S_eigenvalues, S_eigenvectors = linalg.eigh(S)
        X = S_eigenvectors *  S_eigenvalues**(-0.5) 
    else:
        print "No unitary transormation method specified!"

    # COMPUTE CORE (1-electron) HAMILTONIAN MATRIX, H:
    H = T + Vne

    # IF NECESSARY, COMPUTE INTITAL COEFFICIENTS, C, AND ENERGIES, e, FROM H:
    if C == "1el":
        F = H
        F_o = numpy.dot( X.transpose(), numpy.dot( F, X ) )
        e, C_o = linalg.eigh(F_o)
        order = e.argsort()
        e = e[order]
        C_o = C_o[:,order]
        C = numpy.dot( X, C_o)

    #RUN SCF ITERATIONS:
    if print_conv == True:
        print ""
        print "E [hartree]    deltaE [hartree]"
        print ""
    for n in range( maxscf ):

        # COMPUTE DENSITY MATRIX, P:
        P = numpy.zeros( ( nbasis, nbasis ) )
        for v in range( nalpha ):
            for i in range( nbasis ):
                for j in range( nbasis ):
                    P[ i, j ] = P[ i, j ] + 2 * C[ i, v ] * C[ j, v ]

        # COMPUTE FOCK MATRIX, F:
        G = numpy.zeros( ( nbasis, nbasis ) )
        for i in range( nbasis ):
            for j in range( nbasis ):
                gij = 0
                for k in range( nbasis ):
                    for l in range( nbasis ):
                        gij = gij + P[ l, k ] * ( Q[ i, j, k, l ] - 0.5 * Q[ i, l, k, j ] )
                G[ i, j ] = gij
        F = H + G

        # COMPUTE ELECTRONIC ENERGY, E:
        E = 0 
        for i in range( nbasis ):
            for j in range( nbasis ):
                E = E + 0.5 * P[ j, i ] * ( H[ i, j ] + F[ i, j ] )
        if n == 0:
            if print_conv == True:
                print E
                print ""
        else:
            deltaE = E - E_prev
            if print_conv == True:
                print E, deltaE 
                print "" 
            if abs(deltaE) < 10**(-econv):
                break
        E_prev = E  

        # COMPUTE NEW ORTITAL COEFFICIENTS, C:
        F_o = numpy.dot( X.transpose(), numpy.dot( F, X ) )
        e, C_o = linalg.eigh(F_o)
        order = e.argsort()
        e = e[order]
        C_o = C_o[:,order]
        C = numpy.dot( X, C_o)

    #EXIT:
    return E, C, e

def uhf( S, T, Vne, Q, C_alpha, C_beta, run_settings ):
    # run_settings = [ nalpha, nbeta, orth, maxscf, econv, variables.conv ]
    nalpha = run_settings[0]
    nbeta = run_settings[1]
    orth = run_settings[2]
    maxscf = run_settings[3]
    econv = run_settings[4]
    print_conv = run_settings[5]
    nbasis = len(S)

    # COMPUTE TRANSFORMATION MATRIX, X:
    if orth == 'symmetric':
        S_eigenvalues, S_eigenvectors = linalg.eigh(S)
        X = numpy.dot( S_eigenvectors, numpy.dot( S_eigenvalues**(-0.5) * numpy.identity(nbasis), S_eigenvectors.transpose() ) ) 
    elif orth == "canonical":
        S_eigenvalues, S_eigenvectors = linalg.eigh(S)
        X = S_eigenvectors *  S_eigenvalues**(-0.5) 
    else:
        print "No unitary transormation method specified!"

    # COMPUTE CORE (1-electron) HAMILTONIAN MATRIX, H:
    H = T + Vne

    # IF NECESSARY, COMPUTE INTITAL COEFFICIENTS, C_alpha and C_beta, AND ENERGIES, e, FROM H:
    if C_alpha == "1el":
        F = H
        F_o = numpy.dot( X.transpose(), numpy.dot( F, X ) )
        e, C_o = linalg.eigh(F_o)
        order = e.argsort()
        e = e[order]
        C_o = C_o[:,order]
        C_alpha = numpy.dot( X, C_o)
        C_beta = C_alpha 

    #RUN SCF ITERATIONS:
    if print_conv == True:
        print ""
        print "E [hartree]    deltaE [hartree]"
        print ""
    for n in range( maxscf ):

        # COMPUTE DENSITY MATRICES, P_alpha and P_beta:
        P_alpha = numpy.zeros( ( nbasis, nbasis ) )
        for v in range( nalpha ):
            for i in range( nbasis ):
                for j in range( nbasis ):
                    P_alpha[ i, j ] = P_alpha[ i, j ] + C_alpha[ i, v ] * C_alpha[ j, v ]
        P_beta = numpy.zeros( ( nbasis, nbasis ) )
        for v in range( nbeta ):
            for i in range( nbasis ):
                for j in range( nbasis ):
                    P_beta[ i, j ] = P_beta[ i, j ] + C_beta[ i, v ] * C_beta[ j, v ]

        # COMPUTE FOCK MATRICES, F_alpha and F_beta:
        G = numpy.zeros( ( nbasis, nbasis ) )
        for i in range( nbasis ):
            for j in range( nbasis ):
                gij = 0
                for k in range( nbasis ):
                    for l in range( nbasis ):
                        gij = gij + P_alpha[ l, k ] * ( Q[ i, j, k, l ] - Q[ i, l, k, j ] )
                        gij = gij + P_beta[ l, k ] * Q[ i, j, k, l ] 
                G[ i, j ] = gij
        F_alpha = H + G
        G = numpy.zeros( ( nbasis, nbasis ) )
        for i in range( nbasis ):
            for j in range( nbasis ):
                gij = 0
                for k in range( nbasis ):
                    for l in range( nbasis ):
                        gij = gij + P_beta[ l, k ] * ( Q[ i, j, k, l ] - Q[ i, l, k, j ] )
                        gij = gij + P_alpha[ l, k ] * Q[ i, j, k, l ] 
                G[ i, j ] = gij
        F_beta = H + G

        # COMPUTE ELECTRONIC ENERGY, E:
        E = 0 
        for i in range( nbasis ):
            for j in range( nbasis ):
                E = E + 0.5 * ( P_alpha[ j, i ] + P_beta[ j, i ] ) * H[ i, j ] 
                E = E + 0.5 * P_alpha[ j, i ] * F_alpha[ i, j ] 
                E = E + 0.5 * P_beta[ j, i ] * F_beta[ i, j ] 
        if n == 0:
            if print_conv == True:
                print E
                print ""
        else:
            deltaE = E - E_prev
            if print_conv == True:
                print E, deltaE 
                print "" 
            if abs(deltaE) < 10**(-econv):
                break
        E_prev = E  

        # COMPUTE NEW ORTITAL COEFFICIENTS, C_alpha:
        F_o = numpy.dot( X.transpose(), numpy.dot( F_alpha, X ) )
        e, C_o = linalg.eigh(F_o)
        order = e.argsort()
        C_o = C_o[:,order]
        e_alpha = e[order]
        C_alpha = numpy.dot( X, C_o)

        # COMPUTE NEW ORTITAL COEFFICIENTS, C_beta:
        F_o = numpy.dot( X.transpose(), numpy.dot( F_beta, X ) )
        e, C_o = linalg.eigh(F_o)
        order = e.argsort()
        C_o = C_o[:,order]
        e_beta = e[order]
        C_beta = numpy.dot( X, C_o)

    #EXIT:
    return E, C_alpha, C_beta, e_alpha, e_beta


