"""
optimizes basis set to minimize HF ground state energy
"""

import math
import numpy
from numpy import linalg

import integrals
import scf

def eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings ):
    S, T, Vne, Q = integrals.integrals( R, X, Z, basis_settings )
    if basis_settings[0] == basis_settings[1]:
        C_guess = C_alpha
        E, C_RHF, e = scf.rhf( S, T, Vne, Q, C_guess, run_settings )
        C_alpha = C_RHF
        C_beta = C_RHF
    else:
        E, C_alpha, C_beta, e_alpha, e_beta = scf.uhf( S, T, Vne, Q, C_alpha, C_beta, run_settings )
    return E, C_alpha, C_beta 

def basopt( h, scale, check_hessian, maxbasopt, basoptconv, R, X, Z, C_alpha, C_beta, basis_settings, run_settings ):
    # basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
    natoms = basis_settings[0]
    nbasis = basis_settings[1]
    ngtos = basis_settings[2]
    alpha = basis_settings[3]
    coef = basis_settings[4]
    gtos = basis_settings[5]
    igto = basis_settings[6]
    nuclei = basis_settings[7]
    # run_settings = [ nalpha, nbeta, orth, maxscf, econv, variables.conv ]
    nalpha = run_settings[0]
    nbeta = run_settings[1]
    orth = run_settings[2]
    maxscf = run_settings[3]
    econv = run_settings[4]
    print_conv = run_settings[5]

    #COMPUTE GRADIENT, g, AND HESSIAN, H:
    g = numpy.zeros(nbasis)
    H = numpy.zeros((nbasis,nbasis))
    alpha_0 = list(alpha)
    E_0, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
    print alpha, E_0  
    for i in range(nbasis):
        #BACKWARDS MOVE
        alpha = list(alpha_0)
        alpha[i] = alpha[i] - h
        basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
        E_d, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
        print alpha, E_d 
        #FORWARD MOVE
        alpha = list(alpha_0)
        alpha[i] = alpha[i] + h
        basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
        E_u, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
        print alpha, E_u
        g[i] = 0.5 * ( E_u - E_d ) / h
        H[i,i] = ( E_d - 2*E_0 + E_u ) * h**(-2) 
        for j in range(i):
            #UU
            alpha = list(alpha_0)
            alpha[i] = alpha[i] + h
            alpha[j] = alpha[j] + h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_uu, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            print alpha, E_uu
            #UD
            alpha = list(alpha_0)
            alpha[i] = alpha[i] + h
            alpha[j] = alpha[j] - h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_ud, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            print alpha, E_ud
            #DU
            alpha = list(alpha_0)
            alpha[i] = alpha[i] - h
            alpha[j] = alpha[j] + h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_du, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            print alpha, E_du
            #DD
            alpha = list(alpha_0)
            alpha[i] = alpha[i] - h
            alpha[j] = alpha[j] - h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_dd, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            print alpha, E_dd
            H[i,j] = 0.25 * ( E_uu - E_ud - E_du + E_dd ) * h**(-2) 
            H[j,i] = H[i,j]
    print ""
    print "Gradient:"
    print g
    print
    print "Hessian:"
    print H
    print "Eigenvalues:"
    eigval, eigvec = linalg.eigh(H)
    print eigval
    if eigval[0] < 0:
        if check_hessian == "N":
            print "Hessian is not positive definate!"
            print ""
            for i in range(nbasis):
                for j in range(nbasis):
                    if H[i,j] < 0:
                        H[i,j] = 0
            print "Hessian modified. New eigenvalues:"
            print ""
            eigval, eigvec = linalg.eigh(H)
            print eigval
            print ""
        else:
            print "Hessian is not positive definate!"
            print ""
            print eigval
            print ""
            exit()
    # RUN OPTIMIZATION!
    alpha = numpy.asarray(alpha)
    s = numpy.zeros(nbasis)
    y = numpy.zeros(nbasis)
    for n in range(maxbasopt):

        # TAKE STEP, s, AND GET SCF ENERGY, U:
        s = -scale * numpy.dot(linalg.inv(H),g)
        alpha = alpha + s
        basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
        U, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )

        #DO SYMMETRIC FINITE DIFFERENCE TO UPDATE GRADIENT:
        for i in range(nbasis):
            alpha_fd = alpha
            alpha_fd[i] = alpha_fd[i] + h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_u, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            alpha_fd = alpha
            alpha_fd[i] = alpha_fd[i] - h
            basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
            E_d, C_alpha, C_beta = eval( R, X, Z, C_alpha, C_beta, basis_settings, run_settings )
            y[i] = 0.5 * ( E_u - E_d ) / h - g[i]
            g[i] = g[i] + y[i]

        #UPDATE HESSIAN:
        s_m = numpy.asarray([s])
        y_m = numpy.asarray([y])
        H = H + numpy.dot( y_m.transpose(), y_m ) / numpy.dot( y_m, s_m.transpose() ) - numpy.dot( numpy.dot( H, s_m.transpose() ), numpy.dot( s_m, H ) ) / numpy.dot( numpy.dot( s_m, H ), s_m.transpose() )

        #CHECK CONVERGENECE:
        if n > 0:
            deltaU = U - U_prev
            if abs(deltaU) < 10**(-basoptconv):
                basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
                S, T, Vne, Q = integrals.integrals( R, X, Z, basis_settings )
                print("Basis set optimization has converged in " + str(n) + " iterations!")
                print ""
                print alpha
                print ""
                print "Condition number of atomic S matrix:"
                print linalg.cond(S)
                print ""
                print "Coefficients of ground state RHF wavefunction:"
                for i in range(nalpha):
                    print C_alpha[:,i]
                if nalpha != nbeta:
                    print ""
                    for i in range(nbeta):
                        print C_beta[:,i]
                print ""
                print "Final SCF energy:"
                print U
                print ""
                break
        print alpha, U, g
        U_prev = U
    #WRITE CONVERGED BASIS TO OUTPUT
    output = open("output.bas",'w')
    output.write( str(nbasis) + "\n" )
    for i in range(nbasis):
        output.write( "S " + str(alpha[i]) + "\n" )
    print "Final basis vector written to output.bas"
    
    return alpha
