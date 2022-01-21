#!/bin/python
import argparse
import math 
import numpy 
from numpy import linalg 

import config_basis
import distance
import integrals
import scf
import transform
import mp2
import ci
#import cdft
import basopt

# PARSE INPUT
parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help = 'file name for input file')
parser.add_argument('-b', type = str, help = 'file name for basis set')
parser.add_argument('-Z', type = int, help = 'molecular charge', default = 0 )
parser.add_argument('-m', type = int, help = 'multiplicity', default = 1 )
parser.add_argument('-c', type = str, help = 'mp2, ccd, ccsd, fci' )
parser.add_argument('-f', action='store_true', help = 'run frequency calculation' )
parser.add_argument('-o', action='store_true', help = 'run HF basis set optimization')
parser.add_argument('--maxscf', type = int, help = 'maximum SCF cycles', default = 50 )
parser.add_argument('--scfconv', type = float, help = 'requested convergence for SCF energy (use negative log)', default = 14 )
parser.add_argument('--basconv', type = float, help = 'requested convergence for basis set optimization (use negative log)', default = 10 )
parser.add_argument('--orthogonalization', type = str, help = 'method for orthogonalization of basis set (canonical/symmetric)', default = 'symmetric')
parser.add_argument('--orbitals', action='store_true', help = 'print 1- and 2- electron integrals')
parser.add_argument('--integrals', action='store_true', help = 'print optimized HF orbitals')
parser.add_argument('--conv', action='store_true', help = 'print HF convergence')
variables = parser.parse_args()
inputfile = variables.i
basis_def= variables.b
molcharge = variables.Z
mult = variables.m
runtype = variables.c
maxscf = variables.maxscf
econv = variables.scfconv
basoptconv = variables.basconv
orth = variables.orthogonalization
#print variables

############################################
# PARAMETERS NOT YET READ FROM INPUT FILE
############################################

# Self-solvation with cDFT?
do_cdft = "N"
dft_basis_def = "he_dft.bas"
dft_temp = 4.0 

# print options
print_distance = "N"

# finite difference step (0.001 for difficult cases):
h = 0.001

# scale step size (0.0005 for difficult cases):
scale = 0.04

# ensure hessian has positive eigenvalues for initial conditions:
check_hessian = "N"

# maximum number of steps in basis set optimization:
maxbasopt = 10000 

############################################

# READ .INP FILE AND EXTRACT NUCLEAR CHARGES AND POSITIONS:
Z = []
M = []
xang = []
yang = []
zang = []
inputdata = open( inputfile, 'r' )
elementnames = { "H":1, "He":2, "Li":3, "Be":4 }
masses = { "H":1.00782503224, "He":4.00260325415, "Li":7.0160034366, "Be":9.0121822 }
for line in inputdata:
     words = line.split()
     Z.append( elementnames.get( words[ 0 ] ) )
     M.append( masses.get( words[ 0 ] ) )
     xang.append( float( words[ 1 ] ) )
     yang.append( float( words[ 2 ] ) )
     zang.append( float( words[ 3 ] ) )

# DEFINE BASIS SETS: 
#
# basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
#
basis_settings = config_basis.cb( Z, basis_def )
natoms = basis_settings[0]

# CONVERT INPUT TO ATOMIC UNITS:
bohr = 0.52917721092
X = numpy.zeros(( natoms, 3 )) 
for i in range(natoms):
    X[i][0] = xang[i] / bohr
    X[i][1] = yang[i] / bohr
    X[i][2] = zang[i] / bohr

# DETERMINE NUMBER OF ELECTRONS:
nelec = 0
for i in range(natoms):
    nelec = nelec + int(Z[i])
nelec = nelec - molcharge

# DETERMINE NUMBER OF ALPHA AND BETA ELECTRONS
if type(nelec) != type(1) or type(mult) != type(1):
    print("Multiplicity and number of electrons must be integers!")
    exit()
if nelec%2 == 0 and mult%2 ==0:
    print("Invalid combination of multiplicity and number of electrons!")
    exit() 
elif nelec%2 == 1 and mult%2 ==1:
    print("Invalid combination of multiplicity and number of electrons!")
    exit() 
elif nelec%2 == 0 and mult == 1:
    print("scf.py will run RHF calculation.")
    nalpha = nelec / 2
    nbeta = nalpha
elif nelec%2 == 1 or mult != 1:
    print("scf.py will run UHF calculation.")
    nalpha = ( nelec + mult - 1 ) / 2
    nbeta = ( nelec - mult + 1 ) / 2
else:
    print("Number of electrons not specified!")
    exit()

run_settings = [ nalpha, nbeta, orth, maxscf, econv, variables.conv ]

def energy( X, Z, basis_settings, run_settings ):
    natoms = basis_settings[0]
    nbasis = basis_settings[1]

    # COMPUTE DISTANCE MATRIX, R, and NUCLEAR-NUCLEAR REPULSION, Vnn: 
    R, Vnn = distance.distance( X, Z )
    if print_distance == "Y":
        print "Distance Matrix, R:"
        print R
        print ""
        print "Nuclear-Nuclear Repulsion, Vnn:"
        print Vnn

    # COMPUTE INTEGRALS, S, T, V, AND Q:
    S, T, Vne, Q = integrals.integrals( R, X, Z, basis_settings )
    if variables.integrals:
        print "Overlap Matrix, S: "
        print S
        print ""
        print "Kinetic Energy Matrix, T: " 
        print T
        print ""
        print "Nuclear Attraction Potential, V: "
        print V
        print ""
        print "Two-Electron Integrals:"
        print Q
        print ""

    # RUN SCF CALCULATION:
    if nalpha == nbeta:
        E_el, C, e = scf.rhf( S, T, V, Q, "1el", run_settings )
        E = Vnn + E_el
        print "Total RHF Energy:"
        print E
        print ""
        if variables.orbitals:
            print "orbitals:"
            for i in range(nbasis):
                if i < ( nelec /2 ):
                    occ = 2
                else:
                    occ = 0
                print e[i], " ", occ
                for j in range(nbasis):
                    print "    ", basis_settings[7][j], C[j,i]
    else:
        E_el, C_alpha, C_beta, e_alpha, e_beta = scf.uhf( S, T, V, Q, "1el", "1el", run_settings )
        E = Vnn + E_el
        print "Total UHF Energy:"
        print E
        print ""
        if variables.orbitals:
            print "alpha orbitals:"
            for i in range(nbasis):
                if i < nalpha:
                    occ = "a"
                else:
                    occ = "0"
                print e_alpha[i], " ", occ
                for j in range(nbasis):
                    print "    ", basis_settings[7][j], C_alpha[j,i]
            print "beta orbitals:"
            for i in range(nbasis):
                if i < nbeta:
                    occ = "b"
                else:
                    occ = "0"
                print e_beta[i], " ", occ
                for j in range(nbasis):
                    print "    ", basis_settings[7][j], C_beta[j,i]

    # TRANSFORM INTEGRALS:
    if runtype == "cisd" or runtype == "cid" or runtype == "fci" or runtype == "mp2":    
        if nalpha == nbeta:
            H_MO = transform.hamiltonian( nbasis, T, V, C )  
            Q_MO = transform.integrals( nbasis, Q, C, C )
            C_alpha = C
            C_beta = C
            H_AMO = H_MO
            H_BMO = H_MO
            Q_AMO = Q_MO
            Q_BMO = Q_MO
            Q_ABMO = Q_MO
        else:
            H_AMO = transform.hamiltonian( nbasis, T, V, C_alpha )
            H_BMO = transform.hamiltonian( nbasis, T, V, C_beta )
            Q_AMO = transform.integrals( nbasis, Q, C_alpha, C_alpha )
            Q_BMO = transform.integrals( nbasis, Q, C_beta, C_beta )
            Q_ABMO = transform.integrals( nbasis, Q, C_alpha, C_beta ) 

    # RUN MP2 CALCULATION:
    if runtype == "mp2":
        if nalpha != nbeta:
            print "MP2 only available with RHF!"
        E2 = mp2.mp2( nalpha, nbasis, Q_MO, C, e )
        E = E + E2

    # RUN CI CALCULATION:
    if runtype == "cisd" or runtype == "cid" or runtype == "fci":    
        E_CI, C_CI = ci.ci( runtype, nbasis, nalpha, nbeta, H_AMO, H_BMO, Q_AMO, Q_BMO, Q_ABMO, C_alpha, C_beta )  
        E = Vnn + E_CI
        C_CI_0 = abs(C_CI[0])
        print "Correlation Energy:"
        print E_CI - E_el 
        print ""
        print "Total Energy:"
        print E 
        print ""
        print "CI coefficient of HF reference:"
        print C_CI_0 
        print ""
        if runtype == "cisd":
            davidson = (1 - C_CI_0**2) * (E_CI-E_el)
            print "Total Energy with Davidson correction:"
            print Vnn + E_CI + davidson
    return E

E = energy( X, Z, basis_settings, run_settings )

# DO CDFT:
if do_cdft == "Y":
    dft_natoms, dft_nbasis, dft_ngtos, dft_alpha, dft_coef, dft_gtos, dft_igto, dft_nuclei = config_basis.cb( Z, dft_basis_def )
    print "Coefficients:"
    print dft_coef
    print ""
    print "Alpha Values:"
    print dft_alpha
    E = cdft.cdft( Z, dft_nbasis, dft_alpha, dft_coef, dft_temp )

# RUN FREQ CALCULATION:
if variables.f:
    if natoms == 1:
        print "can't do frequency on a single atom!"
        print exit()
    elif natoms == 2:
        throwaway = 5
    else:
        throwaway = 6
    hq = 0.01
    Q = numpy.zeros(( 3 * natoms, 3 * natoms ))
    for i in range( 3 * natoms ):
        ia = i/3
        ix = i%3
        dxi = hq / math.sqrt( M[ia] )  
        for j in range( 3 * natoms ):
            ja = j/3 
            jx = j%3
            dxj = hq / math.sqrt( M[ja] )  
            U = numpy.copy(X) 
            U[ia,ix] = X[ia,ix] + dxi
            U[ja,jx] = X[ja,jx] + dxj
            #print i, j, ia, ix, ja, jx
            #print U
            E_u = energy( U, Z, basis_settings, run_settings )
            D = numpy.copy(X) 
            D[ia,ix] = X[ia,ix] - dxi
            D[ja,jx] = X[ja,jx] - dxj
            #print D
            E_d = energy( D, Z, basis_settings, run_settings )
            Q[i,j] = ( E_u + E_d - 2 * E ) / hq / hq
    #print Q
    vals, vecs = linalg.eigh( Q ) 
    print "frequencies:"
    for i in range( throwaway, 3 * natoms ):
        print math.sqrt( vals[i] )

# OPTIMIZE BASIS SET:
if variables.o:
    R, Vnn = distance.distance( X, Z )
    new_basis = basopt.basopt( h, scale, check_hessian, maxbasopt, basoptconv, R, X, Z, "1el", "1el", basis_settings, run_settings )

