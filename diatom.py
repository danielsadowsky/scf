#!/bin/python
import argparse
import math 
import numpy as np
#from np import linalg as la 

import integrals
#import scf
#import transform
#import mp2
#import ci
#import basopt

# UNITS:
bohr = 0.52917721092

# PARSE INPUT:
parser = argparse.ArgumentParser()
parser.add_argument('-a', type = int, help = 'atomic number of atom a', default = 2 )
parser.add_argument('-b', type = int, help = 'atomic numbe  of atom b', default = 2 )
parser.add_argument('-d', type = float, help = 'distance between nuclei [in angstroms]', default = 1.00 )
parser.add_argument('-n', type = int, help = 'number of basis functions per atom', default = 2 )
parser.add_argument('-q', type = int, help = 'molecular charge', default = 0 )
parser.add_argument('-m', type = int, help = 'multiplicity', default = 1 )
parser.add_argument('-t', type = int, help = 'topology (3 for real space, 4 for curved space)', default = 4 )
parser.add_argument('-c', type = str, help = 'mp2, ccd, ccsd, fci' )
parser.add_argument('-o', action='store_true', help = 'run HF basis set optimization')
parser.add_argument('--maxscf', type = int, help = 'maximum SCF cycles', default = 50 )
parser.add_argument('--scfconv', type = float, help = 'requested convergence for SCF energy (use negative log)', default = 14 )
parser.add_argument('--basconv', type = float, help = 'requested convergence for basis set optimization (use negative log)', default = 10 )
parser.add_argument('--orthogonalization', type = str, help = 'method for orthogonalization of basis set (canonical/symmetric)', default = 'symmetric')
parser.add_argument('--orbitals', action='store_true', help = 'print 1- and 2- electron integrals')
parser.add_argument('--integrals', action='store_true', help = 'print optimized HF orbitals')
parser.add_argument('--conv', action='store_true', help = 'print HF convergence')
variables = parser.parse_args()

Z = [ variables.a, variables.b ]
nbasis = [ variables.n, variables.n ]
nel = variables.a + variables.b - variables.q
mult = variables.m

runtype = variables.c
maxscf = variables.maxscf
econv = variables.scfconv
basoptconv = variables.basconv
orth = variables.orthogonalization

# CONFIGURE NUCLEI: 
R = variables.d / bohr
Unn = Z[0] * Z[1] / R
if variables.t == 4:
    curved = True
elif variables.t == 3:
    curved = False
else:
    print("pick a valid topology!")
    exit()

# DEFINE BASIS:
alpha = []
nuclei = []
for i in range(2):
    for j in range(nbasis[i]):
        nuclei.append(i)
    if Z[i] == 1:  
        if nbasis[i] == 1:
            alpha.append( 0.117740489371 )        
        elif nbasis[i] == 2:
            alpha.append( 1.3218297478 )        
            alpha.append( 0.200264116742 )        
        elif nbasis[i] == 3:
            alpha.append( 4.40473350242 )
            alpha.append( 0.667871415249 )
            alpha.append( 0.149116082126 )
        else:
            print("maximum basis for H is 3.")
            exit()
    elif Z[i] == 2:  
        if nbasis[i] == 2:
            alpha.append( 4.0940459882 )        
            alpha.append( 0.531938227299 )        
        elif nbasis[i] == 3:
            alpha.append( 13.528060945 )
            alpha.append( 1.99193548399 )
            alpha.append( 0.382255870726 )
        else:
            print( "maximum basis for He is 3." )
            exit()
    else:
        print( "only H and He supported!" )
        exit()

# DETERMINE NUMBER OF ALPHA AND BETA ELECTRONS
if type(nel) != int or type(mult) != int:
    print( "Multiplicity and number of electrons must be integers!" )
    exit()
if nel%2 == 0 and mult%2 ==0:
    print( "Invalid combination of multiplicity and number of electrons!" )
    exit() 
elif nel%2 == 1 and mult%2 ==1:
    print( "Invalid combination of multiplicity and number of electrons!" )
    exit() 
elif nel%2 == 0 and mult == 1:
    print( "scf.py will run RHF calculation." )
    nalpha = nel / 2
    nbeta = nalpha
elif nel%2 == 1 or mult != 1:
    print( "scf.py will run UHF calculation." )
    nalpha = ( nel + mult - 1 ) / 2
    nbeta = ( nel - mult + 1 ) / 2
else:
    print( "Number of electrons not specified!" )
    exit()

print( "Nuclear Repulsion Energy, Unn:" )
print(Unn)
print()

S = integrals.curved( R, Z, alpha, nuclei, curved ) 
print( "Overlap Matrix, S:" )
print(S)
exit()
S, T, V, G = integrals.curved( R, Z, basis, nuclei, curved ) 

# TEST REAL SPACE FOR P FUNCTIONS
# DERIVE CURVED SPACE S FUNCTIONS

