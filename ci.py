"""
computes CI correlation energies from RHF or UHF orbitals 
"""

import math
import numpy
from numpy import linalg

def ci( runtype, nbasis, nalpha, nbeta, H_AMO, H_BMO, Q_AMO, Q_BMO, Q_ABMO, C_alpha, C_beta ):

    # INITIALIZE CI OCCUPATION INDEX ARRAY, O_CI:
    O_CI = numpy.zeros((0,2*nbasis))

    # GENERATE REFERENCE DETERMINANT 
    if runtype != "fci":
        occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
        O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))

    # GENERATE DETERMINANTS WITH SINGLE EXCITATIONS (AND SAME MULTIPLICITY AS REFERENCE)
    if runtype == "cisd":
        # ALPHA EXCITATIONS
        for a in range(nalpha):
            for r in range(nalpha,nbasis):
                occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
                occ[a] = 0
                occ[r] = 1
                O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))
        # BETA EXCITATIONS
        for a in range(nbeta):
            for r in range(nbeta,nbasis):
                occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
                occ[a+nbasis] = 0
                occ[r+nbasis] = 1
                O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))
     
    # GENERATE DETERMINANTS WITH DOUBLE EXCITATIONS (AND SAME MULTIPLICITY AS REFERENCE)
    if runtype != "fci":
        # DOUBLE-ALPHA EXCITATIONS
        for a in range(nalpha):
            for b in range(a+1,nalpha):
                for r in range(nalpha,nbasis):
                    for s in range(r+1,nbasis):
                        occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
                        occ[a] = 0
                        occ[b] = 0
                        occ[r] = 1
                        occ[s] = 1
                        O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))
        # ALPHA-BETA EXCITATIONS
        for a in range(nalpha):
            for b in range(nbeta):
                for r in range(nalpha,nbasis):
                    for s in range(nbeta,nbasis):
                        occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
                        occ[a] = 0
                        occ[b+nbasis] = 0
                        occ[r] = 1
                        occ[s+nbasis] = 1
                        O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))
        # DOUBLE-BETA EXCITATIONS
        for a in range(nbeta):
            for b in range(a+1,nbeta):
                for r in range(nbeta,nbasis):
                    for s in range(r+1,nbasis):
                        occ = [1] * nalpha + [0] * ( nbasis - nalpha ) + [1] * nbeta + [0] * ( nbasis - nbeta ) 
                        occ[a+nbasis] = 0
                        occ[b+nbasis] = 0
                        occ[r+nbasis] = 1
                        occ[s+nbasis] = 1
                        O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))
     
    # GENERATE FULL CI DETERMINANT COEFFICIENTS:
    if runtype == "fci":
        ndet = 2**(2*nbasis)
        for n in range(ndet):
            occ = [0] * (2*nbasis)
            for i in range(2*nbasis):
                rem = n%(2**(i+1))
                if rem >= (2**i):
                    occ[i] = 1
            if sum(occ[0:nbasis]) == nalpha:
                if sum(occ[nbasis:(2*nbasis)]) == nbeta:
                    O_CI = numpy.vstack(( O_CI, numpy.asarray(occ) ))

    # COMPUTE HAMILTONIAN, H_CI:
    nconf = O_CI.shape[0]
    H_CI = numpy.zeros((nconf,nconf))
    for i in range(nconf):
        for j in range(i,nconf):
            # DETERMINE NUMBER OF DIFFERING SPIN ORBITAL OCCUPATIONS BETWEEN DETERMINANTS I AND J:
            hole = [] 
            exc = [] 
            for k in range(2*nbasis):
                if O_CI[i,k] != O_CI[j,k]:
                    if O_CI[i,k] == 0:
                        exc.append(k)
                    if O_CI[i,k] == 1:
                        hole.append(k)
            nex = len(hole) 
            # COMPUTE DIAGONAL TERMS:
            if nex == 0:
                hij = 0
                for k in range(nbasis):
                    hij = hij + O_CI[i,k] * H_AMO[k,k]
                    hij = hij + O_CI[i,k+nbasis] * H_BMO[k,k] 
                    for l in range(nbasis):
                        # alpha-alpha (J+K)
                        hij = hij + 0.5 * O_CI[i,k] * O_CI[i,l] * ( Q_AMO[k,k,l,l] - Q_AMO[k,l,l,k] )
                        # beta-beta (J+K)
                        hij = hij + 0.5 * O_CI[i,k+nbasis] * O_CI[i,l+nbasis] * ( Q_BMO[k,k,l,l] - Q_BMO[k,l,l,k] )
                        # alpha-beta (J)
                        hij = hij + 0.5 * O_CI[i,k] * O_CI[i,l+nbasis] * Q_ABMO[k,k,l,l] 
                        hij = hij + 0.5 * O_CI[i,k+nbasis] * O_CI[i,l] * Q_ABMO[l,l,k,k] 
                H_CI[i,j] = hij
            # COMPUTE TERMS FOR DETERMINANTS WITH SINGLE EXCITATIONS:
            if nex == 1:
                m = hole[0]%nbasis
                r = exc[0]%nbasis
                spin = hole[0]/nbasis
                if spin == 0:
                    hij = H_AMO[m,r]
                    for k in range(nbasis):
                        hij = hij + O_CI[i,k] * O_CI[j,k] * ( Q_AMO[m,r,k,k] - Q_AMO[m,k,k,r] )
                        hij = hij + O_CI[i,k+nbasis] * O_CI[j,k+nbasis] * Q_ABMO[m,r,k,k]
                else:
                    hij = H_BMO[m,r]
                    for k in range(nbasis):
                        hij = hij + O_CI[i,k] * O_CI[j,k] * Q_ABMO[k,k,m,r]
                        hij = hij + O_CI[i,k+nbasis] * O_CI[j,k+nbasis] * ( Q_BMO[m,r,k,k] - Q_BMO[m,k,k,r] )
                H_CI[i,j] = hij 
                H_CI[j,i] = hij 
            # COMPUTE TERMS FOR DETERMINANTS WITH DOUBLE EXCITATIONS:
            if nex == 2:
                m = hole[0]%nbasis
                n = hole[1]%nbasis
                r = exc[0]%nbasis
                s = exc[1]%nbasis
                spin1 = hole[0]/nbasis
                spin2 = hole[1]/nbasis
                if spin1 == 0:
                    if spin2 == 0:
                        H_CI[i,j] = Q_AMO[m,r,n,s] - Q_AMO[m,s,n,r] 
                    else:
                        H_CI[i,j] = Q_ABMO[m,r,n,s]
                else:
                    if spin2 == 0:
                        H_CI[i,j] = Q_ABMO[n,s,m,r]
                    else:
                        H_CI[i,j] = Q_BMO[m,r,n,s] - Q_BMO[m,s,n,r]
                H_CI[j,i] = H_CI[i,j]

    #DIAGONALIZE HAMILTONIAN:
    eigenvalues, eigenvectors = linalg.eigh(H_CI)
    E_CI = eigenvalues[0] 
    C_CI = eigenvectors[0]

    return E_CI, C_CI

