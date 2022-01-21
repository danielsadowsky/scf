"""
computes 1- and 2- electron integrals
"""

import math
import numpy as np
from scipy import misc, special 

def integrals( R, X, Z, basis_settings ):
    # basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ]
    natoms = basis_settings[0]
    nbasis = basis_settings[1]
    ngtos = basis_settings[2]
    alpha = basis_settings[3]
    coef = basis_settings[4]
    gtos = basis_settings[5]
    igto = basis_settings[6]
    nuclei = basis_settings[7]

    # COMPUTE OVERLAP MATRIX, S, KINETIC ENERGY MATRIX, T, AND NUCLEAR-ELECTRONIC POTENTIAL, Vne
    S = numpy.zeros( ( nbasis, nbasis ) )
    s = numpy.zeros( ( ngtos,  ngtos  ) )
    T = numpy.zeros( ( nbasis, nbasis ) )
    V = numpy.zeros( ( nbasis, nbasis ) )
    for a in range( nbasis ):
        for b in range( nbasis ):
            n1 = nuclei[ a ]
            n2 = nuclei[ b ]
            if b < a:
                sab = S[ b, a ]
                tab = T[ b, a ]
                vab = V[ b, a ]
            else:
                rab = R[ n1, n2 ]
                sab = 0.0
                tab = 0.0
                vab = 0.0
                for i in range( gtos[a] ):
                    for j in range( gtos[b] ):
                        #print(i,j)
                        p1 = igto[ a ] + i
                        p2 = igto[ b ] + j
                        pab = alpha[p1] + alpha[p2]
                        kp = alpha[p1] * alpha[p2] / pab
                        sij = coef[p1] * coef[p2] * ( 4.0 * kp / pab )**(0.75) * math.exp( -kp * rab * rab )
                        sab = sab + sij
                        tab = tab + kp * ( 3.0 - 2.0 * kp * rab * rab ) * sij
                        for c in range( natoms ):
                            if n1 == n2 and n2 == c:
                                vab = vab - 2.0 * math.sqrt( pab / math.pi ) * Z[c] * sij
                            else:
                                rpc = 0.0
                                for k in range(3):
                                    rpc = rpc + ( X[c,k] - ( alpha[p1] * X[n1,k] + alpha[p2] * X[n2,k] ) / pab ) * ( X[c,k] - ( alpha[p1] * X[n1,k] + alpha[p2] * X[n2,k] ) / pab ) 
                                rpc = math.sqrt( rpc )
                                if rpc == 0.0:
                                    vab = vab - 2.0 * math.sqrt( pab / math.pi ) * Z[c] * sij
                                else:
                                    vab = vab - Z[c] / rpc * math.erf( math.sqrt( pab ) * rpc ) * sij
                        s[ p1, p2 ] = sij
                        s[ p2, p1 ] = sij
            S[ a, b ] = sab
            T[ a, b ] = tab
            V[ a, b ] = vab

    # COMPUTE TWO-ELECTRON INTEGRALS, G (using chemist's notation):
    Q = numpy.zeros( ( ( ( nbasis, nbasis, nbasis, nbasis ) ) ) )
    for a in range( nbasis ):
        for b in range( nbasis ):
            for c in range( nbasis ):
                for d in range( nbasis ):
                    n1 = nuclei[ a ]
                    n2 = nuclei[ b ]
                    n3 = nuclei[ c ]
                    n4 = nuclei[ d ]
                    qab = 0.0
                    if ( n1 == n3 and n2 == n4 ) or ( n1 == n4 and n2 == n3 ):
                        for i in range( gtos[ a ] ):
                            for j in range( gtos [ b ] ):
                                for k in range( gtos [ c ] ):
                                    for l in range( gtos [ d ] ):
                                        p1 = igto[ a ] + i
                                        p2 = igto[ b ] + j
                                        p3 = igto[ c ] + k
                                        p4 = igto[ d ] + l
                                        kpq = ( alpha[p1] + alpha[p2] ) * ( alpha[p3] + alpha[p4] ) / ( alpha[p1] + alpha[p2] + alpha[p3] + alpha[p4] )
                                        rpq = "zero!"
                                        qab = qab + 2.0 * math.sqrt( kpq / math.pi ) * s[ p1, p2 ] * s[ p3, p4 ]
                    else:
                        for i in range ( gtos[ a ] ):
                            for j in range( gtos [ b ] ):
                                for k in range( gtos [ c ] ):
                                    for l in range( gtos [ d ] ):
                                        p1 = igto[ a ] + i
                                        p2 = igto[ b ] + j
                                        p3 = igto[ c ] + k
                                        p4 = igto[ d ] + l
                                        ip = 1 / ( alpha[p1] + alpha[p2] )
                                        iq = 1 / ( alpha[p3] + alpha[p4] )
                                        kpq = ( alpha[p1] + alpha[p2] ) * ( alpha[p3] + alpha[p4] ) / ( alpha[p1] + alpha[p2] + alpha[p3] + alpha[p4] )
                                        rpq = 0.0
                                        for m in range(3):
                                            rpq = ( ip * ( alpha[p1] * X[n1,m] + alpha[p2] * X[n2,m] ) - iq * ( alpha[p3] * X[n3,m] + alpha[p4] * X[n4,m] ) ) * ( ip * ( alpha[p1] * X[n1,m] + alpha[p2] * X[n2,m] ) - iq * ( alpha[p3] * X[n3,m] + alpha[p4] * X[n4,m] ) )  
                                        rpq = math.sqrt( rpq )  
                                        #rpq = math.sqrt( ( ip * ( alpha[p1] * xbohr[n1] + alpha[p2] * xbohr[n2] ) - iq * ( alpha[p3] * xbohr[n3] + alpha[p4] * xbohr[n4] ) )**(2.0) + ( ip * ( alpha[p1] * ybohr[n1] + alpha[p2] * ybohr[n2] ) - iq * ( alpha[p3] * ybohr[n3] + alpha[p4] * ybohr[n4] ) )**(2.0) + ( ip * ( alpha[p1] * zbohr[n1] + alpha[p2] * zbohr[n2] ) - iq * ( alpha[p3] * zbohr[n3] + alpha[p4] * zbohr[n4] ) )**(2.0) )
                                        if rpq == 0.0:
                                            qab = qab + 2.0 * math.sqrt( kpq / math.pi ) * s[ p1, p2 ] * s[ p3, p4 ]
                                        else:
                                            qab = qab + s[ p1, p2 ] * s[ p3, p4 ] * math.erf( kpq * rpq ) / rpq
                    G[ a, b, c, d ] = qab

    return S, T, V, G

def dfact(x):
    if type(x) != int:
        print("dfact only takes integers as arguments!")
        exit()
    elif x < -1:
        print("dfact only takes arguments > -1 !")
        exit()
    elif x == -1:
        y = 1
    else:
        y = misc.factorial( x, exact = True )
    return y

def curved( R, Z, alpha, nuclei, curved ):
    nbasis = len(alpha)
    L = [ [ 0 ] * 3 ] * nbasis
    N = [ math.sqrt( ( 2 * alpha[i] / math.pi )**(1.5) * ( 4*alpha[i] )**(sum(L[i])) / dfact(2*L[i][0]-1) / dfact(2*L[i][1]-1) / dfact(2*L[i][2]-1) ) for i in range(nbasis) ]
    S = np.zeros( ( nbasis, nbasis ) )
    for a in range( nbasis ):
        for b in range( nbasis ):
            if b < a:
                S[ a, b ] = S[ b, a ]
            else: 
                gamma = alpha[a] + alpha[b]
                if nuclei[a] == nuclei[b]:
                    rab = 0.0
                    rp = [ 0.0, 0.0 ]
                else:
                    rab = R 
                    P = ( alpha[a] * nuclei[a] * R + alpha[b] * nuclei[b] * R ) / gamma
                    rp = [ P, P - R ]
                if curved == False: 
                    Sx = [ 0.0 ] * 3
                    for x in range(3):
                        # WHY +1 ?????
                        for k in range( ( L[a][x] + L[b][x] ) // 2 + 1):
                            ck = 0.0
                            # WHY +1 ??????
                            for la in range( L[a][x] + 1 ):
                                for lb in range( L[b][x] + 1 ):
                                    if la + lb == 2*k:
                                        # rp should equal 0 for x and z !!!!
                                        ck = ck + special.binom(L[a][x],la) * special.binom(L[b][x],lb) * rp[0]**(L[a][x]-la) * rp[1]**(L[b][x]-lb)
                            Sx[x] = Sx[x] + ck * math.sqrt( math.pi / gamma ) * dfact(2*k-1) / ( 2 * gamma )**k 
                    S[ a, b ] = N[a] * N[b] * math.exp( -alpha[a] * alpha[b] / gamma * rab * rab ) * Sx[0] * Sx[1] * Sx[2]
                elif curved == True:
                    S[ a, b ] = N[a] * N[b] * math.exp( -alpha[a] * alpha[b] / gamma * rab * rab )
    return S
 
