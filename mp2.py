"""
computes MP2 energies from RHF orbitals 
"""

# this only works with a RHF reference. E2 energies are reasonable, not sure about MP3. 

import math
import numpy
from numpy import linalg

def mp2( nalpha, nbasis, Q_MO, C, e ):

    E2 = 0 
    for a in range(nalpha):
        for b in range(nalpha):
            for r in range(nalpha,nbasis):
                for s in range(nalpha,nbasis):
                    E2 = E2 + ( 2 * Q_MO[a,r,b,s] * Q_MO[r,a,s,b] - Q_MO[a,r,b,s] * Q_MO[r,b,s,a] ) / ( e[a] + e[b] - e[r] - e[s] )
    E3 = 0
    for a in range(nalpha):
        for b in range(nalpha):
            for c in range(nalpha):
                for d in range(nalpha):
                    for r in range(nalpha,nbasis):
                        for s in range(nalpha,nbasis):
                            E3 = E3 + 0.125 * (Q_MO[a,r,b,s]-Q_MO[a,s,b,r]) * (Q_MO[c,a,d,b]-Q_MO[c,b,d,a]) * (Q_MO[r,c,s,d]-Q_MO[r,d,s,c]) / (e[a]+e[b]-e[r]-e[s]) / (e[c]+e[d]-e[r]-e[s])
    for a in range(nalpha):
        for b in range(nalpha):
            for r in range(nalpha,nbasis):
                for s in range(nalpha,nbasis):
                    for t in range(nalpha,nbasis):
                        for u in range(nalpha,nbasis):
                            E3 = E3 + 0.125 * (Q_MO[a,r,b,s]-Q_MO[a,s,b,r]) * (Q_MO[r,t,s,u]-Q_MO[r,u,s,t]) * (Q_MO[t,a,u,b]-Q_MO[t,b,u,a]) / (e[a]+e[b]-e[r]-e[s]) / (e[a]+e[b]-e[t]-e[u])
    for a in range(nalpha):
        for b in range(nalpha):
            for c in range(nalpha):
                for r in range(nalpha,nbasis):
                    for s in range(nalpha,nbasis):
                        for t in range(nalpha,nbasis):
                            E3 = E3 + (Q_MO[a,r,b,s]-Q_MO[a,s,b,r]) * (Q_MO[c,t,s,b]-Q_MO[c,b,s,t]) * (Q_MO[r,a,t,c]-Q_MO[r,c,t,a]) / (e[a]+e[b]-e[r]-e[s]) / (e[a]+e[c]-e[r]-e[t])

    print "MP2:                 ", E2
    print "MP3:                 ", E3
    
    return E2




