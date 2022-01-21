"""
computes first and second derivatives of nuclear repulsion term 
"""
import math
import numpy

def distance( xbohr, ybohr, zbohr, Z ):
    natoms = len(xbohr)
    R = numpy.zeros( ( natoms, natoms ) )
    Vnn = 0.0
    for a in range( natoms ):
        for b in range( natoms ):
            if a == b:
                R[ a, b ] = 0.0
            elif b > a:
                R[ a, b ] = math.sqrt( (xbohr[b]-xbohr[a])*(xbohr[b]-xbohr[a]) + (ybohr[b]-ybohr[a])*(ybohr[b]-ybohr[a]) + (zbohr[b]-zbohr[a])*(zbohr[b]-zbohr[a]) ) 
            else:
                R[ a, b ] = R[ b, a ]
                Vnn = Vnn + Z[a] * Z[b] / R[ a, b ]
    return R, Vnn

# Calculate first-order derivative tensor of Vnn, d1Vnn:
def 1d( xbohr, ybohr, zbohr, Z, R ):
    natoms = len(xbohr)
    d1Vnn = numpy.zeros( ( 3 * natoms ) )
    for a in range( natoms ):
        dx = 0.0
        dy = 0.0
        dz = 0.0
        for b in range( natoms ):
            if a != b:
                #print a, b, R[ a, b ] 
                dx = dx + Z[b] * (xbohr[a]-xbohr[b]) * R[ a, b ]**(-3)
                dy = dy + Z[b] * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-3)
                dz = dz + Z[b] * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-3)
        d1Vnn[ 3*a + 0 ] = -Z[a] * dx
        d1Vnn[ 3*a + 1 ] = -Z[a] * dy
        d1Vnn[ 3*a + 2 ] = -Z[a] * dz
    return d1Vnn

# Calculate second-order derivative tensor of Vnn, d2Vnn:
def 2d( xbohr, ybohr, zbohr, Z, R ):
    natoms = len(xbohr)
    d2Vnn = numpy.zeros( ( 3*natoms, 3*natoms ) )

    # Type I: non-mixed
    for a in range( natoms ):
        dxx = 0.0
        dyy = 0.0
        dzz = 0.0
        for b in range( natoms ):
            if a != b:
                #print a, b, R[ a, b ] 
                dxx = dxx + Z[b] * ( 3 * (xbohr[a]-xbohr[b])**(2) * R[ a, b ]**(-5) - R[ a, b ]**(-3) )
                dyy = dyy + Z[b] * ( 3 * (ybohr[a]-ybohr[b])**(2) * R[ a, b ]**(-5) - R[ a, b ]**(-3) )
                dzz = dzz + Z[b] * ( 3 * (zbohr[a]-zbohr[b])**(2) * R[ a, b ]**(-5) - R[ a, b ]**(-3) )
        d2Vnn[ 3*a + 0, 3*a + 0 ] = Z[a] * dxx
        d2Vnn[ 3*a + 1, 3*a + 1 ] = Z[a] * dyy
        d2Vnn[ 3*a + 2, 3*a + 2 ] = Z[a] * dzz

    # Type II: mixed derivatives that share the same nucleus 
    for a in range( natoms ):
        dxy = 0.0
        dxz = 0.0
        dyz = 0.0
        for b in range( natoms ):
            if a != b:
                #print a, b, R[ a, b ] 
                dxy = dxy + Z[b] * (xbohr[a]-xbohr[b]) * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-5)
                dxz = dxz + Z[b] * (xbohr[a]-xbohr[b]) * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-5)
                dyz = dyz + Z[b] * (ybohr[a]-ybohr[b]) * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-5) 
        d2Vnn[ 3*a + 0, 3*a + 1 ] = 3 * Z[a] * dxy
        d2Vnn[ 3*a + 1, 3*a + 0 ] = 3 * Z[a] * dxy
        d2Vnn[ 3*a + 0, 3*a + 2 ] = 3 * Z[a] * dxz
        d2Vnn[ 3*a + 2, 3*a + 0 ] = 3 * Z[a] * dxz
        d2Vnn[ 3*a + 1, 3*a + 2 ] = 3 * Z[a] * dyz
        d2Vnn[ 3*a + 2, 3*a + 1 ] = 3 * Z[a] * dyz

    # Type III and IV: mixed derivatives that don't share a nucleus
    for a in range( natoms ):
        for b in range( natoms ):
            if a != b:
                d2Vnn[ 3*a + 0, 3*b + 0 ] = Z[a] * Z[b] * ( -R[ a, b ]**(-3) + 3 * (xbohr[a]-xbohr[b]) * (xbohr[a]-ybohr[b]) * R[ a, b ]**(-5) ) 
                d2Vnn[ 3*a + 0, 3*b + 1 ] = 3 * Z[a] * Z[b] * (xbohr[a]-xbohr[b]) * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 0, 3*b + 2 ] = 3 * Z[a] * Z[b] * (xbohr[a]-xbohr[b]) * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 1, 3*b + 0 ] = 3 * Z[a] * Z[b] * (ybohr[a]-ybohr[b]) * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 1, 3*b + 1 ] = Z[a] * Z[b] * ( -R[ a, b ]**(-3) + 3 * (ybohr[a]-ybohr[b]) * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-5) ) 
                d2Vnn[ 3*a + 1, 3*b + 2 ] = 3 * Z[a] * Z[b] * (ybohr[a]-ybohr[b]) * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 2, 3*b + 0 ] = 3 * Z[a] * Z[b] * (zbohr[a]-zbohr[b]) * (xbohr[a]-xbohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 2, 3*b + 1 ] = 3 * Z[a] * Z[b] * (zbohr[a]-zbohr[b]) * (ybohr[a]-ybohr[b]) * R[ a, b ]**(-5)  
                d2Vnn[ 3*a + 2, 3*b + 2 ] = Z[a] * Z[b] * ( -R[ a, b ]**(-3) + 3 * (zbohr[a]-zbohr[b]) * (zbohr[a]-zbohr[b]) * R[ a, b ]**(-5) ) 
    return d2Vnn


