"""
Configures basis sets.
Takes list of atomic numbers and name of basis set as input and
exports all basis set information needed by scf.py
"""
def cb( Z, basis_def ):
    natoms = len(Z)
    nbasis = 0
    alpha  = []
    coef   = []
    gtos   = []
    nuclei = []
    for i in range ( 0, natoms ):
        if ".bas" in basis_def:
            file = open(basis_def,'r')
            nfile = int(file.readline())
            for j in range(nfile):
                line = file.readline()
                words = line.split()
                alpha.append(float(words[1]))
                coef.append(1.0) 
                gtos.append(1)
                nuclei.append(i)
                nbasis = nbasis + 1
            file.close()
        elif Z[i] == 1:
            if basis_def == "STO-2G":
                nbasis = nbasis + 1
                gtos.append(2)
                nuclei.append(i)
                alpha.extend([1.309756377,0.233135974])
                coef.extend([0.430128498,0.678913531])
            elif basis_def == "STO-3G":
                nbasis = nbasis + 1
                gtos.append(3)
                nuclei.append(i)
                alpha.extend([3.42525091,0.62391373,0.16885540])
                coef.extend([0.15432897,0.53532814,0.44463454])
            elif basis_def == "STO-6G":
                nbasis = nbasis + 1
                gtos.append(6)
                nuclei.append(i)
                alpha.extend([35.52322122,6.513143725,1.822142904,0.625955266,0.243076747,0.100112428])
                coef.extend([0.00916359628,0.04936149294,0.16853830490,0.37056279970,0.41649152980,0.13033408410])
            elif basis_def == "6-31G":
                nbasis = nbasis + 2
                gtos.append(3)
                nuclei.append(i)
                alpha.extend([18.7311370,2.8253937,0.6401217])
                coef.extend([0.03349460,0.23472695,0.81375733])
                gtos.append(1)
                nuclei.append(i)
                alpha.extend([0.1612778])
                coef.extend([1.0])
            elif basis_def == "6-31G+":
                nbasis = nbasis + 3
                gtos.append(3)
                nuclei.append(i)
                alpha.extend([18.7311370,2.8253937,0.6401217])
                coef.extend([0.03349460,0.23472695,0.81375733])
                gtos.append(1)
                nuclei.append(i)
                alpha.extend([0.1612778])
                coef.extend([1.0])
                nuclei.append(i)
                gtos.append(1)
                nuclei.append(i)
                alpha.extend([0.036])
                coef.extend([1.0])
            elif basis_def == "6-311G":
                nbasis = nbasis + 3
                gtos.append(3)
                nuclei.append(i)
                alpha.extend([33.8650000,5.0947900,1.1587900])
                coef.extend([0.02549380,0.19037300,0.85216100])
                gtos.append(1)
                nuclei.append(i)
                alpha.extend([0.3258400])
                coef.extend([1.0])
                gtos.append(1)
                nuclei.append(i)
                alpha.extend([0.102741])
                coef.extend([1.0])
            else:
                print "Invalid basis set!"
                print "Currently supported basis sets for H are STO-2G, STO-3G, STO-6G, 6-31G, 6-31G+, and 6-311G."
                exit()
        elif Z[i] == 2:
            if basis_def == "STO-2G":
                nbasis = nbasis + 1
                nuclei.append(i)
                gtos.append(2)
                alpha.extend([2.4328790,0.4330510])
                coef.extend([0.4301280,0.6789140])
            elif basis_def == "STO-3G":
                nbasis = nbasis + 1
                nuclei.append(i)
                gtos.append(3)
                alpha.extend([6.36242139,1.15892300,0.31364979])
                coef.extend([0.15432897,0.53532814,0.44463454])
            elif basis_def == "STO-6G":
                nbasis = nbasis + 1
                nuclei.append(i)
                gtos.append(6)
                alpha.extend([65.98456824,12.09819836,3.384639924,1.162715163,0.451516322,0.185959356])
                coef.extend([0.00916359628,0.04936149294,0.16853830490,0.37056279970,0.41649152980,0.13033408410])
            elif basis_def == "6-311G":
                nbasis = nbasis + 3
                nuclei.append(i)
                gtos.append(3)
                alpha.extend([98.1243000,14.768900,3.3188300])
                coef.extend([0.02874520,0.20806100,0.83763500])
                nuclei.append(i)
                gtos.append(1)
                alpha.extend([0.8740470])
                coef.extend([1.0])
                nuclei.append(i)
                gtos.append(1)
                alpha.extend([0.244564])
                coef.extend([1.0])
            else:
                print "Invalid basis set!"
                print "Currently supported basis sets for He are STO-2G, STO-3G, STO-6G, and 6-311G."
                exit()
        else:
            print "Invalid atomic number! Only H and He supported."
            exit()
    ngtos  = sum( gtos )
    #create list of the identities of the first primitive gtos each basis function 
    igto = []
    count = 0
    for i in range( nbasis ):
        igto.append( count )
        count = count + gtos[i]
    basis_settings = [ natoms, nbasis, ngtos, alpha, coef, gtos, igto, nuclei ] 
    return basis_settings 
