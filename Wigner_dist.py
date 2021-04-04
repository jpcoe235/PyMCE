import re
import numpy as np
from input import initgeom as ig
import Constants


def read_freq(file='freq.out', ndf=18, linear=False):
    if linear:
        nc = ndf - 5
        numcut = 5
    else:
        nc = ndf - 6
        numcut = 6

    var = 1
    str2search = 'Mass Weighted 2nd Derivative Matrix Eigenvalues'
    str2search2 = "Mass Weighted 2nd Derivative Matrix Eigenvectors"
    str2search3 = 'Low Vibration'
    eigenval = []
    eigenvec = []
    ndf_var = 0
    ndf_var2 = 0
    ndf_var3 = 0
    with open(file, 'r') as f:
        for line in f:
            if str2search in line:
                var = 2
                continue
            if var == 2:
                n = re.findall("\d+\.\d+", line)
                if n:
                    eigenval.extend(n)

            if str2search2 in line:
                var = 3
                continue
            if var == 3:
                n = re.findall("\d+\.\d+", line)

                if n:

                    if ndf_var <= ndf - 1:
                        eigenvec.append(n)
                        ndf_var += 1
                    elif ndf_var2 <= ndf - 1:
                        eigenvec[ndf_var2].extend(n)
                        ndf_var2 += 1
                    else:
                        eigenvec[ndf_var3].extend(n)
                        ndf_var3 += 1

            if str2search3 in line:
                var = 1
                continue

    eigenval = np.asarray(eigenval)
    eigenvec = np.asarray(eigenvec)

    return eigenval[numcut:], eigenvec[:, :]


def WignerSampling():
    ph = Constants.physconst()

    geo = ig()
    eigval, eigvec = read_freq('freq.out', geo.ndf, False)

    nc = np.size(eigval)
    inmod=range(6,nc+1)
    alpha=np.zeros(nc)
    for i in range(nc):
        omega = np.sqrt(float(eigval[i]) / ph.amu)
        alpha[i]=0.5*omega


