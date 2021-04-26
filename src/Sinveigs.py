import numpy as np
import numpy.linalg as sla


def Sinv2(A):
    eigvals, eigvecs = sla.eig(A)

    thresh = 0.0001
    CTmpMat = np.zeros_like(A)
    CTmpMat5pi = np.zeros_like(A)
    alpha = np.max(eigvals) / 2.0 * thresh
    for i in range(np.size(eigvals)):
        DStemp = eigvals[i] / (eigvals[i] ** 2 + alpha)
        DStemp5pi = np.sqrt(abs(DStemp))

        for j in range(np.size(eigvals)):
            CTmpMat[i, j] = DStemp * np.conj(eigvecs[j, i])
            CTmpMat5pi[i, j] = DStemp5pi * np.conj(eigvecs[j, i])

    Sinvmat = np.matmul(eigvecs, CTmpMat)

    print(np.matmul(Sinvmat,A))
    return Sinvmat
