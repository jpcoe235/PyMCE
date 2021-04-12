import numpy as np
import src.initialize_traj as it


def coupdotvel(T, i, j):
    if i == j:
        coup = 0.0
    else:
        coup = np.dot(T.getvelocity_traj(), T.getcoupling_traj( i, j))

    return coup
