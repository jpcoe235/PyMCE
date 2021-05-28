
import numpy as np





def output_traj(first, T1,t):
    file = 'output_traj_' + str(T1.gettrajid_traj()) + '.dat'
    if first:
        with open(file, 'w') as f:
            str2w = 'time'+' '*(7)
            for i in range(T1.npart):
                for j in range(3):
                    if j == 0:
                        str2w = str2w + 'pos' + 'X' + str(i + 1)+'  '
                    elif j == 1:
                        str2w = str2w + 'pos' + 'Y' + str(i + 1)+'  '
                    elif j == 2:
                        str2w = str2w + 'pos' + 'Z' + str(i + 1)+'  '
            for i in range(T1.npart):
                for j in range(3):
                    if j == 0:
                        str2w = str2w + 'mom' + 'X' + str(i + 1)+'  '
                    elif j == 1:
                        str2w = str2w + 'mom' + 'Y' + str(i + 1)+'  '
                    elif j == 2:
                        str2w = str2w + 'mom' + 'Z' + str(i + 1)+'  '

            f.write(str2w)
            f.write('\n')

    with open(file, 'a+') as f:
        f.write(str(t)[0:6]+' '*(10-len(str(t))))
        for i in range(T1.ndim):
            if len(str(T1.getposition_traj()[i]))>=6:
                f.write(str(T1.getposition_traj()[i])[0:6])
            else:
                f.write(str(T1.getposition_traj()[i])+'0'*(6-len(str(T1.getposition_traj()[i]))))

            f.write(' ')
        for i in range(T1.ndim):
            if len(str(T1.getmomentum_traj()[i])) >= 6:
                f.write(str(T1.getmomentum_traj()[i])[0:6])
            else:
                f.write(str(T1.getmomentum_traj()[i]) + '0' * (6 - len(str(T1.getmomentum_traj()[i]))))

            f.write(' ')
        f.write('\n')
