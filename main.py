import input
import Wigner_dist
import numpy as np
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
init=input.initgeom()

dynin=input.initdyn()


h,n= Wigner_dist.read_freq('freq.out',ndf=18)
print(np.shape(n))
