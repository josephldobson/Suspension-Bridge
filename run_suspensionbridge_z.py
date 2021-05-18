
import numpy as np
from suspensionbridge import *

N=20

sb = SuspensionBridge_z(N)
sb.plot_spectrum("spectrum_Z_N{}.pdf".format(N))
n_mode = 4
sb.plot_modes(n_mode, "lowest_N{}_z_modes{}.pdf".format(N,n_mode))
print(sb.nu)       
