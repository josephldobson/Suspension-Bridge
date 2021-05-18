import numpy as np
from suspensionbridge import *

N=20

sb = SuspensionBridge_x(N)
sb.plot_spectrum("spectrum_X_N{}.pdf".format(N))
#n_mode = 2
#sb.plot_modes(n_mode, "lowest_N{}_x_modes{}.pdf".format(N,n_mode))
n_mode = 4
sb.plot_modes(n_mode, "lowest_N{}_x_modes{}.pdf".format(N,n_mode))
#n_mode = 10
#sb.plot_modes(n_mode, "lowest_N{}_x_modes{}.pdf".format(N,n_mode))
print(sb.nu)       
