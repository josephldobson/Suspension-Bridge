import numpy as np
from suspensionbridge import *

N=20

sb = SuspensionBridge_xz(N)
sb.prune_nu() # remove null mode
print(sb.Ev)
sb.get_distance_and_depth()
sb.plot_spectrum("spectrum_XZ_N{}.pdf".format(N))
for mode in range(5):
 sb.plot_XZ_modes(mode, "lowest_N{}_XZ_modes{}.pdf".format(N,mode))
print(sb.nu)       