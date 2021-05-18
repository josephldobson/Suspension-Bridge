import numpy as np
import matplotlib.pyplot as plt 
from bridge import Bridge


N = 20
v = np.zeros(N*4) # [X, Z, dX, dZ]
br = Bridge(v, dt=1e-3, t0=0)

br.iterate(100)
# Time profile of node 1/4 of bridge length from the end 
br.plot(br.N+br.N//2,0)
plt.savefig("A")
plt.show()
plt.close()

# Profile of bridge after 0.1 second
br.plot_XZ_profile(0.1)
plt.savefig("B")
plt.show()
# Profile of bridge after 1 second
br.plot_XZ_profile(1)
plt.savefig("C")
plt.show()
# Profile of bridge after 10 seconds
br.plot_XZ_profile(10)
plt.savefig("D")
plt.show()
# Profile of bridge after 100 seconds
br.plot_XZ_profile(100)
plt.savefig("E")
plt.show()