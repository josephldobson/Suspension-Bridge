import bridge
import numpy as np
import matplotlib.pyplot as plt 

N = 20 # Number of nodes on the bridge
pos = 7 # Position on the bridge

nu_vals = np.linspace(0.5,5,2000)# AN array or list of frequency values
max_vals = []

fname = "ResonanceAmp_pos{}.txt".format(pos)
fp = open(fname,"w")

for nu in nu_vals:
  v = np.zeros(N*4) # [X, Z, dX, dZ]
  br = bridge.Bridge_osc(v, dt=1e-3, t0=0, nu=nu, pos=pos)
  br.iterate(60)
  m = br.MaxAmp()
  max_vals.append(m)
  fp.write("{:.3g} {:.8g}\n".format(nu,m))
  print("nu=",nu,"max=",m)

plt.plot(nu_vals,max_vals,"*b")
plt.xlabel(r'$\nu$',fontsize=22)
plt.ylabel(r'$Z_{max}$',fontsize=22)
plt.tight_layout(pad=0.3)
plt.savefig("ResonanceAmp_pos{}.pdf".format(pos))
plt.show()

plt.plot(nu_vals,max_vals,"b")
plt.xlabel(r'$\nu$',fontsize=22)
plt.ylabel(r'$Z_{max}$',fontsize=22)
plt.tight_layout(pad=0.3)
plt.savefig("ResonanceAmp_lpos{}.pdf".format(pos))
plt.show()
