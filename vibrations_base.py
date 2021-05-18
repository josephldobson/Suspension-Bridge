import numpy as np
import math
import matplotlib.pyplot as plt

class VibrationsBase:
    def __init__(self, C):
        """ Solve the eigenvalue problem C v = lambda v.
            C:    Input matrix.
            name: Name for this type of string (used for figures).
        """
        self.Mat= C
        self.n_vec = C.shape[0]
        (Ev , V) = np.linalg.eig(C)
        idx = Ev.argsort()[::-1]   
        self.Ev = Ev[idx]
        self.V = V[:,idx]
        self.nu = np.sqrt(-self.Ev)/(2*math.pi)
        self.pltformats = ["r-", "g-", "b-", "c-", "m-", "y-", "k-"]

    def prune_nu(self):
        pruned_nu = []

        i = 0
        while (abs(self.nu[i]) < 1e-4): i += 1
        self.V = self.V[:,i:]
        self.nu = np.real(self.nu[i:]) # some nu are complex
        
    def index_of_nu(self, nu):
        pass

    def plot_modes(self, Nm, fname=""):
        """ Plot the lowest Nm mode profiles of the string.
            Nm:    Number of modes to display.
            fname: Filename for figure output (if non-empty).
        """
        for i in range(0, Nm):
            # We must add the (fixed) end nodes by hand.
            W = np.pad(self.V[:,i], (1,1), 'constant', constant_values=(0, 0))
            plt.plot(W, self.pltformats[i%len(self.pltformats)],
                    label=r'$\nu={1:.3g}$'.format(i,self.nu[i]))

        plt.legend(loc='lower right', bbox_to_anchor=(1.20, 0.0), prop={'size':10})

        ax = plt.gca()
        ax.set_aspect(1.0/ax.get_data_ratio()*0.5)
        plt.tight_layout(rect=[0, 0, 1.1, 1])
        plt.xlabel("i", fontsize=24)
        plt.ylabel("y", fontsize=24)
        plt.tight_layout(pad=0.1)
        if(fname != ""):
            plt.savefig(fname)
        plt.show()

    def plot_XZ_modes(self, mode, fname=""):
        """ Plot the lowest Nm mode profiles of the string.
            Mode:  Index to display.
            XorZ: Which deformation toplot X or Z?  
            fname: Filename for figure output (if non-empty).
        """
        X = self.V[:self.n_vec//2,:]
        Z = self.V[self.n_vec//2:,:]
        # Add zeros for nodes on the tower
        W = np.pad(X[:,mode], (1,1), 'constant', constant_values=(0, 0))
        plt.title(r'$Mode: {}\,\quad \nu={:.3g}$'.format(mode, self.nu[mode]))
        plt.plot(W, "b*", label="X")
        W = np.pad(Z[:,mode], (1,1), 'constant', constant_values=(0, 0))
        plt.plot(W, "r*", label="Z")

        plt.legend(loc='lower right', bbox_to_anchor=(1.20, 0.0), prop={'size':10})
           # W = np.pad(v[:,i], (1,1), 'constant', constant_values=(0, 0))
           # plt.plot(W, self.pltformats[i%len(self.pltformats)],
           #label=r'$\nu={1:.3g}$'.format(i,self.nu[i]))

        ax = plt.gca()
        ax.set_aspect(1.0/ax.get_data_ratio()*0.5)
        plt.tight_layout(rect=[0, 0, 1.1, 1])
        plt.xlabel("i", fontsize=24)
        plt.ylabel("X/Z", fontsize=24)
        plt.tight_layout(pad=0.1)
        if(fname != ""):
            plt.savefig(fname)
        plt.show()
        
    def plot_spectrum(self, fname=""):   
        """ Plot the vibration spectrum of the bridge in blue.
            Plot multiples of the fundamental frequency in red.
            fname:  Filename for figure output (if non-empty).
        """
        pass
    
    def mk_figs(self, mode, spec=True):
        """ Plot normal modes, spectrum and anharmonicity.
            mode: Index of normal mode profiles to plot.
        """
        if(mode > self.n_vec-1):
            mode = self.n_vec-1
        self.plot_modes(mode, "lowest_N{}_mode{}.pdf".format(self.N,mode))
        if(spec):
           self.plot_spectrum("spectrum_N{}.pdf".format(self.N))


