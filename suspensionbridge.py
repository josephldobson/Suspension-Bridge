import numpy as np
import matplotlib.pyplot as plt
from vibrations_base import VibrationsBase

class SuspensionBridge_x(VibrationsBase):
    """ Class which simulates a guitar string."""
    def __init__(self, N, theta_tower = 30):
        """ Initialise Matrix to calculate X's Normal Modes.
            N : number of nodes including the towers
            theta_tower : the angle at the tower
        """
        self.N=N
        self.Ksc = 5e7 
        self.Kbend = 0
        self.M=15
        self.thetas = np.zeros(self.N, dtype='float64')
        self.set_thetas(theta_tower)
        
        # CALCULATES MATRIX A FOR X OF BRIDGE
        A = np.zeros([N-2,N-2], dtype='float64') # Numpy array ignoring end nodes
        
        # 'If' statements prevent indexing outside of range
        for i in range(1,self.N-1):
            A[i-1,i-1] = -(np.cos(self.thetas[i])**2+np.cos(self.thetas[i+1])**2)
            if i < self.N-2: A[i-1,i] = np.cos(self.thetas[i+1])**2  
            if i > 1: A[i-1,i-2] = np.cos(self.thetas[i])**2
        
        A *= self.Ksc / self.M 
        super().__init__(A)
        
    def get_distance_and_depth(self):
        """ Calculates the distance between the towers and the depth at rest"""
        L = sum([np.cos(self.thetas[i]) for i in range(1,self.N)])
        h = sum([np.sin(self.thetas[i]) for i in range(1,self.N//2)])
        print('The length between the two towers is '+str(L)+'m')
        print('The depth of the bridge is '+str(abs(h))+'m')

    def set_thetas(self, theta_tower = 30):
       """ Initialise the values of theta.
           theta_tower : the angle at the tower.
       """
       self.MoT = 2*np.tan(theta_tower*np.pi/180)/(self.N-2)
       midn = self.N//2
       self.thetas[midn] = 0
       for i in range(1, midn):
           self.thetas[midn-i] = -np.arctan(i*self.MoT)
           self.thetas[midn+i] = np.arctan(i*self.MoT)
      
    def plot_spectrum(self, fname=""):   
        """ Plot the vibration spectrum of the bridge in blue.
            Plot multiples of the fundamental frequency in red.
            fname:  Filename for figure output (if non-empty).
        """
        # Plot actual frequency ration in blue
        Nmax = self.N-2  
        plt.plot(self.nu, "b*")
        # Plot harmonic frequencies in red
        plt.plot(range(1, Nmax+1)*self.nu[0], "r-") 
        plt.xlabel("i", fontsize=24)
        plt.ylabel(r'$\nu$', fontsize=24)
        plt.tight_layout(pad=0.4)
        if(fname != ""):
           plt.savefig(fname)
        plt.show()



class SuspensionBridge_z(SuspensionBridge_x):
    def __init__(self, N, theta_tower=30):
        """ Initialise Matrix to calculate Z's Normal Modes.
            N : number of nodes including towers.
            theta_tower : the angle at the tower.
        """
        self.N=N
        self.Ksc = 5e7
        self.Kbend = 0
        self.M=15
        self.thetas = np.zeros(self.N, dtype='float64')
        self.set_thetas(theta_tower)
        
        # CALCULATES MATRIX A FOR Z OF BRIDGE

        A = np.zeros([N-2,N-2],dtype='float64') # Indexing ignoring end nodes
        
        # 'If' statements prevent indexing outside of range
        for i in range(1,self.N-1):  
            A[i-1,i-1] = -(np.sin(self.thetas[i])**2+np.sin(self.thetas[i+1])**2)
            if i < self.N-2: A[i-1,i] = np.sin(self.thetas[i+1])**2 
            if i > 1: A[i-1,i-2] = np.sin(self.thetas[i])**2
            
        A *= self.Ksc / self.M
        VibrationsBase.__init__(self,A)
        


class SuspensionBridge_xz(SuspensionBridge_x):
    def __init__(self, N, theta_tower=30):
        """ Initialise Matrix to calculate XZ's Normal Modes.
            N : number of nodes including towers.
            theta_tower : the angle at the tower.
        """
        self.N=N
        self.Ksc = 5e7
        self.Kbend = 1e4
        self.M=15
        self.thetas = np.zeros(self.N, dtype='float64')
        self.set_thetas(theta_tower)

        # CALCULATES MATRIX A FOR XZ BRIDGE

        k = self.Kbend/self.Ksc
        n = self.N-2
        A = np.zeros([2*n,2*n], dtype='float64')
        #X[i] = V[i]; Z[i] = V[n+i]
        
        for i in range (1,n+1): # Not computing end nodes
            # (i == j)
            A[i-1,i-1]     = - (np.cos(self.thetas[i])**2+np.cos(self.thetas[i+1])**2)
            A[n+i-1,n+i-1] = - (np.sin(self.thetas[i])**2+np.sin(self.thetas[i+1])**2) - k
            A[i-1,n+i-1]   = - np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                             - np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            A[n+i-1,i-1]   = - np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                             - np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            # (i+1 == j)
            if i < n: 
                A[i-1,i]     = np.cos(self.thetas[i+1])**2
                A[n+i-1,n+i] = np.sin(self.thetas[i+1])**2 +0.5*k 
                A[i-1,n+i]   = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
                A[n+i-1,i]   = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            # (i-1 == j)
            if i > 1:
                A[i-1,i-2]     = np.cos(self.thetas[i])**2 
                A[n+i-1,n+i-2] = np.sin(self.thetas[i])**2  +0.5*k  
                A[n+i-1,i-2]   = np.cos(self.thetas[i])*np.sin(self.thetas[i])
                A[i-1,n+i-2]   = np.cos(self.thetas[i])*np.sin(self.thetas[i])
                
        A = A*self.Ksc/self.M
        VibrationsBase.__init__(self,A)