import numpy as np
import matplotlib.pyplot as plt 
import ode_rk4
from bridge_base import Bridge_base


class Bridge(Bridge_base):
   def __init__(self, V0=[0], dt=0.1, t0=0, theta_tower=30):
       super().__init__(V0, dt, t0, theta_tower)
      
   def F(self, t, V):
      """Solves 4(N-2) ODE's based on derived equation of motion.
         t : time since bridge at rest.
         V : array containing values of X,Z,VX,VZ.
      """
      # Array to store the right hand side of the equation
      eq = np.zeros(self.N*4, dtype='float64')
      # A few array refernces to make the code easier to read
      X = V[:self.N]
      Z = V[self.N:2*self.N]
      VX = V[2*self.N:3*self.N]
      VZ = V[3*self.N:]
      
      FX = eq[:self.N]
      FZ = eq[self.N:2*self.N]
      FVX = eq[2*self.N:3*self.N]
      FVZ = eq[3*self.N:]

      cos_th = np.cos(self.thetas)
      sin_th = np.sin(self.thetas)
      I = np.index_exp[1:self.N-1]
      Im1 = np.index_exp[0:self.N-2]
      Ip1 = np.index_exp[2:self.N]
      L = self.Load(t)
      # IMPROVED SPEED BY USING NO LOOPS 

      FX[I] = VX[I]
      FZ[I] = VZ[I]
      
      FVX[I] = self.Ksc/self.M * (
            (Z[Im1]-Z[I]) * sin_th[I] * cos_th[I]
          + (Z[Ip1]-Z[I]) * sin_th[Ip1] * cos_th[Ip1]
          + (X[Im1]-X[I]) * cos_th[I]**2
          + (X[Ip1]-X[I]) * cos_th[Ip1]**2 )\
          - self.Gamma*VX[I]/self.M    
 
      FVZ[I] = self.Ksc/self.M * (
            (X[Im1]-X[I]) * sin_th[I] * cos_th[I]
          + (X[Ip1]-X[I]) * sin_th[Ip1] * cos_th[Ip1]
          + (Z[Im1]-Z[I]) * sin_th[I]**2
          + (Z[Ip1]-Z[I]) * sin_th[Ip1]**2 )\
          + self.Kbend/self.M * ((Z[Im1]+Z[Ip1])*0.5 - Z[I])\
          - self.Gamma/self.M * VZ[I] + L[I]/self.M
          
      return(eq)

   def Load(self,t):
       """ Return a constant load at the middle of the bridge at time t.
           t : time.
       """
       f_load = np.zeros(self.N, dtype='float64')
       # An 80kg (800N) person at the centre of the bridge
       f_load[self.N//2-1] = f_load[self.N//2] = -200 
       return(f_load)
         

class Bridge_osc(Bridge):
   def __init__(self, V0=[0], dt=0.1, t0=0, nu=1, pos=0):
       super().__init__(V0, dt, t0)
       self.nu = nu # frequency of jumping
       self.pos = pos # index of segment for load
       
   def Load(self, t):
       """ Returns load of man jumping on section pos of the bridge at frequency nu.
           t : time.
       """
       f_load = np.zeros(self.N, dtype='float64')
       # An 80kg (800N) person at position pos 
       if np.sin(2 *np.pi * self.nu * t) > 0:
           f_load[self.pos-1] = f_load[self.pos] = -200 * np.sin(2*np.pi*self.nu*t)
       else:  #No force excerted whilst in the air
           f_load[self.pos-1] = f_load[self.pos] = 0
       return(f_load)

   def MaxAmp(self):
       """ Returns maximum vertical displacement reached by any node at any time
           during simulation.
       """
       MA = [i[self.N:2*self.N] for i in self.V_list]
       MA = map(abs, MA)
       MA = [max(i) for i in MA]
       return(max(MA))