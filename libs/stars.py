from __future__ import print_function

from numpy import pi, sqrt, array, newaxis, ndarray

from libs.const import G, AU

class StarSystem:
    """ Class for creating a single star, or binary system rotating
         contraclockwise.

        Arguments:
            m0  : mass of star A (with units)
            m1  : mass of star B (with units)
            a   : semi-major axis of the binary (with units)
            e   : excentricity of the binary
            cm  : coordinates of the binary's center of mass (with units)
            apo : bool indicating if the binary is at apoastron (True)
                  or at periastron (False)
    """
    def __init__(self, m1=1, m2=1, a=1, e=0, cm=[0,0], rad=0, apo=True):
        
        M = float(m1 + m2)
        T = 2 * pi * sqrt(a**3 / (G * M))
        
        if M == 0:
            print('Binary system not generated.')
            N = pos = vel = mass = rad = 0

        elif (a == 0) | (m1 == 0) | (m2 == 0):
            print('Creating just one sink particle...')
            pos  = array(cm)[newaxis,:]
            vel  = array([0, 0])[newaxis,:]
            mass = M
            N    = 1
            rad  = array([rad])
        
        else:
            print('Creating binary system...')
            # auxiliar values
            if apo:
                eu = 1 + e
                ed = 1 - e
            else:
                eu = 1 - e
                ed = 1 + e
            
            aux_m = max(m1, m2)
            m2    = min(m1, m2)
            m1    = aux_m
            
            aux_p = a * eu / M
            aux_v = 2 * pi / T * a / M * sqrt(ed / eu)
            
            pos   = array([[m2 * aux_p, 0], 
                           [-m1 * aux_p, 0]]) + array(cm)
            vel   = array([[0, m2 * aux_v], 
                           [0, -m1 * aux_v]])
            mass  = array([m1, m2])
            N     = 2
            
            print(' Binary period: {:.3f} [days]'.format(T / 86400.))
            
            if isinstance(rad, (int, float)): rad = array([rad, rad])
            elif isinstance(rad, (ndarray, list)) and (len(rad) != 2):
                raise ValueError('Sinks radii array must be 2D.')
            else: rad = array([0, 0])
            if not all(rad):
                #Extra: Calculating m1 L1 sphere and m2 Hill sphere
                mu = m2 / M
                print(' mu = {:.3f}'.format(mu))
                if mu < 0.5:
                    delta  = (mu / (3 * (1 - mu)))**(1 / 3.)
                    Rhill1 = (1 - delta + delta**2 / 3. - delta**3 / 9. -\
                               23 * delta**4 / 81.) * a
                    Rhill2 = a * (m2 / (3. * m1))**(1 / 3.)
                    rad = array([Rhill1, Rhill2])
            if all(rad):
                print(' Sinks radii:')
                print('  R1 = {:.4f} AU'.format(rad[0] / AU))
                print('  R2 = {:.4f} AU'.format(rad[1] / AU))

             
        
        self.Mbin  = M
        self.Tbin  = T
        self.pos   = pos
        self.vel   = vel
        self.mass  = mass
        self.Nbin  = N
        self.Rbin  = rad
        
