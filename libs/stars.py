from __future__ import print_function

from numpy import pi, sqrt, array, full, newaxis

from libs.const import G

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
    def __init__(self, m1=1., m2=1., a=1., e=0., cm=[0.,0.], apo=True):
        
        M = m1 + m2
        T = 2. * pi * sqrt(a**3 / (G * M))
        
        if M == 0:
            print('Binary system not generated.')
            N = 0
        
        elif (a == 0) | (m1 == 0) | (m2 == 0):
            print('Creating just one star...')
            pos  = array(cm)[newaxis,:]
            vel  = array([0., 0.])[newaxis,:]
            mass = M
            N    = 1
        
        else:
            print('Creating binary system...')
            # auxiliar values
            if apo:
                eu = 1. + e
                ed = 1. - e
            else:
                eu = 1. - e
                ed = 1. + e
                
            aux_p = a * eu / M
            aux_v = 2. * pi / T * a / M * sqrt(ed / eu)
            
            pos   = array([[m2 * aux_p, 0], 
                           [-m1 * aux_p, 0]]) + array(cm)
            vel   = array([[0, m2 * aux_v], 
                           [0, -m1 * aux_v]])
            mass  = array([m1, m2])
            N     = 2
        
            print('Binary period: {:.3f} [days]'.format(T/86400.))
        
        self.Mbin  = M
        self.Tbin  = T
        self.pos   = pos
        self.vel   = vel
        self.mass  = mass
        self.Nbin  = N
        
