from __future__ import print_function
from sys import exit

from numpy import all, any, sum, cos, sin, pi, sqrt
from numpy import array, newaxis, ndarray
from numpy import empty, full, where
from numpy import transpose, unique 

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
            rad : (Accretion) Radius of each star (with units).
            apo : bool indicating if the binary is at apoastron (True)
                  or at periastron (False)
    """
    def __init__(self, m1=1, m2=1, a=1, e=0, cm=[0,0], rad=0, apo=True):
        
        M  = float(m1 + m2)
        T  = 2 * pi * sqrt(a**3 / (G * M))
        cm = array(cm)
        
        if M == 0:
            print('Binary system not generated.')
            N = mass = rad = 0
            pos = vel = None

        elif (a == 0) | (m1 == 0) | (m2 == 0):
            print('Creating just one sink particle...')
            pos  = cm[newaxis,:]
            vel  = array([0, 0])[newaxis,:]
            mass = array([M])
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
                           [-m1 * aux_p, 0]]) + cm
            vel   = array([[0, m2 * aux_v], 
                           [0, -m1 * aux_v]])
            mass  = array([m1, m2])
            N     = 2
            
            print(' Binary period: {:.3f} [days]'.format(T / 86400.))
            
            if isinstance(rad, (int, float)): rad = array([rad, rad])
            elif isinstance(rad, (ndarray, list)) and (len(rad) != 2):
                print('Sinks radii array must be scalar or 2D. Exiting.')
                exit()
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
        self.CM    = cm
        
        
        
    def add_planets(self, m=0, a=0, e=0, f=0, cs=0, rad=0):
        """ Function for creating a planetary system around
            the existing star system.
            
            Arguments:
                m   : mass of the planets (with units).
                a   : semi-major axis of the planets (with units).
                e   : excentricity of the planets.
                f   : true anomaly of the planets in degrees, 
                      considering varpi=0 (argument of the perihelio).
                cs  : planetary orbits's focus.
                       0 : star system's center of mass (~jacobi).
                       1 : most massive (or only) star of the system.
                           if equal stars's mass, the right one is opted.
                       2 : least massive star of the system.
                rad : (Accretion) Radius of each planet (with units).
        """
        
        m = array([m])
        a = array([a])
        if (not all(m)) | (not all(a)) | (self.Nbin == 0):
            print(" No planets created")
            m = rad = array([])
            pos = vel = None

        else:
            print(" Adding planetary system...")
            if any(e >= 1) | any(e < 0):
                print("Planetary eccentricities must be within [0, 1). Exiting.")
                exit()
            if any(a <= 0):
                print("Planetary semi-major axis must be > 0. Exiting.")
                exit()
            
            cs = array([cs])
            Ms = empty(len(cs))
            Rs = empty((len(cs), 2))
            Vs = empty((len(cs), 2))
            
            if any([x not in [0, 1, 2] for x in cs]):
                print("Planetary orbit's focus must be 0, 1 or 2. Exiting.")
                exit()
            
            Ms[where(cs == 0)[0]] = sum(self.mass)
            Rs[where(cs == 0)[0]] = self.CM
            Vs[where(cs == 0)[0]] = array([0,0])            

            Ms[where(cs == 1)[0]] = self.mass[0]
            Rs[where(cs == 1)[0]] = self.pos[0]
            Vs[where(cs == 1)[0]] = self.vel[0]
            
            if any(cs == 2):
                if self.Nbin < 2:
                    print("Missing 2nd star for planets to orbit. Exiting")
                    exit()
            
                Ms[where(cs == 2)[0]] = self.mass[1]
                Rs[where(cs == 2)[0]] = self.pos[1]
                Vs[where(cs == 2)[0]] = self.vel[1]

            if any(m >= Ms):
                print("WARNING: Creating a planet with more mass than")
                print("          it's central star. Make sure this is")
                print("          what you want.") 
        
            f  /= (pi / 180.)
            r   = a * (1 - e**2) / (1 + e * cos(f))
            aux = sqrt(G * (Ms + m) / (a * (1 - e**2)))
            pos = transpose(array([r * cos(f), r * sin(f)])) + Rs
            vel = transpose(array([- aux * sin(f), aux * (e + cos(f))])) + Vs
            
            if len(pos) != len(unique(pos, axis=0)):
                print("Error. Duplicated planets positions found. Exiting.")
                exit()
            
            #Extra: Calculating m2 Hill sphere
            if isinstance(rad, (int, float)): rad = full(len(m), rad)
            elif isinstance(rad, (ndarray, list)) and (len(rad) != len(m)):
                print('Error in planetary radii. Exiting.')
                exit()
            else: rad = full(len(m), 0)
            
            if not all(rad): 
                mu  = m / (m + Ms)
                rad = a * (m / (3. * Ms))**(1 / 3.)
                rad[mu > 0.5] = 0.
                
                
        self.Npla  = len(m)
        self.massP = m
        self.posP  = pos
        self.velP  = vel
        self.Rpla  = rad
        
        
            
        
