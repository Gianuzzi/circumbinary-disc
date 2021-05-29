from __future__ import print_function
from sys import exit

from numpy import sum, max, min, all
from numpy import log, log10, exp, sqrt
from numpy import cos, sin, pi
from numpy import array, arange, diff
from numpy import unique, argsort, zeros
from numpy import full, linspace, logspace
from numpy import mgrid, random, newaxis
from numpy import transpose, append, digitize
from numpy import concatenate, matmul
from numpy.linalg import norm

from functools import partial

from libs.funcs import trapezoidal, simpson, boole, surf_prof
from libs.const import G, KB, mP

class Disk:
    """ Class for creating a distribution of particles in a close-packed
        sphere.

        Arguments:
            n    : total number of desired points to represent the disk.
            cm   : coordinates of the disk's center (with units).
            rmin : disk's outer radius (with units).
            rmax : disk's inner radius (with units).
            mass : disk's mass   (with units).
            h    : disk's aspect ratio.
    """
    def __init__(self, n=10000, cm=[0,0], rmin=0, rmax=1, mass=1, h=0.05):

        if n == 0:
            print("Gaseous disk not created.")
            self.npart = 0

        else:
            print("Creating gaseous disk...")
            self.npart = n
            self.center = array(cm)
            self.rmin   = rmin
            self.rmax   = rmax
            self.mdisk  = mass
            self.h      = h
            
        
    
    def add_profile(self, method=2, thun: bool=True,
                    kappa=0, sigma=1e-4, alpha=1.5, Rgap=0.55, DR=0.055,
                    nbins=250, logb=True, seed=42, AU=1):
        """ Function for setting a radial density profile to a uniform
            disk of particles. Alpha, Rgap and DeltaR are the values to
            introduce for a density profile according to
            Thun et al. (2017) [DOI: 10.1051/0004-6361/201730666]. 

            Arguments:
                method : method to be implemented. It is only availabe
                         for radial power law distribution.
                         1: The profile is created by redistributing the
                            particles positions, considering they have equal
                            mass.
                         2: The profile is created by redistributing the
                            particles masses, without modifying their positions.
                         3: Create particles potitions and masses, following
                            the desired distribution.
                thun   : if the density profile model is the one in Thun et al.
                kappa  : power index of radial density distribution.
                sigma  : reference surface density (with units).
                alpha  : value of alpha (initial slope).
                Rgap   : value of Rgap (estimated size of the gap, with units).
                DR     : value of DeltaR (transition width, with units).
                nbins  : number of bins if method==3 and Rgap!=0.
                logb   : if the bins are separated log-spaced.
                AU     : astronomical unit in input distance units.
                seed   : seed for random numbers generator.
        """
        
        if self.npart == 0:
            self.pos   = transpose(array([[],[]]))
            self.mass  = array([])
            return

        unif = ((not thun) & (kappa == 0)) | ((self.mdisk == 0) & (sigma == 0))
        if (method in [1, 2]) | unif: # we need preset positions
            # first we create a square with uniform distribution, hence we need to sample
            # more particles than the desired N
            side    = 2 * self.rmax
            ratio   = 4 / pi * self.rmax**2 / (self.rmax**2 - self.rmin**2) # ratio between square and disk area
            nsquare = ratio * self.npart
            nside   = int(sqrt(nsquare * 0.25))
            naux    = nside**2
            z3      = mgrid[0:nside, 0:nside].T.reshape(naux, 2)
            grid    = (concatenate((z3,
                                    z3 + [0.0, 0.5], 
                                    z3 + [0.5, 0.0], 
                                    z3 + [0.5, 0.5])) + 0.25)/nside * side
            pos     = grid - self.rmax
            r       = norm(pos, axis=1)
            ig      = (r >= self.rmin) & (r <= self.rmax)
            pos     = pos[ig] + self.center
            npart   = len(pos)
            masses  = full(npart, self.mdisk / float(npart)) # uniform masses

            print(" We placed {:d} gas cells in a close-packed disk.".format(npart))
            self.npart = npart
            
            if unif:
                print(" Disk mass distributed radially uniformingly.")
                self.pos   = pos
                self.mass  = masses
                return

        if not thun:

            if kappa <= -2:   # impossible
                print("kappa must be greater than -2. Exiting.")
                exit()
                
            print(" Setting radial density profile with RHO~r**{}...".format(kappa))

            if method == 1:
                print(" Modifying particles positions...")
                
                # we use centered transpose of pos
                pos  = transpose(pos - self.center)
                
                if self.rmin > 0: 
                    print(" WARNING: This density profile method may create an extra")
                    print("           gap (or modify the existing one) at the center")
                    print("           of the disk.")
                    print("          Please make sure this is what you want.")            

                # redistribute particles uniformingly along radius
                pos *= norm(pos, axis=0)

                # redistribute particles according to kappa
                pos *= norm(pos, axis=0)**(1 / float(2 + kappa) - 1)

                # normalize radii according to radius
                pos *= self.rmax / max(norm(pos, axis=0))
                

                self.pos   = transpose(pos) + self.center
                self.mass  = masses


            elif method == 2:
                
                print(" Modifying particles masses...")
                
                # we use centered transpose of pos
                pos  = transpose(pos - self.center)

                # get radii and sort particles
                radii   = norm(pos, axis=0)
                order   = argsort(radii)
                pos     = pos[:,order]
                radii   = radii[order]

                # normalization cte
                cte     = (kappa + 2) * log(radii[-1] - radii[0]) - log(self.mdisk)

                # radial bins
                eps     = 1e-6
                nbins   = int((radii[-1] - radii[0]) / max(diff(unique(radii))))

                if (self.rmin > 0) & logb:
                    bins = logspace(log10(radii[0] * (1 + eps)), 
                                    log10(radii[-1] * (1 + eps)),
                                    nbins)
                else:
                    bins = linspace(radii[0] * (1 + eps), 
                                    radii[-1] * (1 + eps),
                                    nbins)

                # cumulative mass until each bin line
                CM_b    = bins ** (kappa + 2) * exp(-cte)

                # mass until each bin line
                M_b     = append(CM_b[0], diff(CM_b))

                # particles bin's index
                Pb_ind  = digitize(radii, bins)

                # particles before right bin indexes
                uni, NP_b = unique(Pb_ind, return_counts=True)
                
                nuni = len(uni)
                if not all(arange(nuni) ==  uni):
                    print('  Fixing empty bins...')
                    
                    np_b      = zeros(nbins)
                    np_b[uni] = NP_b
                    
                    naux = zeros(nuni)
                    maux = zeros(nuni)
                    m    = M_b[0]
                    l    = 0
                    for i in range(nbins - 1): # -1
                        if np_b[i] == 0:
                            m += M_b[i+1]
                        else:
                            maux[l] = m
                            naux[l] = np_b[i]
                            l += 1
                            m = M_b[i+1]
                    
                    maux[-1] += M_b[-1]  #-1
                    naux[-1] += NP_b[-1] #-1
                    
                    # particles's mass per bin
                    PM_b   = maux / naux
                    for i in range(nuni): Pb_ind[Pb_ind == uni[i]] = i
                
                else:
                    # particles's mass per bin
                    PM_b    = M_b / NP_b

                # distribute mass
                masses  = PM_b[Pb_ind]


                self.pos  = transpose(pos) + self.center
                self.mass = masses

                
            elif method == 3:
                # STILL TESTING!
                
                print(" Creating particles positions and masses, using {:d} bins...".format(nbins))
                                
                eps   = 1e-7
                if (self.rmin > 0) & logb:
                    bins = logspace(log10(self.rmin * (1 + eps)), 
                                    log10(self.rmax),
                                    nbins)
                else:
                    bins = linspace(self.rmin * (1 + eps), 
                                    self.rmax,
                                    nbins)

                # cumulative mass until each bin line
                CM_b  = bins**(kappa + 2) 
                
                # normalize total mass
                CM_b *= (self.mdisk / CM_b[-1])
                
                # mass until each bin line
                M_b   = append(CM_b[0], diff(CM_b))
                
                # particles mass aproximation
                pmass = sum(M_b) / self.npart
                
                # desired particles until each bin line
                NP_b   = (M_b / pmass).astype(int)
                npart  = sum(NP_b)
                masses = full(npart, pmass)
                
                # now we have to place the particles
                # here we set them randomly, but they
                # could be pre-determined
                random.seed(seed)
                
                # re-define bins edges
                bins  = append(self.rmin, bins)
                
                ## bin centers
                bin_c = bins[:-1] + diff(bins) * 0.5
                                                
                pos   = zeros((npart, 2))
                n     = 0
                
                for i in arange(nbins):
#                     rs = full(NP_b[i], bin_c[i])
                    rs = sqrt(random.uniform(bins[i]**2, 
                                             bins[i+1]**2,
                                             NP_b[i]))
                    th = random.random(NP_b[i]) * 2 * pi
                    
                    pos[n:n+NP_b[i]] = array([rs * cos(th),
                                              rs * sin(th)]).T
                    n += NP_b[i]
                    
                
                self.npart = npart
                self.pos   = pos + self.center
                self.mass  = masses
                
            
            else: 
                print("Method must be 1, 2 or 3. {:d} was given. Exiting.".format(method))
                exit()

        # now the "publication" profile   
        else:
            
            if DR == 0:
                print("Error. dR must be non zero. Exiting")
                exit()

            print(" Setting surface density from Thun et.al 2017,")
            print("  with parameters:") 
            print(" \tAlpha: {}\n \tRgap: {}\n \tdR: {}".format(alpha, Rgap/AU, DR/AU))
            
            # we defined the function we are going to integrate,
            # to get the mass ditribution from r0 to r1.
            # it is surf_prof, in funcs.
            # the three possible integration methods are:
            #   trapezoidal, simpson (even n), boole (n=5*p)
            # default is simpson.
            # the integration result, from r0 to r1, will be the 
            # individual mass in each bin
            
            integ = partial(surf_prof, alpha=alpha, Rgap=Rgap, DR=DR)
 
            # now we offer 2 of the possible methods described before
            
            if method in [1,2]:

                if method == 1:
                    print("WARNING. Method 1 for this distribution isn't implemented.")
                    print("         Switchin to method 2.")
                
                print(" Modifying particles masses...")
                
                # we use centered transpose of pos
                pos   = transpose(pos - self.center)
                
                
                # get radii and sort particles
                radii = norm(pos, axis=0)
                order = argsort(radii)
                pos   = pos[:,order]
                radii = radii[order]

                # radial bins (for integration)
                eps   = 1e-8
                nbins = int((max(radii) - min(radii)) / 
                            max(diff(unique(radii))))
                
                if (self.rmin > 0) & logb:
                    bins = logspace(log10(radii[0] * (1 - eps)), 
                                    log10(radii[-1] * (1 + eps)),
                                    nbins)
                else:
                    bins = linspace(radii[0] * (1 - eps), 
                                    radii[-1] * (1 + eps),
                                    nbins)                
                
                # mass per bin
                M_b = []
                for i in range(nbins - 1): M_b.append(simpson(integ, bins[i], bins[i+1]))

                # normalize, according to Sigma and alpha factor unit
                M_b       = array(M_b) * 2 * pi * sigma * AU**alpha

                # particles bin's index
                Pb_ind    = digitize(radii, bins) - 1
                
                # particles per bin
                uni, NP_b = unique(Pb_ind, return_counts=True)
                
                nuni = len(uni)
                if not all(arange(nuni) ==  uni):
                    print('  Fixing empty bins...')
                    
                    np_b      = zeros(nbins - 1)
                    np_b[uni] = NP_b
                    
                    naux = zeros(nuni)
                    maux = zeros(nuni)
                    m    = M_b[0]
                    l    = 0
                    for i in range(nbins - 2): # -2
                        if np_b[i] == 0:
                            m += M_b[i+1]
                        else:
                            maux[l] = m
                            naux[l] = np_b[i]
                            l += 1
                            m = M_b[i+1]
                            
                    maux[-1] += M_b[-1]  #-1
                    naux[-1] += NP_b[-1] #-1
                    
                    # particles's mass per bin
                    PM_b   = maux / naux
                    for i in range(nuni): Pb_ind[Pb_ind == uni[i]] = i
                
                else:
                    # particles's mass per bin
                    PM_b    = M_b / NP_b
            
                # distribute mass
                masses  = PM_b[Pb_ind]
                
                
                self.pos  = transpose(pos) + self.center
                self.mass = masses
    
    
            elif method == 3:
                
                # this is a bit more complicated, because we can't just
                # shift the positions. we will create new positions.
                # first we calculated mass distribution.
                print(" Creating particles positions and masses, using {:d} bins...".format(nbins))
                
                eps  = 1e-8
                if (self.rmin > 0) & logb:
                    bins = logspace(log10(self.rmin * (1 - eps)), 
                                    log10(self.rmax * (1 + eps)),
                                    nbins)
                else:
                    bins = linspace(self.rmin * (1 - eps), 
                                    self.rmax * (1 + eps),
                                    nbins)
                # mass per bin
                M_b = []
                for i in range(nbins - 1): M_b.append(simpson(integ, bins[i], bins[i+1]))
                
                # normalize, according to Sigma(Rgap)
                M_b    = array(M_b) * 2 * pi * sigma * AU**alpha
                
                # particles mass aproximation
                pmass  = sum(M_b) / self.npart
                
                # desired particles per bin
                NP_b   = (M_b / pmass).astype(int)
                npart  = sum(NP_b)
                masses = full(npart, pmass)
                
                # now we have to place the particles
                # here we set them randomly, but they
                # could be pre-determined
                random.seed(seed)
                
                ## bin centers
                delta = diff(bins)
                bin_c = bins[:-1] + delta * 0.5
                                                
                pos   = zeros((npart, 2))
                n     = 0
                
                for i in arange(nbins - 1):
#                     rs = full(NP_b[i], bin_c[i])
                    rs = sqrt(random.uniform(bins[i]**2, 
                                             bins[i+1]**2,
                                             NP_b[i]))
                    th = random.random(NP_b[i]) * 2 * pi
                    
                    pos[n:n+NP_b[i]] = array([rs * cos(th),
                                              rs * sin(th)]).T
                    n += NP_b[i]
                    
                
                self.npart = npart
                self.pos   = pos + self.center
                self.mass  = masses
                
                                
            else: 
                print("Method must be 1, 2 or 3. {:d} was given. Exiting.".format(method))
                exit()
                
                
                
    def add_velocity(self, Mextra=0, period=0, model=1):
        """ Function for setting the velocities of a disk of particles,
            rotating contraclockwise. Pure circular motion assumed.

            Arguments:
                Mextra : aditional mass (at center of mass) used for calculating
                         velocities.
                period : disk rotation period (with units).
                model  : model to be implemented
                         0: no velocities.
                         1: keplerian velocities from Mextra.
                         2: keplerian velocities from Mextra and disk mass.
                         3: velocities from given period.
        """
    
        if self.npart == 0:
            self.vel = transpose(array([[],[]]))
            return
        
        print(" Adding velocities...")
        
        if model==0: vel = zeros((self.npart, 2))
        
        elif model in [1,2]:
            print("  Setting keplerian velocities...")
            pos        = self.pos - self.center
            radii      = norm(pos, axis=1)
            self.v_kep = sqrt(Mextra * G / radii)
            if model==2: Mextra += sum(self.mass)
            v_kep = sqrt(Mextra * G / radii)
            vel   = matmul(pos / radii[:, newaxis], array([[0, 1], [-1, 0]])) * v_kep[:, newaxis]
            

        elif model==3:
            print("  Setting velocities from binary period...")
            if period==0:
                print("  Incorrect period for setting disk velocities.")
                print("  Disk velocities are set to zero.")
                vel = zeros((self.npart, 2))
                
            else:
                pos   = self.pos - self.center
                v_ang = 1 / float(period)                
                vel   = v_ang * matmul(pos, array([[0, 1], [-1, 0]]))
                
        else:
            print("Model must be 0, 1, 2 or 3.")
            print(" {:d} was given. Exiting.".format(model))
            exit()
            
        
        self.vel = vel
                
    
    
    def add_energy(self, T=0, beta=0, Mkep=0, mu=1, gamma=5/3., Tmin=0., AU=1):
        """ Function for setting the temperature of a disk of particles. This
             function is still being tested and corrected.
            # InternalEnergy = Temperature * kBoltz / (mu * mP * (gamma - 1))

            Arguments:
                T     : temperature at radius = 1 AU. If set to any negative 
                         value, the disk temperature will be set according
                         to the local speed of sound (Mkep must be the value
                         of the central mass generating it).
                beta  : power index of a radial temperature distribution.
                Mkep  : central mass used for calculating the local speed
                         of sound.
                mu    : mu used (mean molecular weight / mP). If set equal
                         to 0, it is calculated dependig the temperature.
                gamma : adiabatic index used. If set equal
                         to 0, it is calculated dependig the temperature.
                Tmin  : Minimum temperature (with temperature units).
                AU    : astronomical unit in input distance units.
        """
        
        if self.npart == 0:
            self.u = array([])
            return

        print(' Adding temperature...')

        radii  = norm(self.pos - self.center, axis=1)
   
        if T < 0:
            print('  Setting Temperature according to c_s(r)**2, with c_s(r) = h * v_kep(r),')
            if not hasattr(self, 'v_kep'): self.v_kep = sqrt(Mkep * G / radii)
            c_s   = self.h * self.v_kep
            if gamma <= 0: gamma_cs = 1
            else: gamma_cs = gamma
            if mu <= 0:
                if gamma == 1.4: mu_cs = 2.01  # Adiabatic molecular
                else: mu_cs = 1.005            # Adiabatic/Isothermal atomic

            print('   using gamma = {:.3f}, and mu = {:.3f}...'.format(gamma_cs, mu_cs))  
            
            T = c_s**2 * mu_cs * mP / (gamma_cs * KB)

        else:
            T = full(self.npart, T * (radii / AU)**beta)


        T[T < Tmin] = Tmin  
        
        # This following settings are still being tested if correct
    
        if mu <= 0:
            print('  Calculating mu used in code...')
            mu = full(self.npart, 0.5) # Ionized
            mu[T < 1.e4] = 1.005       # Atomic
            mu[T < 2.e3] = 2.01        # Molecular (Optional)

        if gamma <= 0:
            print('  Calculating gamma used in code...')
            gamma = full(self.npart, 5/3.) # Adiabatic monoatomic
            gamma[T < 2.e3] = 1.4	   # Adiabatic diatomic

        factor = KB / (mP * mu * (gamma - 1))
        
        if max(T) == min(T): print("  Temperature: {} [K]".format(T[0]))
        else:
            print('  Maximum temperature: {:.2f} [K]'.format(max(T)))
            print('  Minimum temperature: {:.2f} [K]'.format(min(T)))
         

        self.u = T * factor