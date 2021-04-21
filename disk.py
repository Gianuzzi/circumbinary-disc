from __future__ import print_function

import numpy as np
from libs.gaseous_disk import Disk
from libs.stars import StarSystem
from libs.const import G, msol, AU, KB, mP
from libs.utils import save_particles
from libs.options_parser import OptionsParser


if __name__ == "__main__":
    
    op     = OptionsParser()
    args   = op.get_args()

    # where we want to place the center of mass
    r_cm = np.array([0.,0.])
    
    # first we create the binary system 
    m1   = args.m1 * msol
    m2   = args.m2 * msol
    abin = args.a * AU
    ebin = args.e
    
    star = StarSystem(m1=m1, m2=m2, a=abin, e=ebin, cm=r_cm, apo=True)
    
    Nbin = star.Nbin
    Mbin = star.Mbin
    Tbin = star.Tbin
        
    # now we create the disk    
    mdisk = args.mass * msol
    rmax  = args.rmax * AU
    rmin  = args.rmin * AU
    n     = args.num
    h     = args.h

    if n < 1: 
        print("n must be greater than 0. Exiting.")
        exit()
    # set disk positions
    disk = Disk(n=n, cm=r_cm, rmin=rmin, rmax=rmax, mass=mdisk, h=h)
    
    # set disk's mass distribution, velocity, and temperature (energy)
    ## distribution
    kappa = args.kappa
    sigma = args.sigma * msol / AU**2
    alpha = args.alpha
    Rgap  = args.Rgap * abin
    DR    = args.DR * Rgap
    disk.add_profile(kappa=kappa,
                     sigma=sigma, alpha=alpha, Rgap=Rgap, DR=DR,
                     method=3, nbins=250, logb=True, AU=AU)
    
    ## velocity
    model = args.Mrot
    disk.add_velocity(Mextra=Mbin, period=Tbin, model=model)
    
    ## energy
    T    = args.T * KB / mP
    beta = args.beta
    disk.add_energy(T=T, beta=beta, Mkep=Mbin, mu=2.001, gamma=1.0001, AU=AU)

    pos   = disk.pos
    vel   = disk.vel
    ngas  = disk.npart
    mass  = disk.mass
    u     = disk.u
    
    ids   = np.arange(1, ngas + 1)
    types = np.zeros(ngas).astype(int)
    
    print('Total disk mass: {:.2e} [Msol]'.format(sum(disk.mass)/msol))

    if Nbin != 0:
        pos   = np.append(pos, star.pos, axis=0)
        vel   = np.append(vel, star.vel, axis=0)
        mass  = np.append(mass, star.mass)
        ids   = np.append(ids, np.arange(Nbin) + ngas + 1)
        u     = np.append(u, np.zeros(Nbin))
        types = np.append(types, np.ones(Nbin) * 5).astype(int)
        print('Total star mass: {:.3f} [Msol]'.format(Mbin/msol))
    
    # Create z direction
    pos = np.append(pos, np.zeros(len(pos))[:,np.newaxis], axis=1)
    vel = np.append(vel, np.zeros(len(vel))[:,np.newaxis], axis=1)

    print("Writing output file {}...".format(args.outfile))
    save_particles(ids, pos, vel, mass, u, types, args.outfile, args.format, args.units)

    print("done...bye!")
