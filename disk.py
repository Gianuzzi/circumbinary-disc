from __future__ import print_function
from sys import exit

import numpy as np
from libs.gaseous_disk import Disk
from libs.stars import StarSystem
from libs.const import G, msol, AU, rsol 
from libs.utils import save_particles
from libs.options_parser import OptionsParser


if __name__ == "__main__":
    
    op     = OptionsParser()
    args   = op.get_args()

    # where we want to place the center of mass
    r_cm = np.array([0.,0.])
    

    # first we create the star system and the planets
    m1   = args.m1 * msol
    m2   = args.m2 * msol
    abin = args.a * AU
    ebin = args.e
    radb = 0 * rsol
    
    mp   = args.mp * msol
    ap   = args.ap * AU
    ep   = args.ep
    fp   = args.fp
    cs   = args.cs
    radp = 0 * rsol
   
    star = StarSystem(m1=m1, m2=m2, a=abin, e=ebin, cm=r_cm, rad=radb, apo=True)
    star.add_planets(m=mp, a=ap, e=ep, f=fp, cs=cs, rad=radp)
    
    Nbin = star.Nbin
    Mbin = star.Mbin
    Tbin = star.Tbin
    Npla = star.Npla
    
    # now we create the disk 
    thun  = args.thun   
    mdisk = args.mass * msol
    rmax  = args.rmax * AU
    rmin  = args.rmin * AU
    n     = args.num
    h     = args.h

    # set disk positions
    disk = Disk(n=n, cm=r_cm, rmin=rmin, rmax=rmax, mass=mdisk, h=h)
    
    # set disk's mass distribution, velocity, and temperature (energy)
    ## distribution
    thun  = args.thun
    kappa = args.kappa
    sigma = args.sigma * msol / AU**2
    alpha = args.alpha
    Rgap  = args.Rgap * abin
    DR    = args.DR * Rgap
    disk.add_profile(kappa=kappa, thun=thun,
                     sigma=sigma, alpha=alpha, Rgap=Rgap, DR=DR,
                     method=3, nbins=2000, logb=True, AU=AU)
    
    ## velocity
    model = args.Mrot
    disk.add_velocity(Mextra=Mbin, period=Tbin, model=model)
    
    ## energy
    T    = args.T 
    Tmin = 0.
    beta = args.beta
    disk.add_energy(T=T, beta=beta, Mkep=Mbin, mu=0, gamma=1.0001, Tmin=Tmin, AU=AU)

    pos   = disk.pos
    vel   = disk.vel
    ngas  = disk.npart
    mass  = disk.mass
    u     = disk.u 
    types = np.zeros(ngas).astype(int)
    ids   = np.arange(Nbin + Npla + ngas) + 1
    
    if Nbin + Npla + ngas == 0:
        print("Nothing created. Exiting.")
        exit()
    
    u     = np.append(u, np.zeros(Nbin + Npla))
    types = np.append(types, np.ones(Nbin + Npla) * 5).astype(int)
    print('Total disk mass: {:.2e} [Msol]'.format(sum(disk.mass) / msol))

    #dist  = np.diff(np.sort(np.linalg.norm(pos, axis=1)))
    #print('Minimum gas separation: {:.2e} [AU]'.format(np.min(dist) / AU))
    #print('Mean gas separation   : {:.2e} [AU]'.format(np.mean(dist) / AU))

    rad = np.array([])
    if Nbin > 0:
        pos   = np.append(pos, star.pos, axis=0)
        vel   = np.append(vel, star.vel, axis=0)
        mass  = np.append(mass, star.mass)
        rad   = star.Rbin
        print('Total star mass: {:.3f} [Msol]'.format(Mbin / msol))
    if Npla > 0:
        pos   = np.append(pos, star.posP, axis=0)
        vel   = np.append(vel, star.velP, axis=0)
        mass  = np.append(mass, star.massP)
        rad   = np.append(rad, star.Rpla)
        print('Total planetary mass: {:.3f} [Msol]'.format(sum(star.massP) / msol))

    
    # Create z direction
    pos = np.append(pos, np.zeros(len(pos))[:,np.newaxis], axis=1)
    vel = np.append(vel, np.zeros(len(vel))[:,np.newaxis], axis=1)

    print("Writing output file {}...".format(args.outfile))
    save_particles(ids, pos, vel, mass, u, types, rad, args.outfile, args.format, args.units)

    print("done...bye!")
