from __future__ import print_function

from sys import exc_info, exit
from os import path, remove
from numpy import sum, array, zeros, uint32, int32, float32, unique
from logging import warning
from h5py import File

from libs.const import msol, AU

def save_particles(ids, pos, vel, mass, u, types, rad, outfile, format, units):

    nbin = sum(types == 5)
    ngas = len(mass) - nbin
        
    # conversion for different Units
    if units:
        print("[Output Units AU / Msun / km/s]")
        pos  /= AU
        mass /= msol
        vel  /= 1.e5
        u    /= 1.e10
    else:
        print("[Output Units CGS]")

    if path.isfile(outfile):
        print("WARNING: File {} already exist.".format(outfile))
        print("Do yo want to overwrite it? Y/[N]")
        q = input()
        if q.lower() in ['y', 'yes', 's', 'si']:
            remove(outfile)
        else:
            print('Exiting.')
            exit()

    if format == 0:
        # Openning file
        try:
            ofile = open(outfile,'w')
        except IOError as e:
            msg = "IO Error({0}): {1}".format(e.errno, e.strerror)
            warning(msg)
        except:
            print("Unexpected error: {}".format(exc_info()[0]))
            raise

        id_space = len("{}".format(ngas))

        # Preparing every line to print to the file
        for i in range(ngas+nbin):
            # Formatting particle attributes
            ie = '% d'    % ids[i]
            rx = '% 3.8e' % pos[i][0]
            ry = '% 3.8e' % pos[i][1]
            rz = '% 3.8e' % pos[i][2]
            vx = '% 3.8e' % vel[i][0]
            vy = '% 3.8e' % vel[i][1]
            vz = '% 3.8e' % vel[i][2]
            me = '% 3.8e' % mass[i]
            ue = '% 3.8e' % u[i]
            pt = '% d'    % types[i]

            # Right-align the strings
            outstring = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".\
                         format( ie.rjust(id_space + 1),\
                                 rx.rjust(12),\
                                 ry.rjust(12),\
                                 rz.rjust(12),\
                                 vx.rjust(12),\
                                 vy.rjust(12),\
                                 vz.rjust(12),\
                                 me.rjust(12),\
                                 ue.rjust(12),\
                                 pt)
            # Write to file
            ofile.write(outstring)

        # Closing the file
        ofile.close()

    elif format == 1:
        with File(outfile, "w") as f:
            f.create_group("Header")
            if len(unique(mass[:ngas])) > 1: mgas = 0.
            else: mgas = mass[0]
            if bool(nbin) & (len(unique(mass[ngas:])) == 1): mbh = mass[-1]
            else: mbh = 0.
            f["Header"].attrs["NumPart_ThisFile"]     = array([ngas,0,0,0,0,nbin], dtype=uint32)
            f["Header"].attrs["NumPart_Total"]        = array([ngas,0,0,0,0,nbin], dtype=uint32)
            f["Header"].attrs["MassTable"]            = array([mgas,0,0,0,0,mbh], dtype=float32)
            f["Header"].attrs["Time"]                 = 0.0
            f["Header"].attrs["Redshift"]             = 0.0
            f["Header"].attrs["Flag_Sfr"]             = int32(0)
            f["Header"].attrs["Flag_Feedback"]        = int32(0)
            f["Header"].attrs["Flag_DoublePrecision"] = int32(0)
            f.create_group("PartType0")
            f["PartType0"].create_dataset("Masses",         data=mass[:ngas].astype(float32))
            f["PartType0"].create_dataset("Coordinates",    data=pos[:ngas].astype(float32))
            f["PartType0"].create_dataset("Velocities",     data=vel[:ngas].astype(float32))
            f["PartType0"].create_dataset("ParticleIDs",    data=ids[:ngas].astype(int32))
            f["PartType0"].create_dataset("InternalEnergy", data=u[:ngas].astype(float32))
            if bool(nbin):
                f.create_group("PartType5")
                f["PartType5"].create_dataset("Masses",      data=mass[ngas:].astype(float32))
                f["PartType5"].create_dataset("Coordinates", data=pos[ngas:].astype(float32))
                f["PartType5"].create_dataset("Velocities",  data=vel[ngas:].astype(float32))
                f["PartType5"].create_dataset("ParticleIDs", data=ids[ngas:].astype(int32))
                if all(rad): 
                  if nbin == 1: f["PartType5"].create_dataset("SinkRadius", data=float32(rad))
                  else: f["PartType5"].create_dataset("SinkRadius", data=rad.astype(float32))

    else:
        print("Format {} unknown or not implemented. Exiting.".format(format))
        exit()
