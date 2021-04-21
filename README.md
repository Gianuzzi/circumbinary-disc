# circumbinary-disk
A disk-like distribution of particles, with a binary system in it's center.

Te code first creates (can be ommited) a binary system (PartType 5) a the center.
Then it creates a gaseous disk around it, than can have a gap in its center. After
this, a surface denstity profile can be set, such as a radial power law, or following
a relation described in eq(12) of Thun et al. (2017) [DOI:10.1051/0004-6361/201730666].
Then, the code assigns velocities and temperature to the disk particles.

This code is **under development**.

# Author
Emmanuel Gianuzzi

# Usage

The basic usage is
```bash
python disk.py -n NUM
```
where NUM is the total number of particles desired to represent the disk.
This will produce a hdf5 file called 'ics_disk.dat' that works as initial
conditions for GADGET or GIZMO (format 3).

For a complete description of parameters use
```bash
python disk.py -h
```
