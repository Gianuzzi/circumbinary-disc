from __future__ import print_function

from numpy import sqrt

# all constants in cgs

# gravitational constant
G = 6.674e-8                # [ cm^3 g^-1 s^-1 ]

# avogadro constant
NA = 6.0221418e23           # [ ]

# boltzmann constant
KB = 1.3806504e-16          # [ erg K^-1 ]

# planck constant
H = 6.62606896e-27          # [ erg s ]

# speed of light in vacuum
c = 2.99792458e10           # [ cm s^-1 ]

# solar mass
msol = 1.989e33             # [ g ]

# solar radius
rsol = 6.955e10             # [ cm ]

# solar luminosity
lsol = 3.839e33             # [ erg s^-1 ]

# electron charge
qe = 4.80320427e-10         # [ esu ]

# hydrogen mass
mH = 1.6735575e-24          # [ g ]

# proton mass
mP = 1.6726219e-24          # [ g ]

# atomic mass unit
amu = 1.6605390401e-24      # [ g ]

# ev2erg
ev2erg = 1.602177e-12       # [ erg eV^-1 ]

# parsec
parsec = 3.08568025e18      # [ cm ]

# AU
AU = 1.49597871e13          # [ cm ]

# Hydrogen dissociation temperature
Dis_T = 2.e4                # [ K ]

# Hydrogen ionization temperature
Ion_T = 13.6 * ev2erg / KB  # [ K ]

# conversion factor for cosmological magnetic field
bfac = sqrt(1e10 * msol) / sqrt(1e6 * parsec) * 1e5 / (1e6 * parsec)

# golden ratio for image heights
golden_ratio = (sqrt(5)-1)/2

