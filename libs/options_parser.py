from argparse import ArgumentParser, RawTextHelpFormatter

class OptionsParser:
    def __init__(self):
        self.parser = ArgumentParser(description=
                "Generate a disk of particles, with angular velocities, with a\n"+\
                "binary system in it's center. The disk can also have a central gap\n"
                "A radial surface density profile can be setted, following a radial\n"
                "power law, or adjusted to the model shown in equation (12) of\n"+\
                "Thun et al. (2017) [DOI:10.1051/0004-6361/201730666]. If the last\n"+\
                "one is desired, Rgap must be non zero.\n"+\
                "If m1 or m2 equals 0, no binary system is created.",
                formatter_class=RawTextHelpFormatter)

        self.parser.add_argument("-N", "-n",
                                 dest     = "num",
                                 type     = int,
                                 help     = "Number of particles (required).",
                                 required = True)

        self.parser.add_argument("-o",
                                 metavar = "outfile",
                                 dest    = "outfile",
                                 help    = "Name of output file.\n"+\
                                           " [Default = ics_disk.dat]",
                                 default = "ics_disk.dat")

        self.parser.add_argument("-format",
                                 dest    = "format",
                                 type    = int,
                                 help    = "Format of output file:\n"+\
                                           "0 = ASCII,\n"+\
                                           "1 = HDF5 format.\n"+\
                                           " [Default = 1]",
                                 default = 1)
        
        self.parser.add_argument("-rmin", 
                                 dest     = "rmin",
                                 type     = float,
                                 help     = "Inner radius of the disk\n"+\
                                            "(in astronomical units).\n"+\
                                            " [Default = 0]",
                                 default  = 0.)
        
        self.parser.add_argument("-rmax",
                                 dest     = "rmax",
                                 type     = float,
                                 help     = "Outer radius of the disk\n"+\
                                            "(in astronomical units).\n"+\
                                            " [Default = 1]",
                                 default  = 1.)
        
        self.parser.add_argument("-m",
                                 dest     = "mass",
                                 type     = float,
                                 help     = "Total mass of the gaseous disk\n"+\
                                            "(in solar masses).\n"+\
                                            " [Default = 1]",
                                 default  = 1.)
        
        self.parser.add_argument("-H", "-hz",
                                 dest     = "h",
                                 type     = float,
                                 help     = "Razón de aspecto de disco.\n"+\
                                            " [Default = 0.05]",
                                 default  = 0.05)
        
        self.parser.add_argument("-m1",
                                 dest     = "m1",
                                 type     = float,
                                 help     = "Mass of the first binary\n"+\
                                            "component (in solar masses).\n"+\
                                            " [Default = 1]",
                                 default  = 1.)

        self.parser.add_argument("-m2",
                                 dest     = "m2",
                                 type     = float,
                                 help     = "Mass of the second binary\n"+\
                                            "component (in solar masses).\n"+\
                                            " [Default = 0]",
                                 default  = 0.)

        self.parser.add_argument("-abin", "-a",
                                 dest     = "a",
                                 type     = float,
                                 help     = "Semi-major axis of the binary\n"+\
                                            "system (in astronomical units).\n"+\
                                            " [Default = 0]",
                                 default  = 0.)
        
        self.parser.add_argument("-ebin", "-e",
                                 dest     = "e",
                                 type     = float,
                                 help     = "Eccentricity of the binary system.\n"+\
                                            " [Default = 0]",
                                 default  = 0.)

        self.parser.add_argument("-k", "-kappa",
                                 dest     = "kappa",
                                 type     = float,
                                 help     = "Power index of radial surface\n"+\
                                            "density distribution profile.\n"+\
                                            " [Default = 0 (Uniform)]",
                                 default  = 0.)
        
        self.parser.add_argument("-Rgap", "-Rg",
                                 dest     = "Rgap",
                                 type     = float,
                                 help     = "Estimated size of the gap\n"+\
                                            "(in a_bin units).\n"+\
                                            " [Default = 2.5]",
                                 default  = 2.5)
        
        self.parser.add_argument("-s", "-sigma",
                                 dest     = "sigma",
                                 type     = float,
                                 help     = "Reference surface density\n"+\
                                            "(in (solar masses)*(astronomical units)^{-2}).\n"+\
                                            " [Default = 0.0001]",
                                 default  = 1.e-4)

        self.parser.add_argument("-alpha",
                                 dest     = "alpha",
                                 type     = float,
                                 help     = "Initial density profile slope\n"+\
                                            " [Default = 1.5]",
                                 default  = 1.5)

        self.parser.add_argument("-DR", "-dR",
                                 dest     = "DR",
                                 type     = float,
                                 help     = "Transition width (in Rgap units).\n"+
                                            " [Default = 0.1]",
                                 default  = 0.1)
        
        self.parser.add_argument("-Mrot",
                                 dest     = "Mrot",
                                 type     = int,
                                 help     = "Model used to set initial\n"+\
                                            "angular velocities:\n"+\
                                            "0 = No velocities,\n"+\
                                            "1 = Keplerian velocities,\n"+\
                                            "    from binary mass,\n"+\
                                            "2 = Keplerian velocities,\n"+\
                                            "    from disk and binary mass,\n"+\
                                            "3 = Velocities matching\n"+\
                                            "     binary mean motion.\n"
                                            " [Default = 1]",
                                 default  = 1)
        
        self.parser.add_argument("-T", "-Temp",
                                 dest     = "T",
                                 type     = float,
                                 help     = "Temperature at radius (with units) = 1.\n"+
                                            " [Default = 0]",
                                 default  = 0.)
        
        self.parser.add_argument("-b", "-beta",
                                 dest     = "beta",
                                 type     = float,
                                 help     = "Power index of radial temperature\n"+\
                                            "distribution profile.\n"+\
                                            " [Default = 0 (Uniform)]",
                                 default  = 0)
        
        self.parser.add_argument("--units",
                            dest     = "units",
                            help     = "Change units to Msol/AU/km s^{-1} ",
                            action   = "store_true")

    def get_args(self):
        return self.parser.parse_args()

