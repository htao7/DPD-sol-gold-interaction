#!/usr/bin/python
# coding:utf8

# Global parameter definition file
"""global parameter definition.
Including physical parameter and script running flags.
"""


class readin_parameter:
    """readin parameters including temperature,
bond length, angle and etc. from source.txt
    """

    k_BT = 1.0                          # Temperature k_BT
    gamma = None                        # gamma parameter of DPD dissipative force
    beadtypes = 4                       # constant: bead types
    names = ['P', 'B', 'A', 'S']        # type naming
    # pairwise interaction definition
    # BE CAREFUL: In python, array index starts with 0
    rc = [                              # cut-off radius
            [2.0, 1.5, 1.5, 1.5],
            [1.5, 1.0, 1.0, 1.0],
            [1.5, 1.0, 1.0, 1.0],
            [1.5, 1.0, 1.0, 1.0]
        ]
    aij = [                             # DPD interaction parameters
             [25.0, 25.0, 25.0, 150.0],
             [25.0, 25.0, 50.0, 75.0],
             [25.0, 50.0, 25.0, 25.0],
             [150.0, 75.0, 25.0, 25.0],
    ]
    #beads properties
    mass = [8.0, 1.0, 1.0, 1.0]         # mass of beads
    bond_L = [1.29, 0.86, 0.86, 0.86]   # balance bond length: r(P-B) r(B-B) r(B-A) r(A-A)
    k_FENE = 30.0                       # parameters of k_FENE
    kBondBreak = 1.5                    # R0/Rc - bonds breaking criterion
    size = [60.0, 60.0, 60.0]           # system size
    size_half = [30.0, 30.0, 30.0]      # half box size
    rho = 3.0                           # particle density
    concentration = 0.10                # solute bead concentration
    n_Chain = 3                         # number of chain
    n_B = 3                             # number of B type bead in a chain
    n_A = 1                             # number of A type bead in a chain
    n_Amphi = None                      # bead number of a chain: nAmpi = n_B + n_A + ...
    n_DP = None                         # bead number of a polymer: nDP = 1 + nChain*nDP
    n_ABP = None                        # polymer numbers
    n_ABPBead = None                    # total polymer bead numbers
    n_Solvent = None                    # solvent bead number
    n_Bead = None                       # total bead numbers including solvent
    dt = 0.03                           # time step (dt)
    t_run = 1000000                     # total time
    snapshot = 5000                     # making snapshot interval
    sigma = 3.0e0
    k_lambda = 0.65e0

    def __init__(self):
        pass
    pass


class run_variable:
    """Program running flag including debug flag, time passing and etc.
    """

    isDebug = 1
    timestart = None
    timerun = None

    def __init__(self):
        pass
    pass


def boundary(arr):
    """Usage: boundary condition correcting.
    """
    result = []
    for index, item in enumerate(arr):
        if item >= readin_parameter.size_half[index]:
            item -= readin_parameter.size[index]
        elif item < -readin_parameter.size_half[index]:
            item += readin_parameter.size[index]
        result.append(item)
    return result

if __name__ == '__main__':
    print __doc__