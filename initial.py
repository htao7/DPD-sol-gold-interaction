#!/usr/bin/python
# coding:utf8
"""initalize physical environment, model and etc.
"""
from __future__ import division
# from hoomd_script import *
import sys
import numpy
from globalvar import readin_parameter as para
from globalvar import boundary


def readinp(source='source.txt'):
    """readinp(source = 'source.txt') - Reading parameters from source,txt to the program.
    """
    filein = open(source, 'r')
    if filein is None:
        sys.stderr.write("%s not found. Use default parameters\n" % source)
    else:
        lines = filein.readlines()
        filein.close()
        for index, item in enumerate(lines):
            if item.find('Temperature k_BT') >= 0:
                para.k_BT = float(lines[index + 1].strip())
            elif item.find('Rc(Particle-Particle) Rc(Particle-Bead) Rc(Bead-Bead)') >= 0:
                for i in range(para.beadtypes):
                    para.rc[i] = [float(j) for j in lines[index + i + 1].strip().split(' ')]
            elif item.find('System size(Lx, Ly, Lz)') >= 0:
                para.size = [float(j) for j in lines[index + 1].strip().split(' ')]
            elif item.find('rho') >= 0:
                para.rho = float(lines[index + 1].strip())
            elif item.find('ABP concentration') >= 0:
                para.concentration = float(lines[index + 1].strip())
            elif item.find('Number of the chain') >= 0:
                para.n_Chain = int(lines[index + 1].strip())
            elif item.find('Number of the hydrophobic & hydrophilic beads of each chain') >= 0:
                para.n_B, para.n_A = [int(j) for j in lines[index + 1].strip().split(' ')]
            elif item.find('m(P)  m(B)  m(A)  m(S)') >= 0:
                para.mass = [float(j) for j in lines[index + 1].strip().split(' ')]
            elif item.find('r(P-B) r(B-B) r(B-A) r(A-A)') >= 0:
                para.bond_L = [float(j) for j in lines[index + 1].strip().split(' ')]
            elif item.find('time step (dt)') >= 0:
                para.dt = float(lines[index + 1].strip())
            elif item.find('parameters of aij') >= 0:
                for i in range(para.beadtypes):
                    para.aij[i] = [float(j) for j in lines[index + i + 1].strip().split(' ')]
            elif item.find('parameters of sigma') >= 0:
                para.sigma = float(lines[index + 1].strip())
            elif item.find('parameters of k_FENE & R0/Rc') >= 0:
                para.k_FENE, para.kBondBreak = [float(j) for j in lines[index + 1].strip().split(' ')]
            elif item.find('parameters of lambda') >= 0:
                para.k_lambda = float(lines[index + 1].strip())
            elif item.find('total time, snapshot') >= 0:
                para.t_run, para.dt = [float(j) for j in lines[index + 1].strip().split(' ')]
    para.n_Amphi = para.n_B + para.n_A
    para.n_DP = 1 + para.n_Chain * para.n_Amphi
    # equivalent bead number a polymer include, bead sizes considered
    n_eq_DP = int(round(para.mass[0] / para.mass[3] + para.mass[1] / para.mass[3] * para.n_Chain * para.n_B
                        + para.mass[2] / para.mass[3] * para.n_Chain * para.n_A))
    # equivalent total bead number of the system
    n_eq_Bead = int(round(para.rho * para.size[0] * para.size[1] * para.size[2]))
    # polymer number
    para.n_ABP = int(round(n_eq_Bead * para.concentration / n_eq_DP))
    # real polymer bead number
    para.n_ABPBead = para.n_ABP * para.n_DP
    # solvent number
    para.n_Solvent = n_eq_Bead - para.n_ABP * n_eq_DP
    para.n_Bead = para.n_ABPBead + para.n_Solvent
    para.gamma = 0.5 * para.sigma * para.sigma * para.k_BT
    para.size_half = [i / 2 for i in para.size]


def init_polymer():
    """initialize one polymer topo structure.
    """
    r_position = []
    r_type = []
    r_mass = []
    r_bonds = {
        'P_B': [],
        'B_B': [],
        'B_A': [],
        'A_A': []
    }
    # Add nanoparticle
    r_position.append([0.0, 0.0, 0.0])
    r_type.append(para.names[0] + "\n")
    r_mass.append(str(para.mass[0]) + "\n")
    l_P_B = para.bond_L[0]  # distance between P and first A
    for i in xrange(para.n_Chain):
        # Add B particles
        for j in xrange(para.n_B):
            l_to_P = l_P_B + j * para.bond_L[1]
            r_position.append([l_to_P * numpy.cos(2 * numpy.pi * i / para.n_Chain),
                               l_to_P * numpy.sin(2 * numpy.pi * i / para.n_Chain), 0.0])
            r_type.append(para.names[1] + "\n")
            r_mass.append(str(para.mass[1]) + "\n")
        l_P_A = l_P_B + (para.n_B - 1) * para.bond_L[1] + para.bond_L[2]
        # Add A particles
        for k in xrange(para.n_A):
            l_to_P = l_P_A + k * para.bond_L[3]
            r_position.append([l_to_P * numpy.cos(2 * numpy.pi * i / para.n_Chain),
                               l_to_P * numpy.sin(2 * numpy.pi * i / para.n_Chain), 0.0])
            r_type.append(para.names[2] + "\n")
            r_mass.append(str(para.mass[2]) + "\n")
    # Add bond information
    for i in xrange(1, para.n_DP, para.n_Amphi):
        r_bonds['P_B'].append([0, i])
        for j in xrange(1, para.n_B):
            r_bonds['B_B'].append([i + j - 1, i + j])
        r_bonds['B_A'].append([i + para.n_B - 1, i + para.n_B])
        for k in xrange(1, para.n_A):
            r_bonds['A_A'].append([i + para.n_B + k - 1, i + para.n_B + k])
    return r_position, r_type, r_mass, r_bonds


def initial_position():
    """genetaring init_read.xml for hoomd-blue.
    """
    init_xml = open("init_read.xml", "w")
    init_xml.write("""<?xml version="1.0" encoding="UTF-8"?>
<hoomd_xml version="1.5">
<!-- This is the input file of hoomd-blue created automatically by initial.py.
To change system environments, edit source.txt rather. -->
<configuration time_step="0" dimensions="3" vizsigma="1.5">
""")
    # write configuration start
    # Box size
    init_xml.write('<box lx="%f" ly="%f" lz="%f"/>\n' % tuple(para.size))
    r_position, r_type, r_mass, r_bonds = init_polymer()
    # Total lists of bead information
    a_position = []
    a_type = []
    a_mass = []
    a_bonds = []
    for i in xrange(para.n_ABP):
        # initializing nanoparticle position in three dimension
        p_nano = numpy.array([numpy.random.random() * d for d in para.size])
        # Nanoparticle index
        i_P = i * para.n_DP
        # Add bead type info
        a_type += r_type
        # Add bead mass info
        a_mass += r_mass
        # Add bead position info
        for index, item in enumerate(r_position):
            a_position.append("%10.5f %10.5f %10.5f\n" % tuple(boundary(p_nano + numpy.array(item))))
        # Add bead
        for key, values in r_bonds.items():
            for links in values:
                a_bonds.append("%s %d %d\n" % (key, links[0] + i_P, links[1] + i_P))
    a_type += ['S\n'] * para.n_Solvent
    a_mass += [str(para.mass[3]) + "\n"] * para.n_Solvent
    # generating solvent position
    for i in xrange(para.n_ABPBead, para.n_Bead):
        a_position.append("%10.5f %10.5f %10.5f\n" %
                          tuple([numpy.random.random() * para.size[j] - para.size_half[j] for j in range(3)]))
    init_xml.write("<position>\n")
    init_xml.writelines(a_position)
    init_xml.write("</position>\n")
    init_xml.write("<mass>\n")
    init_xml.writelines(a_mass)
    init_xml.write("</mass>\n")
    init_xml.write("<type>\n")
    init_xml.writelines(a_type)
    init_xml.write("</type>\n")
    init_xml.write("<bond>\n")
    init_xml.writelines(a_bonds)
    init_xml.write("</bond>\n")
    # write configuration end
    init_xml.write("""</configuration>
</hoomd_xml>
""")
    init_xml.close()
    pass


if __name__ == '__main__':
    readinp()
    for name, value in vars(para).items():
        if name.find('__') == -1: print('%s=%s' % (name, value))
    initial_position()