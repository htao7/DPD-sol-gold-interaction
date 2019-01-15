#!/usr/bin/python
# coding:utf8

"""Main script of the simulation
"""
from hoomd_script import *
import os
import sys
import initial
from globalvar import readin_parameter as para


if __name__ == '__main__':
    initial.readinp()
    if len(sys.argv)>1 and len(sys.argv[1]) > 0:
        init_file = sys.argv[1]
        if not init_file.endswith('.xml'):
            init_file += '.xml'
    else:
        initial.initial_position()
        init_file = "init_read.xml"
    # read in the file
    init.read_xml(init_file)

    # force field setup
    harmonic = bond.harmonic()
    harmonic.set_coeff('P_B', k=para.k_FENE, r0=para.bond_L[0])
    harmonic.set_coeff('B_B', k=para.k_FENE, r0=para.bond_L[1])
    harmonic.set_coeff('B_A', k=para.k_FENE, r0=para.bond_L[2])
    harmonic.set_coeff('A_A', k=para.k_FENE, r0=para.bond_L[3])

    dpd = pair.dpd(r_cut=para.rc[0][0], T=para.k_BT)
    dpd.pair_coeff.set('P', 'P', A=para.aij[0][0], r_cut=para.rc[0][0], gamma=para.gamma)
    dpd.pair_coeff.set('P', 'B', A=para.aij[0][1], r_cut=para.rc[0][1], gamma=para.gamma)
    dpd.pair_coeff.set('P', 'A', A=para.aij[0][2], r_cut=para.rc[0][2], gamma=para.gamma)
    dpd.pair_coeff.set('P', 'S', A=para.aij[0][3], r_cut=para.rc[0][3], gamma=para.gamma)
    dpd.pair_coeff.set('B', 'B', A=para.aij[1][1], r_cut=para.rc[1][1], gamma=para.gamma)
    dpd.pair_coeff.set('B', 'A', A=para.aij[1][2], r_cut=para.rc[1][2], gamma=para.gamma)
    dpd.pair_coeff.set('B', 'S', A=para.aij[1][3], r_cut=para.rc[1][3], gamma=para.gamma)
    dpd.pair_coeff.set('A', 'A', A=para.aij[2][2], r_cut=para.rc[2][2], gamma=para.gamma)
    dpd.pair_coeff.set('A', 'S', A=para.aij[2][3], r_cut=para.rc[2][3], gamma=para.gamma)
    dpd.pair_coeff.set('S', 'S', A=para.aij[3][3], r_cut=para.rc[3][3], gamma=para.gamma)
    dpd.set_params(T=para.k_BT)


    # dump every few steps
    os.system('mkdir data')
    dump.pdb(filename="data/beads", period=para.snapshot)
    # integrate NVE for a bunch of time steps
    all = group.all()
    integrate.mode_standard(dt=para.dt)
    integrate.nve(group=all, T=para.k_BT, tau=1.0)
    run(para.t_run)
    xml = dump.xml(filename="final.xml", vis=True)
    run(1)
    # convert file format
    os.system('perl filter.pl')