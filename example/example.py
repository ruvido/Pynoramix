#!/usr/bin/python
import sys
sys.path.append('..')

from pyn_cl_set import *
from pyn_cl_unit import *
from pyn_cl_net import *
wbox=cl_set('tip4p_ew.pdb')
wbox.water_analysis(traj_name='md_test.xtc')
net=wbox.get_network()
net_symm=net.symmetrize()
net_symm.gradient_clusters()
