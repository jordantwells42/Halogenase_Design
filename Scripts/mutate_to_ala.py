"""
Mutates the active site of a given residue to all alanines
"""
"""
Adds the modules for PyRosetta into your current iPython environment
"""
import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

from pyrosetta import *
from rosetta import *
from pyrosetta.toolbox.mutants import *
init()

to_mutate = pose_from_pdb("RefPDBs/RebH.pdb")

active_site_res_pdb = [357, 465]
active_site_res_pose = []


for res in active_site_res_pdb:
    for i in range(res - 1, res + 2):
        active_site_res_pose.append(to_mutate.pdb_info().pdb2pose('A', i))
    
for res in active_site_res_pose:
    mutate_residue(to_mutate, res, "A")
    
to_mutate.dump_pdb("PDBs/ala6RebH.pdb")
