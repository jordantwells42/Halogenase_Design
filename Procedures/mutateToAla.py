"""
Mutates the active site of a given residue to all alanines
"""


from pyrosetta import *
from rosetta import *
from pyrosetta.toolbox.mutants import *
init()

to_mutate = pose_from_pdb("PDBs/RebH.pdb")

active_site_res_pdb = [357, 461, 465]
active_site_res_pose = []


for res in active_site_res_pdb:
    for i in range(res - 1, res + 2):
        active_site_res_pose.append(to_mutate.pdb_info().pdb2pose('A', i))
    
for res in active_site_res_pose:
    mutate_residue(to_mutate, res, "A")
    
to_mutate.dump_pdb("PDBs/alaRebH.pdb")
