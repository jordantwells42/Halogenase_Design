import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

from pyrosetta import *
from rosetta import *
from rosetta.protocols.rosetta_scripts import *
init()


pose = pose_from_pdb("RebH")

parser = RosettaScriptsParser()
protocol = parser.generate_mover_and_apply_to_pose(pose, "inputs/min_L1.xml")
protocol.apply(pose)