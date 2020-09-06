"""
Adds the modules for PyRosetta into your current iPython environment
"""
import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

import random
import math
import numpy as np

from pyrosetta import *
from rosetta import *
from pyrosetta.toolbox.generate_resfile import generate_resfile_from_pose
from pyrosetta.toolbox.cleaning import cleanATOM
from rosetta.core.pack.task import parse_resfile,TaskFactory
from pyrosetta.rosetta.protocols.backrub import BackrubMover
from pyrosetta.rosetta.protocols.protein_interface_design import FavorNativeResidue
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.analysis.simple_metrics import SequenceRecoveryMetric
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.docking import setup_foldtree
import rosetta.protocols.rigid as rigid_moves
from pyrosetta.rosetta.core.scoring import *
from rosetta.protocols.rosetta_scripts import *

init()

def design(pose_in, scorefxn, active_site_res_pose = [], ligand_res_pose = [], output = "design", 
jobs = 1, outer_cycles = 4, inner_cycles = 50,
temp_init = 10, temp_final = 0.6, rep_weight_init = 0.05, rep_weight_final = 0.55, 
trans_init = 1.5, trans_final = 0.1, rot_init = 20, rot_final = 2):
	"""
    ===========================
    Design Function:
    Performs high resolution refinement and design of an enzyme active site with an amino acid/protein substrate
    ===========================
        
        pose_in: pose to be designed
        scorefxn: scorefunction to use for design (can have it customly favor native residues, favor certain structural changes etc
        active_site_res_pose: a list containing the pose residue numbers of each residues to be designed/packed
        ligand_res_pose: a list containing the residues of the ligand
        log_output: denotes the name of the log files
        job_output: denotes the names of the decoys
        jobs: how many different decoys to design
        outer_cycles: how many cycles of repulsion ramping
        inner_cycles: how many cycles of refinement per cycle of repulsion ramping
        temp_init: initial temperature for simulated annealing
        temp_final: final temperature for simulated annealing
        rep_weight_init: initial weight for VDW repulsion in the scorefunction
        rep_weight_final: final weight for VDW repulsion in the scorefunction
        trans_init: initial mean of the translation of the ligand
        trans_final: final mean of the translation of the ligand
        rot_init: initial mean of the rotation of the ligand
        rot_init: final mean of the rotation of the ligand
    
    ===========================
		
		Steps:
		1. Create working poses and all necessary log files
		2. Define Movers
		3. Job Cycle
		4. Write all final information
	
	===========================

	===========================
	1. Create working poses and all necessary log files
	===========================
	"""
	pose = Pose()
    pose.assign(pose_in)

    init_seq = pose.sequence()
    
    # csv: contains all of the log information from each job
    csv = open("CSVs/" + log_output + ".csv", "w+")
	
	"""
	===========================
	2. Define Movers
	===========================

		a. Creating necessary movemaps and resfile
    	b. Initializing perturbMover: either translates or rotates the ligand
    	c. Initializing minMovers: changes the conformation to the nearest local minimum
    	d. Initializing packMover: determines which side chain rotamers lead to the lowest energy
    	e. Initializing backrubMover: performs local backbone motion on active site res to 
    		allow for backbone flexibility
    	f. Initializing pyMOLMover: allows us to see the progress in PyMOL
	
	===========================
    """
    
    # a. Creating necessary movemaps and resfile

    # Resfile
    # Using "ALLAAxc" for the active site residues allows them to be mutated to all AAs besides cysteine
	res_to_pack = {res : "ALLAAxc" for res in active_site_res_pose}

	# Using "NATAA" for the ligand residues allows them to switch to different rotamers
    res_to_pack.update({res :  "NATAA" for res in ligand_res_pose})
     
    generate_resfile_from_pose(pose, "design.resfile", False, specific = res_to_pack)
    
    # Creating list of all residues to be worked on
    res_to_work = []
    for res in active_site_res_pose:
        res_to_work.append(res)
        
    for res in ligand_res_pose:
        res_to_work.append(res)
    
    
    #Movemaps
    movemap_bb = MoveMap()
    movemap_chi = MoveMap()
    movemap_rb_bb = MoveMap()
    movemap_bb_chi = MoveMap()
               
    for res in res_to_work:
        movemap_bb.set_bb(res, True)
        movemap_chi.set_chi(res, True)
        
        movemap_bb_chi.set_bb(res, True)
        movemap_bb_chi.set_chi(res,True)
        movemap_rb_bb.set_bb(res, True)
       
    movemap_rb_bb.set_jump(1, True)

    # b. Initializing perturbMover: either translates or rotates the ligand
    # c. Initializing minMovers: changes the conformation to the nearest local minimum
    # d. Initializing packMover: determines which side chain rotamers lead to the lowest energy
    # e. Initializing backrubMover: performs local backbone motion on active site res to 
    # 	   allow for backbone flexibility
    # f. Initializing pyMOLMover: allows us to see the progress in PyMOL




