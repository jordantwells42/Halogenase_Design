"""
Adds the modules for PyRosetta into your current iPython environment
"""
import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

import random
import math
import multiprocessing as mp

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


init("-use_time_as_seed -seed_offset 10")

OUTPUT = mp.Queue()

def design(pose_in, ref, scorefxn, active_site_res_pose = [], ligand_res_pose = [], output = "design", 
    job = 1, outer_cycles = 4, inner_cycles = 60,
    temp_init = 10, temp_final = 0.6, rep_weight_init = 0.15, rep_weight_final = 0.55, 
    trans_init = 1.5, trans_final = 0.1, rot_init = 20, rot_final = 2, backrub_moves = 10):
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
    init()
    pose = Pose()
    pose.assign(pose_in)

    init_seq = pose.sequence()
    
    
    
    """
    ===========================
    2. Define Movers
    ===========================

        a. Creating necessary resfile and movemaps
        b. Initializing perturbMover: either translates or rotates the ligand
        c. Initializing minMovers: changes the conformation to the nearest local minimum
        d. Initializing packMover: determines which side chain rotamers lead to the lowest energy
        e. Initializing backrubMover: performs local backbone motion on active site res to 
            allow for backbone flexibility
        f. Initializing pyMOLMover: allows us to see the progress in PyMOL
        g. (Optional) RosettaScripts protocol
    
    ===========================
    """
    
    # a. Creating necessary resfile and movemaps

    # Using "ALLAAxc" for the active site residues allows them to be mutated to all AAs besides cysteine
    res_to_pack = {res : "ALLAAxc" for res in active_site_res_pose}

    # Using "NATAA" for the ligand residues allows them to switch to different rotamers
    res_to_pack.update({res :  "NATAA" for res in ligand_res_pose})
     
    generate_resfile_from_pose(pose, "design.resfile", False, specific = res_to_pack)
    
    res_to_work = []
    for res in active_site_res_pose:
        res_to_work.append(res)
        
    for res in ligand_res_pose:
        res_to_work.append(res)

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

    perturb_mover = rigid_moves.RigidBodyPerturbMover(1, rot_init, trans_init)

    # c. Initializing minMovers: changes the conformation to the nearest local minimum
    
    minmover_rb_bb = MinMover("dfpmin_armijo_nonmonotone")
    minmover_rb_bb.movemap(movemap_rb_bb)
    minmover_rb_bb.score_function(scorefxn)

    minmover_bb_chi = MinMover("dfpmin_armijo_nonmonotone")
    minmover_bb_chi.movemap(movemap_bb_chi)
    minmover_bb_chi.score_function(scorefxn)

    # d. Initializing pack_mover: determines which side chain rotamers lead to the lowest energy

    pack_design = TaskFactory.create_packer_task(pose)
    parse_resfile(pose, pack_design, "design.resfile")
    pack_mover = PackRotamersMover(scorefxn, pack_design)

    # e. Initializing backrub_mover: performs local backbone motion on active site res to 
    #      allow for backbone flexibility

    backrub_ini_mover = BackrubMover()
    backrub_ini_mover.set_movemap(movemap_bb)
    
    mc_br = MonteCarlo(pose, scorefxn, temp_init)
    backrub_trial_mover = TrialMover(backrub_ini_mover, mc_br)
    backrub_mover = RepeatMover(backrub_trial_mover, backrub_moves)

    # f. Initializing pymol_mover: allows us to see the progress in PyMOL

    pymol_mover = pyrosetta.PyMOLMover()

    """
    ===========================
    3. Job Cycle
    ===========================

        a. Creating incrementation factors for rep_weight, kT, translation, and rotation
        b. Setting up Jobs
        c. Outer Cycle
        d. Inner Cycle
        e. Mover Execution
        f. Outputting Job information
    
    ===========================
    """

    # a. Creating incrementation factors for rep_weight, kT, translation, and rotation

    gamma_kt = math.pow((temp_final / temp_init), (1.0 /(outer_cycles * inner_cycles)))
    slope_rep = (rep_weight_final - rep_weight_init) / (inner_cycles)
    gamma_trans = math.pow((trans_final / trans_init), (1.0 / (outer_cycles * inner_cycles)))
    gamma_rot = math.pow((rot_final / rot_init), (1.0 / (outer_cycles * inner_cycles)))

    # b. Setting up Job

    file_name = str(f"{output}_{str(job)}")
    
    p = Pose()
    p.assign(pose_in)
    p.pdb_info().name(file_name)
    
    score_fa = get_fa_scorefxn()

    kt = temp_init
    rep_weight = rep_weight_init
    trans_mag = trans_init
    rot_mag = rot_init
    
    mc = MonteCarlo(p, scorefxn, kt)

    # c. Outer Cycle
    low_pose = Pose()
    low_pose.assign(pose_in)
    for _ in range(outer_cycles):
        rep_weight = rep_weight_init        

        # d. Inner Cycle

        for _ in range(1, inner_cycles + 1):
            rep_weight += slope_rep
            scorefxn.set_weight(fa_rep, rep_weight) 
            minmover_rb_bb.score_function(scorefxn)
            minmover_bb_chi.score_function(scorefxn)
            mc.score_function(scorefxn)
            mc_br.score_function(scorefxn)

            kt = kt * gamma_kt
            mc.set_temperature(kt)
            mc_br.set_temperature(kt)

            trans_mag *= gamma_trans
            rot_mag *= gamma_rot
            perturb_mover.trans_magnitude(trans_mag)
            perturb_mover.rot_magnitude(rot_mag)
            
            # e. Mover Execution

            perturb_mover.apply(p)
            minmover_rb_bb.apply(p)

            pack_design = TaskFactory.create_packer_task(p)
            parse_resfile(p, pack_design, "design.resfile")
            pack_mover = PackRotamersMover(scorefxn, pack_design)
            pack_mover.apply(p)

            backrub_mover.apply(p)

            minmover_bb_chi.apply(p)
            
            

            mc.boltzmann(p)
            pymol_mover.apply(p)

            low_score = math.inf

            score = score_fa(p)

            if score < low_score:
                low_score = score
                low_pose.assign(p)

            print(score)

        p.assign(low_pose)

    # f. Outputting Job information

    p.assign(low_pose)

    pymol_mover.apply(p)

    p.dump_pdb(f"PDBs/{file_name}.pdb")

    final_seq = p.sequence()
    
    final_mutation = ""
    num_mutations = 0

    for res in active_site_res_pose:
        if init_seq[res - 1] != final_seq[res - 1]:
            final_mutation += str(init_seq[res - 1]) + str(res) + str(final_seq[res - 1]) + " "
            num_mutations += 1
    
    
    scorefxn.set_weight(fa_rep, 0.55)
    final_energy = scorefxn(p)
    p.delete_residue_slow(529)
    final_rmsd = all_atom_rmsd(p, ref)
    
    csv_output = f"{file_name}, {str(final_energy)}, {str(final_rmsd)}, \
    {str(num_mutations)}, {final_mutation}, {final_seq}\n"
    OUTPUT.put((job, csv_output))

    


to_design = pose_from_pdb("RefPDBs/phe-RebH_2.pdb")
ref_enz = pose_from_pdb("RefPDBs/RebH.pdb")

scoreFA = get_fa_scorefxn()

# setup_foldtree(to_design, "A_X", Vector1([1]))
ft = FoldTree()
ft.add_edge(357, 1, -1)
ft.add_edge(357, 528, -1)
ft.add_edge(357, 529, 1)
ft.add_edge(529, 529, -1)
to_design.fold_tree(ft)

active_site_pdb = [357, 461, 465]
active_site_pose = set()

for res in active_site_pdb:
    for i in range(res - 1, res + 2):
        active_site_pose.add(to_design.pdb_info().pdb2pose('A', i))
       
ligand_pdb = [1]
ligand_pose = []
       
for res in ligand_pdb:
    ligand_pose.append(to_design.pdb_info().pdb2pose('X', res))

output = "design_mp_1"


#Running the Code
processes = [mp.Process(target = design, args=(to_design, ref_enz, scoreFA, active_site_pose, ligand_pose, output, x)) for x in range(1, 9)]

for process in processes:
    process.start()

for process in processes:
    process.join()

results = [OUTPUT.get() for p in processes]


# csv: contains all of the log information from each job
csv = open("CSVs/" + output + ".csv", "w+")
csv.write("File Name, Final Energy, Final RMSD, Number of Mutations, Mutations, Final Sequence\n")

for result in results:
    csv.write(result)

csv.close()