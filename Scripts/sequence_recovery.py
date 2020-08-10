"""
Adds the modules for PyRosetta into your current iPython environment
"""
import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

import random
import math

from pyrosetta import *
from rosetta import *
from pyrosetta.toolbox.generate_resfile import generate_resfile_from_pose
from pyrosetta.toolbox.cleaning import cleanATOM
from rosetta.core.pack.task import parse_resfile
from rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.backrub import BackrubMover
from pyrosetta.rosetta.protocols.protein_interface_design import FavorNativeResidue
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.analysis.simple_metrics import SequenceRecoveryMetric
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.docking import setup_foldtree
import rosetta.protocols.rigid as rigid_moves
from pyrosetta.rosetta.core.scoring import *
init()

def recover_sequence(pose_in, ref_pose, scorefxn, active_site_res_pose = [], ligand_res_pose = [], log_output = "designLog", 
jobs = 1, outer_cycles = 4, inner_cycles = 50, linear_temp = False, linear_rep = True, linear_perturb = False,
temp_init = 100, temp_final = 10, rep_weight_init = 0.05, rep_weight_final = 0.55, 
trans_init = 1.5, trans_final = 0.1, rot_init = 20, rot_final = 2, backrub_moves = 10):
    """
    ===========================
    def recover_sequence:
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
        linear_temp: decides whether to use linear temperature/kT decrements from initial to final or a geometric progression
        linear_rep: decides whether to use linear repulsion weight increments from initial to final or a geometric progression 
        linear_perturb: decides whether to use linear perturbation decrements from initial to final or a geometric progression 
        temp_init: initial temperature for simulated annealing
        temp_final: final temperature for simulated annealing
        rep_weight_init: initial weight for VDW repulsion in the scorefunction
        rep_weight_final: final weight for VDW repulsion in the scorefunction
        trans_init: initial mean of the translation of the ligand
        trans_final: final mean of the translation of the ligand
        rot_init: initial mean of the rotation of the ligand
        rot_init: final mean of the rotation of the ligand
        backrub_moves: number of backrub moves to perform per refinement
    ===========================
    """
    #Creating our working pose
    pose = Pose()
    pose.assign(pose_in)
    init_seq = pose.sequence()
    ref_seq = ref_pose.sequence()
    
    log_all = open("Logs/" + log_output + "_all.txt", "w+")
    csv = open("CSVs/" + log_output + ".csv", "w+")
    
    """
    ===========================
    Defining Movers
    ===========================
    1. Creating necessary movemaps and resfile
        Using "ALLAAxc" for the active site residues allows them to be mutated to all AAs besides cysteine
        Using "NATAA" for the ligand residues allows them to switch to different rotamers
    2. Initializing perturbMover: either translates or rotates the ligand
    3. Initializing minMovers: changes the conformation to the nearest local minimum
        minMover_rb_bb: goes to the nearest local minimum by changing rigid body and backbone DoFs
        minMover_bb_chi: goes to the nearest local minimum by changing backbone and chi angle DoFs
    4. Initializing packMover: determines which side chain rotamers lead to the lowest energy
    5. Initializing backrubMover: performs local backbone motion on active site res to allow for backbone flexibility
    6. Initializing pyMOLMover: allows us to see the progress in PyMOL
    ===========================
    """
    
    
    #1. Creating necessary resfile and movemaps
    #Resfile
    specific_design = {res : "ALLAAxc" for res in active_site_res_pose}
    specific_design.update({res :  "NATAA" for res in ligand_res_pose})
    
    print(specific_design)
    
    generate_resfile_from_pose(pose, "design.resfile", False, specific = specific_design)
    
    #Creating list of all residues to be worked on
    to_work_res = []
    for res in active_site_res_pose:
        to_work_res.append(res)
        
    for res in ligand_res_pose:
        to_work_res.append(res)
    
    
    #Movemaps
    movemap_bb = MoveMap()
    movemap_chi = MoveMap()
    movemap_rb_bb = MoveMap()
    movemap_bb_chi = MoveMap()
               
    for res in to_work_res:
        movemap_bb.set_bb(res, True)
        movemap_chi.set_chi(res, True)
        
        movemap_bb_chi.set_bb(res, True)
        movemap_bb_chi.set_chi(res,True)
        movemap_rb_bb.set_bb(res, True)
       
    movemap_rb_bb.set_jump(1, True)
    
    
    #2. Initializing perturbMover
    perturbMover = rigid_moves.RigidBodyPerturbMover(1, rot_init, trans_init)
    
    #3. Initializing minMovers
    #   minMover_rb_bb minimizes rigid body and backbone DoFs
    minMover_rb_bb = MinMover("dfpmin_armijo_nonmonotone")    
    minMover_rb_bb.movemap(movemap_rb_bb)
    minMover_rb_bb.score_function(scorefxn)
   
    #   minMover_bb_chi minimizes backbone and chi angle DoFs
    minMover_bb_chi = MinMover("dfpmin_armijo_nonmonotone")
    minMover_bb_chi.movemap(movemap_bb_chi)
    minMover_bb_chi.score_function(scorefxn)
    
    
    #4. Initializing packMover
    pack_design = TaskFactory.create_packer_task(pose)
    parse_resfile(pose, pack_design, "design.resfile")
    packMover = PackRotamersMover(scorefxn, pack_design)
    
    #5 Initializing backrubMover
    backrubIniMover = BackrubMover()
    backrubIniMover.set_movemap(movemap_bb)
    
    mcBR = MonteCarlo(pose, scorefxn, temp_init)
    backrubTrialMover = TrialMover(backrubIniMover, mcBR)
    backrubMover = RepeatMover(backrubTrialMover, backrub_moves)
    
    #6 Initializing pyMOLMover
    pyMOLMover = pyrosetta.PyMOLMover()
    
    
    """
    Execution: Runs refinement jobs number of times and outputs the poses and log files
    1. Creating incrementation factors for rep_weight, kT, translation, and rotation
    2. Housekeeping (creating log, reestablishing pose, and resetting important values
    3. Outer Cycle
        Keeps tracks of mutations after each cycle, control the repulsion ramping
    4. Inner Cycle
        Meat of the design, applys all the movers and controls simulated annealing
    5. Outputting Job information
    """
 
    p = Pose()
    
    #1. Creating incrementation factors for rep_weight, kT, translation, and rotation
    #   Creating temperature incrementation
    if linear_temp:
        slopekT = (temp_final - temp_init) / (inner_cycles * outer_cycles)
    else:
        gammakT = math.pow((temp_final / temp_init), (1.0 /(outer_cycles * inner_cycles)))
    
    #   Creating repulsion incrementation
    if linear_rep:
        slopeRep = (rep_weight_final - rep_weight_init) / (outer_cycles)
    else:    
        gammaRep = math.pow((rep_weight_final / rep_weight_init), (1.0 / (outer_cycles)))
    
    #   Creating perturbation incrementation
    if linear_perturb:
        slopeTrans = (trans_final - trans_init) / (outer_cycles * inner_cycles)
        slopeRot = (rot_final - rot_init) / (outer_cycles * inner_cycles)
    else:
        gammaTrans = math.pow((trans_final / trans_init), (1.0 / (outer_cycles * inner_cycles)))
        gammaRot = math.pow((rot_final / rot_init), (1.0 / (outer_cycles * inner_cycles)))
    
    
    
    for job in range(1, jobs + 1):
        #2. Housekeeping
        file_name = str(f"{log_output}_{str(job)}")
        
        log = open(f"Logs/{file_name}.txt", "w+")
        
        p.assign(pose_in)
        p.pdb_info().name(file_name)
        
        
        kT = temp_init
        rep_weight = rep_weight_init
        trans_mag = trans_init
        rot_mag = rot_init
        
        mc = MonteCarlo(p, scorefxn, kT)
        
        
        #3. Outer Cycle
        for i in range(1, outer_cycles + 1):
            
            pre_seq = p.sequence()
            
            if linear_rep:
                rep_weight += slopeRep
            else:    
                rep_weight = rep_weight * gammaRep
                
            scorefxn.set_weight(fa_rep, rep_weight)
            
            mc.score_function(scorefxn)
            mcBR.score_function(scorefxn)
            
            #4. Inner Cycle
            for j in range(1, inner_cycles + 1):
                FavorNativeResidue(p, 100000)
                if linear_temp:
                    kT += kT + slopekT
                else:
                    kT = kT * gammakT
                mc.set_temperature(kT)
                mcBR.set_temperature(kT)
                
                if linear_perturb:
                    trans_mag += slopeTrans
                    rot_mag += slopeRot
                    
                else:
                    trans_mag = trans_mag * gammaTrans
                    rot_mag = rot_mag * gammaRot
                  
                perturbMover.trans_magnitude(trans_mag)
                perturbMover.rot_magnitude(rot_mag)
                
                #Mover Execution
                perturbMover.apply(p)
                minMover_rb_bb.apply(p)
                
                pack_design = TaskFactory.create_packer_task(p)
                parse_resfile(p, pack_design, "design.resfile")
                packMover = PackRotamersMover(scorefxn, pack_design)
                packMover.apply(p)
                
                #backrubMover.apply(p)
                minMover_bb_chi.apply(p)
                         
                mc.boltzmann(p)
                
                pyMOLMover.apply(p)
                print(scorefxn(p))
            
            mc.recover_low(p)
            
            pyMOLMover.apply(p)
            post_seq = p.sequence()
            
            #Outputs the mutations caused from each outer cycle
            mutation = ""
            for res in active_site_res_pose:
                if pre_seq[res - 1] != init_seq[res - 1]:
                    mutation += str(pre_seq[res - 1]) + str(res) + str(post_seq[res - 1]) + " "
            
            log.write(f"Outer Cycle {str(i)}:\n")
            log.write(f"Score: {str(scorefxn(p))}\n")
            log.write(f"Mutations : {mutation}\n\n")
        
        #5. Outputting Job Information
        mc.recover_low(p)
        pyMOLMover.apply(p)
        
        p.dump_pdb(f"PDBs/{log_output}_{str(job)}.pdb")
        
        final_seq = p.sequence()
        
        final_mutation = ""
        
        for res in active_site_res_pose:
            if init_seq[res - 1] != final_seq[res - 1]:
                final_mutation += str(init_seq[res - 1]) + str(res) + str(final_seq[res - 1]) + " "
        
        
        final_energy = scorefxn(p)
        final_rmsd = all_atom_rmsd(p, ref_pose)
        
        #Specific to sequence recovering
        away_from_ref = 0
        comparison = ""
        for res in active_site_res_pose:
            if ref_seq[res - 1] != final_seq[res - 1]:
                away_from_ref += 1
            if init_seq[res - 1] != final_seq[res - 1]:
                comparison += f"{ref_seq[res - 1]}:{final_seq[res - 1]} "

        """
        seq_recoverer = SequenceRecoveryMetric()
        seq_recoverer.set_comparison_pose(ref_pose)
        
        selector = ResidueIndexSelector
        selector.append_index(active_site_res_pose)
        for res in active_site_res_pose:
            print(type(res))
            print(type(int(res)))
            selector.append_index(res)
        
        seq_recoverer.set_residue_selector(selector)
        seq_recovery = seq_recoverer.calculate(p)
        """
        seq_recovery = (len(active_site_res_pose) - away_from_ref)/(len(active_site_res_pose))

        #Outputting Final Information
        
        outputLog = ""
        outputLog += f"Final Score : {str(final_energy)}\n"
        outputLog += f"Final Mutations: {final_mutation}\n"          
        outputLog += f"Sequence Identity : {str(seq_recovery)}\n"
        outputLog += f"Comparison (Reference:Designed): {comparison}\n"
        
        log.write("\n\n" + outputLog)
        log_all.write("Job " + str(job) + ": \n" + outputLog + "\n\n")
        log.close()
        
        outputCSV = f"{file_name}, {str(final_energy)}, {str(final_rmsd)}, {str(seq_recovery)}, {str(mutations)}, {final_mutation}, {final_seq}\n"
        csv.write(outputCSV)
        
    log_all.close()
    csv.close()
        
#Running the code!

#Make sure that your protein is chain A and your protein ligand is chain X
cleanATOM("RefPDBs/alaComplex1.pdb")
pose1 = pose_from_pdb("RefPDBs/alaComplex1.clean.pdb")
ref = pose_from_pdb("RefPDBs/complex1.clean.pdb")

#Establishing scorefunction
scoreFA = get_fa_scorefxn()

"""
from scoreDesign import FavorReferenceResidue
fnr = FavorReferenceResidue(pose1).scoreType
scoreFA.set_weight(fnr, 100000000000000000000.0)
"""

#Creating fold tree, A will be rigid and X will be movable
#setup_foldtree(pose1, "A_X", Vector1([1]))

ft = FoldTree()
ft.add_edge(357, 1, -1)
ft.add_edge(357, 528, -1)
ft.add_edge(357, 529, 1)
ft.add_edge(529, 529, -1)
pose1.fold_tree(ft)



active_site_pdb = [357, 461, 465]
active_site_pose = []

#Generating residues to be worked on
for res in active_site_pdb:
    for i in range(res - 1, res + 2):
        active_site_pose.append(pose1.pdb_info().pdb2pose('A', i))
       
ligand_pdb = [1]
ligand_pose = []
       
for res in ligand_pdb:
    ligand_pose.append(pose1.pdb_info().pdb2pose('X', res))
 
#Running the Code
recover_sequence(pose1, ref, scoreFA, active_site_pose, ligand_pose, "scoredesigntest", 8, 3, 5)
            
    
    
    
    