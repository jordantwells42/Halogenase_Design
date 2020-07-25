import random
import math

from pyrosetta import *
from rosetta import *

from pyrosetta.teaching import *
from pyrosetta.toolbox import *
from rosetta.core.pack.task import *
from rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.backrub import *
from pyrosetta.toolbox import generate_resfile_from_pose
from pyrosetta.rosetta.protocols.protein_interface_design import *
from pyrosetta.rosetta.protocols.protein_interface_design.movers import *
import rosetta.protocols.rigid as rigid_moves

init()

"""
===========================
def alaDesign:
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

def alaDesign(pose_in, ref_pose, scorefxn, active_site_res_pose = [], ligand_res_pose = [], log_output = "designLog", 
job_output = "output", jobs = 1, outer_cycles = 4, inner_cycles = 50, linear_temp = False, linear_rep = True, linear_perturb = False,
temp_init = 100, temp_final = 5, rep_weight_init = 0.25, rep_weight_final = 0.55, 
trans_init = 1.2, trans_final = 0.2, rot_init = 20, rot_final = 2, backrub_moves = 15):
    
    #Creating our working pose
    pose = Pose()
    pose.assign(pose_in)
    init_seq = pose.sequence()
    ref_seq = ref_pose.sequence()
    
    
     
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
    specific_design = {}
 
    for res in active_site_res_pose:
        specific_design[res] = "ALLAAxc"
    
    for res in ligand_res_pose:
        specific_design[res] = "NATAA"
    
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
    
    movemap_lig = MoveMap()
    for res in ligand_res_pose:
        movemap_lig.set_chi(res, True)
    
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
    
    minMover_lig = MinMover("dfpmin_armijo_nonmonotone")
    minMover_lig.movemap(movemap_lig)
    minMover_lig.score_function(scorefxn)
    
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
    pyMOLMover.apply(pose)
    
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
    gammakT = math.pow((temp_final / temp_init), (1.0 /(outer_cycles * inner_cycles)))
    gammaRep = math.pow((rep_weight_final / rep_weight_init), (1.0 / (outer_cycles)))
    
    if linear_temp:
        slopekT = (temp_final - temp_init) / (inner_cycles * outer_cycles)
    else:
        gammakT = math.pow((temp_final / temp_init), (1.0 /(outer_cycles * inner_cycles)))
        
    if linear_rep:
        slopeRep = (rep_weight_final - rep_weight_init) / (outer_cycles)
    else:    
        gammaRep = math.pow((rep_weight_final / rep_weight_init), (1.0 / (outer_cycles)))
    
    if linear_perturb:
        slopeTrans = (trans_final - trans_init) / (outer_cycles * inner_cycles)
        slopeRot = (rot_final - rot_init) / (outer_cycles * inner_cycles)
    else:
        gammaTrans = math.pow((trans_final / trans_init), (1.0 / (outer_cycles * inner_cycles)))
        gammaRot = math.pow((rot_final / rot_init), (1.0 / (outer_cycles * inner_cycles)))
    
    for job in range(1, jobs + 1):
    
        #2. Housekeeping
        log = open(log_output + "_" + str(job) + ".txt", "w+")
        
        p.assign(pose_in)
        p.pdb_info().name(job_output + '_' + str(job))
        
        
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
                protein_or_ligand_move = random.randint(0,1)
                
                if protein_or_ligand_move == 0:
                    perturbMover.apply(p)
                    minMover_rb_bb.apply(p)
                    #minMover_lig.apply(p)

                else:
                    backrubMover.apply(p)
                    
                    pack_design = TaskFactory.create_packer_task(p)
                    parse_resfile(p, pack_design, "design.resfile")
                    packMover = PackRotamersMover(scorefxn, pack_design)
                    packMover.apply(p)
                    minMover_bb_chi.apply(p)

                mc.boltzmann(p)
                
                
                
                pyMOLMover.apply(p)
                print(scorefxn(p))
            
            mc.recover_low(p)
            
            minMover_rb_bb.apply(p)
            
            pyMOLMover.apply(p)
            
            post_seq = p.sequence()
            
            #Outputs the mutations caused from each outer cycle
            mutation = ""
            for res in active_site_res_pose:
                if pre_seq[res - 1] != post_seq[res - 1]:
                    mutation += str(pre_seq[res - 1]) + str(res) + str(post_seq[res - 1]) + " "
            
            log.write("Outer Cycle " + str(i) + ": \n")
            log.write("Score: " + str(scorefxn(p)) + "\n")
            log.write("Mutations : " + mutation + " \n\n")
        
        #5. Outputting Job Information
        mc.recover_low(p)
        pyMOLMover.apply(p)
        
        p.dump_pdb(job_output + "_" + str(job) + ".pdb")
        
        designed_seq = p.sequence()
        
        final_mutation = ""
        
        for res in active_site_res_pose:
            if init_seq[res - 1] != designed_seq[res - 1]:
                final_mutation += str(init_seq[res - 1]) + str(res) + str(designed_seq[res - 1]) + " "
            
        log.write("\n\nFinal Score : " + str(scorefxn(p)) + "\n")
        log.write("Final Mutations: " + final_mutation + "\n")
        
        #Specific to alanine sequence recovering
        mismatches = 0
        comparison = ""
        for res in active_site_res_pose:
            if ref_seq[res - 1] != designed_seq[res - 1]:
                mismatches += 1
            comparison += ref_seq[res - 1] + ":" + designed_seq[res - 1] + " "
            
        log.write("Mismatches : " + str(mismatches) + "\n")
        log.write("Comparison (Reference:Designed): " + comparison + "\n")
       
        log.close()
        
#Running the code!

#Make sure that your protein is chain A and your protein ligand is chain X
cleanATOM("alaComplex1.pdb")
pose1 = pose_from_pdb("alaComplex1.clean.pdb")
ref = pose_from_pdb("complex1.clean.pdb")

#Establishing scorefunction
scoreFA = get_fa_scorefxn()

#Creating fold tree, A will be rigid and X will be movable
#setup_foldtree(pose, "A_X", Vector1([1]))
ft = FoldTree()
ft.add_edge(357, 1, -1)
ft.add_edge(357, 528, -1)
ft.add_edge(357, 529, 1)
ft.add_edge(529, 529, -1)
pose1.fold_tree(ft)

active_site_pdb = [357, 461, 465]
active_site_pose = []

for res in active_site_pdb:
    for i in range(res - 1, res + 2):
        active_site_pose.append(pose1.pdb_info().pdb2pose('A', i))

ligand_pdb = [1]
ligand_pose = []
       
for res in ligand_pdb:
    ligand_pose.append(pose1.pdb_info().pdb2pose('X', res))
    
alaDesign(pose1, ref, scoreFA, active_site_pose, ligand_pose, "testLog", "testOutput", 8, 4, 50, linear_perturb = True)
#alaDesign(pose1, ref, scoreFA, active_site_pose, ligand_pose, "test2Log", "default2Output", 4, 4, 50, linear_perturb = True)
            
    
    
    
    
