import rosetta, pyrosetta
from rosetta import *
from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.core.scoring.methods import ContextDependentOneBodyEnergy

@pyrosetta.EnergyMethod()
class FavorReferenceResidue(ContextDependentOneBodyEnergy):

    def __init__(self, reference_pose = Pose()):
        ContextDependentOneBodyEnergy.__init__(self, self.creator())
        self.reference_pose = reference_pose
     
    def residue_energy(self, res, pose, emap):
        if len(self.reference_pose.sequence()) == 0:
            emap.set(self.scoreType, 0)
            return
        elif res == self.reference_pose.residue(res.seqpos()):
            emap.set(self.scoreType, -1.0)
        else:
            emap.set(self.scoreType, 0)
