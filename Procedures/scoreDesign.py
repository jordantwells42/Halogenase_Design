from pyrosetta import *
from rosetta.core.scoring.methods import ContextIndependentOneBodyEnergy

init()

class FavorNativeResidue(ContextIndependentOneBodyEnergy):

    def __init__(self, ref_pose):
        ContextIndependentOneBodyEnergy.__init__(self, self.creator())
        self.ref_pose = ref_pose
     
    def residue_energy(self, res, pose, emap):
        if res == self.ref_pose.residue(res.seqpos()):
            emap.get().set(self.scoreType, -1.0)
