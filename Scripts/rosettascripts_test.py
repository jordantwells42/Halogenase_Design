import sys

#Place your rosetta directory here
sys.path.append('/mnt/c/Users/jorda/Desktop/Software/PyRosetta')

from pyrosetta import *
from rosetta import *
from rosetta.protocols.rosetta_scripts import *
init()



pose = pose_from_pdb("RefPDBs/RebH.pdb")

script = """
<ROSETTASCRIPTS>
	<SCOREFXNS>
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<MOVE_MAP_FACTORIES>
	</MOVE_MAP_FACTORIES>
	<SIMPLE_METRICS>
	</SIMPLE_METRICS>
	<MOVERS>
	<FavorNativeResidue name="favor_native" bonus=
 "1.00"/>
	</MOVERS>
	<PROTOCOLS>
	<Add mover_name="favor_native"/>
	</PROTOCOLS>
</ROSETTASCRIPTS>
"""

xml = XmlObjects.create_from_string(script)
protocol = xml.get_mover("ParsedProtocol")
protocol.apply(pose)