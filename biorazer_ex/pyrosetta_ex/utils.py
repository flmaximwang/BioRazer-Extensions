import site
from pathlib import Path
import pyrosetta as pr
from pyrosetta.rosetta.core.select import get_residue_set_from_subset
from pyrosetta.rosetta.core.pose import PDBInfo, Pose


def print_pyrosetta_dir():
    site_packages_dir = site.getsitepackages()
    pyrosetta_path = Path(site_packages_dir[0]) / "pyrosetta"
    return str(pyrosetta_path)


def generate_pose_from_residue_mask(pose, residue_mask):
    """
    Generate a new pose from the residue mask.
    All water molecules will be lost.
    PDBInfo of the new pose has not been correctlyly updated yet.
    """
    residue_set = get_residue_set_from_subset(residue_mask)
    new_pose = Pose()
    for residue_id in residue_set:
        new_pose.append_residue_by_jump(
            pose.residue(residue_id), len(new_pose.residues)
        )

    pdb_info = pose.pdb_info()
    new_pdb_info = PDBInfo(new_pose.size())

    for new_residue_index, residue_id in enumerate(residue_set, start=1):
        new_pdb_info.set_resinfo(
            new_residue_index, pdb_info.chain(residue_id), pdb_info.number(residue_id)
        )

    new_pose.pdb_info(new_pdb_info)

    return new_pose
