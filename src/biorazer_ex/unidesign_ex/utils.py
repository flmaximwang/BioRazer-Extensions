from pathlib import Path


def fix_unidesign_pdb(pdb_file: str | Path):
    """
    PDBs from UniDesign are not in the correct format for biotite.
    This function replaces ATOM lines[60: 66] with '  0.00'
    """
    pdb_file = Path(pdb_file)
    with open(pdb_file, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            lines[i] = line[:60] + "  0.00" + line[66:]

    fixed_pdb = pdb_file.with_stem(pdb_file.stem + "_fixed")
    with open(fixed_pdb, "w") as file:
        file.writelines(lines)

    return fixed_pdb
