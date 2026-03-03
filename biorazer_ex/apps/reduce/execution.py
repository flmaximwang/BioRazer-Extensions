import tempfile
from pathlib import Path
from biotite.structure import AtomArray
from biorazer.structure.io.protein import PDB2PDB
from ..basic import App


class Reduce(App):

    def __init__(self, app_dir=None, app_bin="reduce", log_name="reduce"):
        super().__init__(app_dir, app_bin, log_name=log_name)


class ReduceFile(Reduce):
    """
    Run reduce in terms of input and output files.
    """

    def run(self, *args, cwd=".", output_file=None, **kwargs):
        result = super().run(*args, cwd=cwd, get_output=True, **kwargs)
        if output_file:
            with open(output_file, "w") as f:
                f.write(result)


class ReduceArray(Reduce):
    """
    Run reduce in terms of Biotite AtomArray.
    """

    def run(self, atom_array: AtomArray, *args, **kwargs):
        input_file = Path(tempfile.gettempdir()) / "reduce_input.pdb"
        output_file = Path(tempfile.gettempdir()) / "reduce_output.pdb"
        PDB2PDB(output_file=input_file).write(atom_array)
        result = super().run(*args, f"{input_file}", get_output=True, **kwargs)
        with open(output_file, "w") as f:
            f.write(result)
        reduced_array = PDB2PDB(input_file=output_file).read()
        return reduced_array
