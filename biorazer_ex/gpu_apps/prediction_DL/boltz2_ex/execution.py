import os, yaml
from pathlib import Path
from ..execution_basic import SingleExecution, BatchExecution


class Boltz2Modification:

    def __init__(self, position: int, ccd: str):
        self.position = position
        self.ccd = ccd

    def output(self):
        return {
            "position": self.position,
            "ccd": self.ccd,
        }


class Boltz2Sequence:

    def __init__(
        self,
        entity_type: str,
        id: list[str],
        sequence: str = None,
        msa: str = None,
        smiles: str = None,
        ccd: str = None,
        modifications: list[Boltz2Modification] = None,
    ):
        self.id = id
        self.entity_type = entity_type

        if entity_type == "protein":
            self.sequence = sequence
            self.msa = msa
            self.modifications = modifications if modifications is not None else []
        elif entity_type in ["dna", "rna"]:
            self.sequence = sequence
            self.modifications = modifications if modifications is not None else []
        elif entity_type == "ligand":
            self.smiles = smiles
            self.ccd = ccd
        else:
            raise ValueError(
                "Invalid entity type. Entity type must be 'protein', 'dna', 'rna', or 'ligand'"
            )

    def seq_from_msa(self, a3m):
        with open(a3m, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    self.sequence = line.strip()
                    break
        self.set_msa(a3m)

    def set_msa(self, a3m: str | Path):
        with open(a3m, "r") as f:
            self.msa = f.read()

    def set_no_msa(self):
        self.msa = f">query\n{self.sequence}\n"

    def output(self):
        output_dict = {
            "id": self.id,
        }
        if self.entity_type == "protein":
            output_dict = {
                self.entity_type: {
                    "id": self.id,
                    "sequence": self.sequence,
                }
            }
            if not self.msa is None:
                output_dict[self.entity_type]["msa"] = self.msa
            if len(self.modifications) > 0:
                output_dict[self.entity_type]["modifications"] = [
                    mod.output() for mod in self.modifications
                ]
            return output_dict
        elif self.entity_type in ["dna", "rna"]:
            output_dict = {
                self.entity_type: {
                    "id": self.id,
                    "sequence": self.sequence,
                }
            }
            if len(self.modifications) > 0:
                output_dict[self.entity_type]["modifications"] = [
                    mod.output() for mod in self.modifications
                ]
            return output_dict
        elif self.entity_type == "ligand":
            # 优先使用 SMILES
            if not self.smiles is None:
                return {
                    self.entity_type: {
                        "id": self.id,
                        "smiles": self.smiles,
                    }
                }
            else:
                return {
                    self.entity_type: {
                        "id": self.id,
                        "ccd": self.ccd,
                    }
                }
        else:
            raise ValueError(
                "Invalid entity type. Entity type must be 'protein', 'dna', 'rna', or 'ligand'"
            )


class Boltz2BondAtom:

    def __init__(self, chain_id, res_id, atom_name):
        self.chain_id = chain_id
        self.res_id = res_id
        self.atom_name = atom_name

    def output(self):
        return [self.chain_id, self.res_id, self.atom_name]


class Boltz2Constraint:

    def __init__(
        self,
        type: str,
        atom1: Boltz2BondAtom = None,
        atom2: Boltz2BondAtom = None,
        binder: str = None,
        contacts: list[list] = None,
    ):
        """

        Parameters
        ----------
        binder: str
            The chain id of the binder molecule.
        contacts: list of list
            A list of contacts, where each contact is represented as a list of [chain_id, res_id].
        """
        self.type = type
        if type == "bond":
            assert (
                atom1 is not None and atom2 is not None
            ), "Bond constraints require atom1 and atom2"
            self.atom1 = atom1
            self.atom2 = atom2
        elif type == "binder":
            assert (
                binder is not None and contacts is not None
            ), "Binder constraints require binder and contacts"
            self.binder = binder
            self.contacts = contacts
        else:
            raise ValueError(
                "Invalid constraint type. Constraint type must be 'bond' or 'binder'"
            )

    def output(self):
        if self.type == "bond":
            return {
                "bond": {
                    "atom1": self.atom1.output(),
                    "atom2": self.atom2.output(),
                }
            }
        elif self.type == "binder":
            return {
                "binder": {
                    "binder": self.binder,
                    "contacts": self.contacts,
                }
            }
        else:
            raise ValueError(
                "Invalid constraint type. Constraint type must be 'bond' or 'binder'"
            )


class SingleBoltz2Ex(SingleExecution):

    def __init__(
        self,
        name: str,
        sequences: list[Boltz2Sequence],
        constraints: list[Boltz2Constraint] = None,
    ):
        self.name = name
        self.sequences = sequences
        self.constraints = constraints if constraints is not None else []

    def output(self):
        output_dict = {
            "sequences": [],
        }
        for seq in self.sequences:
            output_dict["sequences"].append(seq.output())
        if len(self.constraints) > 0:
            output_dict["constraints"] = []
            for con in self.constraints:
                output_dict["constraints"].append(con.output())
        return output_dict


class BatchBoltz2Ex(BatchExecution):

    def __init__(self, executions: list[SingleBoltz2Ex], **extra_info):
        super().__init__(executions)
        self.extra_info = extra_info

    def write(self, target_dir):
        target_dir_path = Path(target_dir)
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True)

        ori_cwd = os.getcwd()
        os.chdir(target_dir_path)
        for execution in self.executions:
            for sequence in execution.sequences:
                if sequence.entity_type == "protein" and not sequence.msa is None:
                    if Path(sequence.msa).suffix == ".a3m":
                        if Path(sequence.msa).exists():
                            continue
                        else:
                            raise FileNotFoundError(
                                f"MSA file {sequence.msa} not found."
                            )
                    else:
                        a3m = f"{execution.name}_{''.join(sequence.id)}.a3m"
                        with open(a3m, "w") as f:
                            f.write(sequence.msa)
                            sequence.msa = f"./{a3m}"
            output_dict: dict = execution.output()
            output_dict.update(self.extra_info)
            with open(f"{execution.name}.yaml", "w") as f:
                yaml.dump(output_dict, f, indent=2)
        os.chdir(ori_cwd)
