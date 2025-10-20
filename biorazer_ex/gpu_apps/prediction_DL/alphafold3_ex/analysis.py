import re, json, shutil, pickle, warnings
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from ..analysis_basic import SinglePrediction, BatchPrediction
from biorazer.structure.io import CIF2PDB, CIF2CIF


class SingleAF3Local(SinglePrediction):
    """
    This class is used to analyze the AF3 local prediction results.
    """

    def format_output(self, target_dir):

        target_marker_dir, _ = super().format_output(target_dir)
        entry_dir = self.get_dir_path()

        for i, sample in enumerate(entry_dir.glob("seed-*_sample-*")):

            if not sample.is_dir():
                continue

            target_sample_dir_path = Path(target_marker_dir) / sample.stem
            if not target_sample_dir_path.exists():
                target_sample_dir_path.mkdir(parents=True)

            cif_src = sample / (entry_dir.stem + "_" + sample.stem + "_model.cif")
            cif_target = target_sample_dir_path / (
                entry_dir.stem + "_" + sample.stem + "_model.cif"
            )
            shutil.copyfile(cif_src, cif_target)

            pdb_target = cif_target.with_suffix(".pdb")
            parser = CIF2PDB(cif_src, pdb_target)
            structure = parser.read(extra_fields=["B_iso_or_equiv"])
            structure.set_annotation("b_factor", structure.B_iso_or_equiv.astype(float))
            parser.write(structure)

            confidences_json = sample / (
                entry_dir.stem + "_" + sample.stem + "_confidences.json"
            )
            confidence_data = json.load(open(confidences_json, "r"))
            summary_confidences_json = sample / (
                entry_dir.stem + "_" + sample.stem + "_summary_confidences.json"
            )
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))

            meta_data = {}
            meta_data.update(summary_confidences_data)
            meta_data.update(confidence_data)
            meta_data_pickle_target = target_sample_dir_path / (
                entry_dir.stem + "_" + sample.stem + "_data.pickle"
            )
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return target_marker_dir, SingleAF3Local(target_marker_dir)


class BatchAF3Local(BatchPrediction):

    single_pred_cls = SingleAF3Local


class SingleAF3Server(SinglePrediction):
    """
    This class is used to analyze the AF3 server prediction results.
    """

    def format_output(self, target_dir):

        warnings.filterwarnings(
            "ignore",
            category=UserWarning,
            module="biotite.structure.io.pdbx",
        )

        marker_dir_path_src = self.get_dir_path()
        marker_src = self.marker
        marker_dir_target, _ = super().format_output(target_dir)
        marker_dir_path_target = Path(marker_dir_target)
        marker_target = marker_dir_path_target.name

        job_request_json = f"fold_{marker_src}_job_request.json"
        job_request = json.load(open(marker_dir_path_src / job_request_json, "r"))
        seed = job_request[0]["modelSeeds"][0]

        cif_list = [f"fold_{marker_src}_model_{i}.cif" for i in range(5)]
        full_data_list = [f"fold_{marker_src}_full_data_{i}.json" for i in range(5)]
        summary_confidences_list = [
            f"fold_{marker_src}_summary_confidences_{i}.json" for i in range(5)
        ]

        for i in range(5):

            sample_seed_dir = f"seed-{seed}_sample-{i+1}"
            sample_seed_dir_path = marker_dir_path_target / sample_seed_dir
            if not sample_seed_dir_path.exists():
                sample_seed_dir_path.mkdir(parents=True)

            job_request_json_target = (
                marker_dir_path_target
                / sample_seed_dir
                / (f"{marker_target}_seed-{seed}_sample-{i+1}_job_request.json")
            )
            json.dump(job_request, open(job_request_json_target, "w"), indent=4)

            cif_src = marker_dir_path_src / cif_list[i]
            cif_target = (
                marker_dir_path_target
                / sample_seed_dir
                / f"{marker_target}_seed-{seed}_sample-{i+1}_model.cif"
            )
            shutil.copyfile(cif_src, cif_target)

            pdb_target = cif_target.with_suffix(".pdb")
            parser = CIF2PDB(cif_src, pdb_target)
            structure = parser.read(extra_fields=["B_iso_or_equiv"])
            structure.set_annotation("b_factor", structure.B_iso_or_equiv.astype(float))
            parser.write(structure)

            full_data_json = marker_dir_path_src / full_data_list[i]
            full_data = json.load(open(full_data_json, "r"))
            summary_confidences_json = marker_dir_path_src / summary_confidences_list[i]
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))

            meta_data = {}
            for key in summary_confidences_data:
                meta_data[key] = summary_confidences_data[key]
            for key in full_data:
                meta_data[key] = full_data[key]
            meta_data_pickle_target = (
                marker_dir_path_target
                / sample_seed_dir
                / (f"{marker_target}_seed-{seed}_sample-{i+1}_data.pickle")
            )
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        warnings.filterwarnings(
            "default",
            category=UserWarning,
            module="biotite.structure.io.pdbx",
        )

        return marker_dir_target, SingleAF3Server(marker_dir_target)


class BatchAF3Server(BatchPrediction):

    single_pred_cls = SingleAF3Server
