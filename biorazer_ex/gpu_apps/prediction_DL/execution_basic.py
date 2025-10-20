from dataclasses import dataclass
from abc import ABC, abstractmethod


@dataclass
class SingleExecution(ABC):

    @abstractmethod
    def output(self):
        """
        Implement in subclasses to return a dict that can easily be parsed by json, yaml, etc.
        """


@dataclass
class BatchExecution(ABC):
    executions: list[SingleExecution]

    @abstractmethod
    def write(self, target_dir, jobs_per_file=None):
        """
        Implement in subclasses to write batch predictions to files.
        This function should call the output() method of each SinglePrediction instance to gather data.

        Parameters
        ----------
        target_dir : str or Path
            The directory where the batch prediction files will be saved.
        jobs_per_file : int, optional
            The number of jobs per file. If None, all predictions are written to a single file.
        """
