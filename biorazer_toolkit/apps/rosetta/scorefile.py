from dataclasses import dataclass
from pathlib import Path
import pandas as pd


@dataclass
class Scorefile:
    data: pd.DataFrame

    @classmethod
    def from_sc(cls, sc_file: str | Path):
        df = pd.read_csv(sc_file, comment="#", sep=r"\s+")
        df = df.iloc[:, 1:]  # drop the first column which is just 'SCORE:'
        return cls(data=df)

    def to_csv(self):
        self.data.to_csv(index=False)
