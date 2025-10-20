from ..basic import App
from pathlib import Path


class MASTERConfig(App):

    def set(self, app_dir):
        self.DIR = app_dir
        self.BIN = f"{self.DIR}/MASTER"

    def check(self):
        assert Path(self.BIN).exists(), f"MASTER binary not found at {self.BIN}"

    def get(self):
        """

        Return
        ------
        dict
            Dictionary with the config parameters.
            - DIR: MASTER directory
            - BIN: MASTER binary path
        """
        return {
            "DIR": self.DIR,
            "BIN": self.BIN,
        }
