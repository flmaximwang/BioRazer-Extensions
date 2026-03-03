from ..basic import App
from pathlib import Path


class MASTERApp(App):
    """
    Do not provide app_bin
    """

    def run_createPDS(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "createPDS"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)

    def run_parsePDS(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "parsePDS"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)

    def run_master(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "master"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)
