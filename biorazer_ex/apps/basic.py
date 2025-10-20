from abc import abstractmethod
import logging


class App:

    def __init__(self, app_dir, app_bin, log_name="abstract_app"):
        self.dir = app_dir
        self.bin = app_bin
        self.logger = logging.getLogger(log_name)

    @abstractmethod
    def from_default_bin(cls, bin_name):
        app = cls(app_dir=None, app_bin=bin_name)
        return app

    @abstractmethod
    def set(self, app_dir):
        """
        Implement in subclass to set the app directory and binary path.
        """

    @abstractmethod
    def check(self):
        """
        Implement in subclass to check if the app is installed correctly.
        """

    @abstractmethod
    def run(self, *args, **kwargs):
        """
        Implement in method to run this app
        """
