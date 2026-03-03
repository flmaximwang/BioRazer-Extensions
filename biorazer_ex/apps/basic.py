from abc import abstractmethod
from pathlib import Path
import logging, subprocess, selectors, pty, os, sys


class App:

    def __init__(self, app_dir=None, app_bin=None, logger=None):
        self.dir = Path(app_dir) if app_dir else None
        self.bin = Path(app_bin) if app_bin else None
        if logger is None:
            self.logger = logging.getLogger(__name__)
            # Check if the logger has handlers already to avoid duplicate logs
            if not self.logger.hasHandlers():
                formatter = logging.Formatter(
                    "[%(asctime)s @ %(name)s] %(levelname)s: %(message)s"
                )
                console_handler = logging.StreamHandler(sys.stdout)
                console_handler.setFormatter(formatter)
                console_handler.setLevel(logging.INFO)
                self.logger.addHandler(console_handler)
                self.logger.setLevel(logging.INFO)
        else:
            self.logger = logger
        self.logger.debug(f"Initialized App with dir: {self.dir}, bin: {self.bin}")

    @abstractmethod
    def from_default_bin(cls, bin_name):
        app = cls(app_dir=None, app_bin=bin_name)
        return app

    @abstractmethod
    def set_dir(self, app_dir):
        self.dir = Path(app_dir)

    @abstractmethod
    def set_bin(self, app_bin):
        self.bin = Path(app_bin)

    @abstractmethod
    def check(self):
        """
        Implement in subclass to check if the app is installed correctly.
        """

    def run(self, *args, cwd=".", get_output=True, verbose=True, mode="pty"):
        """
        Run self.bin with given args and kwargs.


        Parameters
        ----------
        mode : str
            Mode to run the subprocess.
            - pty: Pseudo-terminal mode to ensure line-buffered output.
            - subprocess.run: Using subprocess.run to capture output.
            - subprocess.Popen: Using subprocess.Popen to capture output.
        """
        cmd_args = [f"{self.bin}"]
        for arg in args:
            cmd_args.append(str(arg))
        self.logger.info(f"Running command: {' '.join(cmd_args)}")

        if mode == "pty":
            # 通过 PTY 确保 C++ 的输出采用行缓冲, 而不是全缓冲
            master_stdout, slave_stdout = pty.openpty()
            master_stderr, slave_stderr = pty.openpty()

            p = subprocess.Popen(
                cmd_args,
                cwd=cwd,
                stdout=slave_stdout,
                stderr=slave_stderr,
                close_fds=True,
            )
            os.close(slave_stdout)  # Prevent child from hanging on to the fd
            os.close(slave_stderr)

            output = ""
            error = ""
            finish = False
            sel = selectors.DefaultSelector()
            for master_fd in (master_stdout, master_stderr):
                sel.register(
                    master_fd, selectors.EVENT_READ
                )  # Register for read events

            line_left = ""
            while True:
                self.logger.debug("Process is still running...")

                self.logger.debug(f"Checking if process has finished: {p.poll()}")
                if p.poll() is not None:
                    finish = True

                events = sel.select(timeout=1.0)  # Wait for events with a timeout
                self.logger.debug(f"Events received: {events}")

                if finish and len(events) == 0:
                    break  # Process finished and no more events

                for key, mask in events:
                    if (
                        mask & selectors.EVENT_READ
                    ):  # Bitwise AND to check whether mask contains EVENT_READ
                        data = os.read(key.fileobj, 8192)
                        if not data:
                            sel.unregister(key.fileobj)  # No more data to read
                            continue

                        lines_tmp = data.decode().splitlines(keepends=True)
                        self.logger.debug(f"Line left: {line_left}")
                        self.logger.debug(f"Read data: {lines_tmp}")
                        if line_left:
                            lines_tmp[0] = line_left + lines_tmp[0]
                            line_left = ""
                        if lines_tmp[-1].endswith("\n"):
                            line_left = ""
                        else:
                            line_left = lines_tmp.pop()
                        self.logger.debug(f"Line left: {line_left}")
                        for line in lines_tmp:
                            if key.fileobj == master_stderr:
                                if verbose:
                                    self.logger.error(line.strip())
                                if get_output:
                                    output += line
                                error += line
                            elif key.fileobj == master_stdout:  # master_stdout
                                if verbose:
                                    self.logger.info(line.strip())
                                if get_output:
                                    output += line
                            else:
                                raise RuntimeError(f"Unknown file object {key.fileobj}")
                    else:
                        raise RuntimeError(f"Unknown event mask {mask}")

            if p.returncode != 0:
                if error:
                    raise RuntimeError(f"{self.bin} failed with error:\n{error}")
                else:
                    raise RuntimeError(f"{self.bin} failed with output\n{output}")

            os.close(master_stdout)
            os.close(master_stderr)
            if get_output:
                return output
        elif mode == "subprocess.run":
            result = subprocess.run(
                cmd_args,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode != 0:
                if result.stderr:
                    raise RuntimeError(
                        f"{self.bin} failed with error:\n{result.stderr}"
                    )
                else:
                    raise RuntimeError(
                        f"{self.bin} failed with output\n{result.stdout}"
                    )

            if result.stdout:
                output = result.stdout
            else:
                output = result.stderr
            if verbose:
                self.logger.info("\n" + output)
            if get_output:
                return output
        elif mode == "subprocess.Popen":
            p = subprocess.Popen(
                cmd_args,
                bufsize=1,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stdout = ""
            stderr = ""
            while p.poll() is None:
                self.logger.debug("Process is still running...")
                for line in p.stdout:
                    if verbose:
                        self.logger.info(line.strip())
                    stdout += line
                for line in p.stderr:
                    if verbose:
                        self.logger.error(line.strip())
                    stderr += line

            stderr += p.stderr.read()
            stdout += p.stdout.read()

            if p.returncode != 0:
                if stderr:
                    raise RuntimeError(f"{self.bin} failed with error:\n{stderr}")
                else:
                    raise RuntimeError(f"{self.bin} failed with output\n{stdout}")

            if stdout:
                output = stdout
            else:
                output = stderr
            if verbose:
                self.logger.info("\n" + output)
            if get_output:
                return output

        else:
            raise ValueError(f"Unknown mode: {mode}")
