from ..basic import App
import os
from ...utils.structure_file import atom_array_as_temp_file


class RosettaApp(App):

    def run(self, *args, cwd=".", get_output=False, verbose=True, mode="pty", **kwargs):
        return super().run(
            *args, cwd=cwd, get_output=get_output, verbose=verbose, mode=mode, **kwargs
        )

    def run_with_structure(
        self,
        atom_array,
        *args,
        input_file_format="pdb",
        input_file_flag="-s",
        cwd=".",
        get_output=False,
        verbose=True,
        mode="pty",
        **kwargs,
    ):
        with atom_array_as_temp_file(
            atom_array,
            temp_file_format=input_file_format,
        ) as input_file:
            run_args = list(args)
            if input_file_flag is None:
                run_args.append(str(input_file))
            else:
                run_args.extend([str(input_file_flag), str(input_file)])
            return self.run(
                *run_args,
                cwd=cwd,
                get_output=get_output,
                verbose=verbose,
                mode=mode,
                **kwargs,
            )

    def find_app(self, app_keywords: list[str]):

        results = []
        for parent, dirs, files in self.dir.walk():
            # If the parent directory is not "bin" and does not contain "source", skip it
            if not (
                parent.stem == "bin" and "source" in str(parent).split(os.path.sep)
            ):
                continue
            for file in files:
                flag = True
                for app_keyword in app_keywords:
                    if not app_keyword in file:
                        flag = False
                        break
                if flag:
                    results.append(parent / file)
        return results

    def use_app(self, app_keywords: list[str]):
        results = self.find_app(app_keywords)
        if not results:
            raise RuntimeError(
                f"Failed to find app with keywords {app_keywords} in {self.dir}"
            )
        if len(results) > 1:
            raise RuntimeError(
                f"Found multiple apps with keywords {app_keywords} in {self.dir}: {results}"
            )
        self.bin = results[0]
        self.logger.info(f"Using app {self.bin} for keywords {app_keywords}")

    def find_tool(self, tool_keyword):

        results = []
        for parent, dirs, files in self.dir.walk():
            if not ("tools" in str(parent).split(os.path.sep)):
                continue
            for file in files:
                if tool_keyword in file:
                    results.append(parent / file)
        return results
