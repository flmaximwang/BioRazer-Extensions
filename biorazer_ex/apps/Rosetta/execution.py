from ..basic import App
import os


class RosettaApp(App):

    def run(self, *args, cwd=".", get_output=False, verbose=True, mode="pty", **kwargs):
        return super().run(
            *args, cwd=cwd, get_output=get_output, verbose=verbose, mode=mode, **kwargs
        )

    def find_app(self, app_keyword):

        results = []
        for parent, dirs, files in self.dir.walk():
            if not (
                parent.stem == "bin" and "source" in str(parent).split(os.path.sep)
            ):
                continue
            for file in files:
                if app_keyword in file:
                    results.append(parent / file)
        return results

    def find_tool(self, tool_keyword):

        results = []
        for parent, dirs, files in self.dir.walk():
            if not ("tools" in str(parent).split(os.path.sep)):
                continue
            for file in files:
                if tool_keyword in file:
                    results.append(parent / file)
        return results
