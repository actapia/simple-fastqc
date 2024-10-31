# This trick comes from Paulo Scardine's answer on StackOverflow. It allows
# postponed evaluation of type annotations in older versions of Python.
# https://stackoverflow.com/a/33533514
from __future__ import annotations

import pandas as pd
from dataclasses import dataclass
import subprocess
import io
from pathlib import Path

from typing import Optional, Dict, Union, List

@dataclass
class FastQCModule:
    """Results from a FastQC analysis module.

    This class holds the name, details, and status of one of the analysis
    modules run as part of a FastQC analysis.

    Attributes:
        name (str):    The name of the FastQC module.
        data (str):    The raw results of running the module.
        status (str):  Whether the module passed, failed, or produced a warning.
    """
    name: str = None
    data: str = None
    status: str = None
    
    @classmethod
    def from_module(cls, module: FastQCModule) -> FastQCModule:
        """Construct a module of this class from another module.

        This method may be viewed as a way to copy another FastQCModule.

        Classes inheritting from FastQCModule may override this method. This
        makes it possible for the FastQCResults class to instantiate modules
        with the appropriate class.

        Parameters:
            module:  The FastQCModule from which to construct a FastQCModule.

        Returns:
            A FastQCModule with the same name, data, and status as the argument.
        """
        return cls(module.name, module.data, module.status)
    
def read_csv_or_empty(data: str) -> Optional[pd.DataFrame]:
    """Parse the given string as a CSV or return None if data is empty."""
    if data:
        return pd.read_csv(
            io.StringIO(data),
            sep="\t"
        ).rename(columns=lambda n: n.removeprefix("#"))
    else:
        return None

class OverrepresentedSequencesModule(FastQCModule):
    """A FastQCModule testing for overrepresented sequences.

    See the documentation for FastQC for details on the specifics of how this
    module works.

    Attributes:
        name (str):    The name of the FastQC module
                       (should be "Overrepresented Sequences")
        data (str):    The results of the module, parsed into a dataframe.
        status (str):  Whether the module passed, failed, or produced a warning.
    """
    @classmethod
    def from_module(
            cls,
            module: FastQCModule
    ) -> OverrepresentedSequencesModule:
        """Construct an OverrepresentedSequenceModule from another FastQCModule.

        Unlike the base FastQCModule class, this module parses the results of
        the module to construct a pandas dataframe. If the provided module
        represents the results as a string, this function will parse the string
        to produce a dataframe.
        
        Parameters:
            module:  The FastQCModule from which to construct an
                     OverrepresentedSequencesModule.

        Returns:
            An OverrepresentedSequencesModule with the same name, data, and
            status as the argument.
        """
        return OverrepresentedSequencesModule(
            module.name,
            read_csv_or_empty(module.data),
            module.status
        )

class BasicStatisticsModule(FastQCModule):
    @classmethod
    def from_module(
            cls,
            module: FastQCModule
    ) -> BasicStatisticsModule:
        data_dict = dict(l.split("\t") for l in module.data.splitlines())
        data_dict["FastQC"] = data_dict["##FastQC"]
        keys = list(data_dict)
        for k in keys:
            if k.startswith("#"):
                del data_dict[k]
            else:
                try:
                    as_int = int(data_dict[k])
                    as_float = float(data_dict[k])
                    if as_int == as_float:
                        data_dict[k] = as_int
                except ValueError:
                    pass
        # tb_split = data_dict["Total Bases"].split()
        # data_dict["Total Bases"] = float(tb_split[0]) * si[tb_split[1][0]]
        return BasicStatisticsModule(
            module.name,
            data_dict,
            module.status
        )
    
class FastQCResults:
    """Results of a FastQC analysis.

    This class represents the results of running FastQC on one file containing
    sequence reads. It contains the results of the individual modules run as
    part of FastQC; these results are represented as instances of the
    appropriate subclasses of FastQCModule. Individual modules can be accessed
    by name through the modules attribute (or the equivalent results property).

    Attributes:
        modules (dict):  A dictionary mapping module names to FastQCModules.
    """

    # This maps module names to subclasses of FastQCModule. If no subclass
    # exists for a specific module name, the default FastQCModule class is used.
    module_classes = {
        "Overrepresented sequences": OverrepresentedSequencesModule,
        "Basic Statistics": BasicStatisticsModule
    }
    
    def __init__(self):
        self.modules = {}

    @property
    def results(self) -> Dict[str, FastQCModule]:
        """Return a dictionary mapping module names to FastQCModule results."""
        return self.modules
        
    @classmethod
    def parse_from_file(cls, path: Union[str, Path]) -> FastQCResults:
        """Parse the results of a FastQC analysis from a file.

        The file to be parsed must be in the same format as the fastqc_data.txt
        file produced by FastQC. (Note that this file is normally inside the
        zip file produced by FastQC.)

        Parameters:
            path:  The path to the file to parse.

        Returns:
            A FastQCResults object representing the parsed data from the file.
        """
        results = FastQCResults()
        with open(path, "r") as fastqc_file:
            mod = None
            data_lines = []
            for l in  fastqc_file:
                l = l.rstrip()
                if l.startswith(">>"):
                    l = l.removeprefix(">>")
                    if mod is None:
                        res = l.split("\t")
                        mod = FastQCModule(name=res[0], status=res[1])
                    else:
                        assert l == "END_MODULE"
                        mod.data = "\n".join(data_lines)
                        data_lines = []
                        results.results[mod.name] = cls.module_classes.get(
                            mod.name,
                            FastQCModule
                        ).from_module(mod)
                        mod = None
                else:
                    data_lines.append(l)
        return results
    
    def status_dict(self) -> Dict[str, str]:
        """Return a dictionary mapping modules' names to their statuses."""
        return {m.name: m.status for m in self.results.values()}

class FastQCAnalysis:
    """A read quality analysis to be run with FastQC.

    This class provides a fairly straightforward way of setting up and executing
    a FastQC analysis on one or more files containing sequence reads with
    quality information.

    Parameters of the analysis, including the files to be analyzed, are
    specified as the analysis is constructed. The results are produced lazily
    once they are accessed via the results property.

    This class has no public attributes. The settings of the analysis should be
    considered immutable, and a new analysis should be constructed if the
    need arises to change settings.

    Nevertheless, the settings passed to the constructor may be accessed via
    various properties of this class.
    """
    def __init__(
            self,
            read_paths: List[Union[str, Path]],
            out_dir: Optional[Union[str, Path]] = None,
            threads: int = 1
    ):
        """Construct a FastQCAnalysis with the specified settings.

        This constructor requires paths to the reads to be analyzed.

        Optionally, the caller may provide an alternative output directory and
        may specify the number of concurrent threads to use in the analysis.

        If no output directory is specified, FastQC will default to storing the
        results in the same directory as each input file.

        By default, FastQC will be run with just one thread.

        Parameters:
            read_paths:    Paths to files containing sequence reads.
            out_dir:       An alternative directory in which to store results.
            threads (int): Number of concurrent threads to use in FastQC.
        """
        self._read_paths = [Path(p) for p in read_paths]
        self._out_dir = Path(out_dir)
        self._threads = threads
        self._results = None

    @property
    def read_paths(self) -> List[Path]:
        """Return the list of paths to read files used for the analysis."""
        return list(self._read_paths)

    @property
    def out_dir(self) -> Path:
        """Return the directory in which results will be stored."""
        return self._out_dir

    @property
    def threads(self) -> int:
        """Return the number of concurrent threads to use in FastQC."""
        return self._threads

    @property
    def results(self) -> Dict[str, FastQCResults]:
        """Return the results of the FastQC analysis.

        If necessary, the analysis will be executed first to obtain the results.

        Returns:
            A FastQCResults object representing the FastQC analysis results.
        """
        if self._results is None:
            self._get_results()
        return self._results

    def _build_fastqc_command(self):
        command = [
            "fastqc",
            "--extract",
            "-T",
            str(self._threads),
        ]
        if self._out_dir is not None:
            command = command + ["-o", str(self._out_dir)]
        command = command + [str(p) for p in self._read_paths]
        return command

    def _get_results(self):
        proc = subprocess.Popen(
            self._build_fastqc_command(),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        proc.communicate()
        if proc.returncode:
            raise subprocess.CalledProcessError(proc.returncode, proc.args)
        # Read output files.
        self._results = {
            path : FastQCResults.parse_from_file(
                self._out_dir / (
                    path.name.split(".")[0] + "_fastqc"
                ) / "fastqc_data.txt"
            )
            for path in self._read_paths
        }  
