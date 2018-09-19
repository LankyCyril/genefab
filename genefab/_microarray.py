from ._dataset import GLDS
from ._util import ensure_dir, gunzip, routput
from tqdm import tqdm
from os import walk, getcwd
from os.path import join, isfile, relpath
from re import search, sub
from tarfile import TarFile
from subprocess import call
from glob import iglob
from pandas import read_csv
from functools import lru_cache

AFFY_SCRIPT = """
if (!require("affy")) {{
    source("https://bioconductor.org/biocLite.R")
    biocLite("affy", suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
}}
library("affy")
eset = just.rma(filenames=c({cels}))
write.exprs(eset, file={output})
"""

class MicroarrayExperiment():
    """Implements wrapper of GLDS class that has 'Array Data Files'"""
    glds = None
    accession = None
    design_ref = None
    factors = None
    raw_data = None
    derived_data = None
    _file_list = None
    _storage = None
    _tsv = None
    _R = None
 
    def __init__(self, glds, reextract=False, R="Rscript"):
        """Interpret GLDS, describe; if data has been unpacked before, check and reuse"""
        self.glds = glds
        self.accession = glds.accession
        self._storage = join(glds._storage, "MicroarrayExperiment_source")
        self._R = R
        self.factors = glds.factors(as_fields=False)
        if glds.field_ids("Array Design REF"):
            self.design_ref = glds.property_table("Array Design REF")
        if glds.field_ids("Array Data File"):
            self.raw_data = glds.property_table("Array Data File")
        if glds.field_ids("Derived Array Data File"):
            self.derived_data = glds.property_table("Derived Array Data File")
        if (self.raw_data is None) and (self.derived_data is None):
            raise ValueError("No raw or derived data associated with factors")
        self._tsv = join(getcwd(), self._storage, self.accession+".tsv")
        if (not reextract) and isfile(self._tsv):
            with open(self._tsv, "rt") as tsv_handle:
                putative_file_list = set(map(str.strip, tsv_handle))
                for filename in putative_file_list:
                    if not isfile(filename):
                        break
                else:
                    self._file_list = putative_file_list
 
    @property
    def raw_only(self):
        return (self.raw_data is not None) and (self.derived_data is None)
    @property
    def derived_only(self):
        return (self.raw_data is None) and (self.derived_data is not None)
 
    def _store_file_list(self):
        """List all unpacked files into a TSV file"""
        with open(self._tsv, "wt") as tsv_handle:
            for filename in self._file_list:
                print(filename, file=tsv_handle)
 
    def _unpack(self, force_new_dir=True, update_glds_files=False):
        """Store unpacked contents of associated archive files in experiment subdirectory"""
        self.glds.fetch_files(update=update_glds_files)
        target_dir = join(getcwd(), self._storage, self.accession)
        ensure_dir(target_dir, force_new_dir=force_new_dir)
        file_list_iterator = tqdm(
            self.glds.file_list, desc="Unpacking top-level files", unit="file"
        )
        for filename in file_list_iterator:
            source_file = join(getcwd(), self.glds._storage, filename)
            if search(r'\.tar$|\.tar\.gz$', filename):
                with TarFile(source_file) as tar:
                    tar.extractall(path=target_dir)
            elif search(r'\.gz$', filename):
                gunzip(source_file, target_dir=target_dir, keep_original=True)
            elif search(r'\.zip$', filename):
                call(["unzip", source_file, "-d", target_dir])
            else:
                call(["cp", source_file, target_dir])
        second_level_iterator = tqdm(
            next(walk(target_dir))[2], desc="Unpacking second-level files",
            unit="file"
        )
        for filename in second_level_iterator:
            if search(r'\.gz$', filename):
                gunzip(join(getcwd(), target_dir, filename))
        self._file_list = set(
            relpath(name)
            for name in iglob(join(target_dir, "**"), recursive=True)
            if isfile(name)
        )
        self._store_file_list()
 
    @property
    def file_list(self):
        if self._file_list is None:
            self._unpack()
        if self._file_list:
            return self._file_list
        else:
            raise OSError("No files associated with experiment")
 
    @property
    @lru_cache(maxsize=None)
    def annotation(self):
        if self.raw_data is None:
            raise NotImplementedError("Experiment without raw data")
        _annotation = self.raw_data.copy()
        _property = _annotation.columns[-1]
        _annotation["filename"] = None
        for ix, row in _annotation.iterrows():
            tail = sub(r'^\*', "", row[_property])
            matching_files = set(
                filename
                for filename in self._file_list
                if filename.endswith(tail)
            )
            if len(matching_files) != 1:
                err_msg = "Number of files matching '{}' is not 1".format(tail)
                raise OSError(err_msg)
            _annotation.loc[ix, "filename"] = matching_files.pop()
        del _annotation[_property]
        return _annotation
 
    @property
    def expression_table(self):
        """Run R script with bioconductor-affy on all CELs in experiment"""
        cels = ", ".join(
            repr(cel) # safeguard mainly against single backslashes on Windows
            for cel in self.annotation["filename"]
        )
        expression_table_raw = routput(self._R, AFFY_SCRIPT, cels=cels)
        expression_table = read_csv(expression_table_raw, sep="\t", index_col=0)
        # replace filenames with sample names (trust unique):
        expression_table.columns = [
            self.annotation[self.annotation["filename"]==filename].index[0]
            for filename in expression_table.columns
        ]
        return expression_table
