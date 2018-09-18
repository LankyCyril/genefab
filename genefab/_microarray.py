from ._dataset import GLDS
from ._util import ensure_dir
from os import mkdir, getcwd, walk
from shutil import copyfileobj
from os.path import join
from re import search, sub
from tarfile import TarFile
from gzip import open as gzopen
from subprocess import call

class MicroarrayExperiment():
    """Implements wrapper of GLDS class that has 'Array Data Files'"""
    glds = None
    accession = None
    design_ref = None
    factors = None
    raw_data = None
    derived_data = None
    _file_list = None
 
    def __init__(self, glds):
        self.glds = glds
        self.accession = glds.accession
        self._storage = glds._storage
        self.factors = glds.factors(as_fields=False)
        if glds.field_ids("Array Design REF"):
            self.design_ref = glds.property_table("Array Design REF")
        if glds.field_ids("Array Data File"):
            self.raw_data = glds.property_table("Array Data File")
        if glds.field_ids("Derived Array Data File"):
            self.derived_data = glds.property_table("Derived Array Data File")
        if (self.raw_data is None) and (self.derived_data is None):
            raise ValueError("No raw or derived data associated with factors")
 
    @property
    def raw_only(self):
        return (self.raw_data is not None) and (self.derived_data is None)
    @property
    def derived_only(self):
        return (self.raw_data is None) and (self.derived_data is not None)
 
    def unpack(self, force_new_dir=True):
        target_dir = join(self._storage, self.accession)
        ensure_dir(target_dir, force_new_dir=force_new_dir)
        for filename in self.glds.file_list:
            source_file = join(getcwd(), self._storage, filename)
            if search(r'\.tar$|\.tar\.gz$', filename):
                cmd_a = "untar"
                cmd_b = None
            elif search(r'\.gz$', filename):
                cmd_a = ["gunzip", source_file]
                cmd_b = ["mv", sub(r'\.gz$', "", source_file), "."]
            elif search(r'\.zip$', filename):
                cmd_a = ["unzip", source_file, "-d", "."]
                cmd_b = None
            else:
                cmd_a = ["cp", source_file, "."]
            if cmd_a == "untar": # call(["tar", "xf", ...]) fails on Windows
                with TarFile(source_file) as tar:
                    tar.extractall(path=target_dir)
            else:
                call(cmd_a, cwd=target_dir)
            if cmd_b:
                call(cmd_b, cwd=target_dir)
        self._file_list = set()
        for filename in next(walk(target_dir))[2]:
            if search(r'\.gz$', filename):
                source_file = join(getcwd(), target_dir, filename)
                target_file = sub(r'\.gz$', "", source_file)
                # call(["gunzip", ...]) fails on Windows:
                with gzopen(source_file, "rb") as compressed:
                    with open(target_file, "wb") as uncompressed:
                        copyfileobj(compressed, uncompressed)
                self._file_list.add(sub(r'\.gz$', "", filename))
            else:
                self._file_list.add(filename)
 
    @property
    def file_list(self):
        if self._file_list is None:
            self.unpack()
        if self._file_list:
            return self._file_list
        else:
            raise OSError("No files associated with experiment")
