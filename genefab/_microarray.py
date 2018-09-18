from ._dataset import GLDS
from ._util import ensure_dir, gunzip
from tqdm import tqdm
from os import walk, getcwd
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
        self._storage = join(glds._storage, "MicroarrayExperiment_source")
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
 
    def unpack(self, force_new_dir=True, update_glds_files=False):
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
                gunzip(source_file, target_dir=target_dir)
            elif search(r'\.zip$', filename):
                call(["unzip", source_file, "-d", target_dir])
            else:
                call(["cp", source_file, target_dir])
        self._file_list = set()
        second_level_iterator = tqdm(
            next(walk(target_dir))[2], desc="Unpacking second-level files",
            unit="file"
        )
        for filename in second_level_iterator:
            if search(r'\.gz$', filename):
                gunzip(join(getcwd(), target_dir, filename))
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
