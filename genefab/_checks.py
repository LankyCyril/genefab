from os import walk
from re import search

def safe_file_name(compound_name, as_mask=False, directory="."):
    simplified_name = compound_name.split("/")[-1].lstrip("~")
    if not as_mask:
        return simplified_name
    else:
        regex_filemask = simplified_name.replace("*", ".*")
        candidate_filenames = [
            filename for filename in next(walk(directory))[2]
            if search(regex_filemask, filename)
        ]
        if len(candidate_filenames) == 1:
            return candidate_filenames[0]
        else:
            raise ValueError("Multiple or no files matching " + simplified_name)
