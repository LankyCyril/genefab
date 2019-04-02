from setuptools import setup
from argparse import Namespace

meta = Namespace(
    __name__ = "genefab",
    __version__ = "0.1",
    __author__ = "Kirill Grigorev",
    __git_id__ = "LankyCyril",
    __license__ = "GPLv3",
)

def get_readme(filename):
    with open(filename) as readme_handle:
        return readme_handle.read()

if __name__ == "__main__":
    setup(
        name = meta.__name__,
        version = meta.__version__,
        packages = [meta.__name__],
        url = "https://github.com/LankyCyril/genefab/",
        author = meta.__author__,
        license = meta.__license__,
        zip_safe = True,
        description = "NASA GeneLab API wrapper",
        long_description = get_readme("README.md"),
        long_description_content_type = "text/markdown",
        classifiers = [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 3",
            "Environment :: Console",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        install_requires = ["numpy", "pandas", "tqdm"]
    )
