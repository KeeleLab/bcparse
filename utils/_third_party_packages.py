#! /usr/bin/env python3

import subprocess
import sys
import importlib

# List of required third-party packages
required_packages = [
    'pandas',
    'tqdm',
    'numpy',
    'rapidfuzz',
    'regex',
]

 #'levenshtein', depr in favor of 'rapidfuzz'

def install(package):
    """Install a package using pip."""
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

def check_and_install_packages(packages):
    """Check for each package and install if missing."""
    for package in packages:
        try:
            importlib.import_module(package)
            print(f"{package} is already installed.")
        except ImportError:
            print(f"{package} is not installed. Installing...")
            install(package)

if __name__ == "__main__":
    check_and_install_packages(required_packages)
