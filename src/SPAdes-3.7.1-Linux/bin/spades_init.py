############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

source_dirs = ["", "truspades", "common"]

# developers configuration
spades_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
bin_home = os.path.join(spades_home, 'bin')
python_modules_home = os.path.join(spades_home, 'src')
ext_python_modules_home = os.path.join(spades_home, 'ext', 'src', 'python_libs')
spades_version = ''

def init():
    global spades_home
    global bin_home
    global python_modules_home
    global spades_version
    global ext_python_modules_home

    # users configuration (spades_init.py and spades binary are in the same directory)
    if os.path.isfile(os.path.join(spades_home, 'spades')):
        install_prefix = os.path.dirname(spades_home)
        bin_home = os.path.join(install_prefix, 'bin')
        spades_home = os.path.join(install_prefix, 'share', 'spades')
        python_modules_home = spades_home
        ext_python_modules_home = spades_home

    for dir in source_dirs:
        sys.path.append(os.path.join(python_modules_home, 'spades_pipeline', dir))

    spades_version = open(os.path.join(spades_home, 'VERSION'), 'r').readline().strip()
