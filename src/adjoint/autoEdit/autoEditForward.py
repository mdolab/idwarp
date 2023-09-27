#! /usr/bin/env python
"""
autoEdit - A Python tool to automatically edit a set of files
           according to the specified user rules:
           1) Discard module files
           2) Rename module names in subroutines
"""

# Import modules
import sys
import utils

DIR_ORI = sys.argv[1]
DIR_MOD = sys.argv[2]

utils.checkDir(DIR_ORI)
utils.checkDir(DIR_MOD)

# Check that we actually have files to edit
utils.checkFiles(DIR_ORI)

# Specify file extension
EXT = "_d.f90"

# Specify the list of LINE ID's to find, what to replace and with what
# First set: Find line with 'USE' and replace "_D" with ''
LINE_ID = ["USE"]
STR_OLD = ["_D"]
STR_NEW = [""]

STR_REPLACE_ALL = {"_CD": ""}
# Also, entirely ignore lines with these strings:
LINE_IGNORE = []
FILE_IGNORE = ["gridData_d.f90", "kd_tree_d.f90", "warpMesh_d.f90"]

# Process files
utils.processFiles(EXT, DIR_ORI, DIR_MOD, LINE_ID, STR_OLD, STR_NEW, STR_REPLACE_ALL, LINE_IGNORE, FILE_IGNORE)
