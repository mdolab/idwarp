import os
import sys


def checkDir(directory):
    if not os.path.exists(directory):
        print(f"Directory {directory} does not exist. Check input.")
        sys.exit(1)


def checkFiles(directory):
    if not os.listdir(directory):
        print(f"Directory {directory} is empty. Check input.")
        sys.exit(1)


def processFiles(EXT, DIR_ORI, DIR_MOD, LINE_ID, STR_OLD, STR_NEW, STR_REPLACE_ALL, LINE_IGNORE, FILE_IGNORE):
    print("Directory of input source files  :", DIR_ORI)
    print("Directory of output source files :", DIR_MOD)

    for f in os.listdir(DIR_ORI):
        if f not in FILE_IGNORE:
            if f.endswith(EXT):
                # open original file in read mode
                file_object_ori = open(os.path.join(DIR_ORI, f), "r")
                print(f"\nParsing input file {file_object_ori.name}")

                # open modified file in write mode
                file_object_mod = open(os.path.join(DIR_MOD, f), "w")

                # read the original file, line-by-line
                nEdits = len(LINE_ID)
                nIgnores = len(LINE_IGNORE)

                for line in file_object_ori:
                    # parse original line for relevant identifier
                    # and replace the string
                    line_mod = line.lstrip()  # Strip out Left-hand leading spaces

                    found_ignore = False
                    for i in range(nIgnores):
                        if line_mod.find(LINE_IGNORE[i]) >= 0:
                            found_ignore = True

                    if not found_ignore:
                        for i in range(nEdits):
                            if line_mod[0 : len(LINE_ID[i])] == LINE_ID[i]:
                                line_mod = line_mod.replace(STR_OLD[i], STR_NEW[i])

                    for key in STR_REPLACE_ALL:
                        line_mod = line_mod.replace(key, STR_REPLACE_ALL[key])

                    # write the modified line to new file
                    file_object_mod.write("   " + line_mod)

                # close the files
                file_object_ori.close()
                file_object_mod.close()

                # success message
                print(" Modified file saved", file_object_mod.name)
