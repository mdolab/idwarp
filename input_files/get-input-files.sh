#!/bin/bash
# this file will download the input files for IDWarp regression tests
# and extract them to the right place.

wget http://umich.edu/~mdolaboratory/repo_files/IDWarp/input_files.tar.gz
tar -xzf input_files.tar.gz -C ../
