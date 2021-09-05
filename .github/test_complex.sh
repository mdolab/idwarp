#!/bin/bash
set -e
./input_files/get-input-files.sh
cd tests
testflo -v -n 1 --coverage --coverpkg idwarp
