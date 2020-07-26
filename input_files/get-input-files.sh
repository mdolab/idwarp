#!/bin/bash
# this file will download the input files for IDWarp regression tests
# and extract them to the right place.

wget -O input_files.tar.gz https://umich.box.com/shared/static/rrr6zruxhzbx9eba7wcfzdwmdn0fq52k
tar -xf input_files.tar.gz -C ../
