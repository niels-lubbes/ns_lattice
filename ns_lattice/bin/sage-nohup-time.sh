#! /bin/bash

#
# This script can be used to run this package on a remote server,
# via an ssh session. the standard/error output is written
# in the files err/out respectively.
#

export PYTHONPATH=$PYTHONPATH:../src/
export PYTHONPATH=$PYTHONPATH:../../linear_series/linear_series/src/
export OUTPUT_PATH=../../OUTPUT/

rm err out
#nohup time nice sage -python ../src/ns_lattice/__main__.py  > out 2> err < /dev/null &
time sage -python ../src/ns_lattice/__main__.py  > out 2> err < /dev/null &
cat err
cat out
