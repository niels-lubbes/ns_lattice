#! /bin/bash

export PATH=$PATH:.
export PYTHONPATH=$PYTHONPATH:.
export OUTPUT_PATH=/home/LOCAL/nlubbes/OUTPUT/

rm err out
nohup time nice /home/software/sage-7.3/sage -python __main__.py  > out 2> err < /dev/null &
cat err
cat out

