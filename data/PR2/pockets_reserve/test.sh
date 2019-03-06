#!/bin/bash

for file in *_out/
do
    cd $(file/pockets/)
    echo $*.pdb
    cd ..
done
