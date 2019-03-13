#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/2HPF_pocket}"
done
