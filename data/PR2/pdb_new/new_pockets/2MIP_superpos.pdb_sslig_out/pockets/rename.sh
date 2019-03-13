#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/2MIP_pocket}"
done
