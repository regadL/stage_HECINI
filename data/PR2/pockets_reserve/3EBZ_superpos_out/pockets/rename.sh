#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/3EBZ_pocket}"
done
