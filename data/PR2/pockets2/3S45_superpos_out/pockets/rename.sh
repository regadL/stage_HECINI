#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/3S45_pocket}"
done
