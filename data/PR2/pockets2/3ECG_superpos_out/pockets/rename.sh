#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/3ECG_pocket}"
done
