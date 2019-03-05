#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/1IVQ_pocket}"
done
