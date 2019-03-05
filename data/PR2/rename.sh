#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/_superpos/}"
done
