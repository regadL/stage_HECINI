#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/3EC0_pocket}"
done
