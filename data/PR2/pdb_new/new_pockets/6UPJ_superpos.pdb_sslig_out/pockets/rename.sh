#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/6UPJ_pocket}"
done
