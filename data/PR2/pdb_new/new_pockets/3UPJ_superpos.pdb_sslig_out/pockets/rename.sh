#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/3UPJ_pocket}"
done
