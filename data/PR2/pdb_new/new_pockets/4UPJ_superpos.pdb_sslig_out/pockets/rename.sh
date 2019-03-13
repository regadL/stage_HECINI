#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/4UPJ_pocket}"
done
