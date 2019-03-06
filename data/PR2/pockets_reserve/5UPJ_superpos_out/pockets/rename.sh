#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/5UPJ_pocket}"
done
