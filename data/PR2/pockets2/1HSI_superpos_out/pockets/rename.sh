#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/1HSI_pocket}"
done
