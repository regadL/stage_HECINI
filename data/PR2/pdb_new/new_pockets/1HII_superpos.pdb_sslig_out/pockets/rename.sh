#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/1HII_pocket}"
done
