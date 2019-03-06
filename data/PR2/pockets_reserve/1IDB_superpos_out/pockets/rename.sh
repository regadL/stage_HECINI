#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/1IDB_pocket}"
done
