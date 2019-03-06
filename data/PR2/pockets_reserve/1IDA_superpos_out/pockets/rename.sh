#!/bin/bash

for file in *.pdb
do
    mv "$file" "${file/pocket/1IDA_pocket}"
done
