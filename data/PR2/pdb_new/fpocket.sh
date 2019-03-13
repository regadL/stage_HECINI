#!/bin/bash

for file in *.pdb
do
	fpocket -f $file

done
