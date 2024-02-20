#!/bin/bash

# Assuming the commands are in the same directory as the script
# If not, provide the full paths to the commands

# Run ./slaMEM command
./slaMEM -l 10 Primates_ref.fasta Primates_query.fasta

# Run ./csaCut command
./csaCut

# Run python script
python make_circular.py