#!/bin/bash

# Set the base directory to the directory where the script is located
cd "$(dirname "$0")"

# Clean log and results
rm -f output.log
rm -f core.*
