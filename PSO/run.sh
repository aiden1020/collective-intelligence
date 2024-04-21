#!/bin/bash

# Clean previous build
make clean

# Build the project
make

# Execute the program with parameters
./PSO "$@"
