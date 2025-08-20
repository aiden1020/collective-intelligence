#!/bin/bash

functions=("Ackley" "Rastrigin" "HappyCat" "Rosenbrock" "Zakharov" "Michalewicz" "Schwefel" "BentCigar" "DropWave" "Step")

dimensions=("30")
numParticles=("50")
runs=1
k="1"
c1="1.5"
c2="1.5"

for funct in "${functions[@]}"; do
    for dim in "${dimensions[@]}"; do
        for numP in "${numParticles[@]}"; do
            ./run.sh "$funct" "$runs" "$dim" "$k" "$c1" "$c2" "$numP"
        done
    done
done
