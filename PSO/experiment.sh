#!/bin/bash

# 定义要运行的测试函数
functions=("Ackley" "Rastrigin" "HappyCat" "Rosenbrock" "Zakharov" "Michalewicz" "Schwefel" "BentCigar" "DropWave" "Step")

# 定义维度和粒子数量
dimensions=("30")
numParticles=("50")
runs=1 
k="1"
c1="1.5"
c2="1.5"

# 循环遍历每个测试函数
for funct in "${functions[@]}"; do
    # 循环遍历维度
    for dim in "${dimensions[@]}"; do
        # 循环遍历粒子数量
        for numP in "${numParticles[@]}"; do
            # 运行命令
            ./run.sh "$funct" "$runs" "$dim" "$k" "$c1" "$c2" "$numP"
        done
    done
done
