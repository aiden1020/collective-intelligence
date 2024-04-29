#!/bin/bash

# 定义要运行的测试函数
functions=("Ackley" "Rastrigin" "HappyCat" "Rosenbrock" "Zakharov" "Michalewicz")

# 定义每个函数对应的 k 值
k_values=("0.00007" "0.0002" "0.00005" "0.0003" "0.00007" "0.0000005")

# 定义维度和粒子数量
dimensions=("2" "10" "30")
numParticles=( "200")
runs=30

# 循环遍历每个测试函数
for ((i=0; i<${#functions[@]}; i++)); do
    # 获取当前函数对应的 k 值
    k="${k_values[i]}"
    funct="${functions[$i]}"
        # 循环遍历维度
    for dim in "${dimensions[@]}"; do
        # 循环遍历粒子数量

        for numP in "${numParticles[@]}"; do
            # 运行命令
            ./run.sh "$funct" "30" "$dim" "$k" "1.5" "1.5" "$numP"

        done
    done
done
