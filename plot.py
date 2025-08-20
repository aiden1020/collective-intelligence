import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def read_fitness_data(file_path):
    with open(file_path, 'r') as file:
        fitness = [float(line.strip()) for line in file]
    return fitness


_start = 1000


def plot_fitness(file1, file2, file_name, label1, label2):
    fitness1 = read_fitness_data(file1)[_start:]  # 從指定的 iteration 開始
    fitness2 = read_fitness_data(file2)[_start:]  # 從指定的 iteration 開始
    x = list(range(_start, _start + len(fitness1)))  # x 軸從指定的 iteration 開始

    plt.figure()
    plt.plot(x, fitness1, label=label1)
    plt.plot(x, fitness2, label=label2)
    plt.xlabel('Iteration Time')
    plt.ylabel('Fitness')
    plt.title(f'{file_name[:-4]} Convergence')
    plt.legend()
    plt.xlim(_start, _start + len(fitness1) - 1)  # 設置 x 軸範圍

    # 設置 y 軸為科學記數法和小數點後 3 位
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.savefig(f'plot/{file_name[:-4]}_convergence.png')
    plt.close()


def main(folder_path1, folder_path2):
    files1 = [file_name for file_name in os.listdir(
        folder_path1) if file_name.endswith('.txt')]
    files1.sort()  # 按檔名排序

    files2 = [file_name for file_name in os.listdir(
        folder_path2) if file_name.endswith('.txt')]
    files2.sort()  # 按檔名排序

    common_files = set(files1).intersection(files2)

    for file_name in common_files:
        file1 = os.path.join(folder_path1, file_name)
        file2 = os.path.join(folder_path2, file_name)
        plot_fitness(file1, file2, file_name, 'PSO_improve', 'PSO')


if __name__ == "__main__":
    folder_path1 = 'PSO_improve/result/improve_convergence'
    folder_path2 = 'PSO/result/convergence'
    main(folder_path1, folder_path2)
