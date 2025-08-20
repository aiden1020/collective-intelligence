import os
import csv
from statistics import mean


def calculate_statistics(file_path):
    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file]
    avg = format(mean(numbers), '.3e')  # 科學計數法，小數點後 3 位
    worst = format(min(numbers), '.3e')
    best = format(max(numbers), '.3e')
    return avg, worst, best


def main(folder_path, output_csv):
    files = [file_name for file_name in os.listdir(
        folder_path) if file_name.endswith('.txt')]
    files.sort()  # 按檔名排序

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['file_name', 'avg', 'worst', 'best']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file_name in files:
            file_path = os.path.join(folder_path, file_name)
            avg, best, worst = calculate_statistics(file_path)
            writer.writerow({'file_name': file_name[:-4], 'avg': avg,
                            'worst': worst, 'best': best})


if __name__ == "__main__":
    folder_path = 'PSO/result/50P'
    output_csv = 'PSO/statistics.csv'
    main(folder_path, output_csv)
