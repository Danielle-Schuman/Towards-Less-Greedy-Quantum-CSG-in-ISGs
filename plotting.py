import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import statistics
from tabulate import tabulate


def read_pickle_data(file_path, graph_num_per_size=20):
    all_data = []
    try:
        with open(file_path, 'rb') as file:
            if "total" in file_path:
                all_data = pickle.load(file)
            else:
                while True:
                    group_graph_size = []
                    for _ in range(graph_num_per_size):
                        try:
                            data = pickle.load(file)
                            group_graph_size.append(data)
                        except EOFError:
                            # Reached end of file
                            break
                    if group_graph_size:
                        all_data.append(group_graph_size)
                    else:
                        break
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
    except Exception as e:
        print(f"An error occurred while reading {file_path}: {e}")
    return all_data


def process_folder(folder_path):
    folder_contents = {}
    try:
        for root, dirs, files in os.walk(folder_path):
            for file_name in files:
                if file_name.startswith("data_") and file_name.endswith(".pkl"):
                    # Extracting key from file name
                    parts = file_name.split("_")
                    if "total.pkl" in parts:
                        key = ("_".join(parts[1:-4]), int(parts[-3]))  # Extracting algorithm name and seed for the key
                    else:
                        key = ("_".join(parts[1:-5]), int(parts[-4]))  # Extracting algorithm name and seed for the key
                    # Load data from file
                    file_path = os.path.join(root, file_name)
                    data = read_pickle_data(file_path)

                    # Store data in the dictionary
                    folder_contents[key] = data

    except FileNotFoundError:
        print(f"The folder {folder_path} does not exist.")
    except Exception as e:
        print(f"An error occurred while processing {folder_path}: {e}")
    return folder_contents


def calculate_statistics_over_seeds(folder_dict):
    dict_avgs = {}
    dict_stds = {}
    dict_num_succesful_seeds = {}

    # Group values by algorithm_name
    algorithm_name_groups = {}
    for key, values in folder_dict.items():
        algorithm_name, seed = key
        if algorithm_name not in algorithm_name_groups:
            algorithm_name_groups[algorithm_name] = {}
        algorithm_name_groups[algorithm_name][seed] = values

    # Calculate averages and standard deviations for each algorithm_name group
    for algorithm_name, seed_groups in algorithm_name_groups.items():
        num_graph_sizes = max(len(lst) for lst in seed_groups.values())
        averages = [[] for _ in range(num_graph_sizes)]
        stds = [[] for _ in range(num_graph_sizes)]
        # not every graph could successfully be processed in every run, i.e. for every seed
        # -> we need to count for how many seeds it was "successful" to create average
        num_seeds = [[] for _ in range(num_graph_sizes)]

        for seed, graph_sizes_list in seed_groups.items():
            for i, graph_list in enumerate(graph_sizes_list):
                for j, item in enumerate(graph_list):
                    if len(averages[i]) <= j:
                        averages[i].append((item[1], item[2]))
                        stds[i].append(([item[1]], [item[2]]))
                        num_seeds[i].append(1)
                    else:
                        averages[i][j] = (averages[i][j][0] + item[1], averages[i][j][1] + item[2])
                        num_seeds[i][j] = num_seeds[i][j] + 1
                        stds[i][j][0].append(item[1])
                        stds[i][j][1].append(item[2])

        for i, graph_list in enumerate(averages):
            for j, item in enumerate(graph_list):
                averages[i][j] = (averages[i][j][0] / num_seeds[i][j], averages[i][j][1] / num_seeds[i][j])

        for i, graph_list in enumerate(stds):
            for j, seeds_list in enumerate(graph_list):
                if len(seeds_list[0]) > 1:
                    to_evaluate = [float(x) for x in seeds_list[0]]
                    stds_values = statistics.stdev(to_evaluate)
                else:
                    stds_values = 0.0
                if len(seeds_list[1]) > 1:
                    to_evaluate = [float(x) for x in seeds_list[1]]
                    stds_times = statistics.stdev(to_evaluate)
                else:
                    stds_times = 0.0
                stds[i][j] = (stds_values, stds_times)

        dict_avgs[algorithm_name] = averages
        dict_stds[algorithm_name] = stds
        dict_num_succesful_seeds[algorithm_name] = num_seeds

    return dict_avgs, dict_stds, dict_num_succesful_seeds

# TODO: Decide metric
def relative_solution_quality(solutions_list, optima_list):
    rel_solution_qualities = []
    for i, solution in enumerate(solutions_list):
        rel_solution = solution / optima_list[i]
        #or
        error = abs(solution - optima_list[i]) / optima_list[i]
        rel_solution = 1 - error
        rel_solution_qualities.append(rel_solution)
    return rel_solution_qualities


def plot_line_chart_with_stds(algorithm_names, averages_list, std_devs_list, colors, x_ticks, xlabel, ylabel, title):
    assert len(algorithm_names) == len(averages_list) == len(std_devs_list), "Input lists should have the same length."

    plt.figure(figsize=(10, 6))

    for i in range(len(algorithm_names)):
        plt.errorbar(x_ticks, averages_list[i], yerr=std_devs_list[i], fmt='o-', color=colors[i], label=algorithm_names[i], capsize=5)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend()
    #plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    # plt.show()
    plt.close()

def plot_line_chart(algorithm_names, averages_list, x_ticks, colors, xlabel, ylabel, title):
    assert len(algorithm_names) == len(averages_list), "Input lists should have the same length."

    plt.figure(figsize=(10, 6))

    for i in range(len(algorithm_names)):
        plt.plot(x_ticks, averages_list[i], marker='o', color=colors[i], label=algorithm_names[i])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend()
    #plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    # plt.show()
    plt.close()

def plot_vertical_lines_chart(algorithm_names, values, x_ticks, colors, xlabel, ylabel, title):
    assert len(algorithm_names) == len(values), "Input lists should have the same length."
    line_distance = 0.02
    plt.figure(figsize=(10, 6))

    # Calculate the width of each line
    line_distance_total = len(algorithm_names) * line_distance
    line_distance_shift = np.linspace(-line_distance_total / 2, line_distance_total / 2, len(algorithm_names))

    for i in range(len(algorithm_names)):
        assert len(values[i]) == len(
            x_ticks), f"Values for category {algorithm_names[i]} should have the same length as x_ticks."
        plt.vlines(x_ticks + line_distance_shift[i], 0, values[i], linestyle='-', linewidth=3, color=colors[i], label=algorithm_names[i])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend(loc='upper left')
    #plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    #plt.show()
    plt.close()

def generate_latex_table(algorithm_names, values_list, std_devs_list, row_names):
    # Ensure input lists have the same length
    assert len(algorithm_names) == len(values_list) == len(std_devs_list), "Input lists should have the same length."

    # Initialize the table data
    table_data = []

    # Add headers as the first row
    headers = [""] + algorithm_names
    table_data.append(headers)

    # Iterate through row names and add corresponding data
    for i, row_name in enumerate(row_names):
        row = [row_name]
        for j in range(len(algorithm_names)):
            value = f"{values_list[j][i]} $\\pm$ {std_devs_list[j][i]}"
            row.append(value)
        table_data.append(row)

    # Generate LaTeX table
    latex_table = tabulate(table_data, headers="firstrow", tablefmt="latex_raw")

    return latex_table


def generate_latex_table_double(algorithm_names, value_sums, stds_value_sums, rel_values, std_rel_values, row_names):
    # Ensure input lists have the same length
    assert len(algorithm_names) == len(value_sums) == len(stds_value_sums), "Input lists should have the same length."
    assert len(algorithm_names) == len(rel_values) == len(std_rel_values), "Input lists should have the same length."

    # Initialize the table data
    table_data = []

    # Add header for algorithm names
    header = [""]
    for name in algorithm_names:
        header.append(f"\\multicolumn{{2}}{{c}}{{{name}}}")
    table_data.append(header)

    # Add subheader for averages and standard deviations
    subheader = [""]
    for _ in algorithm_names:
        subheader.extend(["V $\pm$ std", "R $\pm$ std"])
    table_data.append(subheader)

    # Add data rows
    for i, row_name in enumerate(row_names):
        row = [row_name]
        for j in range(len(algorithm_names)):
            value1 = f"{value_sums[j][i]} $\pm$ {stds_value_sums[j][i]}"
            value2 = f"{rel_values[j][i]} $\pm$ {std_rel_values[j][i]}"
            row.extend([value1, value2])
        table_data.append(row)

    # Generate LaTeX table
    latex_table = tabulate(table_data, headers="firstrow", tablefmt="latex_raw")

    # Header lines need to be edited by hand a bit afterwards (remove unnecessary columns in header, add \hline below subheaders)
    return latex_table


def sum_over_same_sized_graphs(dict_avg):
    dict_sums = {}

    for algorithm_name, graph_size_lists in dict_avg.items():
        summed_graph_size_lists = []
        for graph_list in graph_size_lists:
            value_sum = sum(item[0] for item in graph_list)
            time_sum = sum(item[1] for item in graph_list)
            summed_graph_size_lists.append((value_sum, time_sum))
        dict_sums[algorithm_name] = summed_graph_size_lists

    return dict_sums


def filter_dict_by_substring(original_dict, substring):
    filtered_dict = {}

    for key, value in original_dict.items():
        if substring in key:
            filtered_dict[key] = value

    return filtered_dict


def remove_substring_from_keys(input_dict, substring):
    new_dict = {}

    for key, value in input_dict.items():
        new_key = key.replace(substring, "")
        new_dict[new_key] = value

    return new_dict


def split_dicts_by_graphsizes(original_dict, graph_sizes):
    graphsize_dicts = {graph_size: {} for graph_size in graph_sizes}

    for algorithm_name, all_results_for_algo in original_dict.items():
        start_name = algorithm_name.split("_", 1)[0]
        if start_name.isnumeric() and start_name != '2':
            k = int(start_name)
            index_k = graph_sizes.index(k) if k in graph_sizes else graph_sizes.index(k + 1)
            # Insert measured values and times into graphsize_dicts
            for i, graph_size in enumerate(graph_sizes[index_k:]):
                graphsize_dicts[graph_size][algorithm_name] = all_results_for_algo[i]
        elif start_name == 'GCS-Q' or start_name == '2':
            # k = 2
            index_k = 0
            # Insert measured values and times into graphsize_dicts
            for i, graph_size in enumerate(graph_sizes[index_k:]):
                graphsize_dicts[graph_size][algorithm_name] = all_results_for_algo[i]
        elif start_name == 'n' or "ours" or "R-QUBO":
            # Insert measured values and times into graphsize_dicts
            for i, graph_size in enumerate(graph_sizes):
                graphsize_dicts[graph_size][algorithm_name] = all_results_for_algo[i]
        else:
            raise Exception("k for split not known")

    return graphsize_dicts


def compare_ks_GCSQ_exactly(dict_sums):
    # get results from k=2 (GCS-Q) to k=n (n_split), and combine them into one dict
    dict_GCSQ = filter_dict_by_substring(dict_sums, "GCS-Q_")
    dict_k_split_GCSQ = filter_dict_by_substring(dict_sums, "_split_GCSQ_exactly")
    dict_n_split_GCSQ = filter_dict_by_substring(dict_sums, "n_split_GCSQ")
    dict_all = {**dict_GCSQ, **dict_k_split_GCSQ, **dict_n_split_GCSQ}
    plot_values_and_times(dict_all, "GCS-Q exactly")


def compare_ks_GCSQ_at_most(dict_sums):
    # get results from k=2 (GCS-Q) to k=n-1, and combine them into one dict
    dict_GCSQ = filter_dict_by_substring(dict_sums, "GCS-Q_")
    dict_k_split_GCSQ = filter_dict_by_substring(dict_sums, "_split_GCSQ_at_most")
    dict_all = {**dict_GCSQ, **dict_k_split_GCSQ}
    plot_values_and_times(dict_all, "GCS-Q at most")


def compare_ours_iterative_exactly(dict_sums):
    # get results from k=2 to k=n (ours_n), and combine them into one dict
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_ours_iterative_exactly")
    dict_n_split = filter_dict_by_substring(dict_sums, "ours_n_q")  # added the q to exclude the ours_n_half
    dict_all = {**dict_k_split, **dict_n_split}
    plot_values_and_times(dict_all, "Ours exactly")


def compare_ours_iterative_at_most(dict_sums):
    # get results from k=2 to k=n (ours_n_half), and combine them into one dict
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_ours_iterative_at_most")
    dict_n_split = filter_dict_by_substring(dict_sums, "ours_n_half")
    dict_all = {**dict_k_split, **dict_n_split}
    plot_values_and_times(dict_all, "Ours at most")


def compare_r_qubo_iterative(dict_sums):
    # get results from k=2 to k=n (r_qubo), and combine them into one dict
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_R_QUBO-iterative")
    dict_n_split = filter_dict_by_substring(dict_sums, "R-QUBO")
    dict_all = {**dict_k_split, **dict_n_split}
    plot_values_and_times(dict_all, "R-QUBO")


def get_labels_x_axis(input_dict, title):
    if "GCS-Q exactly" in title:
        return ['GCS-Q'] + sorted(int(key) for key in input_dict.keys() if key.isdigit()) + ['n']
    elif "GCS-Q at most" in title:
        return ['GCS-Q'] + sorted(int(key) for key in input_dict.keys() if key.isdigit())
    elif "Ours exactly" in title:
        return sorted(int(key) for key in input_dict.keys() if key.isdigit()) + ['ours_n']  # TODO: Or Glover?
    elif "Ours at most" in title:
        # TODO: Put ours_n_half in middle, next to k = n/2
        return sorted(int(key) for key in input_dict.keys() if key.isdigit()) + ['ours_n_half']  # TODO: Or just ours?
    elif "R-QUBO" in title:
        return sorted(int(key) for key in input_dict.keys() if key.isdigit()) + ['R-QUBO']
    else:
        raise Exception("Algorithm to plot unknown -> Don't know how to label x-axis")

def plot_values_and_times(dict_all, title):
    # split dict for different solvers
    dict_qbsolv = remove_substring_from_keys(filter_dict_by_substring(dict_all, "_qbsolv"), "_qbsolv")
    dict_qaoa = remove_substring_from_keys(filter_dict_by_substring(dict_all, "_qaoa"), "_qaoa")
    # split dicts according to graph sizes
    graph_sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
    dicts_qbsolv = split_dicts_by_graphsizes(dict_qbsolv, graph_sizes)
    dicts_qaoa = split_dicts_by_graphsizes(dict_qaoa, graph_sizes)
    # plot value-sums and time-sums for every graph size
    for graph_size in graph_sizes:
        # split dicts into coalition values and times
        dict_qbsolv_values = {key: value[0] for key, value in dicts_qbsolv[graph_size].items()}
        dict_qbsolv_times = {key: value[1] for key, value in dicts_qbsolv[graph_size].items()}
        dict_qaoa_values = {key: value[0] for key, value in dicts_qaoa[graph_size].items()}
        dict_qaoa_times = {key: value[1] for key, value in dicts_qaoa[graph_size].items()}
        plot_bar_chart(dict_qbsolv_values, dict_qaoa_values, "Sums of coalition values", f"Sums of coalition values of 20 graphs of size {graph_size} for {title}", "k for k-split")
        plot_bar_chart(dict_qbsolv_times, dict_qaoa_times, "Sums of times in seconds",
                   f"Time taken for 20 graphs of size {graph_size} for {title}", "k for k-split")


def split_keys_at_underscore(original_dict):
    new_dict = {}

    for key, value in original_dict.items():
        if "split" in key:
            new_key = key.split('_', 1)[0]  # Split at the first occurrence of '_', before "split"
        elif "parallel" in key:
            new_key = key.split('_p', 1)[0]  # Split of the '_parallel'
        else:
            new_key = key
        new_dict[new_key] = value

    return new_dict


def plot_bar_chart(dict_qbsolv, dict_qaoa, y_label, title, label_x_axis):
    dict_qbsolv_short_keys = split_keys_at_underscore(dict_qbsolv)
    dict_qaoa_short_keys = split_keys_at_underscore(dict_qaoa)
    labels_x_axis = get_labels_x_axis(dict_qbsolv_short_keys, title)
    values_qbsolv = [dict_qbsolv_short_keys[str(key)] for key in labels_x_axis]
    #values_qaoa = [dict_qaoa_short_keys[str(key)] for key in labels_x_axis]

    x = np.arange(len(labels_x_axis))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))  # Adjust the figure size as needed

    bars_qbsolv = ax.bar(x - width/2, values_qbsolv, width, label="qbsolv")
    #bars_qaoa = ax.bar(x + width/2, values_qaoa, width, label="qaoa")

    ax.set_xlabel(label_x_axis)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_x_axis, rotation=45, ha='right')  # Rotate x-axis labels by 45 degrees
    ax.legend()

    # Calculate the maximum value of the data and adjust the y-axis limit
    max_value = max(values_qbsolv)
    y_limit = max_value * 1.25  # Set the y-axis limit to 25% higher than the maximum value

    ax.set_ylim(0, y_limit)  # Set the y-axis limit

    # Annotate values on top of bars
    for bar in bars_qbsolv:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',  # Format value rounded to three decimal places
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', rotation=90, fontsize=8)  # Rotate and reduce font size

    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    #plt.show()
    plt.close()


# TODO: Update this method so that it works
def get_barchart_win_lose_draw(algo1_name, algo2_name, results_dict, seed, note=""):
    win_distribution = {
        algo1_name: [],
        "Draw": [],
        algo2_name: [],
    }
    for test in results_dict:
        w1, d, w2 = 0, 0, 0
        for v1, v2 in zip(results_dict[test][algo1_name], results_dict[test][algo2_name]):
            if v1 == v2:
                d += 1
            elif v1 > v2:
                w1 += 1
            else:
                w2 += 1
        win_distribution[algo1_name].append(w1)
        win_distribution["Draw"].append(d)
        win_distribution[algo2_name].append(w2)
    width = 0.5
    fig, ax = plt.subplots()
    bottom = np.zeros(len(results_dict))

    species = tuple(["num_agents=" + str(i) for i in results_dict])
    for boolean, weight_count in win_distribution.items():
        p = ax.bar(species, weight_count, width, label=boolean, bottom=bottom)
        bottom += weight_count

    ax.legend(loc="upper right")

    if note:
        plotname = f"barchart_win_lose_draw_{algo1_name}_{algo2_name}_{note}_seed_{seed}"
    else:
        plotname = f"barchart_win_lose_draw_{algo1_name}_{algo2_name}_seed_{seed}"

    plt.savefig(f"plots/{plotname}.pdf", format='pdf')
    #plt.show()
    plt.close()

if __name__ == "__main__":
     data_path = "results/eon_data/classical"
     #data_path = "results/eon_data/quantum/qbsolv/parallel/k=12/data_12_split_GCSQ_exactly_qbsolv_parallel__3949468976__2024-01-26_13-12-53.747804.pkl"
     data = process_folder(data_path)
     #data = read_pickle_data(data_path)
     #print(data)

     data_avgs, data_stds, data_num_successful_seeds = calculate_statistics_over_seeds(data)
     print(data_avgs, "\n", data_stds)
     avg_data_summed = sum_over_same_sized_graphs(data_avgs)
     stds_summed = sum_over_same_sized_graphs(data_stds)
     sorted_avg_summed = sorted(avg_data_summed.items())
     print(sorted_avg_summed)
     print("\n")
     sorted_stds_summed = sorted(stds_summed.items())
     print(sorted_stds_summed)
     sorted_num_seeds = sorted(data_num_successful_seeds.items())
     print(sorted_num_seeds)
     '''
     compare_ks_GCSQ_exactly(avg_data_summed)
     compare_ks_GCSQ_at_most(avg_data_summed)
     compare_ours_iterative_exactly(avg_data_summed)
     compare_ours_iterative_at_most(avg_data_summed)
     compare_r_qubo_iterative(avg_data_summed)

    # Example usage:
    algorithm_names = ['Algorithm A', 'Algorithm B', 'Algorithm C']
    averages_list1 = [[10, 20, 30], [15, 25, 35], [12, 22, 32]]
    std_devs_list1 = [[1, 2, 3], [1.5, 2.5, 3.5], [1.2, 2.2, 3.2]]
    averages_list2 = [[11, 21, 31], [16, 26, 36], [13, 23, 33]]
    std_devs_list2 = [[1.1, 2.1, 3.1], [1.6, 2.6, 3.6], [1.3, 2.3, 3.3]]
    row_names = [4, 6, 28]

    latex_table = generate_latex_table_double(algorithm_names, averages_list1, std_devs_list1, averages_list2, std_devs_list2, row_names)
    print(latex_table)
    '''

