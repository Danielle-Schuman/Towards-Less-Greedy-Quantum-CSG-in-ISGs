import matplotlib.pyplot as plt
import numpy as np
import bisect
import numbers
import math
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


def coalitions_feasible(coalition_structure):
    agents_used = []
    coalition_structure_feasible = True
    for coalition in coalition_structure:
        for agent in coalition:
            if agent in agents_used:
                coalition_structure_feasible = False
            else:
                agents_used.append(agent)
    return coalition_structure_feasible


def filter_infeasible_coalitions(folder_dict):
    for key, graph_size_list in folder_dict.items():
        for i, graph_list in enumerate(graph_size_list):
            for j, tuple in enumerate(graph_list):
                if not coalitions_feasible(tuple[0]):
                    folder_dict[key][i][j] = (tuple[0], "infeasible", tuple[2])
    return folder_dict


def calculate_statistics_over_seeds(folder_dict):
    dict_avgs_values = {}
    dict_stds_values = {}
    dict_avgs_times = {}
    dict_stds_times = {}
    dict_num_succesful_seeds = {}
    dict_num_feasible_seeds = {}

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
        averages_values = [[] for _ in range(num_graph_sizes)]
        stds_values = [[] for _ in range(num_graph_sizes)]
        averages_times = [[] for _ in range(num_graph_sizes)]
        stds_times = [[] for _ in range(num_graph_sizes)]
        # not every graph could successfully be processed in every run, i.e. for every seed
        # -> we need to count for how many seeds it was "successful" to create average
        num_seeds = [[] for _ in range(num_graph_sizes)]
        num_seeds_feasible = [[] for _ in range(num_graph_sizes)]

        for seed, graph_sizes_list in seed_groups.items():
            for i, graph_list in enumerate(graph_sizes_list):
                for j, item in enumerate(graph_list):
                    if len(averages_values[i]) <= j:
                        averages_values[i].append(item[1])
                        stds_values[i].append([item[1]])
                        # if solution is feasible
                        if isinstance(item[1], (numbers.Number, np.number)):
                            num_seeds_feasible[i].append(1)
                        else:
                            num_seeds_feasible[i].append(0)
                        averages_times[i].append(item[2])
                        stds_times[i].append([item[2]])
                        num_seeds[i].append(1)
                    else:
                        # if solution is feasible
                        if isinstance(item[1], (numbers.Number, np.number)):
                            if isinstance(averages_values[i][j], (numbers.Number, np.number)):
                                averages_values[i][j] = averages_values[i][j] + item[1]
                                stds_values[i][j].append(item[1])
                            else:
                                # current average value / last value in stds says "infeasible" -> needs to be replaced
                                averages_values[i][j] = item[1]
                                stds_values[i][j][-1] = item[1]
                            num_seeds_feasible[i][j] = num_seeds_feasible[i][j] + 1
                        averages_times[i][j] = averages_times[i][j] + item[2]
                        stds_times[i][j].append(item[2])
                        num_seeds[i][j] = num_seeds[i][j] + 1

        for i, graph_list in enumerate(averages_values):
            for j, item in enumerate(graph_list):
                if isinstance(averages_values[i][j], (numbers.Number, np.number)):
                    averages_values[i][j] = averages_values[i][j] / num_seeds_feasible[i][j]
                averages_times[i][j] = averages_times[i][j] / num_seeds[i][j]

        for i, graph_list in enumerate(stds_values):
            for j, seeds_list in enumerate(graph_list):
                if len(seeds_list) > 1:
                    to_evaluate = [float(x) for x in seeds_list]
                    stds_value = statistics.stdev(to_evaluate)
                else:
                    # for only infeasible solutions, we will also use 0 as std
                    stds_value = 0.0
                stds_values[i][j] = stds_value

        for i, graph_list in enumerate(stds_times):
            for j, seeds_list in enumerate(graph_list):
                if len(seeds_list) > 1:
                    to_evaluate = [float(x) for x in seeds_list]
                    stds_time = statistics.stdev(to_evaluate)
                else:
                    stds_time = 0.0
                stds_times[i][j] = stds_time

        dict_avgs_values[algorithm_name] = averages_values
        dict_stds_values[algorithm_name] = stds_values
        dict_num_feasible_seeds[algorithm_name] = num_seeds_feasible
        dict_avgs_times[algorithm_name] = averages_times
        dict_stds_times[algorithm_name] = stds_times
        dict_num_succesful_seeds[algorithm_name] = num_seeds

    return dict_avgs_values, dict_stds_values, dict_num_feasible_seeds, dict_avgs_times, dict_stds_times, dict_num_succesful_seeds


def relative_solution_quality(solutions_dict, optima_list):
    rel_values_dict = {}
    for algorithm_name, solutions_list in solutions_dict.items():
        rel_solution_qualities = []
        # check if is k-split with k > 4 -> solution_list starts at later graph size
        parts = algorithm_name.split("_")
        if parts[0].isnumeric():
            k = int(parts[0])
            graph_sizes = [4,6,8,10,12,14,16,18,20,22,24,26,28]
            start = bisect.bisect_left(graph_sizes, k)
            optima_list_relevant = optima_list[start:]
        else:
            optima_list_relevant = optima_list
        # shorten in case algorithm couldn't solve the large sizes
        num_graph_sizes = len(solutions_list)
        optima_list_relevant = optima_list_relevant[:num_graph_sizes]
        for i, solution_sublist in enumerate(solutions_list):
            rel_solution_qualities.append([])
            for j, solution in enumerate(solution_sublist):
                if solution != "infeasible":
                    rel_solution = round(solution, 6) / optima_list_relevant[i][j]
                    #or
                    #error = abs(round(solution, 6) - optima_list_relevant[i][j]) / optima_list_relevant[i][j]
                    #rel_solution = 1 - error
                else:
                    # We assign a relative solution quality of 0 to infeasible solutions
                    rel_solution = 0.0
                rel_solution_qualities[i].append(rel_solution)
        rel_values_dict[algorithm_name] = rel_solution_qualities
    return rel_values_dict


def optimum_found(solutions_dict, optima_list):
    opt_found_dict = {}
    for (algorithm_name, seed), solutions_list in solutions_dict.items():
        opt_found = []
        # check if is k-split with k > 4 -> solution_list starts at later graph size
        parts = algorithm_name.split("_")
        if parts[0].isnumeric():
            k = int(parts[0])
            graph_sizes = [4,6,8,10,12,14,16,18,20,22,24,26,28]
            start = bisect.bisect_left(graph_sizes, k)
            optima_list_relevant = optima_list[start:]
        else:
            optima_list_relevant = optima_list
        # shorten in case algorithm couldn't solve the large sizes
        num_graph_sizes = len(solutions_list)
        optima_list_relevant = optima_list_relevant[:num_graph_sizes]
        for i, solution_sublist in enumerate(solutions_list):
            opt_found.append([])
            for j, (coalitions, solution, time) in enumerate(solution_sublist):
                if solution != "infeasible":
                    if round(solution, 6) == optima_list_relevant[i][j]:
                        found = 1
                    else:
                        found = 0
                else:
                    # We assign a found value 0 to infeasible solutions
                    found = 0
                opt_found[i].append((coalitions, found, time))
        opt_found_dict[(algorithm_name, seed)] = opt_found
    return opt_found_dict


def plot_line_chart_with_stds(algorithm_names, values, std_devs, colors, x_ticks, xlabel, ylabel, title):
    assert len(algorithm_names) == len(values) == len(std_devs), "Input lists should have the same length."

    plt.figure(figsize=(10, 6))

    # Calculate the width of each line
    line_distance = 0.05
    line_distance_total = len(algorithm_names) * line_distance
    line_distance_shift = np.linspace(-line_distance_total / 2, line_distance_total / 2, len(algorithm_names))

    for i in range(len(algorithm_names)):
        plt.errorbar(x_ticks + line_distance_shift[i], values[i], yerr=std_devs[i], fmt='o-', color=colors[i], label=algorithm_names[i], capsize=5)

    if "Time" in title and "qaoa" in title:
        plt.yscale('log')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    # plt.show()
    plt.close()


def plot_line_chart(algorithm_names, values, x_ticks, colors, xlabel, ylabel, title, std_devs=None):
    assert len(algorithm_names) == len(values), "Input lists should have the same length."

    plt.figure(figsize=(10, 6))

    for i in range(len(algorithm_names)):
        plt.plot(x_ticks, values[i], marker='o', color=colors[i], label=algorithm_names[i])

    if "Time" in title and "qaoa" in title:
        plt.yscale('log')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    # plt.show()
    plt.close()


def plot_vertical_lines_chart(algorithm_names, values, x_ticks, colors, xlabel, ylabel, title, std_devs=None):
    assert len(algorithm_names) == len(values), "Input lists should have the same length."
    line_distance = 0.1
    plt.figure(figsize=(10, 6))

    # Calculate the width of each line
    line_distance_total = len(algorithm_names) * line_distance
    line_distance_shift = np.linspace(-line_distance_total / 2, line_distance_total / 2, len(algorithm_names))

    for i in range(len(algorithm_names)):
        assert len(values[i]) == len(
            x_ticks), f"Values for category {algorithm_names[i]} should have the same length as x_ticks."
        plt.vlines(x_ticks + line_distance_shift[i], 0, values[i], linestyle='-', linewidth=1.5, color=colors[i], label=algorithm_names[i])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    #plt.show()
    plt.close()


def plot_vertical_lines_chart_with_stds(algorithm_names, values, std_devs, x_ticks, colors, xlabel, ylabel, title):
    assert len(algorithm_names) == len(values), "Input lists should have the same length."
    line_distance = 0.2
    plt.figure(figsize=(10, 6))

    # Calculate the width of each line
    line_distance_total = len(algorithm_names) * line_distance
    line_distance_shift = np.linspace(-line_distance_total / 2, line_distance_total / 2, len(algorithm_names))

    for i in range(len(algorithm_names)):
        assert len(values[i]) == len(
            x_ticks), f"Values for category {algorithm_names[i]} should have the same length as x_ticks."
        plt.vlines(x_ticks + line_distance_shift[i], 0, values[i], linestyle='-', linewidth=1, color=colors[i], label=algorithm_names[i])
        # Add error bars for standard deviations
        plt.errorbar(x_ticks + line_distance_shift[i], values[i], yerr=std_devs[i], fmt='o', color='black')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x_ticks)
    plt.legend(loc='upper right')
    plt.grid(True)
    #plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    plt.show()
    #plt.close()


def get_max(values_list_padded, num_rows, num_columns):
    max = []
    for i in range(num_rows):
        max_row = float('-inf')
        for j in range(num_columns):
            if isinstance(values_list_padded[j][i], (numbers.Number, np.number)):
                if values_list_padded[j][i] > max_row:
                    max_row = values_list_padded[j][i]
        max.append(max_row)
    return max


# pad with -
def pad(values_list):
    values_list_padded = []
    for i, column in enumerate(values_list):
        if (i + 2) > 4:
            padded_column = (((math.ceil((i + 2) / 2)) - 2) * ["-"]) + column
        else:
            padded_column = column
        while len(padded_column) < 13:
            padded_column.append("-")
        values_list_padded.append(padded_column)
    return values_list_padded

def generate_latex_table(column_names, values_list, std_devs_list, row_names):
    # Ensure input lists have the same length
    assert len(row_names) == len(values_list) == len(std_devs_list), "Input lists should have the same length."

    # Initialize the table data
    table_data = []

    # Add headers as the first row
    headers = ["k"] + [f"n={c}" for c in column_names]
    table_data.append(headers)

    # pad with -
    values_list_padded = pad(values_list)
    std_devs_list_padded = pad(std_devs_list)

    # Find the maximum value in each column
    max_values = get_max(values_list_padded, len(column_names), len(row_names))

    # Iterate through row names and add corresponding data
    for i, row_name in enumerate(row_names):
        row = [row_name]
        for j in range(len(column_names)):
            if isinstance(values_list_padded[i][j], (numbers.Number, np.number)):
                value = f"{round(values_list_padded[i][j], 6)} $\\pm$ {round(std_devs_list_padded[i][j], 6)}"
                if values_list_padded[i][j] == max_values[j]:  # If cell value equals maximum in the column
                    value = "\\textbf{" + value + "}"  # Make it bold
            else:
                value = f"{values_list_padded[i][j]}"
            row.append(value)
        table_data.append(row)

    # Generate LaTeX table
    latex_table = tabulate(split_table_horizontally(table_data), headers="firstrow", tablefmt="latex_raw")

    return latex_table


def split_table_horizontally(table):
    # Calculate the index to split each row
    split_index = math.ceil(len(table[0]) / 2) + 1

    # Split each row horizontally
    left_half = []
    right_half = []
    for row in table:
        if row[1:split_index] != len(row[1:split_index]) * ["-"]:
            left_half.append(row[:split_index])
        else:
            pass
        right_half.append([row[0]] + row[split_index:])

    # Concatenate the right half to the bottom of the left half
    left_half.extend(right_half)

    return left_half


def generate_latex_table_horizontal(column_names, values_list, std_devs_list, row_names):
    # Ensure input lists have the same length
    assert len(column_names) == len(values_list) == len(std_devs_list), "Input lists should have the same length."

    # Initialize the table data
    table_data = []

    # Add headers as the first row
    headers = ["n"] + [f"k={c}" for c in column_names]
    table_data.append(headers)

    # pad with -
    values_list_padded = pad(values_list)
    std_devs_list_padded = pad(std_devs_list)

    # Find the maximum value in each row
    max_values = get_max(values_list_padded, len(row_names), len(column_names))

    # Iterate through row names and add corresponding data
    for i, row_name in enumerate(row_names):
        row = [row_name]
        for j in range(len(column_names)):
            if isinstance(values_list_padded[j][i], (numbers.Number, np.number)):
                value = f"{round(values_list_padded[j][i], 3)} $\\pm$ {round(std_devs_list_padded[j][i], 3)}"
                if values_list_padded[j][i] == max_values[i]:  # If cell value equals maximum in the row
                    value = "\\textbf{" + value + "}"  # Make it bold
            else:
                value = f"{values_list_padded[j][i]}"
            row.append(value)
        table_data.append(row)

    # Generate LaTeX table
    latex_table = tabulate(split_table_horizontally(table_data), headers="firstrow", tablefmt="latex_raw")

    return latex_table


def generate_latex_table_double(column_names, value_sums, stds_value_sums, rel_values, std_rel_values, row_names):
    # Ensure input lists have the same length
    assert len(column_names) == len(value_sums) == len(stds_value_sums), "Input lists should have the same length."
    assert len(column_names) == len(rel_values) == len(std_rel_values), "Input lists should have the same length."

    # Initialize the table data
    table_data = []

    # Add header for algorithm names
    header = [""]
    for name in column_names:
        header.append(f"\\multicolumn{{2}}{{c}}{{{name}}}")
    table_data.append(header)

    # Add subheader for averages and standard deviations
    subheader = [""]
    for _ in column_names:
        subheader.extend(["V $\\pm$ std", "R $\\pm$ std"])
    table_data.append(subheader)

    # Find the maximum value in each row for both V and R
    max_value_v = [max(row) for row in value_sums]
    max_value_r = [max(row) for row in rel_values]

    # Add data rows
    for i, row_name in enumerate(row_names):
        row = [row_name]
        for j in range(len(column_names)):
            value1 = f"{value_sums[j][i]} $\\pm$ {stds_value_sums[j][i]}"
            value2 = f"{rel_values[j][i]} $\\pm$ {std_rel_values[j][i]}"
            if value_sums[j][i] == max_value_v[j] or rel_values[j][i] == max_value_r[j]:  # If cell value equals maximum in the row
                value1 = "\\textbf{" + value1 + "}"  # Make it bold
                value2 = "\\textbf{" + value2 + "}"  # Make it bold
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
            feasible_graphlist = [item for item in graph_list if isinstance(item, (numbers.Number, np.number))]
            if len(feasible_graphlist) > 0:
                sum_graphs = sum(feasible_graphlist)
            else:
                sum_graphs = "infeasible"
            summed_graph_size_lists.append(sum_graphs)
        dict_sums[algorithm_name] = summed_graph_size_lists

    return dict_sums


def average_over_same_sized_graphs(dict_avg):
    dict_avgs = {}
    dict_stds = {}

    for algorithm_name, graph_size_lists in dict_avg.items():
        averaged_graph_size_lists = []
        stds_graph_size_lists = []
        for graph_list in graph_size_lists:
            feasible_graphlist = [item for item in graph_list if isinstance(item, (numbers.Number, np.number))]
            if len(feasible_graphlist) > 0:
                avg = sum(feasible_graphlist) / len(feasible_graphlist)
                if len(feasible_graphlist) > 1:
                    std = statistics.stdev(feasible_graphlist)
                else:
                    std = 0.0
            else:
                avg = "infeasible"
                std = 0.0
            averaged_graph_size_lists.append(avg)
            stds_graph_size_lists.append(std)
        dict_avgs[algorithm_name] = averaged_graph_size_lists
        dict_stds[algorithm_name] = stds_graph_size_lists

    return dict_avgs, dict_stds


def average_over_graph_sizes(dict_avg):
    dict_avgs = {}
    for algorithm_name, graph_size_lists in dict_avg.items():
        feasible_graph_size_list = [item for item in graph_size_lists if isinstance(item, (numbers.Number, np.number))]
        if len(feasible_graph_size_list) > 0:
            avg = np.mean(np.array(feasible_graph_size_list))
        else:
            # if all infeasible, we set rel. solution quality to zero
            avg = 0.0
        dict_avgs[algorithm_name] = avg
    return dict_avgs


def std_over_graph_sizes(dict_avg):
    dict_stds = {}
    for algorithm_name, graph_size_lists in dict_avg.items():
        feasible_graph_size_list = [item for item in graph_size_lists if isinstance(item, (numbers.Number, np.number))]
        if len(feasible_graph_size_list) > 0:
            std = np.std(np.array(feasible_graph_size_list))
        else:
            # if all infeasible, we set rel. solution quality to zero
            std = 0.0
        dict_stds[algorithm_name] = std
    return dict_stds


def unzip_dictionaries(zipped_dictionaries):
    values_dict = {}
    times_dict = {}

    for algorithm_name, graph_size_lists in zipped_dictionaries.items():
        # unzip
        if len(graph_size_lists[0]) == 2:
            value_list_total = [[tup[0] for tup in graph_list] for graph_list in graph_size_lists]
            times_list_total = [[tup[1] for tup in graph_list] for graph_list in graph_size_lists]
        elif len(graph_size_lists[0]) == 3:
            value_list_total = [[tup[1] for tup in graph_list] for graph_list in graph_size_lists]
            times_list_total = [[tup[2] for tup in graph_list] for graph_list in graph_size_lists]
        else:
            raise Exception("Behaviour of unzip undefined for this tuple length")
        values_dict[algorithm_name] = value_list_total
        times_dict[algorithm_name] = times_list_total

    return values_dict, times_dict


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
    # get results from k=2 (GCS-Q) to k=n-1, and combine them into one dict
    dict_GCSQ = filter_dict_by_substring(dict_sums, "GCS-Q_")
    dict_k_split_GCSQ = filter_dict_by_substring(dict_sums, "_split_GCSQ_exactly")
    #dict_n_split_GCSQ = filter_dict_by_substring(dict_sums, "n_split_GCSQ")
    dict_all = {**dict_GCSQ, **dict_k_split_GCSQ} #, **dict_n_split_GCSQ}
    #plot_values_and_times(dict_all, "GCS-Q exactly")
    list_all = sort_by_k(dict_all, "GCS-Q exactly")
    return list_all


def compare_ks_GCSQ_at_most(dict_sums):
    # get results from k=2 (GCS-Q) to k=n-1, and combine them into one dict
    dict_GCSQ = filter_dict_by_substring(dict_sums, "GCS-Q_")
    dict_k_split_GCSQ = filter_dict_by_substring(dict_sums, "_split_GCSQ_at_most")
    dict_all = {**dict_GCSQ, **dict_k_split_GCSQ}
    #plot_values_and_times(dict_all, "GCS-Q at most")
    list_all = sort_by_k(dict_all, "GCS-Q at most")
    return list_all



def compare_ours_iterative_exactly(dict_sums):
    # get results from k=2 to k=n-1
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_ours_iterative_exactly")
    #dict_n_split = filter_dict_by_substring(dict_sums, "ours_n_q")  # added the q to exclude the ours_n_half
    #dict_all = {**dict_k_split, **dict_n_split}
    #plot_values_and_times(dict_all, "Ours exactly")
    list_all = sort_by_k(dict_k_split, "Ours exactly")
    return list_all


def compare_ours_iterative_at_most(dict_sums):
    # get results from k=2 to k=n-1
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_ours_iterative_at_most")
    #dict_n_split = filter_dict_by_substring(dict_sums, "ours_n_half")
    #dict_all = {**dict_k_split, **dict_n_split}
    #plot_values_and_times(dict_all, "Ours at most")
    list_all = sort_by_k(dict_k_split, "Ours at most")
    return list_all


def compare_r_qubo_iterative(dict_sums):
    # get results from k=2 to k=n-1
    dict_k_split = filter_dict_by_substring(dict_sums, "_split_R_QUBO-iterative")
    #dict_n_split = filter_dict_by_substring(dict_sums, "R-QUBO")
    #dict_all = {**dict_k_split, **dict_n_split}
    #plot_values_and_times(dict_all, "R-QUBO")
    list_all = sort_by_k(dict_k_split, "R-QUBO")
    return list_all


def sort_by_k(input_dict, title):
    if "GCS-Q exactly" in title:
        labels_x_axis = ['GCS-Q'] + list(range(3,28))
    elif "GCS-Q at most" in title:
        labels_x_axis = ['GCS-Q'] + list(range(3,28))
    elif "Ours exactly" in title:
        labels_x_axis = list(range(2,28))
    elif "Ours at most" in title:
        labels_x_axis = list(range(2,28))
    elif "R-QUBO" in title:
        labels_x_axis = list(range(2,28))
    else:
        raise Exception("Algorithm to plot unknown -> can't sort by k")
    input = split_keys_at_underscore(input_dict)
    values_qbsolv = [input[str(key)] for key in labels_x_axis if (str(key)) in input]
    return values_qbsolv


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
        else:
            if '_q' in key:
                new_key = key.split('_q', 1)[0]  # Split of the '_qbsolv_parallel' or '_qaoa_parallel'
            elif '_d' in key:
                new_key = key.split('_d', 1)[0]  # Split of the '_dwave_parallel'
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


def pad_with_nan(value_list, num_graphsizes):
    for sub_list in value_list:
        while len(sub_list) < num_graphsizes:
            sub_list.append(np.nan)
    return value_list


def adjust_lists_with_existing(values_dict, keys_list, names_list, colors_list):
    values_list = []
    these_names = names_list.copy()
    these_colors = colors_list.copy()
    for (i, key) in enumerate(keys_list):
        if key in values_dict:
            values_list.append(values_dict[key])
        else:
            these_names.remove(names_list[i])
            these_colors.remove(colors_list[i])
    if '5' in keys_list[0]:
        num_graphsizes = 12
    else:
        num_graphsizes = 13
    return pad_with_nan(values_list, num_graphsizes), these_names, these_colors


def plot_over_graph_sizes(averages_dict, stds_dict, solver, ylabel, title, function):
    # plot graph sizes
    graph_sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
    # non-iterative algorithms
    algorithm_names_non_iterative = [f"GCS-Q", "Kochenberger", "our approach", "R-QUBO", "n-split GCS-Q"]
    algorithm_keys_non_iterative = [f"GCS-Q_{solver}_parallel", f"ours_n_{solver}", f"ours_n_half_{solver}",
                                    f"R-QUBO_{solver}", f"n_split_GCSQ_{solver}"]
    colors_non_iterative = ["C6", "C0", "C1", "C2", "C3"]

    avgs_list_non_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_non_iterative, algorithm_names_non_iterative, colors_non_iterative)
    stds_list_non_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_non_iterative, algorithm_names_non_iterative, colors_non_iterative)
    function(algorithm_names=names, values=avgs_list_non_it,
                              std_devs=stds_list_non_it, colors=colors, x_ticks=graph_sizes,
                              xlabel="n", ylabel=ylabel,
                              title=f"{title} non-iterative approaches using {solver} with respect to n")

    # without n_split GCS-Q and GCS-Q
    function(algorithm_names=names[1:-1], values=avgs_list_non_it[1:-1],
                              std_devs=stds_list_non_it[1:-1], colors=colors[1:-1], x_ticks=graph_sizes,
                              xlabel="n", ylabel=ylabel,
                              title=f"{title} non-iterative approaches using {solver} with respect to n (2)")

    print("Plotting ks")
    for k in [2, 3, 4, 5]:
        #print(f"k={k}")
        if k == 2:
            algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                         "GCS-Q"]
            algorithm_keys_iterative = [f"{k}_split_ours_iterative_exactly_{solver}_parallel",
                                        f"{k}_split_ours_iterative_at_most_{solver}_parallel",
                                        f"{k}_split_R_QUBO-iterative_{solver}_parallel",
                                        f"GCS-Q_{solver}_parallel"]
            colors_iterative = ["C0", "C1", "C2", "C6"]
        else:
            algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                         "k-split GCS-Q (exactly)", "k-split GCS-Q (at most)"]
            algorithm_keys_iterative = [f"{k}_split_ours_iterative_exactly_{solver}_parallel",
                                        f"{k}_split_ours_iterative_at_most_{solver}_parallel",
                                        f"{k}_split_R_QUBO-iterative_{solver}_parallel",
                                        f"{k}_split_GCSQ_exactly_{solver}_parallel",
                                        f"{k}_split_GCSQ_at_most_{solver}_parallel"]
            colors_iterative = ["C0", "C1", "C2", "C3", "C4"]
        avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_iterative,
                                                                     algorithm_names_iterative,
                                                                     colors_iterative)
        stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_iterative,
                                                            algorithm_names_iterative, colors_iterative)
        # if runs with this k exist for this solver
        if len(names) > 0:
            if k == 5:
                x_ticks = graph_sizes[1:]
            else:
                x_ticks = graph_sizes
            function(algorithm_names=names, values=avgs_list_it,
                                      std_devs=stds_list_it, colors=colors, x_ticks=x_ticks,
                                      xlabel="n", ylabel=ylabel,
                                      title=f"{title} {k}-split approaches using {solver} with respect to n")
            if k > 2:
                function(algorithm_names=names[:4], values=avgs_list_it[:4],
                                          std_devs=stds_list_it[:4], colors=colors[:4], x_ticks=x_ticks,
                                          xlabel="n", ylabel=ylabel,
                                          title=f"{title} {k}-split approaches using {solver} with respect to n (2)")
                function(algorithm_names=names[:3], values=avgs_list_it[:3],
                                          std_devs=stds_list_it[:3], colors=colors[:3], x_ticks=x_ticks,
                                          xlabel="n", ylabel=ylabel,
                                          title=f"{title} {k}-split approaches using {solver} with respect to n (3)")
                function(algorithm_names=names[1:3], values=avgs_list_it[1:3],
                                          std_devs=stds_list_it[1:3], colors=colors[1:3], x_ticks=x_ticks,
                                          xlabel="n", ylabel=ylabel,
                                          title=f"{title} {k}-split approaches using {solver} with respect to n (4)")

    algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                 "k-split GCS-Q (exactly)", "k-split GCS-Q (at most)"]
    # colors_iterative = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0), (1.0, 0.4980392156862745, 0.054901960784313725, 1.0), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0), (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0), (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0)]
    # C6: (0.8901960784313725, 0.4666666666666667, 0.7607843137254902, 1.0)
    colors_iterative = ["C5", "C7", "C8", "C9"]
    algorithm_keys_iterative = [f"_split_ours_iterative_exactly_{solver}_parallel",
                                f"_split_ours_iterative_at_most_{solver}_parallel",
                                f"_split_R_QUBO-iterative_{solver}_parallel",
                                f"_split_GCSQ_exactly_{solver}_parallel",
                                f"_split_GCSQ_at_most_{solver}_parallel"]

    print("Plotting iterative algorithms")
    ks = [2, 3, 4, 5]
    k_names = ["k=2", "k=3", "k=4", "k=5"]
    for i, a in enumerate(algorithm_keys_iterative):
        #print(f"algorithm={a}")
        if not "GCS" in a:
            algorithm_keys_with_k = [(str(k) + a) for k in ks]
            colors_now = colors_iterative
        else:
            k_names = ["GCS-Q", "k=3", "k=4", "k=5"]
            algorithm_keys_with_k = [f"GCS-Q_{solver}_parallel"] + [(str(k) + a) for k in ks[1:]]
            colors_now = ["C6"] + colors_iterative[1:]
        avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_with_k, k_names,
                                                                 colors_now)
        stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_with_k, k_names, colors_now)
        # make the k=5 lists long enough (if they exist)
        if "k=5" in names:
            avgs_list_it[-1].insert(0, np.nan)
            avgs_list_it[-1] = avgs_list_it[-1][:-1]
            stds_list_it[-1].insert(0, np.nan)
            stds_list_it[-1] = stds_list_it[-1][:-1]
        function(algorithm_names=names, values=avgs_list_it,
                                  std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} {algorithm_names_iterative[i]} using {solver} with respect to n")
        function(algorithm_names=names[1:], values=avgs_list_it[1:],
                                  std_devs=stds_list_it[1:], colors=colors[1:], x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} {algorithm_names_iterative[i]} using {solver} with respect to n (2)")
        function(algorithm_names=names[2:], values=avgs_list_it[2:],
                                  std_devs=stds_list_it[2:], colors=colors[2:], x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} {algorithm_names_iterative[i]} using {solver} with respect to n (3)")

    # total best
    if solver == "qbsolv":
        # total best (at least for QBSolv)
        algorithm_names_all_1 = [f"GCS-Q", "Kochenberger", "our approach", "R-QUBO", "iterative R-QUBO (4-split)",
                               "our iterative approach (4-split)"]
        algorithm_keys_all_1 = [f"GCS-Q_{solver}_parallel", f"ours_n_{solver}", f"ours_n_half_{solver}",
                              f"R-QUBO_{solver}", f"4_split_R_QUBO-iterative_{solver}_parallel",
                              f"4_split_ours_iterative_at_most_{solver}_parallel"]
        colors_all_1 = ["C6", "C0", "C1", "C2", "C8", "C5"]
        avgs_list_all, algorithm_names_all, colors_all = adjust_lists_with_existing(averages_dict, algorithm_keys_all_1,
                                                                     algorithm_names_all_1, colors_all_1)
        stds_list_all, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_all_1,
                                                            algorithm_names_all_1, colors_all_1)
        function(algorithm_names=algorithm_names_all[3:], values=avgs_list_all[3:],
                                  std_devs=stds_list_all[3:], colors=colors_all[3:],
                                  x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} the best approaches using {solver} with respect to n")
        function(algorithm_names=algorithm_names_all, values=avgs_list_all,
                                  std_devs=stds_list_all, colors=colors_all,
                                  x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} the best approaches using {solver} with respect to n (+ GCS-Q)")
    elif solver == "qaoa":
        pass
        # iterative approaches are alredy best.


def plot_over_graph_sizes_with_classical_baseline(averages_dict, stds_dict, solver, ylabel, title, function, classical):
    classical_stds = 13 * [0.0]
    averages_dict["belyi"] = classical
    stds_dict["belyi"] = classical_stds
    # plot graph sizes
    graph_sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]

    if solver == "qaoa":
        k = 2
    elif solver == "qbsolv":
        k = 4

    # non-iterative algorithms
    algorithm_names_non_iterative = ["Belyi", f"GCS-Q", "Kochenberger", "our approach", "R-QUBO", "n-split GCS-Q"]
    algorithm_keys_non_iterative = ["belyi", f"GCS-Q_{solver}_parallel", f"ours_n_{solver}", f"ours_n_half_{solver}",
                                    f"R-QUBO_{solver}", f"n_split_GCSQ_{solver}"]
    colors_non_iterative = ["0.8", "C6", "C0", "C1", "C2", "C3"]

    avgs_list_non_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_non_iterative,
                                                                 algorithm_names_non_iterative, colors_non_iterative)
    stds_list_non_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_non_iterative,
                                                        algorithm_names_non_iterative, colors_non_iterative)
    function(algorithm_names=names, values=avgs_list_non_it,
             std_devs=stds_list_non_it, colors=colors, x_ticks=graph_sizes,
             xlabel="n", ylabel=ylabel,
             title=f"{title} non-iterative approaches using {solver} with respect to n")

    algorithm_names_iterative = ["Belyi", "GCS-Q", "iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                 "k-split GCS-Q (exactly)", "k-split GCS-Q (at most)"]
    algorithm_keys_iterative = ["belyi", f"GCS-Q_{solver}_parallel", f"{k}_split_ours_iterative_exactly_{solver}_parallel",
                                f"{k}_split_ours_iterative_at_most_{solver}_parallel",
                                f"{k}_split_R_QUBO-iterative_{solver}_parallel",
                                f"{k}_split_GCSQ_exactly_{solver}_parallel",
                                f"{k}_split_GCSQ_at_most_{solver}_parallel"]
    colors_iterative = ["0.8", "C6", "C0", "C1", "C2", "C3", "C4"]
    avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_iterative,
                                                             algorithm_names_iterative,
                                                             colors_iterative)
    stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_iterative,
                                                    algorithm_names_iterative, colors_iterative)
    x_ticks = graph_sizes
    function(algorithm_names=names, values=avgs_list_it,
             std_devs=stds_list_it, colors=colors, x_ticks=x_ticks,
             xlabel="n", ylabel=ylabel,
             title=f"{title} {k}-split approaches using {solver} with respect to n")

    # total best (at least for QBSolv)
    if solver == "qbsolv":
        algorithm_names_all_1 = ["Belyi", f"GCS-Q", "our approach", "R-QUBO", "iterative R-QUBO (4-split)",
                               "our iterative approach (4-split)"]
        algorithm_keys_all_1 = ["belyi", f"GCS-Q_{solver}_parallel", f"ours_n_half_{solver}",
                              f"R-QUBO_{solver}", f"4_split_R_QUBO-iterative_{solver}_parallel",
                              f"4_split_ours_iterative_at_most_{solver}_parallel"]
        colors_all_1 = ["0.8", "C6", "C1", "C2", "C8", "C5"]
        avgs_list_all, algorithm_names_all, colors_all = adjust_lists_with_existing(averages_dict, algorithm_keys_all_1,
                                                                                    algorithm_names_all_1, colors_all_1)
        stds_list_all, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_all_1,
                                                         algorithm_names_all_1, colors_all_1)
        function(algorithm_names=algorithm_names_all, values=avgs_list_all,
                                  std_devs=stds_list_all, colors=colors_all,
                                  x_ticks=graph_sizes,
                                  xlabel="n", ylabel=ylabel,
                                  title=f"{title} the best approaches using {solver} with respect to n (+ GCS-Q)")




if __name__ == "__main__":
    # data_path = "results/eon_data/quantum/qbsolv/parallel/k=12/data_12_split_GCSQ_exactly_qbsolv_parallel__3949468976__2024-01-26_13-12-53.747804.pkl"
    # data = read_pickle_data(data_path)
    # print(data)
    data_path = "results/eon_data/quantum/qaoa"
    solver = [solver for solver in ["qbsolv", "qaoa", "dwave", "classical"] if solver in data_path][0]
    data_path_optima = "results/eon_data/classical"
    data_different_seeds = process_folder(data_path)
    feasible_data_different_seeds = filter_infeasible_coalitions(data_different_seeds)

    # get baseline
    data_classical = process_folder(data_path_optima)
    optima_dict, _, _, times_classical_dict, _, _ = calculate_statistics_over_seeds(data_classical)
    optima_list, times_classical_list = optima_dict['belyi'], times_classical_dict['belyi']
    times_classical_summed = [sum(l) for l in times_classical_list]
    optima_summed = [sum(l) for l in optima_list]

    data_values, stds_over_seeds, num_feasible_seeds, times, stds_times, num_successful_seeds = calculate_statistics_over_seeds(feasible_data_different_seeds)
    #print(data_values, "\n", stds_over_seeds)

    # Absolute values
    data_summed = sum_over_same_sized_graphs(data_values)
    stds_summed = sum_over_same_sized_graphs(stds_over_seeds)
    times_summed = sum_over_same_sized_graphs(times)
    stds_times_summed = sum_over_same_sized_graphs(stds_times)
    '''
    sorted_data_summed = sorted(data_summed.items())
    print(sorted_data_summed)
    print("\n")
    sorted_stds_summed = sorted(stds_summed.items())
    print(sorted_stds_summed)
    sorted_num_seeds = sorted(num_successful_seeds.items())
    print(sorted_num_seeds)
    '''

    # Relative values
    relative_values = relative_solution_quality(data_values, optima_list)
    for algorithm_name, graph_size_list in relative_values.items():
        for i, graph_list in enumerate(graph_size_list):
            for j, n in enumerate(graph_list):
                if n > 1.0:
                    print(algorithm_name, i, j, n)
    rel_values_over_same_sized_graphs, stds_rel_values_over_same_sized_graphs = average_over_same_sized_graphs(relative_values)
    nums_successful_over_same_sized_graphs, std_num_over_same_sized_graphs = average_over_same_sized_graphs(num_successful_seeds)
    nums_feasible_over_same_sized_graphs, std_num_feasible_over_same_sized_graphs = average_over_same_sized_graphs(num_feasible_seeds)
    rel_values_over_graph_sizes = average_over_graph_sizes(rel_values_over_same_sized_graphs)
    rel_values_over_graph_sizes_std = std_over_graph_sizes(rel_values_over_same_sized_graphs)

    # Num optima found values
    optima_found_dict = optimum_found(feasible_data_different_seeds, optima_list)
    avg_optima_found_over_seeds, std_optima_found_over_seeds, _, _, _, _ = calculate_statistics_over_seeds(optima_found_dict)
    optima_found_over_same_sized_graphs = sum_over_same_sized_graphs(avg_optima_found_over_seeds)
    stds_optima_found_over_same_sized_graphs = sum_over_same_sized_graphs(std_optima_found_over_seeds)
    num_opt_found_avg_over_graph_sizes = average_over_graph_sizes(optima_found_over_same_sized_graphs)
    stds_num_opt_found_avg_over_graph_sizes = std_over_graph_sizes(optima_found_over_same_sized_graphs)


    # Plots k's
    # relative values

    algorithm_names = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO", "k-split GCS-Q (exactly)", "k-split GCS-Q (at most)"]
    colors = ["C0", "C1", "C2", "C3", "C4"]

    averages_list_rel = [compare_ours_iterative_exactly(rel_values_over_graph_sizes), compare_ours_iterative_at_most(rel_values_over_graph_sizes), compare_r_qubo_iterative(rel_values_over_graph_sizes), compare_ks_GCSQ_exactly(rel_values_over_graph_sizes), compare_ks_GCSQ_at_most(rel_values_over_graph_sizes)]
    stds_list_rel = [compare_ours_iterative_exactly(rel_values_over_graph_sizes_std),
                     compare_ours_iterative_at_most(rel_values_over_graph_sizes_std),
                     compare_r_qubo_iterative(rel_values_over_graph_sizes_std),
                     compare_ks_GCSQ_exactly(rel_values_over_graph_sizes_std),
                     compare_ks_GCSQ_at_most(rel_values_over_graph_sizes_std)]
    x_ticks = list(range(2, 28))[:len(averages_list_rel[0])]
    plot_line_chart(algorithm_names=algorithm_names, values=averages_list_rel, x_ticks=x_ticks, colors=colors, xlabel="k", ylabel="Relative solution quality", title=f"Relative solution quality of k-split approaches using {solver}")
    # num optima found
    averages_list_opt = [compare_ours_iterative_exactly(num_opt_found_avg_over_graph_sizes),
                         compare_ours_iterative_at_most(num_opt_found_avg_over_graph_sizes),
                         compare_r_qubo_iterative(num_opt_found_avg_over_graph_sizes),
                         compare_ks_GCSQ_exactly(num_opt_found_avg_over_graph_sizes),
                         compare_ks_GCSQ_at_most(num_opt_found_avg_over_graph_sizes)]
    stds_list_opt = [compare_ours_iterative_exactly(stds_num_opt_found_avg_over_graph_sizes),
                     compare_ours_iterative_at_most(stds_num_opt_found_avg_over_graph_sizes),
                     compare_r_qubo_iterative(stds_num_opt_found_avg_over_graph_sizes),
                     compare_ks_GCSQ_exactly(stds_num_opt_found_avg_over_graph_sizes),
                     compare_ks_GCSQ_at_most(stds_num_opt_found_avg_over_graph_sizes)]
    plot_vertical_lines_chart(algorithm_names=algorithm_names, values=averages_list_opt, x_ticks=x_ticks,colors=colors,xlabel="k",ylabel="Average number of optima found per graph size",title=f"Average number of optima found of k-split approaches using {solver}")

    # absolute values
    averages_list_sums = [compare_ours_iterative_exactly(data_summed),
                         compare_ours_iterative_at_most(data_summed),
                         compare_r_qubo_iterative(data_summed),
                         compare_ks_GCSQ_exactly(data_summed),
                         compare_ks_GCSQ_at_most(data_summed)]
    stds_list_sums = [compare_ours_iterative_exactly(stds_summed),
                     compare_ours_iterative_at_most(stds_summed),
                     compare_r_qubo_iterative(stds_summed),
                     compare_ks_GCSQ_exactly(stds_summed),
                     compare_ks_GCSQ_at_most(stds_summed)]
    for a, algorithm in enumerate(algorithm_names):
        print(algorithm, ": ")
        row_names = list(range(2, 28))[:len(averages_list_sums[a])]
        table = generate_latex_table(column_names=[4,6,8,10,12,14,16,18,20,22,24,26,28], values_list=averages_list_sums[a], std_devs_list=stds_list_sums[a], row_names=row_names)
        print(table)
        print("\n\n")

    # plot over graph sizes
    # relative values
    print("Plotting relative values")
    plot_over_graph_sizes(rel_values_over_same_sized_graphs, stds_rel_values_over_same_sized_graphs, solver, "Relative solution quality", "Relative solution quality of", plot_line_chart_with_stds)
    # optima
    print("Plotting optima")
    plot_over_graph_sizes(optima_found_over_same_sized_graphs, stds_optima_found_over_same_sized_graphs, solver, "Number of optima found", "Number of optima found by", plot_vertical_lines_chart)
    # num solutions found
    print("Plotting num solutions")
    plot_over_graph_sizes(nums_successful_over_same_sized_graphs, std_num_over_same_sized_graphs, solver,
                          "Number of solutions found", "Number of solutions found by",
                          plot_vertical_lines_chart)
    # num feasible solutions found
    print("Plotting num feasible solutions")
    plot_over_graph_sizes(nums_feasible_over_same_sized_graphs, std_num_feasible_over_same_sized_graphs, solver,
                          "Number of feasible solutions found", "Number of feasible solutions found by", plot_vertical_lines_chart)
    # time taken
    print("Plotting time")
    plot_over_graph_sizes(times_summed, stds_times_summed, solver,
                          "Time (in seconds)", "Time taken by", plot_line_chart_with_stds)
    plot_over_graph_sizes_with_classical_baseline(times_summed, stds_times_summed, solver,
                          "Time (in seconds)", "Time taken by", plot_line_chart_with_stds, times_classical_summed)
    print("Done.")


