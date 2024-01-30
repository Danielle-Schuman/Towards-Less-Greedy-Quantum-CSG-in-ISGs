import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import statistics


def read_pickle_data(file_path, graph_num_per_size=20):
    all_data = []
    try:
        with open(file_path, 'rb') as file:
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

    # Group values by algorithm_name
    algorithm_name_groups = {}
    for key, values in folder_dict.items():
        algorithm_name, seed = key
        if algorithm_name not in algorithm_name_groups:
            algorithm_name_groups[algorithm_name] = {}
        algorithm_name_groups[algorithm_name][seed] = values

    # Calculate averages and standard deviations for each algorithm_name group
    for algorithm_name, seed_groups in algorithm_name_groups.items():
        num_sublists = len(next(iter(seed_groups.values())))
        averages = [[] for _ in range(num_sublists)]
        stds = [[] for _ in range(num_sublists)]

        for seed, graph_sizes_list in seed_groups.items():
            for i, graph_list in enumerate(graph_sizes_list):
                for j, item in enumerate(graph_list):
                    if len(averages[i]) <= j:
                        averages[i].append((item[1], item[2]))
                        stds[i].append(([item[1]], [item[2]]))
                    else:
                        averages[i][j] = (averages[i][j][0] + item[1], averages[i][j][1] + item[2])
                        stds[i][j][0].append(item[1])
                        stds[i][j][1].append(item[2])

        for i, graph_list in enumerate(averages):
            averages[i] = [(value_sum / len(seed_groups), time_sum / len(seed_groups)) for (value_sum, time_sum) in graph_list]

        for i, graph_list in enumerate(stds):
            for j, seeds_list in enumerate(graph_list):
                if len(seeds_list[0]) > 1:
                    stds_values = statistics.stdev(seeds_list[0])
                else:
                    stds_values = 0.0
                if len(seeds_list[1]) > 1:
                    stds_times = statistics.stdev(seeds_list[1])
                else:
                    stds_times = 0.0
                stds[i][j] = (stds_values, stds_times)

        dict_avgs[algorithm_name] = averages
        dict_stds[algorithm_name] = stds

    return dict_avgs, dict_stds


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


def compare_ks_GCSQ_exactly(dict_sums):
    # get results from k=2 (GCS-Q) to k=n (n_split), and combine them into one dict
    dict_GCSQ = filter_dict_by_substring(dict_sums, "_GCS-Q_")
    dict_k_split_GCSQ = filter_dict_by_substring(dict_sums, "_split_GCSQ_exactly")
    dict_n_split_GCSQ = filter_dict_by_substring(dict_sums, "n_split_GCSQ")
    dict_all = {**dict_GCSQ, **dict_k_split_GCSQ, **dict_n_split_GCSQ}
    plot_values_and_times(dict_all, "GCS_Q exactly")


def plot_values_and_times(dict_all, title):
    # split dict for different solvers
    dict_qbsolv = remove_substring_from_keys(filter_dict_by_substring(dict_all, "_qbsolv"), "_qbsolv")
    dict_qaoa = remove_substring_from_keys(filter_dict_by_substring(dict_all, "_qaoa"), "_qaoa")
    # split dicts into coalition values and times
    dict_qbsolv_values = {key: value[0] for key, value in dict_qbsolv.items()}
    dict_qbsolv_times = {key: value[1] for key, value in dict_qbsolv.items()}
    dict_qaoa_values = {key: value[0] for key, value in dict_qaoa.items()}
    dict_qaoa_times = {key: value[1] for key, value in dict_qaoa.items()}
    # TODO: split into different graph sizes!!! -> One plot per graph size
    plot_bar_chart(dict_qbsolv_values, dict_qaoa_values, "Sums of coalition values", f"Sums of coalition values of 20 graphs for {title}")
    plot_bar_chart(dict_qbsolv_times, dict_qaoa_times, "Sums of times in seconds",
                   f"Time taken for 20 graphs for {title}")


def plot_bar_chart(dict_qbsolv, dict_qaoa, y_label, title):
    keys = dict_qbsolv.keys()
    values_qbsolv = dict_qbsolv.values()
    values_qaoa = dict_qaoa.values()

    x = np.arange(len(keys))
    width = 0.35

    fig, ax = plt.subplots()
    bars_qbsolv = ax.bar(x - width/2, values_qbsolv, width, label="qbsolv")
    bars_qaoa = ax.bar(x + width/2, values_qaoa, width, label="qaoa")

    #ax.set_xlabel('Categories')
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(keys)
    ax.legend()

    plt.savefig(f"plots/{title.replace(' ', '_')}.pdf", format='pdf')
    plt.show()
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
    plt.show()
    plt.close()

if __name__ == "__main__":
     data_path = "results/eon_data/quantum/qbsolv/test"
     data = process_folder(data_path)
     #data = read_pickle_data(data_path)
     data_avgs, data_stds = calculate_statistics_over_seeds(data)
     print(data_avgs, "\n", data_stds)
     avg_data_summed = sum_over_same_sized_graphs(data_avgs)
     print(avg_data_summed)
     ours_n_half = filter_dict_by_substring(avg_data_summed, "half")
     print(ours_n_half)