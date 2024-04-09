from plotting import *
import copy


def process_folder_physical_qubits(folder_path):
    folder_contents = {}
    try:
        for root, dirs, files in os.walk(folder_path):
            for file_name in files:
                if file_name.startswith("physical_qubits_"):
                    # Extracting key from file name
                    parts = file_name.split(".")[0].split("_")
                    key = "_".join(parts[2:])
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


def prune_zeros(data):
    pruned_data = copy.deepcopy(data)
    for algorithm_name in data:
        for g, graph_size_list in enumerate(data[algorithm_name]):
            all_zeros = True
            for item in graph_size_list:
                if not(item == 0):
                    all_zeros = False
                    break
            if all_zeros:
                pruned_data[algorithm_name].remove(graph_size_list)
        for g, graph_size_list in enumerate(pruned_data[algorithm_name]):
            for i, item in enumerate(graph_size_list):
                if item == 0:
                    pruned_data[algorithm_name][g].pop(i)
    return pruned_data




if __name__ == "__main__":
    data_path = "embeddings/"
    data = process_folder_physical_qubits(data_path)
    data = prune_zeros(data)
    averages, stds = average_over_same_sized_graphs(data)
    plot_over_graph_sizes(averages, stds, "dwave", "Number of qubits needed", "Number of qubits needed by", plot_line_chart_with_stds)
