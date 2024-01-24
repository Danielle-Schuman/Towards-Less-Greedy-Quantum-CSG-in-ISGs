import matplotlib.pyplot as plt
import numpy as np
import pickle


def read_pickle_data(file_path):
    all_data = []
    try:
        with open(file_path, 'rb') as file:
            while True:
                data = pickle.load(file)
                all_data.append(data)
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
    except EOFError:
        # Reached end of file
        pass
    except Exception as e:
        print(f"An error occurred while reading {file_path}: {e}")
    return all_data

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
     data_path = "results/eon_data/quantum/qbsolv/data_ours_n_half_qbsolv__1334442952__2024-01-24_00-57-15.964398.pkl"
     data = read_pickle_data(data_path)
     print(data)