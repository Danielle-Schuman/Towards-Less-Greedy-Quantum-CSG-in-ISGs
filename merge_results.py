import os
import fnmatch
import pickle
import copy

from plotting import read_pickle_data


def get_matching_done_file(data_file, done_files):
    algorithm_name = data_file.split('data_')[1].split('__')[0]
    seed = data_file.split('data_')[1].split('__')[1].split('__')[0]
    timestamp, _ = os.path.splitext(data_file.split('data_')[1].split('__')[2])
    # should contain only one file
    matching_done_files = [f for f in done_files if f == f'done_{algorithm_name}__{seed}__{timestamp}.pkl']
    return matching_done_files[0]

def merge_and_sort_files(directory):
    for root, dirs, all_files in os.walk(directory):
        # Filter files based on their names
        data_files = fnmatch.filter(all_files, 'data_*.pkl')
        data_files_to_do = copy.copy(data_files)
        done_files = fnmatch.filter(all_files, 'done_*.pkl')

        # Load and process the files
        for data_file in data_files:
            if data_file in data_files_to_do:
                algorithm_name = data_file.split('data_')[1].split('__')[0]
                seed = data_file.split('data_')[1].split('__')[1].split('__')[0]
                data_files_same_algo_and_seed = fnmatch.filter(data_files_to_do, f'data_{algorithm_name}__{seed}__*.pkl')
                matching_done_files = [get_matching_done_file(a_data_file, done_files) for a_data_file in data_files_same_algo_and_seed]

                # Load data and done files
                data_same_algo_and_seed = []
                for a_data_file in data_files_same_algo_and_seed:
                    data = read_pickle_data(os.path.join(root, a_data_file))
                    data_same_algo_and_seed.extend(data)
                done_same_algo_and_seed = []
                for a_done_file in matching_done_files:
                    done = read_pickle_data(os.path.join(root, a_done_file))
                    done_same_algo_and_seed.extend(done)

                # Merge and sort data according to done
                zipped_lists = [[(x, y) for x, y in zip(sublist1, sublist2)] for sublist1, sublist2 in zip(done_same_algo_and_seed, data_same_algo_and_seed)]
                merged_list = [item for sublist in zipped_lists for item in sublist]
                sorted_merged_list = sorted(merged_list, key=lambda tup: (tup[0][0], tup[0][1]))
                # redivide into sublists by graph_size
                sublists = []
                current_sublist = []
                current_key = None
                for tup in sorted_merged_list:
                    key = tup[0][0]
                    if key != current_key:
                        if current_sublist:
                            sublists.append(current_sublist)
                        current_sublist = []
                        current_key = key
                    current_sublist.append(tup)
                if current_sublist:
                    sublists.append(current_sublist)
                # unzip
                done_list_total = [[tup[0] for tup in sublist] for sublist in sublists]
                data_list_total = [[tup[1] for tup in sublist] for sublist in sublists]

                # Save the merged and sorted data
                with open(os.path.join(root, f'data_{algorithm_name}__{seed}__total.pkl'), 'wb') as f:
                    pickle.dump(data_list_total, f)
                with open(os.path.join(root, f'done_{algorithm_name}__{seed}__total.pkl'), 'wb') as f:
                    pickle.dump(done_list_total, f)

                data_files_to_do = [x for x in data_files_to_do if x not in data_files_same_algo_and_seed]



# Example usage
merge_and_sort_files('results/eon_data/quantum/qaoa/parallel/k=2')
