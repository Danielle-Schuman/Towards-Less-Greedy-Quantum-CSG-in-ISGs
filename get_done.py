import os
from plotting import read_pickle_data

def process_folder(folder_path):
    folder_contents = {}
    try:
        for root, dirs, files in os.walk(folder_path):
            for file_name in files:
                if file_name.startswith("done_") and file_name.endswith(".pkl"):
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

if __name__ == "__main__":
    '''
    data_path = "results/eon_data/quantum/unused_for_plotting/other_unused_qaoa/k=2_originals/data_2_split_ours_iterative_at_most_qaoa_parallel__3949468976__2024-03-21_23-39-40.139813.pkl"
    data = read_pickle_data(data_path)
    print(data)
    '''
    data_path = "results/eon_data/quantum/qaoa/parallel/k=2"
    done_different_seeds = process_folder(data_path)
    sorted = sorted(done_different_seeds.items())
    for algorith_name, data in sorted:
        print(algorith_name, ": ", sorted)