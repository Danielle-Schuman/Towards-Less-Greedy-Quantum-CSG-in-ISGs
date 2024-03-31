import pickle

# Define the filename
path = "../results/eon_data/classical/"
filename = path + "belyi.txt"

# Initialize an empty dictionary to store data
data = {}

# Open the file and read line by line
with open(filename, 'r') as file:
    for line in file:
        # Split the line by spaces
        parts = line.split()

        # Check if the line contains network size information
        if len(parts) >= 5 and parts[0] == "Network" and parts[1] == "size":
            # Extract graph / network size
            graph_size = int(parts[3].strip('.'))

            # Extract optimum and "Batch" Time (aka total wall clock time for this graph)
            optimum = float(parts[28].strip(','))
            time = float(parts[6].strip(','))

            # If graph size not in the dictionary, add it with an empty list
            if graph_size not in data:
                data[graph_size] = []

            # Append a tuple containing optimum and time to the list
            data[graph_size].append(('_', optimum, time))

# Convert the dictionary to a list of sublists
data_list = [[values for values in data[graph_size]] for graph_size in
                sorted(data.keys())]

# Dump the list into a pickle file
pickle_filename = path + "data_belyi__0__total.pkl"
with open(pickle_filename, 'wb') as pickle_file:
    pickle.dump(data_list, pickle_file)

print("Data has been dumped into", pickle_filename)