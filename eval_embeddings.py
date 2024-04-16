from qubits import qubits, couplers
import pickle
import utils
from jonas import Jonas

embeddings_path = "embeddings/embedding_ours_n_half_dwave_2024-03-27_02-59-05.840991.pkl"
#embeddings_path = "embeddings/embedding_ours_n_half_dwave_2024-03-27_02-10-43.579582.pkl"
with open(embeddings_path, 'rb') as file:
    embedding = pickle.load(file)

for logical_qubit, physical_qubits in embedding.items():
    for physical_qubit in physical_qubits:
        if physical_qubit not in qubits:
            print(f"Qubit {physical_qubit} missing from machine!")

# loading E.ON data
data = pickle.load(open('data/data_new_20samples_4_28.pkl', 'rb'))
edges = utils.transform(data[22][0])
algorithm = Jonas(seed=0, num_graph_sizes=20, solver="dwave")
qubo = algorithm.get_qubo(22, edges)
couplers_missing = 0
for (q1,q2) in qubo:
    physical_qubits1 = embedding[q1]
    physical_qubits2 = embedding[q2]
    coupler_found = False
    for p1 in physical_qubits1:
        for p2 in physical_qubits2:
            if ([p1, p2] in couplers or [p2, p1] in couplers):
                coupler_found = True
    if not coupler_found:
        print(f"Coupler ({q1}, {q2}) missing from machine!")
        couplers_missing += 1
print(f"{couplers_missing} Couplers missing in total")

# Result: 0 Couplers missing in total -> ?