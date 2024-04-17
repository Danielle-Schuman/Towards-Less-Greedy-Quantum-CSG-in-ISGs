/*                                                                            
    Copyright 2022
    Alexander Belyi <alexander.belyi@gmail.com>,
    Stanislav Sobolevsky <stanly@mit.edu>                                               
                                                                            
    This file is part of Combo algorithm.

    Combo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Combo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Combo.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Graph.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <locale>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
using namespace std;

Graph ReadFromEdgelist(const string& file_name, double mod_resolution, bool treat_as_modularity)
{
	ifstream file(file_name.c_str());
	if (!file.is_open()) {
        cerr << "File " << file_name << " can not be opened." << endl;
		return Graph();
    }
	vector<int> sources, destinations;
	vector<double> weights;
	int min_vertex_number = INT_MAX;
	int max_vertex_number = 0;
	while (file.good()) {
		string line;
		getline(file, line);
		int src = -1, dst = -1;
		double weight = 1.0;
		stringstream str_stream(line);
		str_stream >> src >> dst;
		if (!str_stream.eof())
			str_stream >> weight;
		if (!str_stream.fail() && src > -1 && dst > -1)
		{
			min_vertex_number = min(min_vertex_number, min(src, dst));
			max_vertex_number = max(max_vertex_number, max(src, dst));
			sources.push_back(src);
			destinations.push_back(dst);
			weights.push_back(weight);
		}
	}
	file.close();
	bool is_directed;

    // This if-statement has been added by Daniëlle Schuman (translated from Python with ChatGPT)
	if (file_name.find("eon") != std::string::npos) {
    	is_directed = false;
	} else {
    	is_directed = true;
	}
	size_t size = size_t(max(0, 1 + max_vertex_number - min_vertex_number));
	for (size_t i = 0; i < sources.size(); ++i) {
		sources[i] -= min_vertex_number;
		destinations[i] -= min_vertex_number;
	}
	return Graph(size, sources, destinations, weights, is_directed, mod_resolution, treat_as_modularity);
}

Graph ReadFromPajek(const string& file_name, double mod_resolution, bool treat_as_modularity)
{
	ifstream file(file_name.c_str());
	if (!file.is_open()) {
		cerr << "File " << file_name << " can not be opened." << endl;
		return Graph();
    }
    std::locale locale;
	vector<int> sources, destinations;
	vector<double> weights;
	int min_vertex_number = INT_MAX;
	int max_vertex_number = 0;
	bool is_directed = false;
	bool skip = true;
	while (file.good()) {
		string line;
		getline(file, line);
		stringstream str_stream(line);
        string trimmed_string;
        str_stream >> trimmed_string; // Strips carriage return on Windows
        transform(trimmed_string.begin(), trimmed_string.end(), trimmed_string.begin(), [&locale](char elem){return std::tolower(elem, locale);});
		if (trimmed_string == "*edges") {
			skip = false;
			is_directed = false;
		} else if (trimmed_string == "*arcs") {
			skip = false;
			is_directed = true;
		} else if (skip) {
			str_stream.str(line);
			int v = -1;
			str_stream >> v;
			if (!str_stream.fail() && v > -1) {
				min_vertex_number = min(min_vertex_number, v);
				max_vertex_number = max(max_vertex_number, v);
			}
		} else if (!skip) {
			int src = -1, dst = -1;
			double weight = 1.0;
			str_stream.str(line);
			str_stream >> src >> dst;
			if (!str_stream.eof())
				str_stream >> weight;
			if (!str_stream.fail() && src > -1 && dst > -1) {
				min_vertex_number = min(min_vertex_number, min(src, dst));
				max_vertex_number = max(max_vertex_number, max(src, dst));
				sources.push_back(src);
				destinations.push_back(dst);
				weights.push_back(weight);
			}
		}
	}
	file.close();
	size_t size = size_t(max(0, 1 + max_vertex_number - min_vertex_number));
	for (size_t i = 0; i < sources.size(); ++i) {
		sources[i] -= min_vertex_number;
		destinations[i] -= min_vertex_number;
	}
	return Graph(size, sources, destinations, weights, is_directed, mod_resolution, treat_as_modularity);
}

Graph ReadFromCSV(const string& file_name, double mod_resolution, bool treat_as_modularity)
{
	ifstream file(file_name.c_str());
	if (!file.is_open()) {
		cerr << "File " << file_name << " can not be opened." << endl;
		return Graph();
	}
    vector<vector<double>> matrix;
    size_t src = 0;
    while (file.good()) {
		string line;
		getline(file, line);
		stringstream str_stream(line);
		string trimmed_line;
		str_stream >> trimmed_line;
		if (trimmed_line != "")
			matrix.push_back(vector<double>());
		else
			break;
		str_stream = stringstream(line);
		double weight = 1.0;
		char separator;
		while (str_stream.good()) {
			str_stream >> weight;
			str_stream >> separator;
			matrix[src].push_back(weight);
		}
		if (src > 0 && matrix[src].size() != matrix[src-1].size()) {
			cerr << "Error in ReadFromCSV: matrix is not square." << endl;
			return Graph();
		}
		++src;
	}
	file.close();
	if (matrix.size() != matrix.back().size()) {
		cerr << "Error in ReadFromCSV: matrix is not square." << endl;
		return Graph();
	}
    return Graph(matrix, mod_resolution, treat_as_modularity);
}

Graph ReadGraphFromFile(const string& file_name, double mod_resolution, bool treat_as_modularity)
{
	size_t dot_position = file_name.rfind('.');
	if (dot_position != string::npos) {
		string ext = file_name.substr(file_name.rfind('.'), file_name.length() - file_name.rfind('.'));
		if (ext == ".edgelist")
			return ReadFromEdgelist(file_name, mod_resolution, treat_as_modularity);
		else if (ext == ".net")
			return ReadFromPajek(file_name, mod_resolution, treat_as_modularity);
		else if (ext == ".csv")
			return ReadFromCSV(file_name, mod_resolution, treat_as_modularity);
	}
	cerr << "Error in ReadGraphFromFile: unsupported file format. Must be Pajek .net, .edgelist or .csv." << endl;
	return Graph();
}

void swap(Graph& left, Graph& right)
{
	using std::swap;
	swap(left.m_number_of_communities, right.m_number_of_communities);
	swap(left.m_is_directed, right.m_is_directed);
	swap(left.m_modularity_resolution, right.m_modularity_resolution);
	swap(left.m_modularity_matrix, right.m_modularity_matrix);
	swap(left.m_communities, right.m_communities);
}

Graph::Graph(bool is_directed, double modularity_resolution)
{
	m_is_directed = is_directed;
	m_number_of_communities = 0;
	m_modularity_resolution = modularity_resolution;
}

Graph::Graph(size_t size, const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights,
	bool is_directed, double modularity_resolution, bool treat_as_modularity) :
	Graph(is_directed, modularity_resolution)
{
	if (treat_as_modularity)
		FillModMatrix(size, sources, destinations, weights);
	else
		CalcModMatrix(size, sources, destinations, weights);
}

Graph::Graph(size_t size, const std::vector<std::tuple<int, int, double>>& edges, bool is_directed, double modularity_resolution, bool treat_as_modularity) :
	Graph(is_directed, modularity_resolution)
{
	if (treat_as_modularity)
		FillModMatrix(size, edges);
	else
		CalcModMatrix(size, edges);
}

bool IsMatrixSymmetric(const std::vector<std::vector<double>>& matrix)
{
	for (size_t i = 0; i < matrix.size(); ++i)
		for (size_t j = i+1; j < matrix.size(); ++j)
			if (matrix[i][j] != matrix[j][i])
				return false;
	return true;
}

Graph::Graph(const vector<vector<double>>& matrix, double modularity_resolution, bool treat_as_modularity) :
	Graph(false, modularity_resolution)
{
	for (size_t i = 0; i < matrix.size(); ++i) {
		if (matrix.size() != matrix[i].size()) {
			cerr << "Error in Graph(matrix): matrix must be a square matrix." << endl;
			return;
		}
	}
	m_is_directed = !IsMatrixSymmetric(m_modularity_matrix);
	if (treat_as_modularity)
		FillModMatrix(matrix);
	else
		CalcModMatrix(matrix);
}

Graph::Graph(vector<vector<double>>&& matrix, double modularity_resolution, bool treat_as_modularity) :
	Graph(false, modularity_resolution)
{
	for (size_t i = 0; i < matrix.size(); ++i) {
		if (matrix.size() != matrix[i].size()) {
			cerr << "Error in Graph(matrix): matrix must be a square matrix." << endl;
			return;
		}
	}
	m_is_directed = !IsMatrixSymmetric(m_modularity_matrix);
	if (treat_as_modularity)
		FillModMatrix(std::move(matrix));
	else
		CalcModMatrix(matrix);
}

Graph::Graph(const Graph& graph) :
	m_modularity_matrix(graph.m_modularity_matrix),
	m_communities(graph.m_communities)
{
	m_is_directed = graph.IsDirected();
	m_modularity_resolution = graph.ModularityResolution();
	m_number_of_communities = graph.NumberOfCommunities();
}

Graph::Graph(Graph&& graph) noexcept :
	Graph()
{
	swap(*this, graph);
}

Graph& Graph::operator=(Graph graph)
{
	swap(*this, graph);
	return *this;
}

void SymmetrizeMatrix(vector<vector<double>>& matrix)
{
	for (size_t i = 0; i < matrix.size(); ++i)
		for (size_t j = i+1; j < matrix.size(); ++j)
			matrix[i][j] = matrix[j][i] = (matrix[i][j] + matrix[j][i]) / 2;
}

void Graph::CalcModMatrix(size_t size, const vector<int>& sources, const vector<int>& destinations, const vector<double>& weights)
{
	double total_weight = 0.0;
	for (double weight : weights)
		total_weight += weight;
	if (!m_is_directed)
		total_weight *= 2;
	int max_index = max(*max_element(sources.begin(), sources.end()), *max_element(destinations.begin(), destinations.end()));
	if (max_index < 0) {
		cerr << "Error in CalcModMatrix: vertices' index cannot be negative" << endl;
		return;
	}
	if (size_t(max_index) >= size) {
		cerr << "Error in CalcModMatrix: vertices' index cannot be greater than size-1" << endl;
		return;
	}
	m_modularity_matrix.assign(size, vector<double>(size, 0));
	vector<double> sumQ2(size, 0.0);
	vector<double> sumQ1(size, 0.0);
	for (size_t i = 0; i < sources.size(); ++i) {
		size_t source = size_t(sources[i]);
		size_t destination = size_t(destinations[i]);
		m_modularity_matrix[source][destination] += weights[i] / total_weight;
		if (!m_is_directed)
			m_modularity_matrix[destination][source] += weights[i] / total_weight;
		sumQ1[source] += weights[i] / total_weight;
		sumQ2[destination] += weights[i] / total_weight;
		if (!m_is_directed) {
			sumQ1[destination] += weights[i] / total_weight;
			sumQ2[source] += weights[i] / total_weight;
		}
	}
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < size; ++j)
			m_modularity_matrix[i][j] -= m_modularity_resolution * sumQ1[i]*sumQ2[j];
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::FillModMatrix(size_t size, const vector<int>& sources, const vector<int>& destinations, const vector<double>& weights)
{
	int max_index = max(*max_element(sources.begin(), sources.end()), *max_element(destinations.begin(), destinations.end()));
	if (max_index < 0) {
		cerr << "Error in FillModMatrix: vertices' index cannot be negative" << endl;
		return;
	}
	if (size_t(max_index) >= size) {
		cerr << "Error in FillModMatrix: vertices' index cannot be greater than size-1" << endl;
		return;
	}
	m_modularity_matrix.assign(size, vector<double>(size, 0));
	for (size_t i = 0; i < sources.size(); ++i) {
		size_t source = size_t(sources[i]);
		size_t destination = size_t(destinations[i]);
		if (m_is_directed)
			m_modularity_matrix[source][destination] += weights[i];
		else {
			m_modularity_matrix[source][destination] += weights[i] / 2;
			m_modularity_matrix[destination][source] += weights[i] / 2;
		}
	}
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::CalcModMatrix(size_t size, const std::vector<std::tuple<int, int, double>>& edges)
{
	double total_weight = 0.0;
	int max_index = 0;
	for (const tuple<int, int, double>& edge : edges) {
		total_weight += get<2>(edge);
		max_index = max(max_index, max(get<0>(edge), get<1>(edge)));
	}
	if (max_index < 0) {
		cerr << "Error in CalcModMatrix: vertices' index cannot be negative" << endl;
		return;
	}
	if (size_t(max_index) >= size) {
		cerr << "Error in CalcModMatrix: vertices' index cannot be greater than size-1" << endl;
		return;
	}
	if (!m_is_directed)
		total_weight *= 2;
	m_modularity_matrix.assign(size, vector<double>(size, 0));
	vector<double> sumQ2(size, 0.0);
	vector<double> sumQ1(size, 0.0);
	for (const tuple<int, int, double>& edge : edges) {
		size_t source = size_t(get<0>(edge));
		size_t destination = size_t(get<1>(edge));
		double weight = get<2>(edge);
		m_modularity_matrix[source][destination] += weight / total_weight;
		if (!m_is_directed)
			m_modularity_matrix[destination][source] += weight / total_weight;
		sumQ1[source] += weight / total_weight;
		sumQ2[destination] += weight / total_weight;
		if (!m_is_directed) {
			sumQ1[destination] += weight / total_weight;
			sumQ2[source] += weight / total_weight;
		}
	}
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < size; ++j)
			m_modularity_matrix[i][j] -= m_modularity_resolution * sumQ1[i]*sumQ2[j];
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::FillModMatrix(size_t size, const std::vector<std::tuple<int, int, double>>& edges)
{
	int max_index = 0;
	for (const tuple<int, int, double>& edge : edges)
		max_index = max(max_index, max(get<0>(edge), get<1>(edge)));
	if (max_index < 0) {
		cerr << "Error in FillModMatrix: vertices' index cannot be negative" << endl;
		return;
	}
	if (size_t(max_index) >= size) {
		cerr << "Error in FillModMatrix: vertices' index cannot be greater than size-1" << endl;
		return;
	}
	m_modularity_matrix.assign(size, vector<double>(size, 0));
	for (const tuple<int, int, double>& edge : edges) {
		size_t source = size_t(get<0>(edge));
		size_t destination = size_t(get<1>(edge));
		double weight = get<2>(edge);
		if (m_is_directed)
			m_modularity_matrix[source][destination] += weight;
		else {
			m_modularity_matrix[source][destination] += weight / 2;
			m_modularity_matrix[destination][source] += weight / 2;
		}
	}
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::CalcModMatrix(const std::vector<std::vector<double>>& matrix)
{
	size_t size = matrix.size();
	double total_weight = 0.0;
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < size; ++j)
			total_weight += matrix[i][j];
	// double loop edges for undirected graphs
	if (!m_is_directed)
		for (size_t i = 0; i < size; ++i)
			total_weight += matrix[i][i];
	m_modularity_matrix.assign(size, vector<double>(size, 0.0));
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j)
			m_modularity_matrix[i][j] = matrix[i][j] / total_weight;
		if (!m_is_directed)
			m_modularity_matrix[i][i] *= 2;
	}
	vector<double> sumQ2(size, 0.0);
	vector<double> sumQ1(size, 0.0);
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < size; ++j) {
			sumQ1[i] += m_modularity_matrix[i][j];
			sumQ2[j] += m_modularity_matrix[i][j];
		}
	for (size_t i = 0; i < size; ++i)
		for (size_t j = 0; j < size; ++j)
			m_modularity_matrix[i][j] -= m_modularity_resolution * sumQ1[i]*sumQ2[j];
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::FillModMatrix(const std::vector<std::vector<double>>& matrix)
{
	m_modularity_matrix = matrix;
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::FillModMatrix(std::vector<std::vector<double>>&& matrix)
{
	std::swap(m_modularity_matrix, matrix);
	if (m_is_directed)
		SymmetrizeMatrix(m_modularity_matrix);
}

void Graph::Print() const
{
	cout << "Modularity matrix:" << endl;
	for (size_t i = 0; i < Size(); ++i) {
		for (size_t j = 0; j < Size(); ++j) {
			cout << m_modularity_matrix[i][j] << '\t';
		}
		cout << endl;
	}
}

void Graph::PrintCommunity(const string& file_name) const
{
    if (file_name.empty()) {
        cout << "Nodes\' communities:\n";
        for (size_t i = 0; i < Size(); ++i)
		    cout << m_communities[i] + 1 << ' ';
        cout << endl;
    } else {
        ofstream file(file_name.c_str());
        if (!file.is_open()) {
            cerr << "File " << file_name << " can not be opened." << endl;
            return;
        }
        for (size_t i = 0; i < Size(); ++i)
            file << m_communities[i] << endl;
        file.close();
    }
}

void Graph::SetCommunities(const vector<size_t>& new_communities, size_t number)
{
	if (Size() != new_communities.size()) {
		cerr << "Error in SetCommunities: number of elements in new_communities must be equal to graph size." << endl;
		return;
	}
	m_communities = new_communities;
	if (number == 0)
		m_number_of_communities = set<size_t>(m_communities.begin(), m_communities.end()).size();
	else
		m_number_of_communities = number;
}

double Graph::Modularity() const
{
	double modularity = 0;
	for (size_t i = 0; i < m_modularity_matrix.size(); ++i)
		for (size_t j = 0; j < m_modularity_matrix.size(); ++j)
			if (m_communities[i] == m_communities[j])
				modularity += m_modularity_matrix[i][j];
	return modularity;
}

void Graph::PerformSplit(size_t origin, size_t destination, const vector<bool>& to_be_moved)
{
	if (destination > m_number_of_communities) {
		cerr << "WARNING: in PerformSplit, destination community is greater than number of communities." << endl;
		destination = m_number_of_communities;
	}
	if (destination == m_number_of_communities)
		++m_number_of_communities;
	for (size_t i = 0; i < Size(); ++i)
		if (m_communities[i] == origin && to_be_moved[i])
			m_communities[i] = destination;
}

bool Graph::IsCommunityEmpty(size_t community) const
{
	for (size_t i = 0; i < Size(); ++i)
		if (m_communities[i] == community)
			return false;
	return true;
}

bool Graph::DeleteCommunityIfEmpty(size_t community)
{
	if (IsCommunityEmpty(community)) {
		set<size_t> community_labels;
        for (size_t i = 0; i < Size(); ++i) {
			if (m_communities[i] > community)
				--m_communities[i];
			community_labels.insert(m_communities[i]);
		}
		m_number_of_communities = community_labels.size();
        return true;
	}
	return false;
}

vector<size_t> Graph::CommunityIndices(size_t community) const
{
	vector<size_t> res;
	for (size_t i = 0; i < Size(); ++i)
		if (m_communities[i] == community)
			res.push_back(i);
	return res;
}

vector< vector<double> > Graph::GetModularitySubmatrix(const vector<size_t>& indices) const
{
    return Submatrix(m_modularity_matrix, indices);
}

vector<double> Graph::GetCorrectionVector(const vector<size_t>& orig_comm_ind, const vector<size_t>& dest_comm_ind) const
{
	vector<double> res(orig_comm_ind.size(), 0.0);
	for (size_t i = 0; i < orig_comm_ind.size(); ++i)
		for (size_t j = 0; j < dest_comm_ind.size(); ++j)
			res[i] += m_modularity_matrix[dest_comm_ind[j]][orig_comm_ind[i]];
	return res;
}

std::vector<size_t> Graph::MergeIdenticalNodes()
{
	// first, find and mark identical connected nodes
	// for each node in a group of identical nodes, their representative is minimal node from the group
    vector<size_t> representatives(m_modularity_matrix.size());
    iota(representatives.begin(), representatives.end(), 0);
    for (size_t i = 0; i < m_modularity_matrix.size(); ++i)
		if (representatives[i] == i)
			for (size_t j = i+1; j < m_modularity_matrix.size(); ++j)
				if (m_modularity_matrix[i][j] >= 0) {
					bool same = true;
					for (size_t k = 0; k < m_modularity_matrix.size(); ++k)
						if (k != i && k != j && fabs(m_modularity_matrix[i][k] - m_modularity_matrix[j][k]) > EPS) {
							same = false;
							break;
						}
					if (same)
						representatives[j] = i;
				}
	// then merge identical nodes
    int cnt_merged = 0;
	for (size_t i = 1; i < representatives.size(); ++i)
        if (representatives[i] != i - cnt_merged) {
            m_modularity_matrix = MergeTwoNodes(m_modularity_matrix, representatives[i], i - cnt_merged);
            for (size_t &r : representatives)
                if (r > i - cnt_merged)
                    --r;
            ++cnt_merged;
        }
    return representatives;
}

std::vector<size_t> Graph::MergeStronglyConnected()
{
    vector<size_t> representatives(m_modularity_matrix.size());
    iota(representatives.begin(), representatives.end(), 0);
	for (size_t i = 0; i < m_modularity_matrix.size(); /*empty*/) {
		double sum_pos = 0;
		double sum_neg = 0;
		double max_weight = 0;
		size_t positive_cnt = 0;
		size_t max_ind = 0;
		for (size_t j = 0; j < m_modularity_matrix.size(); ++j)
			if (i != j) {
				if (m_modularity_matrix[i][j] > 0) {
					sum_pos += m_modularity_matrix[i][j];
					if (m_modularity_matrix[i][j] > max_weight) {
						max_weight = m_modularity_matrix[i][j];
						max_ind = j;
					}
					++positive_cnt;
				} else {
					sum_neg += m_modularity_matrix[i][j];
				}
			}
		if (positive_cnt > 0 && max_weight >= sum_pos - max_weight - sum_neg) {
			size_t src = i;
			size_t dst = max_ind;
			if (i < max_ind) {
				++i;
				swap(src, dst);
			}
			m_modularity_matrix = MergeTwoNodes(m_modularity_matrix, dst, src);
			for (size_t r = 1; r < representatives.size(); ++r) {
				if (representatives[r] == src)
					representatives[r] = dst;
				else if (representatives[r] > src)
					--representatives[r];
			}
		} else {
			if (positive_cnt == 0)
				cout << "Disconnected node i=" << i << endl;
			++i;
		}
	}
    return representatives;
}

std::vector<size_t> Graph::ReduceSize()
{
    vector<size_t> old_to_new(m_modularity_matrix.size());
    iota(old_to_new.begin(), old_to_new.end(), 0);
	bool merged_identical = true;
	bool merged_strong = true;
	size_t current_size = m_modularity_matrix.size();
	while (true) {
    	vector<size_t> representatives = MergeIdenticalNodes();
		merged_identical = (current_size > m_modularity_matrix.size());
		if (!merged_identical && !merged_strong)
			break;
		else if (merged_identical)
			for (size_t i = 0; i < old_to_new.size(); ++i)
				old_to_new[i] = representatives[old_to_new[i]];
		current_size = m_modularity_matrix.size();
    	representatives = MergeStronglyConnected();
		merged_strong = (current_size > m_modularity_matrix.size());
		if (!merged_identical && !merged_strong)
			break;
		else if (merged_strong)
			for (size_t i = 0; i < old_to_new.size(); ++i)
				old_to_new[i] = representatives[old_to_new[i]];
	}
	return old_to_new;
}
