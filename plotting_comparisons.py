from plotting import *


def plotting_comparisons(averages_dict, stds_dict, ylabel, title, function):
    labels_template = ["QBSolv", "D-Wave", "QAOA"]
    solver_keys = ["qbsolv", "dwave", "qaoa"]
    colors_template = ["C0", "C1", "C2"]
    graph_sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]

    # non-iterative algorithms and GCS-Q
    algorithm_names_non_iterative = ["GCS-Q", "Kochenberger", "our approach", "R-QUBO", "n-split GCS-Q"]
    algorithm_keys_non_iterative = [[f"GCS-Q_{solver}_parallel" for solver in solver_keys], [f"ours_n_{solver}" for solver in solver_keys], [f"ours_n_half_{solver}" for solver in solver_keys],
                                    [f"R-QUBO_{solver}" for solver in solver_keys], [f"n_split_GCSQ_{solver}" for solver in solver_keys]]

    for a, algorithm_name in enumerate(algorithm_names_non_iterative):
        avgs_list_non_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_non_iterative[a], labels_template, colors_template)
        stds_list_non_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_non_iterative[a], labels_template, colors_template)
        function(algorithm_names=names, values=avgs_list_non_it,
                                  std_devs=stds_list_non_it, colors=colors, x_ticks=graph_sizes,
                                  xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
                                  title=f"{title} {algorithm_name}")

    # iterative algorithms
    for k in [2, 3, 4, 5]:
        if k == 2:
            algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO"]
            algorithm_keys_iterative = [
                [f"{k}_split_ours_iterative_exactly_{solver}_parallel" for solver in solver_keys],
                [f"{k}_split_ours_iterative_at_most_{solver}_parallel" for solver in solver_keys],
                [f"{k}_split_R_QUBO-iterative_{solver}_parallel" for solver in solver_keys]]
        else:
            algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                         "GCS-Q (exactly)", "GCS-Q (at most)"]
            algorithm_keys_iterative = [[f"{k}_split_ours_iterative_exactly_{solver}_parallel" for solver in solver_keys],
                                        [f"{k}_split_ours_iterative_at_most_{solver}_parallel" for solver in solver_keys],
                                        [f"{k}_split_R_QUBO-iterative_{solver}_parallel" for solver in solver_keys],
                                        [f"{k}_split_GCSQ_exactly_{solver}_parallel" for solver in solver_keys],
                                        [f"{k}_split_GCSQ_at_most_{solver}_parallel" for solver in solver_keys]]
        for a, algorithm_name in enumerate(algorithm_names_iterative):
            avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_iterative[a],
                                                                         labels_template,
                                                                         colors_template)
            stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_iterative[a],
                                                                labels_template, colors_template)
            if k == 5:
                x_ticks = graph_sizes[1:]
            else:
                x_ticks = graph_sizes
            function(algorithm_names=names, values=avgs_list_it,
                                          std_devs=stds_list_it, colors=colors, x_ticks=x_ticks,
                                          xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
                                          title=f"{title} {k}-split {algorithm_name}")

    # with best k's for solvers
    algorithm_names_iterative = ["iterative Kochenberger", "our iterative approach", "iterative R-QUBO",
                                 "GCS-Q (exactly) & GCS-Q"]
    # (this list comprehension has been coded with the help of ChatGPT)
    algorithm_keys_iterative = [[f"{4 if solver == 'qbsolv' else 2}_split_ours_iterative_exactly_{solver}_parallel" for solver in solver_keys],
                                [f"{4 if solver == 'qbsolv' else 2}_split_ours_iterative_at_most_{solver}_parallel" for solver in solver_keys],
                                ["4_split_GCSQ_exactly_qbsolv_parallel", "GCS-Q_dwave_parallel", "GCS-Q_qaoa_parallel"],
                                [f"{4 if solver == 'qbsolv' else 2}_split_R_QUBO-iterative_{solver}_parallel" for solver in solver_keys]]
    for a, algorithm_name in enumerate(algorithm_names_iterative):
        avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_iterative[a],
                                                                 labels_template,
                                                                 colors_template)
        stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_iterative[a],
                                                        labels_template, colors_template)
        function(algorithm_names=names, values=avgs_list_it,
                 std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
                 xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
                 title=f"{title} the best {algorithm_name} variant")

    # total best
    # TODO: Possibly adapt after receiving final it. R-QUBO results for QAOA
    algorithm_names_best = ["QBSolv: 4-split iterative R-QUBO", "D-Wave: GCS-Q", "D-Wave: 2-split iterative R-QUBO", "QAOA: GCS-Q", "QAOA: 2-split iterative R-QUBO"]
    algorithm_keys_best = ["4_split_R_QUBO-iterative_qbsolv_parallel", "GCS-Q_dwave_parallel", "2_split_R_QUBO-iterative_dwave_parallel", "GCS-Q_qaoa_parallel", "2_split_R_QUBO-iterative_qaoa_parallel"]
    colors_best = ["C0", "C1", "C5", "C2", "C8"]
    avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_best,
                                                             algorithm_names_best,
                                                             colors_best)
    stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_best,
                                                    algorithm_names_best, colors_best)
    function(algorithm_names=names, values=avgs_list_it,
             std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
             xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
             title=f"{title} the best algorithm variants")


def plotting_comparisons_with_classical_baseline(averages_dict, stds_dict, ylabel, title, function, classical):
    classical_stds = 13 * [0.0]
    averages_dict["belyi"] = classical
    stds_dict["belyi"] = classical_stds

    labels_template = ["QBSolv", "D-Wave", "QAOA", "Belyi"]
    solver_keys = ["qbsolv", "dwave", "qaoa"]
    colors_template = ["C0", "C1", "C2", "0.8"]
    graph_sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]

    # GCS-Q
    algorithm_keys_non_iterative = [f"GCS-Q_{solver}_parallel" for solver in solver_keys]
    algorithm_keys_non_iterative.append("belyi")
    avgs_list_non_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_non_iterative,
                                                                         labels_template, colors_template)
    stds_list_non_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_non_iterative, labels_template,
                                                        colors_template)
    function(algorithm_names=names, values=avgs_list_non_it,
             std_devs=stds_list_non_it, colors=colors, x_ticks=graph_sizes,
             xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
             title=f"{title} GCS-Q")

    # with best k's for solvers
    # (this list comprehension has been coded with the help of ChatGPT)
    algorithm_keys_iterative = [f"{4 if solver == 'qbsolv' else 2}_split_R_QUBO-iterative_{solver}_parallel" for solver in solver_keys]
    algorithm_keys_iterative.append("belyi")
    avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_iterative,
                                                             labels_template,
                                                             colors_template)
    stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_iterative,
                                                    labels_template, colors_template)
    function(algorithm_names=names, values=avgs_list_it,
             std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
             xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
             title=f"{title} the best iterative R-QUBO variant")

    # total best
    # TODO: Possibly adapt after receiving final it. R-QUBO results for QAOA
    algorithm_names_best = ["QBSolv: 4-split iterative R-QUBO", "D-Wave: GCS-Q", "D-Wave: 2-split iterative R-QUBO",
                            "QAOA: GCS-Q", "QAOA: 2-split iterative R-QUBO", "Belyi"]
    algorithm_keys_best = ["4_split_R_QUBO-iterative_qbsolv_parallel", "GCS-Q_dwave_parallel",
                           "2_split_R_QUBO-iterative_dwave_parallel", "GCS-Q_qaoa_parallel",
                           "2_split_R_QUBO-iterative_qaoa_parallel", "belyi"]
    colors_best = ["C0", "C1", "C5", "C2", "C8", "0.8"]
    avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_best,
                                                             algorithm_names_best,
                                                             colors_best)
    stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_best,
                                                    algorithm_names_best, colors_best)
    function(algorithm_names=names, values=avgs_list_it,
             std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
             xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
             title=f"{title} the best algorithm variants")

    # for paper
    algorithm_names_best = ["GCS-Q (qb)", "4-split it. Kochenb. (qb)", "ZEnS (qb)", "4-split it. ZEnS (qb)", "R-QUBO (qb)", "4-split it. R-QUBO (qb)",
                            "GCS-Q (dw)", "2-split it. Kochenb. (dw)", "2-split it. ZEnS (dw)", "2-split it. R-QUBO (dw)", "Belyi"]
    algorithm_keys_best = ["GCS-Q_qbsolv_parallel", "4_split_ours_iterative_exactly_qbsolv_parallel", "ours_n_half_qbsolv", "4_split_ours_iterative_at_most_qbsolv_parallel", "R-QUBO_qbsolv", "4_split_R_QUBO-iterative_qbsolv_parallel",
                           "GCS-Q_dwave_parallel", "2_split_ours_iterative_exactly_dwave_parallel", "2_split_ours_iterative_at_most_dwave_parallel", "2_split_R_QUBO-iterative_dwave_parallel", "belyi"]
    colors_best = ["C4", "C9", "C1", "C5", "C2", "C8", "C6", "C0", "C7", "C3", "0.8"]
    avgs_list_it, names, colors = adjust_lists_with_existing(averages_dict, algorithm_keys_best,
                                                             algorithm_names_best,
                                                             colors_best)
    stds_list_it, _, _ = adjust_lists_with_existing(stds_dict, algorithm_keys_best,
                                                    algorithm_names_best, colors_best)
    function(algorithm_names=names, values=avgs_list_it,
             std_devs=stds_list_it, colors=colors, x_ticks=graph_sizes,
             xlabel="$n \in \mathbb{N}$", ylabel=ylabel,
             title=f"{title} our best approaches (& GCS-Q)")


def plotting_ks(averages_list, ylabel, title):
    # for paper
    names = ["it. Kochenb. (qb)", "it. ZEnS (qb)", "it. R-QUBO (qb)", "k-split GCS-Q (qb)",
                            "it. Kochenb. (dw)", "it. ZEnS (dw)", "it. R-QUBO (dw)", "k-split GCS-Q (dw)"]
    colors = ["C9", "C5", "C8", "C4", "C0", "C7", "C3", "C6"]
    x_ticks = list(range(2, 28))[:len(averages_list[0])]
    averages_list = pad_with_nan(averages_list, len(x_ticks))
    plot_line_chart(algorithm_names=names, values=averages_list,
             std_devs=None, colors=colors, x_ticks=x_ticks,
             xlabel="k", ylabel=ylabel,
             title=f"{title} (& GCS-Q)")


if __name__ == "__main__":
    # get classical baseline
    data_path_optima = "results/eon_data/classical"
    data_classical = process_folder(data_path_optima)
    optima_dict, _, _, times_classical_dict, _, _ = calculate_statistics_over_seeds(data_classical)
    optima_list, times_classical_list = optima_dict['belyi'], times_classical_dict['belyi']
    times_classical_summed = [sum(l) for l in times_classical_list]
    optima_summed = [sum(l) for l in optima_list]

    data_path_root = "results/eon_data/quantum/"
    solvers = ["qbsolv", "dwave"] #, "qaoa"]
    times_summed = {}
    stds_times_summed = {}
    rel_values_over_same_sized_graphs, stds_rel_values_over_same_sized_graphs = {}, {}
    nums_successful_over_same_sized_graphs, std_num_over_same_sized_graphs = {}, {}
    nums_feasible_over_same_sized_graphs, std_num_feasible_over_same_sized_graphs = {}, {}
    percent_feasible_over_same_sized_graphs, stds_percent_feasible_over_same_sized_graphs = {}, {}
    optima_found_over_same_sized_graphs = {}
    stds_optima_found_over_same_sized_graphs = {}
    averages_list_rel = []
    averages_list_opt = []

    for solver in solvers:
        data_path = data_path_root + solver
        data_different_seeds = process_folder(data_path)
        feasible_data_different_seeds = filter_infeasible_coalitions(data_different_seeds)
        data_values, stds_over_seeds, num_feasible_seeds, times, stds_times, num_successful_seeds = calculate_statistics_over_seeds(
            feasible_data_different_seeds)

        # Absolute values times
        times_summed.update(sum_over_same_sized_graphs(times))
        stds_times_summed.update(sum_over_same_sized_graphs(stds_times))

        # Relative values
        relative_values = relative_solution_quality(data_values, optima_list)
        for algorithm_name, graph_size_list in relative_values.items():
            for i, graph_list in enumerate(graph_size_list):
                for j, n in enumerate(graph_list):
                    if n > 1.0:
                        print(algorithm_name, i, j, n)
        rel_values_over_same_sized_graphs_this_solver, stds_rel_values_over_same_sized_graphs_this_solver = average_over_same_sized_graphs(
            relative_values)
        rel_values_over_same_sized_graphs.update(rel_values_over_same_sized_graphs_this_solver)
        stds_rel_values_over_same_sized_graphs.update(stds_rel_values_over_same_sized_graphs_this_solver)
        rel_values_over_graph_sizes = average_over_graph_sizes(rel_values_over_same_sized_graphs_this_solver)
        averages_list_rel.extend([compare_ours_iterative_exactly(rel_values_over_graph_sizes),
                             compare_ours_iterative_at_most(rel_values_over_graph_sizes),
                             compare_r_qubo_iterative(rel_values_over_graph_sizes),
                             compare_ks_GCSQ_exactly(rel_values_over_graph_sizes)])

        # Num solutions
        nums_successful_over_same_sized_graphs_this_solver, std_num_over_same_sized_graphs_this_solver = average_over_same_sized_graphs(
            num_successful_seeds)
        nums_successful_over_same_sized_graphs.update(nums_successful_over_same_sized_graphs_this_solver)
        std_num_over_same_sized_graphs.update(std_num_over_same_sized_graphs_this_solver)
        nums_feasible_over_same_sized_graphs_this_solver, std_num_feasible_over_same_sized_graphs_this_solver = average_over_same_sized_graphs(
            num_feasible_seeds)
        nums_feasible_over_same_sized_graphs.update(nums_feasible_over_same_sized_graphs_this_solver)
        std_num_feasible_over_same_sized_graphs.update(std_num_feasible_over_same_sized_graphs_this_solver)
        percent_feasible = get_percent_feasible(num_successful_seeds, num_feasible_seeds)
        percent_feasible_over_same_sized_graphs_this_solver, stds_percent_feasible_over_same_sized_graphs_this_solver = average_over_same_sized_graphs(
            percent_feasible)
        percent_feasible_over_same_sized_graphs.update(percent_feasible_over_same_sized_graphs_this_solver)
        stds_percent_feasible_over_same_sized_graphs.update(stds_percent_feasible_over_same_sized_graphs_this_solver)

        # Num optima found values
        optima_found_dict = optimum_found(feasible_data_different_seeds, optima_list)
        avg_optima_found_over_seeds, std_optima_found_over_seeds, _, _, _, _ = calculate_statistics_over_seeds(
            optima_found_dict)
        num_optima_this_solver = sum_over_same_sized_graphs(avg_optima_found_over_seeds)
        optima_found_over_same_sized_graphs.update(num_optima_this_solver)
        stds_optima_found_over_same_sized_graphs.update(sum_over_same_sized_graphs(std_optima_found_over_seeds))
        num_opt_found_avg_over_graph_sizes = average_over_graph_sizes(num_optima_this_solver)
        averages_list_opt.extend([compare_ours_iterative_exactly(num_opt_found_avg_over_graph_sizes),
                             compare_ours_iterative_at_most(num_opt_found_avg_over_graph_sizes),
                             compare_r_qubo_iterative(num_opt_found_avg_over_graph_sizes),
                             compare_ks_GCSQ_exactly(num_opt_found_avg_over_graph_sizes)])

    # plot over graph sizes
    # num solutions found
    print("Plotting num solutions")
    plotting_comparisons(nums_successful_over_same_sized_graphs, std_num_over_same_sized_graphs,
                          "Number of solutions found", "Number of solutions found by",
                          plot_vertical_lines_chart)
    # num feasible solutions found
    print("Plotting num feasible solutions")
    plotting_comparisons(nums_feasible_over_same_sized_graphs, std_num_feasible_over_same_sized_graphs,
                          "Number of feasible solutions found", "Number of feasible solutions found by",
                          plot_vertical_lines_chart)
    # percentage feasible solutions found
    print("Plotting percentage feasible solutions")
    plotting_comparisons(percent_feasible_over_same_sized_graphs, stds_percent_feasible_over_same_sized_graphs,
                         "Proportion of feasible solutions found",
                          "Proportion of feasible solutions found by", plot_vertical_lines_chart)

    # relative values
    print("Plotting relative values")
    plotting_comparisons(rel_values_over_same_sized_graphs, stds_rel_values_over_same_sized_graphs,
                         "Approximation ratio", "Approximation ratio of", plot_line_chart_with_stds)
    # optima
    print("Plotting optima")
    plotting_comparisons(optima_found_over_same_sized_graphs, stds_optima_found_over_same_sized_graphs,
                         "Number of optima found", "Number of optima found by", plot_vertical_lines_chart)

    # time taken
    print("Plotting time")
    plotting_comparisons(times_summed, stds_times_summed,
                          "Time (in seconds)", "Runtime of", plot_line_chart_with_stds)
    plotting_comparisons_with_classical_baseline(times_summed, stds_times_summed,
                                                  "Wall clock time (in seconds)", "Runtime of", plot_line_chart_with_stds,
                                                  times_classical_summed)

    # plotting ks
    print("Plotting ks")
    plotting_ks(averages_list_rel, "Approximation ratio", "AR of our k-split approaches")
    plotting_ks(averages_list_opt, "Average number of optima found per graph size", "Average NO by our k-split approaches")

    print("Done.")
