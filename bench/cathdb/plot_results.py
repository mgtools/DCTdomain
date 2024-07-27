"""Plots AUC1 values of all benchmarked methods for each query in the CATH dataset.

__author__ = "Ben Iovino"
__date__ = "4/23/24"
"""

import matplotlib.pyplot as plt


def get_fams(file: str) -> dict[str, list]:
    """Returns a dictionary of family names and the sequences that belong to them.

    Args:
        file (str): Path to queries file

    Returns:
        dict[str, list]: key: family name, value: list of sequences
    """

    with open(file, 'r', encoding='utf8') as file:
        fams: dict[str, list] = {}
        for line in file:
            if line.startswith('>'):
                line = line.split('|')
                dom, fam = line[0], line[1].strip()
                fams[fam] = fams.get(fam, []) + [dom]

    return fams


def read_results(path: str, query_ind: int, result_ind: int) -> dict[str, set]:
    """Returns a dictionary of query PID's and their top hits until the first FP. Only one hit per
    unique sequence is counted, hence the use of a set. This allows for the calculation of both
    AUC1 and top1 scores.

    Args:
        path (str): Path to results file
        query (int): Column index of query PID
        result (int): Column index of result PID

    Returns:
        dict[str, set]: key: query PID, value: set of top hits until the first FP
    """

    with open(f'{path}', 'r', encoding='utf8') as file:
        results, curr_query = {}, ''
        for line in file:
            line = line.split()

            # If query had FP, continue until new query
            query = line[query_ind]
            result = line[result_ind]
            if query == curr_query:
                continue
            results[query] = results.get(query, set())

            # Ignore self-hits
            if query == result:
                continue

            # Stop counting hits for current query if domains are different
            query_dom = query.split('|')[1]
            result_dom = result.split('|')[1]
            if query_dom != result_dom:
                curr_query = query
                continue

            # Add hit to results
            results[query].add(result)

    return results


def eval_scores(fams: dict[str, list], results: dict[str, set]) -> dict[str, float]:
    """Returns dict of AUC1 scores for each query. AUC1 is calculated as the number of TP's
    up to the 1st FP divided by the number of sequences in the family.

    Args:
        path (str): Path to results file
        results (dict[str, set]): Dictionary of query PID's and a set of TP's up to the 1st FP.

    Returns:
        dict[str, float]: key: query PID, value: AUC1 score
    """
    
    # Read results and calculate AUC1 score for each query
    auc_scores: dict[str, float] = {}
    for query, tps in results.items():
        fam = query.split('|')[1]  # i.e. 16vpA00|3.30.930.10

        # Calculate AUC1 score
        try:
            score = len(tps) / (len(fams[fam])-1)  # ignore self hit
        except ZeroDivisionError:
            score = 0
        auc_scores[query] = auc_scores.get(query, 0) + score  # ignore self hit
    
    return auc_scores


def graph_results(scores: list[dict[str, float]], methods: list[str]):
    """Graphs AUC1 scores for each query.

    Args:
        scores (list[dict[str, float]]): List of dictionaries containing AUC1 scores.
        methods (list[str]): List of methods being evaluated
    """

    averages = [sum(sco.values()) / len(sco) for sco in scores]
    labels = [f'{m} (mean: {a:.2f})' for m, a in zip(methods, averages)]
    colors = ['blue', 'orange', 'red', 'green']
    _, ax = plt.subplots()
    for i, sco in enumerate(scores):
        y = range(len(sco))
        x = sorted(list(sco.values()), reverse=True)
        ax.plot(x, y, label=methods[i], color=colors[i])
    ax.set_xlabel('AUC1')
    ax.set_ylabel('Query')
    ax.legend(title='Search Tool', labels=labels)
    ax.set_title('AUC1 Scores for CATH20 Queries')
    plt.savefig('bench/cathdb/cath_auc1.pdf')


def main():
    """
    """

    path = 'bench/cathdb/data'
    fams = get_fams(f'{path}/cath20_queries.fa')

    # Read results after running each tool
    dct_res = read_results(f'{path}/results_dct.txt', 1, 5)
    mean_pt5_res = read_results(f'{path}/results_mean_pt5.txt', 1, 5)
    mean_esm2_res = read_results(f'{path}/results_mean_esm2.txt', 1, 5)
    mmseqs_res = read_results(f'{path}/results_mmseqs.txt', 0, 1)

    # Plot AUC1 scores for each query
    scores = []
    scores.append(eval_scores(fams, dct_res))
    scores.append(eval_scores(fams, mean_pt5_res))
    scores.append(eval_scores(fams, mmseqs_res))
    scores.append(eval_scores(fams, mean_esm2_res))
    methods = ['DCTdomain', 'ProtT5-Mean', 'MMseqs2-sens', 'ESM2-Mean']
    graph_results(scores, methods)


if __name__ == '__main__':
    main()
