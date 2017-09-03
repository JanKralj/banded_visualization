import argparse
import re

import matplotlib.pyplot as plt
import numpy as np

import core
try:
    from .helpers import load_data, load_rules, load_hierarchy
    from .visualization import plot_matrix_clusters, plot_matrix_rules
except ValueError:
    from helpers import load_data, load_rules, load_hierarchy
    from visualization import plot_matrix_clusters, plot_matrix_rules


def run():
    parser = argparse.ArgumentParser(description='Banded matrix visualization')
    parser.add_argument('mode', metavar='MODE', type=str, choices=['clusters', 'rules'], help='Choose mode (coloring clusters or displaying rules)')
    parser.add_argument('folder', metavar='FOLDER', type=str, help='Folder containing the data')
    parser.add_argument('-s', '--steps',
                        metavar='STEPS',
                        type=int,
                        help='Number of iterations of the algorithm',
                        default=100)
    parser.add_argument('-a', '--algorithm',
                        metavar='ALGORITHM',
                        type=str,
                        choices=['biMBA', 'barycentric', 'MBA'],
                        help='Number of iterations of the algorithm',
                        default='biMBA')

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Name of output file (leave blank to display in window)',
                        default=None)

    args = parser.parse_args()
    step = args.steps
    method = args.algorithm
    folder = args.folder
    y_label = folder.capitalize()
    save_to = args.output
    dna_f, band_f, clus_f, rules_f, hierarchy_f = [folder + '/data.txt',
                                                   folder + '/column_names.txt',
                                                   folder + '/row_labels.txt',
                                                   folder + '/rules.txt',
                                                   folder + '/hierarchy.txt']
    matrix, column_names, clusters = load_data(dna_f, band_f, clus_f)
    if args.mode == 'clusters':
        rows(matrix, column_names, clusters, y_label, method, step, save_to)
    else:
        hierarchy = load_hierarchy(hierarchy_f)
        new_hierarchy = {}
        for band in column_names:
            new_hierarchy[band] = hierarchy[band] if band in hierarchy else []
        rules = load_rules(rules_f, new_hierarchy, column_names, top_k_rules=5)
        columns(matrix, column_names, clusters, rules, save_to, y_label, step)


def rows(matrix, bands, clusters, ylabel, method, step, save_to):
    """run banded matrices on only chromosome 1 without permuting the columns
    (chromosome bands). Show clusters (each row is in one) with colors."""
    if method == 'biMBA':
        find_permutation = core.alternatingBi
    elif method == 'barycentric':
        find_permutation = core.baryCentric
    elif method == 'MBA':
        find_permutation = core.alternating
    else:
        raise Exception('Unknown method string: %s' % method)

    new_matrix, row_per, col_per = find_permutation(matrix, step)
    assert (np.linalg.norm(new_matrix - matrix[row_per, :][:, col_per]) == 0)
    new_bands = [bands[i] for i in col_per]

    new_clusters = [0 for _ in range(len(row_per))]
    for i in range(len(row_per)):
        new_clusters[i] = clusters[row_per[i]]
    cc = [[] for _ in range(len(set(clusters)))]
    for k in range(len(new_clusters)):
        cc[new_clusters[k] - 1].append(k)
    plot_matrix_clusters(new_matrix, new_bands, new_clusters, cc, ylabel, save_to=save_to)
    if save_to is not None:
        np.save(save_to + '.npy', new_matrix)


def columns(matrix, bands, clusters, rules, save_to, y_label, step):
    # method = core.baryCentric
    method = core.alternatingBi
    # method = core.alternating

    new_matrix, row_per, col_per = method(matrix, step)
    new_bands = [bands[i] for i in col_per]
    assert (np.linalg.norm(new_matrix - matrix[row_per, :][:, col_per]) == 0)

    # apply row permutation to clusters:----------------------------------------------------------------------------
    new_clusters = [0 for _ in range(len(row_per))]
    for i in range(len(row_per)):
        new_clusters[i] = clusters[row_per[i]]
    cc = [[] for _ in range(len(set(clusters)))]
    for key in range(len(new_clusters)):
        cc[new_clusters[key] - 1].append(key)

    # apply column permutation to rules:------------------------------------
    n = len(set(clusters))
    new_rules = dict([(key, []) for key in rules])
    for key in rules.keys():
        new_rules[key] = [(x[0], col_per.index(x[1]), x[2]) for x in rules[key]]
    rules = new_rules
    rr = [[] for _ in range(n)]
    for key in rules.keys():
        iiiii = int(re.match(r'.*(?P<nn>\d+).*', key).group('nn')) - 1
        rr[iiiii] = [(x[0], x[1], x[2]) for x in rules[key]]
    plot_matrix_rules(new_matrix, new_bands, n, rr, cc, True, save_to, y_label)


def chromosome_all(matrix, bands):
    hotspots = []
    with open('hotspot_sorted_calc.txt') as hf:
        for line in hf:
            hotspots.append(line.strip())
    virus_sites = []
    with open('virus_int_sites_calc.txt') as vf:
        for line in vf:
            virus_sites.append(line.strip())

    bar_res = core.baryCentric(matrix, 200, clusters=None)
    band_perm = bar_res[2]
    hotspot_indices = [i for i in range(len(bands)) if bands[band_perm[i]] in hotspots]
    virus_indices = [i for i in range(len(bands)) if bands[band_perm[i]] in virus_sites]
    for item in virus_indices:
        plt.axvspan(item, item + 1, alpha=0.5, lw=0, color='red')
    plt.pcolormesh(bar_res[0], cmap='Greys')
    plt.gca().invert_yaxis()
    plt.show()


def perm_multiply(perm1, perm2):
    """returns perm1(perm2), perm1 applied to the result of perm2"""
    return [perm1[x] for x in perm2]


def inverse(perm):
    inv = [0 for i in range(len(perm))]
    for i in range(len(perm)):
        inv[perm[i]] = i
    return inv


if __name__ == '__main__':
    run()