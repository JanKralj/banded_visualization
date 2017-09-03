import re
import numpy as np
from operator import itemgetter
from collections import defaultdict

NORMAL = 'normal'
COMMON = 'common'
RARE = 'rare'


def load_data(amps_f, bands_f, clusters_f, **kwargs):
    matrix, bands, clusters = load_data_basic(amps_f, bands_f, clusters_f)
    matrix, bands, clusters = filter_data(matrix, bands, clusters, **kwargs)
    return matrix, bands, clusters


def filter_data(matrix, bands, clusters):
    matrix, bands, clusters = filter_rows(matrix, bands, clusters, 0)
    matrix, bands, clusters = filter_rows(matrix, bands, clusters, 1)
    return matrix, bands, clusters


def load_data_basic(amps_f, bands_f, clusters_f):
    matrix = np.loadtxt(amps_f, dtype=int)
    bands = []
    clusters = []
    if bands_f:
        with open(bands_f) as bands_file:
            for line in bands_file:
                bands.append(line.split()[1])
    if clusters_f:
        with open(clusters_f) as clusters_file:
            for line in clusters_file:
                d = [int(x) for x in line.split()]
                clusters.append(d[1])
    return matrix, bands, clusters


def expand_hierarchy(parent_dict):
    to_return = {}
    for child in parent_dict.keys():
        # child = 'http://dbpedia.org/resource/Category:Sociology'
        tocheck = [(child, None)]
        currparents = []
        i=0
        while len(tocheck) > 0:
            i+=1
            # print tocheck[0]
            currchild = tocheck.pop(0)[0]
            currparents = currparents + parent_dict[currchild]
            tocheck = tocheck + [(x, currchild) for x in parent_dict[currchild]]
        to_return[child] = currparents
    return to_return



def load_hierarchy(hierarchy_f):
    with open(hierarchy_f) as f:
        parent_dict = defaultdict(list)
        for line in f:
            l = line.split('\t')
            parent = l[0]
            l = l[1].strip('\n').split(';')
            for child in l:
                parent_dict[child].append(parent)
    return expand_hierarchy(parent_dict)


def load_rules(rules_f, hierarchy, bands, top_k_rules=5):
    rules = {}
    with open(rules_f) as af:
        state = 0
        # go through text between ----
        for line in af:
            if line[:3] == '---':
                state += 1
                if state == 2:
                    break
        curr_key = None
        rule_number = 0
        for line in af:
            if line.split()[1] == '<--':
                rule_number = 0
                curr_key = line.split()[0]
                rules[curr_key] = []
            else:
                if rule_number < top_k_rules:
                    add_rule(rules[curr_key], hierarchy, line, bands, rule_number)
                rule_number += 1
    return rules


# -------AGGREGATE FUNCTIONS--------------------------------------------------------------------------------------------


def dist_types(types):
    return list(set(types))


def class_frequencies(types):
    r_dict = {}
    for tt in types:
        if tt in r_dict:
            r_dict[tt] += 1
        else:
            r_dict[tt] = 1
    return r_dict


def cumulate(types, freqdict=None):
    if freqdict is None:
        freqdict = class_frequencies(types)
    frequencies = sorted(freqdict.items(), key=itemgetter(1), reverse=True)
    cumulative_frequencies = []
    s = 0
    for i in range(len(frequencies)):
        s = s + frequencies[i][1]
        cumulative_frequencies.append((frequencies[i][0], s))
    return cumulative_frequencies


def calc_weights(types):
    freqdict = class_frequencies(types)
    for key in freqdict.keys():
        freqdict[key] = 1.0 / freqdict[key]
    return freqdict

# -------FILTER FUNCTIONS----------------------------------------------------------------------------------------


def filter_chromosome(matrix, bands, chromosome):
    """Filters data to keep only chromosome regions from "bands" which are part of chromosome "chr".
    Returns a tuple of 2 objects. The first is a filtered "matrix" with only relevant columns left. The second
    is a list of chromosome regions represented by the new matrix columns."""
    regex = str(chromosome) + r'(p|q).*'
    to_keep = [i for i in range(len(bands)) if re.match(regex, bands[i])]
    return matrix[:, to_keep], [bands[i] for i in to_keep]


def filter_rows(matrix, bands, clusters, value):
    """Filters out all matrix rows in which all elements equal "value".
    Returns a tuple of 2 objects. The first is a filtered "matrix".
    The second is a list of cancer types represented by the new matrix rows."""
    z_rows = np.all(matrix == value, axis=1)
    to_keep = [i for i in range(matrix.shape[0]) if not z_rows[i]]
    return_matrix = matrix[to_keep, :]
    return_clusters = clusters
    return return_matrix, bands, return_clusters


def fragiles(bands_list, sites_f):
    bands = {}
    for band in bands_list:
        bands[band] = NORMAL
    with open(sites_f) as f:
        for line in f:
            bands[line.split()[1]] = line.split()[0]
    return bands


def add_rule(rule_list, hierarchy, rule, bands, rule_number):
    full_rule = re.sub(r'\[.*\]', '', rule.strip())
    for rule in full_rule.split('), '):
        pos_rule = not rule.startswith('~')
        rightlist = rule.split('annotated_with(X, ')[1].rstrip(' )')
        for right in rightlist.split(','):
            if right in bands:
                rule_list.append((pos_rule, bands.index(right), rule_number))
            else:
                rule_list += [(pos_rule, i, rule_number) for i in range(len(bands)) if right in hierarchy[bands[i]]]
        #
        # split = re.match(r'(?P<neg>~?)(?P<arm>\d*[pq])(?P<from>\d*)(-(?P<to>\d*))?', rule)
        # split = re.match(r'(?P<neg>~?)(?P<arm>\d*[pq])(?P<from>\d*)(-(?P<to>\d*))?', rule)
        # if split is None:
        #     continue
        # pos_rule = len(split.group('neg')) == 0
        # if split.group('to') is None:
        #     to_add = split.group('arm') + split.group('from')
        #     if to_add in bands:
        #         rule_list.append((pos_rule, bands.index(to_add), rule_number))
        #     else:
        #         rule_list += [(pos_rule, i, rule_number) for i in range(len(bands)) if bands[i].startswith(to_add)]
        #     continue
        # fr = int(split.group('from'))
        # to = int(split.group('to'))
        # for i in range(fr, to + 1):
        #     rule_list.append((pos_rule, bands.index(split.group('arm') + str(i)), rule_number))


def add_rule_2014(rule_list, rule, bands, rule_number):
    full_rule = re.sub(r'\[.*\]', '', rule)
    for rule in full_rule.split():
        split = re.match(r'(?P<neg>~?)(?P<arm>\d*[pq])(?P<from>\d*)(-(?P<to>\d*))?', rule)
        if split is None:
            continue
        pos_rule = len(split.group('neg')) == 0
        if split.group('to') is None:
            to_add = split.group('arm') + split.group('from')
            if to_add in bands:
                rule_list.append((pos_rule, bands.index(to_add), rule_number))
            else:
                rule_list += [(pos_rule, i, rule_number) for i in range(len(bands)) if bands[i].startswith(to_add)]
            continue
        fr = int(split.group('from'))
        to = int(split.group('to'))
        for i in range(fr, to + 1):
            rule_list.append((pos_rule, bands.index(split.group('arm') + str(i)), rule_number))