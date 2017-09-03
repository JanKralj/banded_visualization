__author__ = 'jank'
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm
import numpy as np
mpl.rcParams['xtick.labelsize'] = 8

# def plot_matrix(matrix, bands):
#     msh = plt.pcolormesh(matrix, cmap='Greys')
#     plt.axis([0, matrix.shape[1], 0, matrix.shape[0]])
#     plt.subplots_adjust(top=0.95, bottom=0.15)
#     plt.gca().set_ylabel('Cancer patients')
#     plt.gca().set_xlabel('Chromosome regions')
#     if False:
#         msh.axes.set_xticks([j + 0.5 for j in range(28)], minor=False)
#         msh.axes.set_xticks(range(28), minor=True)
#         msh.axes.set_xticklabels([item for item in bands], rotation=90)
#         plt.tick_params(axis='x', which='major', bottom='off', top='off', labelbottom='on')
#     plt.gca().invert_yaxis()
#     plt.show()


def plot_matrix(matrix, bands, types, s):
    plt.figure(figsize=(7, 15))
    use_map = matplotlib.cm.get_cmap('hsv')
    msh = plt.pcolormesh(matrix, cmap='Greys')
    plt.subplots_adjust(top=0.92, bottom=0.15)
    plt.axis([0, matrix.shape[1], 0, matrix.shape[0]])
    msh.axes.invert_yaxis()
    # blue red yellow green purple black
    plt.gca().set_ylabel('Cancer patients')
    plt.gca().set_xlabel('Chromosome regions')

    msh.axes.set_xticks([j + 0.5 for j in range(len(bands))], minor=False)
    if types:
        msh.axes.set_yticks([j + 0.5 for j in range(len(types))], minor=False)
        msh.axes.set_yticklabels([item for item in types])
    msh.axes.set_xticks(range(len(bands)), minor=True)
    msh.axes.set_xticklabels([item for item in bands], rotation=45)
    plt.tick_params(axis='x', which='major', bottom='off', top='off', labelbottom='on')
    # plt.show()
    # plt.close("all")
    plt.savefig('c:/SLO-pizza/basic_MBA_%i.png' % s, dpi=300)


def plot_matrix_clusters(matrix, bands, new_clusters, cc, ylabel, save_to=None):
    fig = plt.figure(figsize=(7, 8))
    n = len(set(new_clusters))
    use_map = matplotlib.cm.get_cmap('hsv')
    colour_map = colors.ListedColormap(
        [(1, 1, 1, 0)] + [tuple([(use_map(i * 1.0 / n)[j]) * 0.4 for j in range(3)] + [1]) for i in range(1, n + 1)])
    alt_map = colors.ListedColormap(
        [[(use_map(i * 1.0 / n)[j]) * 0.7 for j in range(3)] + [1] for i in range(1, n + 1)])

    draw_matrix = matrix * np.array(new_clusters)[:, np.newaxis]
    msh = plt.pcolormesh(draw_matrix, cmap=colour_map)
    plt.subplots_adjust(top=0.904)
    plt.axis([0, draw_matrix.shape[1], 0, draw_matrix.shape[0]])
    msh.axes.invert_yaxis()

    c_bar = plt.colorbar(msh, shrink=0.5)
    c_bar.set_cmap(alt_map)
    c_bar.draw_all()
    c_bar.set_ticks([i + 0.5 for i in range(n)])
    c_bar.ax.set_yticklabels(range(1, n + 1))
    c_bar.set_label('Cluster number')
    al = 0.3
    for i in range(n):
        for item in cc[i]:
            plt.axhspan(item, item + 1, alpha=al, lw=0, color=use_map((i + 1.0) / n))
    plt.gca().set_ylabel(ylabel)
    plt.gca().set_xlabel('Annotations')

    for i, x in enumerate(bands):
        if len(x.split('/')[-1].rstrip('0123456789')) > 20:
            bands[i] = 'LocatedEntity'
            # print x.split('/')[-1].rstrip('0123456789')

    msh.axes.set_xticks([j + 0.5 for j in range(len(bands))], minor=False)
    msh.axes.set_xticks(range(len(bands)), minor=True)
    msh.axes.set_xticklabels([item.split('/')[-1].rstrip('0123456789') for item in bands], rotation=90)
    plt.tick_params(axis='x', which='major', bottom='off', top='off', labelbottom='on')
    fig.autofmt_xdate(bottom=0.19, ha='right', rotation=60)
    # plt.show()
    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to, dpi=400)


def plot_matrix_clusters_transpose(matrix, bands, new_clusters, cc, s):
    fig = plt.figure(figsize=(20, 15))
    n = len(set(new_clusters))
    use_map = matplotlib.cm.get_cmap('hsv')
    colour_map = colors.ListedColormap(
        [(1, 1, 1, 0)] + [tuple([(use_map(i * 1.0 / n)[j]) * 0.4 for j in range(3)] + [1]) for i in range(1, n + 1)])
    alt_map = colors.ListedColormap(
        [[(use_map(i * 1.0 / n)[j]) * 0.7 for j in range(3)] + [1] for i in range(1, n + 1)])

    draw_matrix = matrix * np.array(new_clusters)[:, np.newaxis]
    msh = plt.pcolormesh(draw_matrix.transpose(), cmap=colour_map)
    plt.subplots_adjust(top=0.95, bottom=0.15)
    plt.axis([0, draw_matrix.shape[0], 0, draw_matrix.shape[1]])
    msh.axes.invert_yaxis()
    c_bar = plt.colorbar(msh, shrink=0.5)
    c_bar.set_cmap(alt_map)
    c_bar.draw_all()
    c_bar.set_ticks([i + 0.5 for i in range(n)])
    c_bar.ax.set_yticklabels(range(1, n + 1))
    c_bar.set_label('Cluster number')
    al = 0.3

    plt.draw()

    for i in range(n):
        for item in cc[i]:
            plt.axvspan(item, item + 1, alpha=al, lw=0, color=use_map((i + 1.0) / n))
    # blue red yellow green purple black
    plt.gca().set_xlabel('Cancer patients')
    plt.gca().set_ylabel('Chromosome regions')

    msh.axes.set_yticks([j + 0.5 for j in range(len(bands))], minor=False)
    msh.axes.set_yticks(range(len(bands)), minor=True)
    msh.axes.set_yticklabels([item for item in bands])
    plt.tick_params(axis='y', which='major', bottom='off', top='off', labelbottom='on')
    # plt.show()
    plt.savefig('c:/SLO-headlines/transposed_biMBA_%i.png' % s, dpi=300)


def plot_matrix_rules(matrix, bands, n, rr, cc, all_rules, save_to, ylabel):
    al = 0.2
    for i in range(n):
        fig = plt.figure(figsize=(7, 8))
        # fig.set_size_inches([8, 8], forward=True)

        msh = plt.pcolormesh(matrix, cmap='Greys')
        plt.subplots_adjust(top=0.904)
        plt.axis([0, matrix.shape[1], 0, matrix.shape[0]])
        msh.axes.invert_yaxis()
        bands_in_rules = list(set([item[1] for item in rr[i]]))
        d = dict([(x, []) for x in bands_in_rules])

        plt.draw()
        fc = msh.get_facecolors()
        fc_grid = fc.reshape(matrix.shape[0], matrix.shape[1], -1)

        for item in rr[i]:
            fc_grid[:, item[1]:item[1] + 1] = (1 - al) * fc_grid[:, item[1]:item[1] + 1] + al * np.array([0, 0, 1, 1])
            d[item[1]].append(str(item[2] + 1))
        if all_rules:
            for key in d.keys():
                plt.text(key + 0.4, -0.01 * matrix.shape[0], ','.join(sorted(list(set(d[key])))), rotation=90, va='bottom')
        for item in cc[i]:
            fc_grid[item:item + 1, :] = (1 - al) * fc_grid[item:item + 1, :] + al * np.array([0, 0, 0, 1])

        for iii, x in enumerate(bands):
            if len(x.split('/')[-1].rstrip('0123456789')) > 20:
                bands[iii] = 'LocatedEntity'
                # print x.split('/')[-1].rstrip('0123456789')
        msh.axes.set_ylabel(ylabel)
        msh.axes.set_xlabel('Annotations')
        msh.axes.set_xticks([j + 0.5 for j in range(len(bands))], minor=False)
        msh.axes.set_xticks(range(len(bands)), minor=True)
        msh.axes.set_xticklabels([item.split('/')[-1].rstrip('0123456789') for item in bands], rotation=90)
        msh.axes.set_xticklabels([item for item in bands], rotation=90)
        plt.tick_params(axis='x', which='major', bottom='off', top='off', labelbottom='on')
        fig.autofmt_xdate(bottom=0.19, ha='right', rotation=60)
        if save_to is None:
            plt.show()
        else:
            plt.savefig(save_to + '_cluster_%i.png' % i, dpi=400)
        plt.show()
        plt.close("all")


def plot_matrix_rules_transpose(matrix, bands, n, rr, cc, all_rules, s):
    al = 0.2
    for i in range(n):
        fig = plt.figure(figsize=(20, 15))
        # fig = plt.figure()
        fig.set_size_inches([8, 8], forward=True)
        msh = plt.pcolormesh(matrix.transpose(), cmap='Greys')
        plt.subplots_adjust(top=0.95, bottom=0.15)
        plt.axis([0, matrix.shape[0], 0, matrix.shape[1]])
        msh.axes.invert_yaxis()
        bands_in_rules = list(set([item[1] for item in rr[i]]))
        d = dict([(x, []) for x in bands_in_rules])

        plt.draw()
        fc = msh.get_facecolors()
        fc_grid = fc.reshape(matrix.shape[1], matrix.shape[0], -1)

        for item in rr[i]:
            fc_grid[item[1]:item[1] + 1, :] = (1 - al) * fc_grid[item[1]:item[1] + 1, :] + al * np.array([0, 0, 1, 1])
            d[item[1]].append(str(item[2] + 1))
        if all_rules:
            for key in d.keys():
                plt.text(4120, key + 0.4, ','.join(d[key]), va='bottom')
        for item in cc[i]:
            fc_grid[:, item:item + 1] = (1 - al) * fc_grid[:, item:item + 1] + al * np.array([0, 0, 0, 1])

        msh.axes.set_xlabel('Cancer patients')
        msh.axes.set_ylabel('Chromosome regions')
        # fig.autofmt_ydate(bottom=0.12, ha='center')

        msh.axes.set_yticks([j + 0.5 for j in range(len(bands))], minor=False)
        msh.axes.set_yticks(range(len(bands)), minor=True)
        msh.axes.set_yticklabels([item for item in bands])
        plt.tick_params(axis='y', which='major', bottom='off', top='off', labelbottom='on')
        plt.savefig('c:/SLO_LEUK/rules_transpose_biMBA_%i_%i.png' % (s, i + 1), dpi=300)
        plt.close("all")
        # plt.show()