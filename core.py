import numpy as np
import matplotlib.pyplot as plt
import random
#import pdb

def MBAStep(matrix, nRows, nCols):
    As =[]
    Bs = []
    final = []
    for row in matrix:
        first = 0; last = nCols
        # find first one:
        while first < nCols and row[first] == 0:
            first += 1
        if first == nCols:
            # matrices with all-zero rows:
            As.append(-1)
            Bs.append(-1)
            continue
        # find last one:
        while row[last - 1] == 0:
            last -= 1
        As.append(first)
        Bs.append(last)
    # print "Second for loop:"
    for i in range(nRows):
        # print i
        a = As[i]
        b = Bs[i]
        C = []
        # check each row if this one is included in it:
        for j in range(nRows):
            if a > As[j] and b < Bs[j]:
                # j-th row is in C:
                C.append(j)
        if len(C) == 0:
            best = [a, b]
        else:
            # A-Case candidate: [min{a_j, [a_j, b_j] in C}, b]
            best = [min([As[j] for j in C]), b]
            # B-Case candidate: [a, max{b_j, [a_j, b_j] in C}]
            candidate = [a, max([Bs[j] for j in C])]
            if candidate[1] - candidate[0] < best[1] - best[0]:
                best = candidate
            # C-Case candidate:
            for j in C:
                # for each [a_j, b_j] in C, candidate is: [a_j, max{b_k, [a_k, b_k] in C and a_k < a_j}]
                cc = [Bs[k] for k in C if As[k] < As[j]]
                if len(cc) > 0:
                    candidate = [As[j], max(cc)]
                    if candidate[1] - candidate[0] < best[1] - best[0]:
                        best = candidate
        final.append(best)
    return sorted(range(len(final)), key=lambda _i: final[_i])

def BiMBAStep(matrix, nRows, nCols):
    As =[]; Bs = []; final = [];
    for row in matrix:
        ms = maxSubArray(row)
        As.append(ms[0])
        Bs.append(ms[1])
    for i in range(nRows):
        a = As[i]; b = Bs[i];
        if a == -1:
            pass
        # looking at (i,j) pairs: (random order here will find multiple paths)
        for j in range(nRows):
            # since b != -1, this also checks if Bs[j] != -1:
            if As[j] < a and Bs[j] > b:
                # print "%i is subinterval of %i" % (i, j)
                # [a, b] < [aj, bj]. Solve consecutive ones on array A of tipe [zeros], ones, zeros, ones, [zeros]
                # some randomization here might be a good idea?
                firstOnes = a - As[j]
                secondOnes = Bs[j] - b
                midZeros = b - a + 1
                m = min(firstOnes, secondOnes, midZeros)
                if m == firstOnes:
                    # delete 1. row of ones in A which becomes [b+1, bj]. This means adding ones to [a,b] which becomes [aj, b]
                    a = As[j]
                elif m == secondOnes:
                    # delete 2. row of ones in A which becomes [aj, a-1]. This measn adding ones to [a,b] which becomes [a, bj]
                    b = Bs[j]
                else:
                    a = -1; b = -1;
                    break
        final.append([a, b])
    return sorted(range(nRows), key = lambda i: final[i]) # [[As[i], Bs[i]] for i in range(len(As))]

def alternating(matrix, steps):
    rowper = range(matrix.shape[0])
    colper = range(matrix.shape[1])
    stateDict = {False: rowper, True: colper}
    state = False
    for i in range(steps):
        # print "step %i" % i
        newPerm = MBAStep(matrix, matrix.shape[0], matrix.shape[1])
        if newPerm == range(len(stateDict[state])):
            print "optimal permutation achieved at i = %i" % i
            break
        stateDict[state] = [stateDict[state][x] for x in newPerm]
        newmatrix = matrix[newPerm, :]
        matrix = newmatrix.transpose()
        state = not state
    if state:
        matrix = matrix.transpose()
    return [matrix, stateDict[False], stateDict[True]]

def alternatingBi(matrix, steps):
    rowper = range(matrix.shape[0])
    colper = range(matrix.shape[1])
    stateDict = {False: rowper, True: colper}
    state = False
    for i in range(steps):
        print "step %i" % i
        newPerm = BiMBAStep(matrix, matrix.shape[0], matrix.shape[1])
        if newPerm == range(len(stateDict[state])):
            print "optimal permutation achieved at i = %i" % i
            break
        stateDict[state] = [stateDict[state][x] for x in newPerm]
        newmatrix = matrix[newPerm, :]
        matrix = newmatrix.transpose()
        state = not state
    if state:
        matrix = matrix.transpose()
    return matrix, stateDict[False], stateDict[True]


def baryCentric(mat, steps, clusters = None):
    #mat.shape = (nRows, nColumns)
    ones0 = np.ones(mat.shape[1])
    ones1 = np.ones(mat.shape[0])
    range0 = np.arange(mat.shape[1])
    range1 = np.arange(mat.shape[0])
    rowper = range(mat.shape[0])
    colper = range(mat.shape[1])

    stateDict = {False: [ones0, range0, rowper], True: [ones1, range1, colper]}
    state = False # tells if matrix is transposed
    for i in range(steps):
        num = np.dot(mat, stateDict[state][1])
        den = np.dot(mat, stateDict[state][0])
        baryCenters = num / (den + 0.00001)

        if (clusters is None) or (not state):
            permutation = sorted(range(baryCenters.shape[0]), key = lambda i: [baryCenters[i], den[i]])
        else:
            permutation = sorted(range(baryCenters.shape[0]), key = lambda i: [baryCenters[i], den[i], clusters[stateDict[state][2][i]]])
        #multiply stateDict[state][2] with permutation (from left)
        stateDict[state][2] = [stateDict[state][2][x] for x in permutation]
        if permutation == range(baryCenters.shape[0]):
            print "barycenters alined at i = %i" % i
            break
        mat = mat[permutation, :].transpose()
        state = not state
    if state:
        mat = mat.transpose()
    return [mat, stateDict[False][2], stateDict[True][2]]


def maxSubArray(list):
    currStart = 0; currBest = 0; optStart = 0; optEnd = 0; optimal = 0
    for i in range(len(list)):
        if list[i] == 1:
            currBest = currBest + 1
            if currBest > optimal:
                optStart = currStart
                optEnd = i
                optimal = currBest
        elif currBest > 1:
            currBest = currBest - 1
        else:
            currBest = 0
            currStart = i + 1
    return (optStart, optEnd + 1)


def matrixFromIntervals(intervals, nRows, nCols, permuteRows = False):
    if permuteRows:
        intervals.sort()
    matrix = [[0 for i in range(nCols)] for i in range(nRows)]
    for i, interval in enumerate(intervals):
        if interval[0] >= 0:
            for j in range(interval[0], interval[1]):
                matrix[i][j] = 1
    return matrix