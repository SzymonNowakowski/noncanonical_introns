import sys  
sys.path.insert(0, '/scratch/introns/noncanonical_introns')

import kmeans
import random_forest
from Bio.SeqIO.FastaIO import SimpleFastaParser
from sklearn.cluster import KMeans
import numpy as np

#head -n 20000 selected_introns.fasta > subset.fasta
file = '/scratch/introns/subset.fasta'
with open(file, "r") as handle:
    sequences = list(SimpleFastaParser(handle))

    types = np.zeros((len(sequences),), dtype=int)
    for i, s in enumerate(sequences):
        #left_anchor = s[1][3:5]
        #right_anchor = s[1][-5:-3]
        #if left_anchor == 'GT' and right_anchor == 'AG' or left_anchor == 'CT' and right_anchor == 'AC':
        class_signature=s[0][-2:]
        if class_signature == "KX":    
                # 0  for conventional intron
            types[i] = (0)
        else:
                                                                    # 1 for nonconventional
            types[i] = (1)
print(sum(types), len(types))

_,_,seq4 = kmeans.preprocess(0, 4, sequences)
clusters4 = np.zeros((7, len(sequences)))
for x in range(1,8):
    n_clusters = 2**x
    kmeans_model = KMeans(n_clusters=n_clusters, init='k-means++', max_iter=100, n_init=10, random_state=0, n_jobs=12)
    kmeans_model.fit(seq4)
    clusters4[x - 1] = kmeans_model.labels_
print(clusters4)

_,_,seq7 = kmeans.preprocess(0, 7, sequences)
clusters7 = np.zeros((7, len(sequences)))
for x in range(1,8):
    n_clusters = 2**x
    kmeans_model = KMeans(n_clusters=n_clusters, init='k-means++', max_iter=100, n_init=10, random_state=0, n_jobs=12)
    kmeans_model.fit(seq4)
    clusters7[x - 1] = kmeans_model.labels_
print(clusters7)

def assess_cluster(clusters, n_of_clusters, true_classes):
    clus_dict = dict([(x, [0, 0]) for x in range(n_of_clusters)])
    for i, x in enumerate(clusters):
        true_type = types[i]
        clus_dict[int(x)][true_type] += 1    
        
    results = [[], []]
    for i, x in enumerate(types):
        cluster = clusters[i]
        res = clus_dict[cluster][x] / sum(clus_dict[cluster])
        results[int(x)].append(res)
    print('conv: ')
    print(sum(results[0]) / len(results[0]))
    print('nonconv: ')
    print(sum(results[1]) / len(results[1]))
    return clus_dict

for x in range(1,8):
    n_of_clusters =  2**x
    print('Word length: 4, number of clusters: %d' % (n_of_clusters))
    assess_cluster(clusters4[x-1], n_of_clusters, types)
    print('\n')
for x in range(1,8):
    n_of_clusters =  2**x
    print('Word length: 7, number of clusters: %d' % (n_of_clusters))
    assess_cluster(clusters7[x-1], n_of_clusters, types)
    print('\n') 

for n in range(4, 12):
    print('\nlen of ogligonucleotides: ', n)
    for cut in range(2):
        if cut == 0:
            print('uncut sequences: ')
        if cut == 1:
            print('conventional splices cut: ')
        data = random_forest.preprocess(cut, n, sequences)
        acc = random_forest.forest(10, data, types)
        print(str(acc))

###VISUALISATION OF WHAT IS HAPPENNING
sequences_cut, split_seqs, seq4 = kmeans.preprocess(0, 4, sequences)
print(split_seqs[1])
print(sequences_cut[1])
print(seq4)


