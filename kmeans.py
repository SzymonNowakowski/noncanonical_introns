from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans


def cut_junctions(cut, sequences):
    '''0 - dont cut, 1 - only conventional, 2 - both'''
    new_seqs = []
    for i in sequences:
        if cut == 0:
            new_seqs.append(i)
            continue
        left_anchor = i[1][3:5]
        right_anchor = i[1][-5:-3]
        if left_anchor == 'GT' and right_anchor == 'AG' or left_anchor == 'CT' and right_anchor == 'AC':
            new_seqs.append(i[1][5:-5])
        else:
            if cut == 1:
                new_seqs.append(i)
            else:
                new_seqs.append(i[1][8:-9])
    return new_seqs


def kmers(sequence, k):
    return [sequence[i: i + k] for i in range(len(sequence)-1)]


def split(k, seqs):
    return [" ".join(kmers(sequence, k)) for sequence in seqs]


def preprocess(cut, k, sequences, file_out):
    sequences_cut = cut_junctions(cut, [x[1] for x in sequences])
    split_seqs = split(k, sequences_cut)
    tfidf_vectorizer = TfidfVectorizer(use_idf=True, lowercase=False)
    tfidf_vectorizer_vectors = tfidf_vectorizer.fit_transform(split_seqs)
    print(type(tfidf_vectorizer_vectors))
    np.save(file_out, tfidf_vectorizer_vectors)

#
# def kmeans():
#

def main():
    file = '/home/rozycka/sequences_euglena/good_introns50+3.fasta'
    # file = '/home/julia/Documents/licencjat/good_introns50+3.fasta'
    with open(file, "r") as handle:
        sequences = list(SimpleFastaParser(handle))[:100]

    types = []
    for s in sequences:
        left_anchor = s[1][3:5]
        right_anchor = s[1][-5:-3]
        if left_anchor == 'GT' and right_anchor == 'AG' or left_anchor == 'CT' and right_anchor == 'AC':
            types.append(0)
        else:
            types.append(1)

    preprocess(0, 4, sequences, './seq_4_no_cut.npy')
    seq4 = np.load('./seq_4_no_cut.npy')
    for x in [2**x for x in range(2,3)]:
        kmeans = KMeans(n_clusters=x, init='k-means++', max_iter=50, n_init=1, random_state=0, n_jobs=8)
        kmeans.fit(seq4)
        print(kmeans.labels_)
        print(type(kmeans.labels_))


    # preprocess(0, 7, sequences, './seq_7_no_cut.npy')


if __name__ == "__main__":
    main()
