from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.model_selection import train_test_split


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


def preprocess(cut, k, sequences):
    sequences_cut = cut_junctions(cut, [x[1] for x in sequences])
    split_seqs = split(k, sequences_cut)
    tfidf_vectorizer = TfidfVectorizer(use_idf=True, lowercase=False)
    tfidf_vectorizer_vectors = tfidf_vectorizer.fit_transform(split_seqs)
    return tfidf_vectorizer_vectors


def forest(validation, x, types):
    y = np.array(types)
    acc = []
    for i in range(validation):
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3)
        clf = RandomForestClassifier(n_estimators=10,n_jobs=8)
        clf.fit(x_train, y_train)
        y_pred = clf.predict(x_test)
        acc.append(metrics.accuracy_score(y_test, y_pred))
    return sum(acc) / len(acc)
