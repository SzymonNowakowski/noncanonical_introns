from introns import Intron


def process_file(file_path):
    try:
        with open(file_path) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                yield [int(x) if x.isnumeric() else x for x in line.split()]
    except (IOError, OSError):
        print("Error opening / processing file")


def intron_from_line(line):
    try:
        if isinstance(line[3], str) and line[3] in '+-':
            return Intron(line[0], line[1], line[2], line[3])
        else:
            return Intron(line[0], line[1], line[2])
    except IndexError:
        return Intron(line[0], line[1], line[2])


def extract_introns_from_gtf(file, file_out):
    introns_p = []
    unique = set()
    for line in process_file(file):
        if line[2] == 'transcript':
            new_gene = True
        elif line[2] == 'exon':
            if new_gene:
                new_gene = False
                start = line[4]
            else:
                end = line[3] - 1
                scaffold = line[0]
                sign = line[6]
                i = Intron(scaffold, start, end, sign)
                if ' '.join([scaffold, str(start), str(end)]) not in unique:
                    unique.add(' '.join([scaffold, str(start), str(end)]))
                    introns_p.append(i)
                start = line[4]

    with open(file_out, 'w') as f_out:
        for intron in introns_p:
            f_out.write(str(intron))
            f_out.write('\n')


def intron_dict(file):
    my_introns = {}
    for line in process_file(file):
        intron = intron_from_line(line)
        if intron.scaffold not in my_introns.keys():
            my_introns[intron.scaffold] = [intron]
        else:
            my_introns[intron.scaffold].append(intron)
    return my_introns


def compare_introns(file, introns2, file_out):
    '''to be improved'''
    all_my_introns = 0
    intron_pairs = []
    no_pair = []
    diff_lengths = []
    for line in process_file(file):
        all_my_introns += 1
        intron = intron_from_line(line)
        pairing = False
        try:
            for i in introns2[intron.scaffold]:
                if intron.intersect(i):
                    intron_pairs.append([intron, i])
                    pairing = True

            if not pairing:
                no_pair.append(intron)
        except KeyError:
            no_pair.append(intron)

    print('All my introns: ', all_my_introns)
    print('n of pairs: ', len(intron_pairs))
    print('no match: ', len(no_pair))

    exact_match = []
    one_side_match = []
    with open(file_out, 'w') as f_out:
        for i1, i2 in intron_pairs:
            intersection = i1.classify_intersection(i2)
            if intersection == 'same':
                exact_match.append([i1, i2])
                f_out.write('\t'.join([i1.scaffold, str(i1.start - 3), str(i1.end + 3)]))
                f_out.write('\n')
            elif intersection == 'one end':
                one_side_match.append([i1, i2])
                diff_lengths.append(abs(i1.length() - i2.length()))
    return exact_match, one_side_match, diff_lengths


def main():
    path = '/home/julia/Documents/licencjat/'
    extract_introns_from_gtf(path + 'longa_stringtie_strand_informed.gtf', path + 'introns_pawel.bed')
    introns_pawel = path + 'introns_pawel.bed'
    p_introns = intron_dict(introns_pawel)

    file = path + 'best_introns_hisat50.bed'
    file_out = path + 'introns_hisat_50_cross_checked+3.bed'
    exact_match, one_side_match, diff_lengths = compare_introns(file, p_introns, file_out)

    print('exact: ', len(exact_match))
    # print(exact_match)

    print('one side: ', len(one_side_match))
    # print(one_side_match)

    print([(x, diff_lengths.count(x)) for x in range(10)])


if __name__ == '__main__':
    main()
