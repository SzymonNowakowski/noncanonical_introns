from intron_comparison import process_file
from introns import Intron
from math import floor
from itertools import chain
from collections import defaultdict


def junctions_to_introns(file_in, file_out):
    """
    Take junctions created by RegTools and save to .bed file only the location and support of introns.

    :param file_in: (str) Path to the RegTools in file.
    :param file_out: (str) Path to the out file.
    """
    introns_dict = {}
    for line in process_file(file_in):
        scaffold = line[0]
        start = line[1] + int(line[-2].split(',')[0])
        end = line[2] - int(line[-2].split(',')[1])
        score = line[4]
        i = Intron(scaffold, start, end, support=score)

        if scaffold in introns_dict.keys():
            introns_dict[scaffold].append(i)
        else:
            introns_dict[scaffold] = [i]

    with open(file_out, 'w') as f_out:
        for scaff in introns_dict.keys():
            # introns in every scaffold are sorted by start for easier access
            int_list = sorted(introns_dict[scaff], key=lambda intron: intron.start)
            for intron in int_list:
                f_out.write(str(intron))
                f_out.write('\n')


def intron_stats(file):
    """Calculate basic statistics of introns from a file."""
    introns = []
    for line in process_file(file):
        i = Intron(line[0], line[1], line[2], support=line[3])
        introns.append(i)
    print(f'Number of introns: {len(introns)}')
    print(f'Mean support: {sum([i.support for i in introns]) / len(introns)}')

    def median(l):
        l.sort()
        return l[floor(len(l) / 2)]
    print(f'Median support: {median([intron.support for intron in introns])}')


def junctions_to_stats(program, m, path):
    print('Program: ', program)
    file_in = f'{path}junctions_{program}_m{m}.bed'
    file_out = f'{path}good_junctions_{program}.bed'
    junctions_to_introns(file_in, file_out)
    intron_stats(file_out)


def choose_best_introns(file_in, file_out, cutoff):
    """
    Choose one best intron over every position.

    :param file_in: (str) Path to the .bed file with introns in format: scaffold start end support. All introns from
    a scaffold must come one after another in the file, and within one scaffold introns have to be sorted by start.
    :param file_out: (str) Path to the out file with best introns.
    :param cutoff: (int) Minimum support of the best intron.
    :return: Two dictionaries where key is scaffold and value is the list of introns on the scaffold:
    one containing all the introns from the input file and one with the best introns.
    """
    with open(file_out, 'w') as f_out:

        best_introns = defaultdict(list)
        all_introns = defaultdict(list)

        chrom_old = 'scaffold_0'
        start_old = 0
        end_old = 0
        score_old = 0

        def write_junction():
            junction = '\t'.join([str(x) for x in [chrom, start_old, end_old, score_old]])
            f_out.write(junction)
            f_out.write('\n')
            best_introns[chrom].append(i)

        for line in process_file(file_in):
            chrom, start, end, score = line
            i = Intron(chrom, start, end, support=score)
            all_introns[chrom].append(i)
            # we only consider introns with high enough support
            if score < cutoff:
                continue

            if chrom == chrom_old:
                if start < end_old:
                    # we are still in the same intron
                    if score > score_old:
                        # we want one best intron in each position
                        start_old, end_old, score_old = start, end, score
                else:
                    # we are in a new intron, so we need to write down the old one
                    if not start_old - end_old == 0:
                        write_junction()
                    start_old, end_old, score_old = start, end, score

            else:
                # new scaffold, so we need to write down the last intron
                write_junction()
                chrom_old, start_old, end_old, score_old = chrom, start, end, score

        # now we need to write the last one
        write_junction()

        return all_introns, best_introns


def compare_best_introns(all_i, best_i):
    mean_supports = []
    for scaffold in best_i.keys():
        for best_intron in best_i[scaffold]:
            supports = []
            for intron in all_i[scaffold]:
                if best_intron.intersect(intron):
                    supports.append(intron.support)
            mean_supports.append(best_intron.support / sum(supports))
    print('Mean share of the best intron: ', sum(mean_supports) / len(mean_supports))

    sup_of_best = [intron.support for intron in chain.from_iterable(best_i.values())]
    sup_of_all = [intron.support for intron in chain.from_iterable(all_i.values())]
    print('Share of best introns in all: ', sum(sup_of_best) / sum(sup_of_all))
    print('Mean support of the best intron: ', sum(sup_of_all) / len(sup_of_best))
    print('\n')


def introns_for_seq(file_in, file_out, margin):
    """Create a file with introns elongated by a margin to extract sequences with parts of neighbouring exons."""
    with open(file_out, 'w') as f_out:
        for line in process_file(file_in):
            new_line = '\t'.join([str(x) for x in [line[0], line[1] - margin, line[2] + margin]])
            f_out.write(new_line)
            f_out.write('\n')


def main():
    m = 20
    path = '/home/julia/Documents/licencjat/'

    programs = ['hisat', 'star']
    # for program in programs:
    #     print('Program: ', program)
    #     file_in = f'{path}junctions_{program}_m{m}.bed'
    #     file_out = f'{path}good_junctions_{program}.bed'
    #     # junctions_to_introns(file_in, file_out)
    #     intron_stats(file_out)
    #     print('\n')

    # for program in programs:
    #     print('Program: ', program)
    #     cutoff = 1
    #
    #     file_in = f'{path}good_junctions_{program}.bed'
    #     file_out = f'{path}best_introns_{program}{cutoff}.bed'
    #     all_i, best_i = choose_best_introns(file_in, file_out, cutoff)
    #     intron_stats(file_out)
    #
    #     compare_best_introns(all_i, best_i)

    introns_for_seq(path + 'introns_hisat_50_cross_checked+3.bed', path + 'introns_hisat_50_cross_checked.bed', -3)
    # introns_for_seq(path + 'best_introns_hisat50.bed', path + 'best_introns_hisat_50+0.bed', 0)


if __name__ == "__main__":
    main()
