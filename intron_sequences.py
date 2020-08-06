from introns import Intron
import re

path = '/home/julia/Documents/licencjat/'
introns_seq_file = path + 'good_introns50+3.fasta'


def file_to_seq_introns(file, margin):
    introns = []
    with open(file) as f_in:
        for line in f_in.readlines():
            if line[0] == '>':
                line = re.split('>|:|-', line.strip())
                scaffold, start, end = line[1:]
            else:
                i = Intron(scaffold, start, end, margin=margin, sequence=line.strip())
                introns.append(i)
    return introns


def seq_statistics(introns, k):
    length = sum([intron.length() - 2 * k for intron in introns]) / len(introns)
    print('Numberof introns: ', len(introns))
    print('Mean length: ', length)

    def gc_content(intron):
        gc = 0
        for nukleotide in intron.sequence:
            if nukleotide in 'GC':
                gc += 1
        return gc / (intron.length() - 2 * k)
    gc = sum([gc_content(intron) for intron in introns]) / len(introns)
    print('Mean gc: ', gc)

    def extract_junction(intron):
        return intron.sequence[k:k+2] + intron.sequence[-k-2:-k]
    junctions = {}
    for intron in introns:
        junction = extract_junction(intron)
        if junction in junctions.keys():
            junctions[junction] += 1
        else:
            junctions[junction] = 1
    top10 = [(j, v) for j, v in sorted(junctions.items(), key=lambda item: item[1], reverse=True)][:10]
    print('Top junctions: ', top10)


def main():
    introns = file_to_seq_introns(introns_seq_file, 3)
    conv = []
    non_conv = []
    for intron in introns:
        if intron.check_conventional():
            conv.append(intron)
        else:
            non_conv.append(intron)

    print('conv stats: ')
    seq_statistics(conv, 3)
    print('nonconv stats: ')
    seq_statistics(non_conv, 3)


if __name__ == "__main__":
    main()
