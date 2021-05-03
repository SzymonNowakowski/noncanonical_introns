from collections import defaultdict


def getter_setter_gen(name, type_):
    def getter(self):
        return getattr(self, "__" + name)

    def setter(self, value):
        if not isinstance(value, type_) and value is not None:
            raise TypeError("%s attribute must be set to an instance of %s" % (name, type_))
        setattr(self, "__" + name, value)
    return property(getter, setter)


def auto_attr_check(cls):
    new_dct = {}
    for key, value in cls.__dict__.items():
        if isinstance(value, type):
            value = getter_setter_gen(key, value)
        new_dct[key] = value
    # Creates a new class, using the modified dictionary as the class dict:
    return type(cls)(cls.__name__, cls.__bases__, new_dct)


class GenomicSequence:
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None):
        self.scaffold_name = scaffold_name
        self.scaffold_start = scaffold_start
        self.scaffold_end = scaffold_end
        if sequence and len(sequence) != self.length():
            raise ValueError('Incorrect sequence length.')
        self.sequence = sequence
        if strand and strand not in ['+', '-', '.']:
            raise ValueError('Strand can only be +, - or .')
        else:
            self.strand = strand

    def length(self):
        return self.scaffold_end - self.scaffold_start
    
    def __repr__(self):
        return ' '.join([self.scaffold_name, str(self.scaffold_start), str(self.scaffold_end)])
        
    def __str__(self):
        return self.__repr__()


class Gene(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='', transcript=None, exons=None, introns=None,
                 name=''):
        # if strand == '-':
        #     start, end = end, start
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.transcript = transcript
        self.exons = [] if exons is None else exons
        self.introns = [] if introns is None else introns
        self.name = name
        self.expanded_sequence = ''
        self.expansion_left = 0
        self.expansion_right = 0

    def append_exons(self, exon):
        self.exons.append(exon)
    
    def append_introns(self, intron):
        self.introns.append(intron)
    
    def extract_sequence(self, genome):
        # elif self.strand == '-':
        #     sequence = genome[self.scaffold_name][self.end:self.start]
        # else:
        #     raise Exception('cojest')
        scaffold_seq = genome[self.scaffold_name]
        sequence = scaffold_seq[self.scaffold_start:self.scaffold_end]
        expanded_sequence, expansion_left, expansion_right = self.get_expanded_sequence(scaffold_seq)
        if self.strand == '-':
            self.sequence = reverse_complement(sequence)
            self.expanded_sequence = reverse_complement(expanded_sequence)
            self.expansion_left = expansion_right
            self.expansion_right = expansion_left
        else:
            self.sequence = sequence
            self.expanded_sequence = expanded_sequence
            self.expansion_right = expansion_right
            self.expansion_left = expansion_left
        for exon in self.exons:
            if self.strand == '+':
                exon.sequence = self.sequence[exon.scaffold_start - self.scaffold_start:exon.scaffold_end - self.scaffold_start]
            elif self.strand == '-':
                exon.sequence = self.sequence[- exon.scaffold_end + self.scaffold_end:- exon.scaffold_start + self.scaffold_end]
            else:
                print(self.name, exon.scaffold_name, exon.scaffold_start, exon.scaffold_end)
                # raise Exception('co jest')
        transcript_sequence = self.get_transcript_sequence()
        self.transcript = Transcript(self.scaffold_name, self.scaffold_start, self.scaffold_end, strand=self.strand,
                                     sequence=transcript_sequence)

    def get_expanded_sequence(self, scaffold_seq):
        start = max(0, self.scaffold_start - 500)
        end = min(len(scaffold_seq), self.scaffold_end + 500)
        expanded_sequence = scaffold_seq[start:end]
        expansion_left = self.scaffold_start - start
        expansion_right = end - self.scaffold_end
        return expanded_sequence, expansion_left, expansion_right

    def get_transcript_sequence(self):
        exons_seqs = []
        for exon in self.exons:
            exons_seqs.append((exon.sequence, exon.scaffold_start))
        if self.strand == '+':
            exons_seqs.sort(key=lambda tup: tup[1])
        elif self.strand == '-':
            exons_seqs.sort(key=lambda tup: tup[1], reverse=True)
        sequence = ''.join([exon_seq[0] for exon_seq in exons_seqs])
        return sequence

    def get_transcript_with_gaps_sequence(self, expanded=False):
        to_be_joined = []
        start, end = None, None
        reverse = True if self.strand == '-' else False
        exons_sorted = sorted(self.exons, key=lambda obj: obj.start, reverse=reverse)
        if not reverse:
            for exon in exons_sorted:
                start = exon.start
                if end:
                    intron_length = start - end
                    assert intron_length > -1
                    to_be_joined.append((''.join(['-' for i in range(intron_length)])))
                to_be_joined.append(exon.sequence)
                end = exon.end
        else:
            for exon in exons_sorted:
                start = exon.end
                if end:
                    intron_length = end - start
                    assert intron_length > -1
                    to_be_joined.append((''.join(['-' for i in range(intron_length)])))
                to_be_joined.append(exon.sequence)
                end = exon.start
        sequence = ''.join(to_be_joined)
        if expanded:
            sequence = self.expansion_left * '-' + sequence + self.expansion_right * '-'
        return sequence

    def create_introns(self):
        if len(self.exons) < 1:
            raise ValueError('No exons specified.')
            return
        else:
            self.introns = []
            start, end = 0, 0
            for exon in self.exons:
                prev_e = exon.prev_exon
                end = exon.start - 1
                if start:
                    if self.strand == '+':
                        sequence = self.sequence[start - self.scaffold_start:end - self.scaffold_start]
                    elif self.strand == '-':
                        sequence = self.sequence[- end + self.scaffold_end: - start + self.scaffold_end]
                    int=Intron(self.scaffold_name, start, end, strand=self.strand, sequence=sequence, gene=self, prev_exon=prev_e, next_exon=exon)
                    self.append_introns(int)
                    
                    prev_e.next_intron = int
                    exon.prev_intron=int

                    prev_e=exon
                    # elif self.strand == '-':
                    #     self.append_introns(Intron(self.scaffold_name, end, start))
                    # else:
                    #     raise Exception('co jest')
                start = exon.end
        for intron in self.introns:
            intron.movable_boundary_no_margins()
        
        # for intron in self.introns:
        #     if self.strand == '-':
        #         print(str(intron), intron.start, self.start, intron.end, self.start, self.strand, intron.strand)
        #         intron.sequence = self.sequence[intron.start - self.start:intron.end - self.start]
        #         print(self.sequence)
        #         print(intron.sequence)
        #     elif self.strand == '+':
        #         print(str(intron), intron.end, self.end, intron.start, self.end, self.strand, intron.strand)
        #         intron.sequence = self.sequence[- intron.end + self.end:- intron.start + self.end]
        #         print(self.sequence)
        #         print(intron.sequence)


@auto_attr_check
class Intron(GenomicSequence):
    """
    This is a class for working with biological introns.

    :param scaffold_name: (str) Name of scaffold or chromosome on which the intron is located.
    :param start: (int) Location of the first base of the intron (within scaffold).
    :param end: (int) Location of the first base of the next exon (within scaffold).
    :param start: (int) Location of the first base of the intron (within gene).
    :param end: (int) Location of the first base of the next exon (within gene).
    :param gene: (str) Gene in which the intron is located.
    :param strand (str): Optional, defines the strand on which intron is located. Must be either + or -.
    :param support: (int) Optional, how many reads support the intron.
    :param margin_left: (int) Optional, how many nucleotides from the preceding exon are included.
    :param margin_right: (int) Optional, how many nucleotides from the following exon are included.
    :param margin_left: (str) Optional, end sequence from the preceding exon.
    :param margin_right: (str) Optional, beginning sequence of the following exon.
    :param sequence: (str) Optional, genomic sequence of the intron.
    """
    scaffold_name = str
    scaffold_start = int
    scaffold_end = int
    gene_end = int
    scaffold_start = int
    gene_start = int
    scaffold_end = int
    gene = Gene
    strand = str
    support = int
    margin_left = int
    margin_right = int
    sequence = str
    margin_left_seq = str
    margin_right_seq = str
    is_conventional = int
    is_nonconventional = int

    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None, gene=None, support=None, margin_left=0,
                 margin_right=0, margin_left_seq='', margin_right_seq='', prev_exon=None, next_exon=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.gene = gene
        self.support = support
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_left_seq = margin_left_seq
        self.margin_right_seq = margin_right_seq
        self.variations = []
        self.is_conventional = 0
        self.is_nonconventional = 0
        self.prev_exon = prev_exon
        self.next_exon = next_exon
        
        if self.strand=='-':
            self.gene_start, self.gene_end = -self.scaffold_end + self.gene.scaffold_end, -self.scaffold_start + self.gene.scaffold_end
        elif self.strand=='+':
            self.gene_start, self.gene_end = self.scaffold_start - self.gene.scaffold_start, self.scaffold_end - self.gene.scaffold_start

    # def intersect(self, other_intron):
    #     """
    #     Checks if two introns cover partially the same area.
    #
    #     :param other_intron: (Intron) Intron with which intersection of self is checked.
    #     :return: Bool: True if introns share at least one base, False otherwise.
    #     """
    #     if not isinstance(other_intron, Intron):
    #         raise ValueError('Both introns must be instances of Intron class')
    #
    #     else:
    #         if self.scaffold_name != other_intron.scaffold_name:
    #             raise ValueError('Introns on different scaffold cannot intersect.')
    #         if (self.scaffold_start in range(other_intron.scaffold_start, other_intron.scaffold_end))\
    #                 or (self.scaffold_end in range(other_intron.start, other_intron.end))\
    #                 or (self.scaffold_start < other_intron.scaffold_start and self.scaffold_end > other_intron.scaffold_end):
    #             return True
    #         else:
    #             return False

    # def classify_intersection(self, intron):
    #     """
    #     For two intersecting introns return the type of intersection
    #
    #     :param intron: (Intron) Intron with which type of intersection with self is checked.
    #     :return: 2 if the introns are identical, 1 if only either start or end is the same, 0 otherwise.
    #     """
    #     if self.scaffold_name != intron.scaffold_name:
    #         raise ValueError('Introns do not intersect at all.')
    #     if self.start == intron.start and self.end == intron.end:
    #         return 2
    #     elif self.start == intron.start or self.end == intron.end:
    #         return 1
    #     else:
            return 0

    def movable_boundary(self):
        # moglem zepsuc
        # ja zaraz jeszcze bardziej popsuje
        """
        Check if there are repeats on the intron junctions so the intron position could be shifted without changing
        transcript sequence. If there are, add new possible introns to self.variations.
        """
        i = 1
        # start checking for repeats left from the junction
        check = 'left'

        while True:
            left_base_index = self.margin_left - i
            right_base_index = -self.margin_right - i
            if left_base_index < 0 or right_base_index > -1:
                # index out of boundary, change direction
                if check == 'left':
                    # end of going left, time to go right from the junction
                    check = 'right'
                    i = 0
                    continue
                else:
                    # end of going right, both directions checked
                    break

            left_base = self.sequence[left_base_index]
            right_base = self.sequence[right_base_index]
            if left_base == right_base:
                # there is a repeat on the junction
                if check == 'left':
                    new_left_margin, new_right_margin = self.margin_left - i, self.margin_right + i
                else:
                    new_left_margin, new_right_margin = self.margin_left - (i - 1), self.margin_right + (i - 1)
                # TODO tu musi byc zmieniona sekwencja nowej wariacji
                new_variation = Intron(self.scaffold_name, self.scaffold_start, self.scaffold_end, margin_left=new_left_margin,
                                       margin_right=new_right_margin, sequence=self.sequence)
                self.variations.append(new_variation)
            else:
                if check == 'left':
                    # end of going left, time to go right from the junction
                    check = 'right'
                    i = 0
                    continue
                else:
                    # end of going right, both directions checked
                    break

            if check == 'left':
                # going further left
                i += 1
            else:
                # going further right
                i -= 1
        
    def movable_boundary_no_margins(self):
        # moglem zepsuc
        # ja zaraz jeszcze bardziej popsuje
        """
        Check if there are repeats on the intron junctions so the intron position could be shifted without changing
        transcript sequence. If there are, add new possible introns to self.variations.
        """
        
        if not self.prev_exon or not self.next_exon:
            return
        mls, mrs = self.prev_exon.sequence, self.next_exon.sequence
        
        i = 1
        # start checking for repeats left from the junction
        check = 'left'
        
        while True:
            if check == 'left':
                if i > len(mls) or i > len(self.sequence) or mls[-i] != self.sequence[-i]:
                    check = 'right'
                    i = 0
                    continue
                new_seq=mls[-i:]+self.sequence[:-i]
                new_variation = Intron(self.scaffold_name, scaffold_start=self.scaffold_start - i, scaffold_end=self.scaffold_end - i, gene=self.gene, sequence=new_seq)
                
            else:
                if i + 1 > len(mrs) or i + 1 > len(self.sequence) or self.sequence[i] != mrs[i]:
                    break
                new_seq=self.sequence[i:]+mrs[:i]
                new_variation = Intron(self.scaffold_name, scaffold_start=self.scaffold_start + i, scaffold_end=self.scaffold_end + i, gene=self.gene, sequence=new_seq)
            self.variations.append(new_variation)
            i += 1

    def check_conventional(self):
        """ Check if the intron junctions suggest the intron is conventional."""
        left_anchor = self.sequence[0:2]
        right_anchor = self.sequence[-2:]
        
        if left_anchor in ['GT', 'GC'] and right_anchor == 'AG':
            return True
        elif left_anchor == 'CT' and right_anchor in ['AC', 'GC']:
            return True
        else:
            return False

    def check_unconventional(self):
        """Check if intron may be unconventional according to our current knowledge, meaning it can form
        secondary structure in specific positions. Also check if variations with shifted junctions may be
        unconventional."""
        def complimentary(n1, n2):
            if {n1, n2} in [{'A', 'T'}, {'C', 'G'}, {'C', 'T'}]:
                return True
            else:
                return False

        left_anchor = self.sequence[3: 5]
        right_anchor = self.sequence[-7:-5]

        # try:
        if complimentary(left_anchor[0], right_anchor[1]) and complimentary(left_anchor[1], right_anchor[0]):
            return True
        else:
            # checking variations
            if len(self.variations) > 0:
                for son in self.variations:
                    if son.check_unconventional():
                        return True
            return False
    
    def conventional_version(self):
        if self.sequence[0:2] in ['GT', 'GC'] and self.sequence[-2:] == 'AG':
            i = 4
        else:
            return
        
        if not self.prev_exon or not self.next_exon:
            return
        mls, mrs = self.prev_exon.sequence[-3:], self.next_exon.sequence[:3]
        
        #mls, mrs = self.margin_left_seq, self.margin_right_seq
        #if check == 0:
        if self.sequence[-3] == 'C':
            i = 3
            if mls and mls[-1] == 'G' and \
                mrs and mrs[0] == 'G':
                i = 2
                if len(mls) > 1 and len(mrs) > 1 and \
                    mls[-2] == 'A' and mrs[1] == 'T':
                    i = 1
        if self.sequence[1] == 'C': i += 4
        self.is_conventional = i
        # else:  # if check = 1
        #     if self.sequence[2] == 'G':
        #         i = 3
        #         if self.margin_left_seq and self.margin_right_seq and \
        #            self.margin_left_seq[-1] == 'C' and self.margin_right_seq[0] == 'C':
        #             i = 2
        #             if len(self.margin_left_seq) > 1 and len(self.margin_right_seq) > 1 and \
        #                self.margin_left_seq[-2] == 'A' and self.margin_right_seq[1] == 'T':
        #                 i = 1
        #     if self.sequence[-2] == 'G': i += 4
        #     self.is_conventional = i


class Exon(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='', prev_exon=None, next_exon=None, prev_intron=None, next_intron=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.prev_exon = prev_exon
        self.next_exon = next_exon
        self.prev_intron = prev_intron
        self.next_intron = next_intron


class Transcript():
    def __init__(self, scaffold_name, start, end, sequence='', strand=''):
        self.scaffold_name = scaffold_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.strand = strand


def process_file(file_path):
    try:
        with open(file_path) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                yield [int(x) if x.isnumeric() else x for x in line.split()]
    except (IOError, OSError) as exc:
        print("Error opening / processing file")
        raise exc


def read_genome(file):
    genome = defaultdict(str)
    with open(file) as f:
        for line in f.readlines():
            if line[0] == '>':
                gene = line.strip()[1:]
            else:
                genome[gene] += line.strip()
    return genome
    

def read_genes(file):
    genes = {}  # slownik genow
    gene, exon = None, None
    prev = None
    for line in process_file(file):
        if line[0] == '#':
            continue
        if line[2] == 'transcript':
            if gene:
                genes[gene.name] = gene
            gene = Gene(line[0], line[3] - 1, line[4], name=line[11].strip('";'), strand=line[6], exons=[])
        elif line[2] == 'exon':
            exon = Exon(line[0], line[3] - 1, line[4], strand=line[6], prev_exon = prev)
            gene.append_exons(exon)
            if prev: prev.next_exon=exon
            prev=exon
    if gene:
        genes[gene.name] = gene
    return genes


def complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    letters = [complement_dict[base] for base in seq]
    return ''.join(letters)


def reverse_complement(seq):
    return complement(seq[::-1])

