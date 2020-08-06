class Intron:
    """
    This is a class for working with biological introns.

    :param scaffold (str): Scaffold or chromosome on which the intron is located.
    :param start (str or int): : Location of the first base of the intron. If str must be convertible to int.
    :param end (str or int): Location of the last base of the intron. If str must be convertible to int.
    :param strand (str): Optional, defines the strand on which intron is located. Must be either + or -.
    :param support (str or int): Optional, how many reads support the intron. If str must be convertible to int.
    :param margin (str or int): Optional, how much of exons sequences is included on both ends. If str must be
    convertible to int.
    :param sequence (str): Optional, genomic sequence of the intron.
    """

    def __init__(self, scaffold, start, end, strand=None, support=None, margin=0, sequence=None):
        self.scaffold = scaffold
        if start < end:
            self.start = int(start)
            self.end = int(end)
        else:
            self.start = int(end)
            self.end = int(start)
        if strand:
            self.strand = strand
        if support:
            self.support = int(support)
        if margin:
            self.margin = int(margin)
        else:
            self.margin = 0
        if sequence:
            if len(sequence) != self.length():
                raise ValueError('The length of the sequence does not correspond to the intron length.')
            self.sequence = sequence

    def __repr__(self):
        if self.support:
            return ' '.join([self.scaffold, str(self.start), str(self.end), str(self.support)])
        else:
            return ' '.join([self.scaffold, str(self.start), str(self.end)])

    def __str__(self):
        return self.__repr__()

    def intersect(self, intron):
        """
        Checks if two introns cover partially the same area.

        :param intron (Intron): Intron with which intersection of self is checked.
        :return: Bool: True if introns share at least one base, False otherwise.
        """
        if not isinstance(intron, Intron):
            raise ValueError('Both introns must be instances of Intron class')

        else:
            if self.scaffold != intron.scaffold:
                raise ValueError('Introns on different scaffold cannot intersect.')
            if (self.start in range(intron.start, intron.end + 1))\
                    or (self.end in range(intron.start, intron.end + 1))\
                    or (self.start < intron.start and self.end > intron.end):
                return True
            else:
                return False

    def classify_intersection(self, intron):
        """
        For two intersecting introns return the type of intersection

        :param intron (Intron): Intron with which type of intersection with self is checked.
        :return: 'same' if the introns are identical, 'one end' if only start or end is the same,
        'rest' otherwise.
        """
        if self.scaffold != intron.scaffold:
            raise ValueError('Introns do not intersect at all.')
        if self.start == intron.start and self.end == intron.end:
            return 'same'
        elif self.start == intron.start or self.end == intron.end:
            return 'one end'
        else:
            return 'rest'

    def length(self):
        """Return length of the intron."""
        return self.end - self.start + 1

    def check_conventional(self):
        """ Check if the intron junctions suggest the intron is conventional."""
        left_anchor = self.sequence[self.margin:self.margin + 2]
        right_anchor = self.sequence[-self.margin - 2:-self.margin]
        if left_anchor in ['GT', 'GC'] and right_anchor == 'AG':
            return True
        elif left_anchor == 'CT' and right_anchor in ['AC', 'GC']:
            return True
        else:
            return False

    def check_unconventional(self):
        """Check if intron may be nonconventional according to our current knowledge, meaning it can form
        secondary structure in specific positions."""
        left_anchor = self.sequence[self.margin + 3:self.margin + 5]
        right_anchor = self.sequence[-self.margin - 7: -self.margin - 5]
        if left_anchor == 'CA' and right_anchor == 'TG':
            return True
        else:
            return False
