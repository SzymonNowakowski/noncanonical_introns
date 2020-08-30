class Intron:
    """
    This is a class for working with biological introns.

    :param scaffold (str): Scaffold or chromosome on which the intron is located.
    :param start (str or int): : Location of the first base of the intron. If str must be convertible to int.
    :param end (str or int): Location of the fist base of the next exon. If str must be convertible to int.
    :param strand (str): Optional, defines the strand on which intron is located. Must be either + or -.
    :param support (str or int): Optional, how many reads support the intron. If str must be convertible to int.
    :param margin (str or int): Optional, how much of exons sequences is included on both ends. If str must be
    convertible to int.
    :param sequence (str): Optional, genomic sequence of the intron.
    """

    def __init__(self, scaffold, start, end, strand=None, support=None, margin_left=0, margin_right=0, sequence=None):
        self.scaffold = scaffold
        if int(start) < int(end):
            self.start = int(start)
            self.end = int(end)
        else:
            self.start = int(end)
            self.end = int(start)
        self.strand = strand
        try:
            self.support = int(support)
        except TypeError:
            self.support = None
        if margin_left:
            self.margin_left = int(margin_left)
        else:
            self.margin_left = 0
        if margin_right:
            self.margin_right = int(margin_right)
        else:
            self.margin_right = 0
        if sequence:
            # if len(sequence) != self.length():
            #     raise ValueError('The length of the sequence does not correspond to the intron length.')
            self.sequence = sequence
        else:
            self.sequence = None
        self.sons = []

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

        :param intron: (Intron) Intron with which intersection of self is checked.
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

        :param intron: (Intron) Intron with which type of intersection with self is checked.
        :return: 2 if the introns are identical, 1 if only either start or end is the same, 0 otherwise.
        """
        if self.scaffold != intron.scaffold:
            raise ValueError('Introns do not intersect at all.')
        if self.start == intron.start and self.end == intron.end:
            return 2
        elif self.start == intron.start or self.end == intron.end:
            return 1
        else:
            return 0

    def length(self):
        """Return length of the intron."""
        if self.end - self.start < 0:
            print(self, self.margin_left, self.margin_right)
        return self.end - self.start

    def movable_boundary(self):
        i = 1
        # start checking left from the junction
        check = 'left'

        while True:
            left_base_index = self.margin_left - i
            right_base_index = -self.margin_right - i
            if left_base_index < 0 or right_base_index > -1:
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
                if check == 'left':
                    new_left_margin, new_right_margin = self.margin_left - i, self.margin_right + i
                else:
                    new_left_margin, new_right_margin = self.margin_left - (i - 1), self.margin_right + (i - 1)
                self.sons.append(Intron(self.scaffold, self.start, self.end, margin_left=new_left_margin,\
                                        margin_right=new_right_margin, sequence=self.sequence))
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

    def check_conventional(self):
        """ Check if the intron junctions suggest the intron is conventional."""
        left_anchor = self.sequence[self.margin_left:self.margin_left + 2]
        right_anchor = self.sequence[-self.margin_right - 2:-self.margin_right]
        if left_anchor == 'GT' and right_anchor == 'AG':
            return True
        elif left_anchor == 'CT' and right_anchor == 'AC':
            return True
        else:
            return False

    def check_unconventional(self):
        """Check if intron may be unconventional according to our current knowledge, meaning it can form
        secondary structure in specific positions."""
        def complimentary(n1, n2):
            if {n1, n2} in [{'A', 'T'}, {'C', 'G'}, {'C', 'T'}]:
                return True
            else:
                return False

        left_anchor = self.sequence[self.margin_left + 3:self.margin_left + 5]
        right_anchor = self.sequence[-self.margin_right - 7: -self.margin_right - 5]

        # try:
        if complimentary(left_anchor[0], right_anchor[1]) and complimentary(left_anchor[1], right_anchor[0]):
            return True
        else:
            if len(self.sons) > 0:
                for son in self.sons:
                    if son.check_unconventional():
                        if son.margin_left < 0 or son.margin_right < 0:
                            print('somethings fucked', son, son.margin_left, son.margin_right,\
                                  son.sequence[:10], son.sequence[-10:])
                        return True
            return False
        # except IndexError:
        #     print('fuckedup: ', self, self.margin_left, self.margin_right, self.sequence)
        #     return False
