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


@auto_attr_check
class Intron:
    """
    This is a class for working with biological introns.

    :param scaffold: (str) Scaffold or chromosome on which the intron is located.
    :param start: (int) Location of the first base of the intron.
    :param end: (int) Location of the fist base of the next exon.
    :param strand (str): Optional, defines the strand on which intron is located. Must be either + or -.
    :param support: (int) Optional, how many reads support the intron.
    :param margin: (int) Optional, how much of exons sequences is included on both ends.
    :param sequence: (str) Optional, genomic sequence of the intron.
    """
    scaffold = str
    start = int
    end = int
    strand = str
    support = int
    margin = int
    sequence = str

    def __init__(self, scaffold, start, end, strand=None, support=None, margin_left=0, margin_right=0, sequence=None):
        self.scaffold = scaffold
        if start < end:
            self.start = start
            self.end = end
        else:
            self.end = start
            self.start = end
        if strand and strand not in ['+', '-']:
            raise ValueError('Strand can only be + or -')
        else:
            self.strand = strand
        self.support = support
        self.margin_left = margin_left
        self.margin_right = margin_right
        if sequence and len(sequence) != self.length():
            raise ValueError('The length of the sequence does not correspond to the intron length.')
        else:
            self.sequence = sequence
        self.variations = []

    def __repr__(self):
        return ' '.join([self.scaffold, str(self.start), str(self.end), str(self.support)])

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
            if (self.start in range(intron.start, intron.end))\
                    or (self.end in range(intron.start, intron.end))\
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
        return self.end - self.start

    def movable_boundary(self):
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
                new_variation = Intron(self.scaffold, self.start, self.end, margin_left=new_left_margin,\
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
        secondary structure in specific positions. Also check if variations with shifted junctions may be
        unconventional."""
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
            # checking variations
            if len(self.variations) > 0:
                for son in self.variations:
                    if son.check_unconventional():
                        return True
            return False
