import pandas
import numpy

class WrappingException(Exception):
    def __init__(self, error_message):
        super(WrappingException, self).__init__(error_message)
        self.error_message = error_message

class RepresentationWrapper:
    def __init__(self, sequences, background_corpus, type_of_input, alphabet: [], space_character, classes):
        #type_of_input is one of:
        #"list_of_FASTA_pairs" - this is how Bio reads FASTA format, with the first pair element being the gene/sequence name/annotation, and the second pair element is the sequence string
        if type_of_input == "list_of_FASTA_pairs":
            self.sequences = sequences
            self.background_corpus = background_corpus
            self.alphabet = alphabet
            self.space_character = space_character
            self.classes = classes
        else:
            raise WrappingException("Unknown input type %s"%type_of_input)
    
    def sanity_check(self):
        #checks data consistency.
        #returns True on success and False with comment on failure
        #check if only alphabet sequences or space is present: 
        extended_alphabet = self.alphabet
        for ind, pair in enumerate(self.sequences):
            letters = set(list(pair[1]))
            if len(letters - {self.space_character} - set(self.alphabet)) > 0:
                return False, "Sequence %d has some letters outside of the alhpabet (and space): %s"%(ind, letters - {self.space_character} - set(self.alphabet))
    
        # check if len of classes is the same as len of sequences:
        if len(self.sequences) != len(self.classes):
            return False, "Length of classes list doesn't match sequence count"
        
        return True, "Everything OK"
    
    
    def to_TfIdf(self, n=4) -> []:
        pass
    
    def to_pandas_dataframe(self):
        #first, constructing a list of pandas.Series objects
        pandas_objects = []
        lengths = []
        for _, pair in enumerate(self.sequences):
            pandas_objects.append(pandas.Series(data=list(pair[1]), name=pair[0]))
            lengths.append(len(pair[1]))
        pandas_dataframe = pandas.DataFrame(pandas_objects, dtype = pandas.SparseDtype("object", numpy.nan))
        pandas_dataframe["class"] = self.classes
        pandas_dataframe["length"] = lengths
        return pandas_dataframe
         