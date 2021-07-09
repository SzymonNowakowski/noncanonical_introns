import pandas
import numpy
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer, HashingVectorizer

class WrappingException(Exception):
    def __init__(self, error_message):
        super(WrappingException, self).__init__(error_message)
        self.error_message = error_message

class RepresentationWrapper:
    def __init__(self, sequences, type_of_input, alphabet: [], space_character, classes):
        #type_of_input is one of:
        #"list_of_FASTA_pairs" - this is how Bio reads FASTA format, with the first pair element being the gene/sequence name/annotation, and the second pair element is the sequence string
        if type_of_input == "list_of_FASTA_pairs":
            self.sequences = []
            self.names = []
            for ind, pair in enumerate(sequences):
                self.names.append(pair[0])
                self.sequences.append(pair[1]) 
            self.alphabet = alphabet
            self.space_character = space_character
            self.classes = classes
        else:
            raise WrappingException("Unknown input type: %s"%type_of_input)
    
    def sanity_check(self):
        #checks data consistency.
        #returns True on success and False with comment on failure
        #check if only alphabet sequences or space is present: 
        
        comments = ""
        
        extended_alphabet = set(self.alphabet).union({self.space_character})
        for ind, seq in enumerate(self.sequences):
            letters = set(seq)
            if len(letters - extended_alphabet) > 0:
                return False, "Sequence %d has some letters outside of the alhpabet (and space): %s"%(ind, letters - {self.space_character} - set(self.alphabet))
    
        # check if len of classes is the same as len of sequences:
        if self.classes is None:
            comments = comments + ", but classes are not provided"
        else:
            if len(self.sequences) != len(self.classes):
                return False, "Length of classes list doesn't match sequence count"

        
          # check if len of names is the same as len of sequences:
        if len(self.sequences) != len(self.names):
            return False, "Length of names list doesn't match sequence count"
        

        return True, "Everything OK" + comments
    
    
    def to_TfIdf(self, ngram_length=4, space_treatment="exclude") -> []:
        #space_treatment is one of
        #"include" - it will treat space as a regular alphabet characters e.g. for "atcg_attcg" will be decomposed into 9 n-grams with length n=4: atcg, tcg_, cg_a, g_at, _att, attc, ttcg
        #"exclude" (default) - it will exclude all n-grams with spaces, e.g. string "atcg_attcg" will be decomposed into 3 ng-rams with length n=4: atcg, attc, ttcg
        #"average" (not implemented)- treat space as any possible letter from alhpabet, "atcg_attcg" should be treated as a sequence representing 25% chances of being "atcgAattcg", 25% chances of being "atcgTattcg", 25% chances of being "atcgCattcg", 25% chances of being "atcgGattcg".
        
        if space_treatment == "include":
            alphabet_str = self.alphabet+self.space_character
        elif space_treatment == "exclude":
            alphabet_str = self.alphabet
        elif space_treatment == "average":
            raise NotImplementedError
        else:
            raise WrappingException("Unknown space treatment type: %s"%space_treatment)
            
        vectorizer = TfidfVectorizer(lowercase=False, analyzer = "word", dtype=numpy.float32, token_pattern="(?=([%s]{%d}))"%(alphabet_str, ngram_length))
        #positivelookahed used for matching overlapping expressions: https://exceptionshub.com/how-to-find-overlapping-matches-with-a-regexp.html
        self.vector = vectorizer.fit_transform(self.sequences)
        return self.vector

    def to_fast_bag_of_words(self, ngram_length=4, space_treatment="exclude") -> []:
        #space_treatment is one of
        #"include" - it will treat space as a regular alphabet characters e.g. for "atcg_attcg" will be decomposed into 9 n-grams with length n=4: atcg, tcg_, cg_a, g_at, _att, attc, ttcg
        #"exclude" (default) - it will exclude all n-grams with spaces, e.g. string "atcg_attcg" will be decomposed into 3 ng-rams with length n=4: atcg, attc, ttcg
        #"average" (not implemented)- treat space as any possible letter from alhpabet, "atcg_attcg" should be treated as a sequence representing 25% chances of being "atcgAattcg", 25% chances of being "atcgTattcg", 25% chances of being "atcgCattcg", 25% chances of being "atcgGattcg".
        
        if space_treatment == "include":
            alphabet_str = self.alphabet+self.space_character
        elif space_treatment == "exclude":
            alphabet_str = self.alphabet
        elif space_treatment == "average":
            raise NotImplementedError
        else:
            raise WrappingException("Unknown space treatment type: %s"%space_treatment)
            
        vectorizer = HashingVectorizer(n_features = 2**len(alphabet_str), lowercase=False, analyzer = "word", token_pattern="(?=([%s]{%d}))"%(alphabet_str, ngram_length))
        #positivelookahed used for matching overlapping expressions: https://exceptionshub.com/how-to-find-overlapping-matches-with-a-regexp.html
        self.vector = vectorizer.fit_transform(self.sequences)
        return self.vector

    def to_bag_of_words(self, ngram_length=4, space_treatment="exclude") -> []:
        #space_treatment is one of
        #"include" - it will treat space as a regular alphabet characters e.g. for "atcg_attcg" will be decomposed into 9 n-grams with length n=4: atcg, tcg_, cg_a, g_at, _att, attc, ttcg
        #"exclude" (default) - it will exclude all n-grams with spaces, e.g. string "atcg_attcg" will be decomposed into 3 ng-rams with length n=4: atcg, attc, ttcg
        #"average" (not implemented)- treat space as any possible letter from alhpabet, "atcg_attcg" should be treated as a sequence representing 25% chances of being "atcgAattcg", 25% chances of being "atcgTattcg", 25% chances of being "atcgCattcg", 25% chances of being "atcgGattcg".
        
        if space_treatment == "include":
            alphabet_str = self.alphabet+self.space_character
        elif space_treatment == "exclude":
            alphabet_str = self.alphabet
        elif space_treatment == "average":
            raise NotImplementedError
        else:
            raise WrappingException("Unknown space treatment type: %s"%space_treatment)
            
        vectorizer = CountVectorizer(lowercase=False, analyzer = "word", dtype=numpy.int16, token_pattern="(?=([%s]{%d}))"%(alphabet_str, ngram_length))
        #positivelookahed used for matching overlapping expressions: https://exceptionshub.com/how-to-find-overlapping-matches-with-a-regexp.html
        self.vector = vectorizer.fit_transform(self.sequences)
        return self.vector

    def to_kmer(self, k, filename):
        # Creates a new file
        with open(filename, 'w') as handle:
            handle.write("sequence\tlabel\n")
            for i, seq in enumerate(self.sequences):
                if i%100000 == 0:
                    print(i)
                handle.write(seq2kmer(seq, k)+"\t%d\n"%self.classes[i])
                
            
        
        
    def last_vector_representation(self, n=4) -> []:
        return self.vector
    
    def to_pandas_dataframe(self) -> pandas.DataFrame:
        #first, constructing a list of pandas.Series objects
        pandas_objects = []
        lengths = []
        for ind in range(len(self.sequences)):
            pandas_objects.append(pandas.Series(data=list(self.sequences[ind]), name=self.names[ind]))
            lengths.append(len(self.sequences[ind]))
        pandas_dataframe = pandas.DataFrame(pandas_objects, dtype = pandas.SparseDtype("object", numpy.nan))
        pandas_dataframe["class"] = self.classes
        pandas_dataframe["length"] = lengths
        return pandas_dataframe
         
        
def seq2kmer(seq, k):
    """
    from DNABBERT: https://github.com/jerryji1993/DNABERT
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space

    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers