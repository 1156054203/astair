
import sys

if sys.version[0] == '3':
    
    def complementary(sequence):
        """Takes an input DNA string and gives its complementary."""
        final_string = sequence.translate({ord("T"): "A", ord("A"): "T", ord("G"): "C", ord("C"): "G", ord("t"): "a", ord("a"): "t", ord("g"): "c", ord("c"): "g"})
        return final_string
    
    def reverse_complementary(sequence):
        """Takes an input DNA string and gives its reverse complementary."""
        reverse_string = list(sequence)
        reverse_string.reverse()
        final_string = "".join(reverse_string).translate({ord("T"): "A", ord("A"): "T", ord("G"): "C", ord("C"): "G", ord("t"): "a", ord("a"): "t", ord("g"): "c", ord("c"): "g"})
        return final_string


elif sys.version[0] == '2':
    
    import string 

    def complementary(sequence):
        """Takes an input DNA string and gives its complementary."""
        translate_DNA = string.maketrans("TtCcGgAa", "AaGgCcTt")
        final_string = sequence.translate(translate_DNA)
        return final_string
     
    def reverse_complementary(sequence):
        """Takes an input DNA string and gives its reverse complementary."""
        reverse_string = list(sequence)
        reverse_string.reverse()
        translate_DNA = string.maketrans("TtCcGgAa", "AaGgCcTt")
        final_string = "".join(reverse_string).translate(translate_DNA)
        return final_string
    
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

def reverse(sequence):
    """Takes an input DNA string and gives its reverse."""
    reverse_string = list(sequence)
    reverse_string.reverse()
    final_string = "".join(reverse_string)
    return final_string

