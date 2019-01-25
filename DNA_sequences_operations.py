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

