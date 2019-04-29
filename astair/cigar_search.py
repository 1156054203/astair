import re

def cigar_search(read_data):
    """Looks whether there are indels, soft clipping or pads the CIGAR string"""
    changes = [int(s) for s in re.findall(r'\d+', read_data)]
    non_overlap = [x + 1 if x == 0 else x for x in changes]
    names = list(re.findall(r'[^\W\d_]+', read_data))
    positions = [x for x in [sum(non_overlap[0:i]) for i in range(1, len(non_overlap)+1)]]
    return names, positions, changes
