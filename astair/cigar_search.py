import re


def cigar_search(read_data):
    """Looks whether there are indels, soft clipping or pads the CIGAR string"""
    changes = [int(s) for s in re.findall(r'\d+', read_data)]
    non_overlap = [x + 1 if x == 0 else x for x in changes]
    names = list(re.findall(r'[^\W\d_]+', read_data))
    positions = [x for x in [sum(non_overlap[0:i]) for i in range(1, len(non_overlap)+1)]]
    return names, positions, changes


def position_correction_cigar(read, method, random_sample, positions, reverse_modification):
    """Uses the CIGAR string information to correct the expected cytosine positions."""
    names, positions_cigar, changes = cigar_search(read.cigarstring)
    index = 0
    for change in names:
        if len(positions) != 0:
            if change == 'D':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(read.qstart - read.reference_start) - changes[index] for x in positions if x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(
                            read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(
                            read.qstart - read.reference_start) - changes[index] for x in subsample]
                else:
                    corrected_positions = [x if x < positions_cigar[index] else x - changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'I':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart-read.reference_start) if (x[1] - abs(read.qstart-read.reference_start)) < positions_cigar[index] else x[1] - abs(read.qstart-read.reference_start) + changes[index] for x in positions if x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(
                            read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(
                            read.qstart - read.reference_start) + changes[index] for x in subsample]
                else:
                    corrected_positions = [x if x < positions_cigar[index] else x + changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'S' or change == 'H':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) for x in positions if
                                               x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart-read.reference_start) for x in subsample]
                else:
                    if index == 0:
                        corrected_positions = [x for x in positions if x > positions_cigar[index]]
                    else:
                        corrected_positions = [x for x in positions if x < positions_cigar[index]]
                index += 1
                positions = corrected_positions
            else:
                index += 1
    return positions
