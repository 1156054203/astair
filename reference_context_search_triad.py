import re
import pdb
import itertools
import ahocorasick


from DNA_sequences_operations import complementary


def sequence_context_set_creation(desired_sequence, user_defined_context):
    """Prepares sets of possible cytosine contexts."""
    letters_top = ['A', 'C', 'T', 'a', 'c', 't']
    if user_defined_context:
        user = list(map(''.join, itertools.product(*zip(user_defined_context.upper(), user_defined_context.lower()))))
    if desired_sequence == 'all':
        CHG = [y + x + z for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]
        CHGb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in ['G', 'g']]
        CHH = [y + x + z for x in letters_top for y in ['C', 'c'] for z in letters_top]
        CHHb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in letters_top]
        CG = [y + z + x for y in ['C', 'c'] for z in ['G', 'g'] for x in ['A', 'C', 'T', 'a', 'c', 't', 'g', 'G']]
        CGb = [x+z + y for y in ['C', 'c'] for z in ['G', 'g'] for x in ['A', 'C', 'T', 'a', 'c', 't', 'g', 'G']]
        CN = [y + x + z for x in ['N', 'n', 'A', 'C', 'T', 'a', 'c', 't'] for y in ['C', 'c'] for z in ['N', 'n']]
        CNb = [y + x + z for x in ['N', 'n', 'A', 'C', 'T', 'a', 'c', 't'] for z in ['C', 'c'] for y in ['N', 'n']]
        if user_defined_context:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'CHH': list(CHH), 'CHHb': list(CHHb), 'CG': list(CG), 'CGb': list(CGb),'CN': list(CN), 'CNb': list(CNb), 'user': list(user)}
            all_keys = list(('CHG', 'CHGb', 'CHH', 'CHHb', 'CG', 'CGb', 'CN', 'CNb', 'user'))
        else:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'CHH': list(CHH), 'CHHb': list(CHHb), 'CG': list(CG), 'CGb': list(CGb), 'CN': list(CN), 'CNb': list(CNb)}
            all_keys = list(('CHG', 'CHGb', 'CHH', 'CHHb', 'CG', 'CGb', 'CN', 'CNb'))
    elif desired_sequence == 'CpG':
        CG = [y + z + x for y in ['C', 'c'] for z in ['G', 'g'] for x in ['A', 'C', 'T', 'a', 'c', 't', 'g', 'G']]
        CGb = [x+z + y  for y in ['C', 'c'] for z in ['G', 'g'] for x in ['A', 'C', 'T', 'a', 'c', 't', 'g', 'G']]
        if user_defined_context:
            contexts = {'CG': list(CG), 'CGb': list(CGb), 'user': list(user)}
            all_keys = list(( 'CG', 'CGb', 'user'))
        else:
            contexts = {'CG': list(CG), 'CGb': list(CGb)}
            all_keys = list(('CG', 'CGb'))
    elif desired_sequence == 'CHG':
        CHG = [y + x + z for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]
        CHGb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in ['G', 'g']]
        if user_defined_context:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'user': list(user)}
            all_keys = list(( 'CHG', 'CHGb', 'user'))
        else:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb)}
            all_keys = list(('CHG', 'CHGb'))
    elif desired_sequence == 'CHH':
        CHH = [y + x + z for x in letters_top for y in ['C', 'c'] for z in letters_top]
        CHHb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in letters_top]
        if user_defined_context:
            contexts = {'CHH': list(CHH), 'CHHb': list(CHHb), 'user': list(user)}
            all_keys = list(( 'CHH', 'CHHb', 'user'))
        else:
            contexts = {'CHH': list(CHH), 'CHHb': list(CHHb)}
            all_keys = list(('CHH', 'CHHb'))
    return contexts, all_keys


def ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts):
    """Looks for cytosine contexts in the reference fasta file."""
    auto = ahocorasick.Automaton()
    for pattern in context[objects]:
        auto.add_word(pattern, pattern)
    auto.make_automaton()
    if objects[-1] == 'b':
        for end_ind, found in auto.iter(complementary(string)):
            context_total_counts[objects] += 1
            reversed = list(found.upper())
            reversed.reverse()
            context_total_counts["".join(reversed)] += 1
            if objects != 'CGb':
                data_context[(string_name, end_ind, end_ind + 1)] = tuple(("".join(reversed), objects[0:-1], 'A', 'G'))
            else:
                data_context[(string_name, end_ind, end_ind + 1)] = tuple(("".join(reversed), 'CpG', 'A', 'G'))
    elif objects == 'CG':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            data_context[(string_name, end_ind - 2, end_ind - 1)] = tuple((found.upper(), 'CpG', 'T', 'C'))
    elif objects == 'CHG' or objects == 'CHH':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            data_context[(string_name, end_ind - 2, end_ind - 1)] = tuple((found.upper(), objects, 'T', 'C'))
    elif objects == 'CN':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            data_context[(string_name, end_ind - 1, end_ind)] = tuple((found.upper(), 'CN', 'T', 'C'))
    elif objects == 'user':
        index_c = [user_defined_context.find('C') if user_defined_context.find('C') >= 0 else user_defined_context.find('c') if user_defined_context.find('c') >= 0 else print('The user defined context does not contain cytosines.')][0]
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            data_context[(string_name, end_ind - index_c - 1, end_ind - index_c)] = tuple(
                (found.upper(), 'user defined context', 'T', 'C'))
        for end_ind, found in auto.iter(complementary(string)):
            context_total_counts[objects] += 1
            reversed = list(found.upper())
            reversed.reverse()
            context_total_counts["".join(reversed)] += 1
            data_context[(string_name, end_ind - index_c, end_ind - index_c + 1)] = tuple(
                ("".join(reversed), 'user defined context', 'A', 'G'))
    if context in ['all', 'CG', 'CGb']:
        if re.match(r'CG', string[0:3]):
            context_total_counts['CG'] += 2
            data_context[(string_name, re.match(r'CG', string[0:3]).start(), re.match(r'CG', string[0:3]).start() + 1)] = tuple(('CGN', 'CpG', 'T', 'C'))
            data_context[(string_name, re.match(r'CG', string[0:3]).start() + 1, re.match(r'CG', string[0:3]).start() + 2)] = tuple(('CGN', 'CpG', 'A', 'G'))
        elif re.match(r'CG', string[-2:]):
            context_total_counts['CG'] += 2
            data_context[(string_name, re.match(r'CG', string[-2:]).start(), re.match(r'CG', string[-2:]).start() + 1)] = tuple(('CGN', 'CpG', 'T', 'C'))
            data_context[(string_name, re.match(r'CG', string[-2:]).start() + 1, re.match(r'CG', string[-2:]).start() + 2)] = tuple(('CGN', 'CpG', 'A', 'G'))


def context_sequence_search(context, key, fastas, string_name, user_defined_context, context_total_counts):
    """Starts the search for cytosine contexts in the reference fasta file."""
    data_context = {}
    string = fastas[string_name]
    if key.count('C') == 0:
        for objects in key:
            ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts)
    else:
        objects = "".join(key)
        ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts)
    return data_context
