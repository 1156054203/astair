from __future__ import print_function

import re
import pdb
import logging
import itertools
import ahocorasick

from DNA_sequences_operations import complementary

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

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
        CN = [y + z + x for x in ['N', 'n', 'A', 'C', 'T', 'a', 'c', 't'] for y in ['C', 'c'] for z in ['N', 'n']]
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


def ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts, region):
    """Looks for cytosine contexts in the reference fasta file."""
    auto = ahocorasick.Automaton()
    for pattern in context[objects]:
        if not isinstance(pattern, str):
            pattern = str(pattern)
        auto.add_word(pattern, pattern)
    auto.make_automaton()
    if objects[-1] == 'b':
        for end_ind, found in auto.iter(complementary(string)):
            context_total_counts[objects] += 1
            reversed = list(found.upper())
            reversed.reverse()
            context_total_counts["".join(reversed)] += 1
            if region == None or (end_ind >= region[1] and end_ind+1 <=region[2]):
                if objects != 'CGb':
                    data_context[(string_name, end_ind, end_ind + 1)] = tuple(("".join(reversed), objects[0:-1], 'A', 'G'))
                else:
                    data_context[(string_name, end_ind, end_ind + 1)] = tuple(("".join(reversed), 'CpG', 'A', 'G'))
    elif objects == 'CG':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            if region == None or (end_ind -2  >= region[1] and end_ind - 1 <= region[2]):
                data_context[(string_name, end_ind - 2, end_ind - 1)] = tuple((found.upper(), 'CpG', 'T', 'C'))
    elif objects == 'CHG' or objects == 'CHH':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            if region == None or (end_ind - 2 >= region[1] and end_ind - 1 <= region[2]):
                data_context[(string_name, end_ind - 2, end_ind - 1)] = tuple((found.upper(), objects, 'T', 'C'))
    elif objects == 'CN':
        for end_ind, found in auto.iter(string):
            context_total_counts[objects] += 1
            context_total_counts[found.upper()] += 1
            if region == None or (end_ind - 2 >= region[1] and end_ind - 1 <= region[2]):
                data_context[(string_name, end_ind - 2, end_ind -1 )] = tuple((found.upper(), 'CN', 'T', 'C'))
    elif objects == 'user':
        index_c = user_defined_context.upper().find('C')
        if index_c == -1:
            index_c = None
        try:
            for end_ind, found in auto.iter(string):
                context_total_counts['user defined context'] += 1
                if region == None or (end_ind - len(user_defined_context) + index_c + 1 >= region[1] and end_ind - len(user_defined_context) + index_c + 2<= region[2]):
                    data_context[(string_name, end_ind - len(user_defined_context) + index_c + 1, end_ind - len(user_defined_context) + index_c + 2)] = tuple(
                        (user_defined_context, 'user defined context', 'T', 'C'))
            for end_ind, found in auto.iter(complementary(string)):
                reversed = list(found.upper())
                reversed.reverse()
                index_c = "".join(reversed).upper().find('C')
                context_total_counts['user defined context'] += 1
                if region == None or (end_ind - index_c >= region[1] and end_ind - index_c + 1 <= region[2]):
                    data_context[(string_name, end_ind - index_c, end_ind - index_c + 1)] = tuple((user_defined_context, 'user defined context', 'A', 'G'))
        except TypeError:
            logs.error('The user-provided context does not contain cytosines. asTair will not output any user-provided context summary.', exc_info=True)
            raise
    if context in ['all', 'CG', 'CGb']:
        if re.match(r'CG', string[0:3]) and region == None or (re.match(r'CG', string[0:3]).start() >= region[1] and re.match(r'CG', string[0:3]).start() + 2 <= region[2]):
            context_total_counts['CG'] += 2
            data_context[(string_name, re.match(r'CG', string[0:3]).start(), re.match(r'CG', string[0:3]).start() + 1)] = tuple(('CGN', 'CpG', 'T', 'C'))
            data_context[(string_name, re.match(r'CG', string[0:3]).start() + 1, re.match(r'CG', string[0:3]).start() + 2)] = tuple(('CGN', 'CpG', 'A', 'G'))
        elif re.match(r'CG', string[-2:]) and region == None or (re.match(r'CG', string[-2:]).start() >= region[1] and re.match(r'CG', string[-2:]).start() + 2 <= region[2]):
            context_total_counts['CG'] += 2
            data_context[(string_name, re.match(r'CG', string[-2:]).start(), re.match(r'CG', string[-2:]).start() + 1)] = tuple(('CGN', 'CpG', 'T', 'C'))
            data_context[(string_name, re.match(r'CG', string[-2:]).start() + 1, re.match(r'CG', string[-2:]).start() + 2)] = tuple(('CGN', 'CpG', 'A', 'G'))


def context_sequence_search(context, key, fastas, string_name, user_defined_context, context_total_counts, region):
    """Starts the search for cytosine contexts in the reference fasta file."""
    data_context = {}
    string = fastas[string_name]
    if key.count('C') == 0:
        for objects in key:
            try:
                ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts, region)
            except TypeError:
                pass
    else:
        objects = "".join(key)
        try:
            ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context, context_total_counts, region)
        except TypeError:
            pass
    return data_context
