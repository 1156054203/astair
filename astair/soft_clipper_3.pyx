#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import logging


logs = logging.getLogger(__name__)


def soft_clipper(reads, start_clip, end_clip, add_indels):
    """Ensures that only bases whitin the desired read positions are used for modification calling."""
    new_pileups = [reads.pileups[i] for i in range(reads.n) if reads.get_query_positions()[i]>start_clip and  reads.get_query_positions()[i]<(reads.pileups[i].alignment.rlen - end_clip)]
    try:
        sequences = [reads.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=add_indels)[i]  for i in range(reads.n) if reads.get_query_positions()[i]>start_clip and  reads.get_query_positions()[i]<(reads.pileups[i].alignment.rlen - end_clip)] 
    except AssertionError:
        logs.exception("Failed getting query sequences (AssertionError, pysam). Please decrease the max_depth parameter.")
    return new_pileups, sequences
