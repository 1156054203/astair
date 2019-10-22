
from __future__ import division
from __future__ import print_function

import os
import sys
import pdb
import pysam
import logging


from astair.DNA_sequences_operations import reverse_complementary

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)


def read_vcf(vcf_file, chromosome, fasta, threads):
    """Opens neatly and separately the vcf file and takes all SNPs of interest from the current chromosomes."""
    try:
        if os.path.isfile(vcf_file) and not os.path.isfile(os.path.join(str(vcf_file) + '.tbi')):
            logs.info("Creating a TABIX index for the provided VCF file. This may take a while...")
            pysam.tabix_index(vcf_file, preset='vcf')
        true_variants, possible_mods, all_contexts = set(), dict(), dict()
        CHG = ['C' + x + 'G' for x in ['A', 'C', 'T']]
        CHH = ['C' + x + z for x in ['A', 'C', 'T'] for z in ['A', 'C', 'T']]
        CG = ['C' + 'G' + x for x in ['A', 'C', 'T', 'G', 'N']]
        CN = ['C' + 'N' + z for z in ['N', 'A', 'C', 'T', 'G']]
        CN.extend(['C' + z + 'N' for z in ['N', 'A', 'C', 'T', 'G']])
        for item_c in CHH:
            all_contexts[item_c] = 'CHH'
        for item_c in CHG:
            all_contexts[item_c] = 'CHG'
        for item_c in CG:
            all_contexts[item_c] = 'CpG' 
        for item_c in CN:
            all_contexts[item_c] = 'CN'
        file_to_open = pysam.VariantFile(vcf_file, 'r', threads=threads)
        for variant in file_to_open:
            if variant.chrom == chromosome:
                if (variant.ref in ["C", "G"]) or ((variant.ref in ["C", "G"] and variant.filter["PASS"])):
                    true_variants.add(tuple((variant.chrom, variant.start, variant.start+1)))
                elif (set(variant.alts).intersection({'C'}) or  set(variant.alts).intersection({'G'})) or (set(variant.alts).intersection({'C'}) or  set(variant.alts).intersection({'G'}) and variant.filter["PASS"]):
                    if 'C' in variant.alts:
                        if variant.start+3 < len(fasta):
                            subcontext = ('C'+fasta[variant.start+1:variant.start+3].upper())
                            context = all_contexts[subcontext]
                        else:
                            subcontext, context = 'CNN', 'CN'
                        possible_mods[tuple((variant.chrom, variant.start, variant.start+1))] = (subcontext, context, 'C', variant.ref)
                    else:
                        if variant.start-2 >= 0:
                            subcontext = reverse_complementary(fasta[variant.start-2:variant.start].upper()+'G')
                            context = all_contexts[subcontext]
                        else:
                            subcontext, context = 'CNN', 'CN'
                        possible_mods[tuple((variant.chrom, variant.start, variant.start+1))] = (subcontext, context, 'G', variant.ref)
        file_to_open.close()
        return true_variants, possible_mods 
    except Exception:
        logs.error('The input VCF file does not exist or is truncated.', exc_info=True)
        raise
