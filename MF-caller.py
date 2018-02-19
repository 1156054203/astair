import click
import pysam
import tempfile
import os
import sys
import scipy
import allel
import h5py
import itertools
import csv



@click.command()
@click.argument('fasta_file')
@click.argument('input_file')
def caller_exec(fasta_file, input_file):
    caller(fasta_file, input_file)

def caller(fasta_file, input_file):
    '''input_file is file name without extension
    fasta_file is reference.fasta
    '''
    try:
        tempdir = tempfile.TemporaryDirectory(suffix='tempVCF', dir='./')
        def piler(NAME,DIR,FASTA,FILE):
            temp = tempfile.NamedTemporaryFile(mode='w+', prefix=NAME, suffix='.vcf', dir=DIR)
            temp.write(pysam.mpileup("--excl-flags", "260", "-t" "DP,AD,ADF,ADR", "-u", "-v","-f", FASTA, FILE))
            return temp

        temp1 = piler('OT',tempdir.name, fasta_file, input_file+"OT.bam")
        temp2 = piler('OB',tempdir.name,  fasta_file, input_file+"OB.bam")

        def caller(NAME, alele,  REF_base, ALT_base):  
            allel.vcf_to_hdf5(NAME,'temp.h5', fields='*', overwrite=True, log=sys.stdout, compression_opts=9)
            callset = h5py.File('temp.h5', mode='r')
            chrom = callset['variants/CHROM'][:]
            pos = callset['variants/POS'][:]
            ref = callset['variants/REF'][:]
            alt = callset['variants/ALT'][:]
            mods = callset[alele][:,0,1:]
            org = callset[alele][:,0,0]
            depth = callset['calldata/DP'][:,0]
            results = [(i,*j) for i,j in itertools.zip_longest(ref,alt, fillvalue='BLANK')]
            indexed = [(i,*j) for i, j in enumerate(results)]
            index = [(x[0], x.index(ALT_base)-2) for x in indexed if ALT_base in x[:] and   REF_base in x[1]] # gives row of c and t and t position + 1
            data = [(chrom[i], pos[i]-1 , pos[i], mods[i,k]/(org[i]+mods[i,k]), mods[i,k], org[i], ALT_base, REF_base, depth[i]) for n in range(int(len(index)/1000)+1)  for i, k in index[n*1000:n*1000+1000] if mods[i,k] > 0]
            os.system("rm -f temp.h5")
            return data

        data_OT = caller(temp1.name, 'calldata/AD', 'C', 'T')
        data_OB = caller(temp2.name, 'calldata/AD', 'G', 'A')
        

        with open(input_file+"OT.mods", 'w', newline='') as myfile:
            wr=csv.writer(myfile, delimiter='\t', lineterminator='\n')
            wr.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "NMOD", "ALT", "REF", "DEPTH"])
            for row in data_OT:
                wr.writerow(row)

        with open(input_file+"OB.mods",  'w', newline='') as myfile:
            wr=csv.writer(myfile, delimiter='\t', lineterminator='\n')
            wr.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "NMOD", "ALT", "REF", "DEPTH"])
            for row in data_OB:
                wr.writerow(row)
                
        data_OT.extend(data_OB)
        data_OT.sort()
        
        with open(input_file+"_total.mods",  'w', newline='') as myfile:
            wr=csv.writer(myfile, delimiter='\t', lineterminator='\n')
            wr.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "NMOD", "ALT", "REF", "DEPTH"])
            for row in data_OT:
                    wr.writerow(row)
                
                
        with open(input_file+"_total.bedGraph",  'w', newline='') as myfile:
            wr=csv.writer(myfile, delimiter='\t', lineterminator='\n')
            for row in data_OT:
                    wr.writerow(row[0:4])
    finally:       
        temp1.close()
        temp2.close() 


if __name__ == '__main__':
    caller_exec()
