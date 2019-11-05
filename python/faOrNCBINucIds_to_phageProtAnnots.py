import os
import time
import argparse
import urllib.error

from Bio import Entrez, SeqIO
from collections import OrderedDict
from itertools import groupby

# HARD CODED PARAMS AND PATHS:

def check_paths(*args):
    path = next((arg for arg in args if os.path.exists(arg)), None)
    if not path:
        print('Warning: no path found in: {}'.format(*args))
    return path

# Tool paths
METAGENEMARK_PATH = ''  # not needed use Prodigal
PRODIGAL_PATH = check_paths('/path/prodigal')
DIAMOND_PATH = check_paths('/path/diamond')
HMMSEARCH_PATH = check_paths('/path/hmmer3/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch')

# Data paths
METAGENEMARK_MODEL_PATH = ''  # not needed use Prodigal 
NR_PATH = check_paths('/path/nr.dmnd') # The protein database to scan against, can use full NR path (slow) or custom phage database of proteins
PVOGS_HMMS_PATH = check_paths('/path/allpvoghmms')
PVOGS_TABLE_PATH = check_paths('/path/VOGTable.txt')

# OTHER
NR_ANNOTS_TO_IGNORE = ('hypothetical', 'MULTISPECIES: hypothetical', 'putative uncharacterized', 'unknown', 'uncharacterized',
                            'Uncharacterized', 'Uncharacterised', 'predicted', 'conserved hypothetical')

PVOG_ANNOTS_TO_IGNORE = ('hypothetical', 'No annotation provided')

# BioEntrez Params
Entrez.email = 'example@hotmail.com'
#Entrez.api_key = 'key here'

def parse_args():
    parser = argparse.ArgumentParser(description='Get phage protein annotations from initial nucleotide fasta file or NCBI nucleotide ids @Chris Rands 2018')
    parser.add_argument('--input', required=True, help='Input, either nucleotide fasta (ending .fa) or NCBI nucleotide ids (.ids), 1 per line')
    parser.add_argument('--evalue', required=True, type=float, help='Evalue cutoff for diamond and for hmmer3, recommend 0.1 relaxed paramter')
    parser.add_argument('--threads', required=True, type=int, help='Threads for "diamond" and "hmmer3"')
    parser.add_argument('--genepred', required=True, help='Gene prediction tool "metagenemark" or "prodigal"')
    args = parser.parse_args()
    return args

def get_fasta_from_ids(ids_lst, prefix):
    print('ids list: {}'.format(ids_lst))
    handle = Entrez.efetch(db='nucleotide', id=','.join(ids_lst), rettype='fasta', retmode='text')
    SeqIO.write(SeqIO.parse(handle, 'fasta'), '{}.fa'.format(prefix), 'fasta')

def print_and_run(cmd):
    print(cmd +'\n')
    os.system(cmd)

def run_geneprediction(in_file, out_file, genepred):
    if genepred == 'metagenemark':
        cmd = '{} {} -m {} -A {}'.format(METAGENEMARK_PATH, in_file, METAGENEMARK_MODEL_PATH, out_file)
        print_and_run(cmd)
    elif genepred == 'prodigal':
        cmd = '{} -i {} -a {}.tmp > {}.prodical_log'.format(PRODIGAL_PATH, in_file, out_file, out_file)
        print_and_run(cmd)
        print('reformating prodigal faa file')
        with open(out_file+'.tmp') as in_f, open(out_file, 'w') as out_f:
            for line in in_f:
                if line.startswith('>'):
                    line = '{} {}\n'.format(line.split(';')[0].replace(' ', '_'), line.split(' #')[0])
                out_f.write(line)
        cmd = 'rm {}.tmp'.format(out_file)
        print_and_run(cmd)
    else:
        raise ValueError('How did you get here?')

def run_diamond(in_file, evalue, out_file, threads):
    NR_PATH = NR_PATH.rsplit('.dmnd', 1)[0]
    cmd = '{} blastp --query {} --db {} --evalue {} --outfmt 6 --out {} --threads {}'.format(DIAMOND_PATH, in_file, NR_PATH, evalue, out_file, threads)
    print_and_run(cmd)

def run_hmmer3(in_file, evalue, out_file, threads):
    cmd = '{} --tblout {} -E {} --cpu {} {} {}'.format(HMMSEARCH_PATH, out_file, evalue, threads, PVOGS_HMMS_PATH, in_file)
    print_and_run(cmd)

def get_record_info(s, line):
    count = 0
    while True:
        try:
            handle = Entrez.efetch(db='protein', id=s, retmode='xml')
            record = Entrez.read(handle)
            try:
                return record[0]['GBSeq_definition']
            except IndexError:
                print('Warning: index error when trying to parse Entrez handle, returning NA')
                return 'NA'
        except urllib.error.HTTPError:
            time.sleep(5)  # wait a bit
            count += 1
            if count == 10:  # don't try more than 10 attempts
                print('Warning: urllib.error.HTTPError for id {} at line {}, returning NA'.format(s, line))
                return 'NA'
        except Exception as e:
            print('Warning: {}: {}: for id "{}" at line "{}", returning NA'.format(e.__class__.__name__, e, s, line))
            return 'NA'

def parse_diamond_hmmer3_outputs(diamond_file, hmmer3_file, faa_file, out_file):

    print('Get gene num to scaf info...\n')
    gene2scaf_d = OrderedDict()
    with open(faa_file) as in_f:
        for line in in_f:
            if line.startswith('>'):
                gene, scaf = line.rstrip().split(None, 1)
                gene = gene.lstrip('>')
                scaf = scaf.lstrip('>')
                gene2scaf_d[gene] = scaf

    print('Parsing diamond output...\n')
    gene2diamondAnnot_d = {}
    with open(diamond_file) as in_f:
        groups = groupby((line for line in in_f), key = lambda x: x.split()[0])
        for gene, group in groups:
            first_line = next(group)
            NCBI_id = first_line.split('\t')[1].split('|')[3]
            first_record = get_record_info(NCBI_id, first_line)
            if not first_record.startswith(NR_ANNOTS_TO_IGNORE):
                gene2diamondAnnot_d[gene] = first_record
            else:
                for line in group:
                    NCBI_id = line.split('\t')[1].split('|')[3]
                    record = get_record_info(NCBI_id, line)
                    if not record.startswith(NR_ANNOTS_TO_IGNORE):
                        gene2diamondAnnot_d[gene] = record
                        break
                else:  # If all hits starts with forbiden terms, take the first one
                    gene2diamondAnnot_d[gene] = first_record

    print('Get pVOG table info...\n')
    vog2Annot_d = {}
    with open(PVOGS_TABLE_PATH) as in_f:
        for line in in_f:
            line = line.split(None, 6)
            vog = line[0]
            annots = line[6].split(';')
            annots.pop()  # remove POG mapping
            annot = next((annot for annot in annots if not annot.startswith(PVOG_ANNOTS_TO_IGNORE)), annots[0])
            vog2Annot_d[vog] = annot

    print('Parsing HMMER3 output...\n')
    gene2Hmmer3Annot_d = {}
    with open(hmmer3_file) as in_f:
         for line in in_f:
             if line.startswith('#'):
                 continue
             line = line.split()
             gene, vog, evalue = line[0], line[2], line[4]
             evalue = float(evalue)
             if gene in gene2Hmmer3Annot_d:
                 if evalue > gene2Hmmer3Annot_d[gene][1]:
                     continue
             gene2Hmmer3Annot_d[gene] = vog, evalue
       
    print('Merge annotations...\n')
    with open(out_file, 'w') as out_f:
        out_f.write('gene_id\tscaffold_info\tNR_hit_info\tVOG_hit\tVOG_desc\n')
        for gene, scaf in gene2scaf_d.items():
            NR_hit = gene2diamondAnnot_d.get(gene, 'NA')
            VOG_hit = '_'.join(map(str, gene2Hmmer3Annot_d.get(gene, ('NA', 'NA'))))
            VOG_info = vog2Annot_d.get(VOG_hit.split('_')[0], 'NA')
            gene = gene.replace('GeneMark.hmm|','')

            out_f.write('\t'.join([gene, scaf, NR_hit, VOG_hit, VOG_info]) + '\n')
    

def main():
    args = parse_args()

    if not args.input.endswith(('.ids', '.fa')):
        raise ValueError('"in" file must end with ".fa" or ".ids"')
    prefix, suffix = args.input.rsplit('.', 1)  

    if not args.genepred in ('metagenemark', 'prodigal'):
        raise ValueError('Gene prediction tool must be "metagenemark" or "prodigal"')

    print('Prefix: {}, Suffix: {}\n'.format(prefix, suffix))
    print('Evalue cutoff: {}\n'.format(args.evalue))
    print('Threads: {}\n'.format(args.threads))

    if suffix == 'ids':
        print('Getting nucleotide fasta sequence from ids...\n')
        with open(args.input) as f:
            ids_lst = [line.strip() for line in f]
            get_fasta_from_ids(ids_lst, prefix)
    elif suffix == 'fa':
        print('Input is nucleotide fasta sequence...\n')
    else:
        raise ValueError('How did you get here??')

    params = OrderedDict()
    params['fa_file'] = '{}.fa'.format(prefix)
    params['faa_file'] = '{}.faa'.format(prefix)
    params['diamond_file'] = '{}.diamondVsNR.out'.format(prefix)
    params['hmmer3_file'] =  '{}.HMMER3VsPVOGs.out'.format(prefix)
    params['merged_file'] = '{}.merged_final.out'.format(prefix)

    print('File names..')
    for name, param in params.items():
        print('{}: {}'.format(name, param))
    print()

    done_hmmer3 = False

    if all(os.path.isfile(f_name) for f_name in [params['fa_file'], params['faa_file'], params['diamond_file']]):
        print('Diamond output file already exists, skipping these expensive steps\n')
        if os.path.isfile(params['hmmer3_file']):
            done_hmmer3 = True

    else:
        print('Running gene prediction...\n')
        run_geneprediction(params['fa_file'], params['faa_file'], args.genepred)

        print('Running diamond vs NR...\n')
        run_diamond(params['faa_file'], args.evalue, params['diamond_file'], args.threads)

    if done_hmmer3:
        print('Hmmer3 output file already exists, skipping step\n')
    else:
        print('Running HMMER3 vs pVOGs...\n')
        run_hmmer3(params['faa_file'], args.evalue, params['hmmer3_file'], args.threads)

    print('Parsing diamond and HMMER3 outputs...\n')
    parse_diamond_hmmer3_outputs(params['diamond_file'], params['hmmer3_file'],
                                 params['faa_file'], params['merged_file'])

if __name__ == '__main__':
    print('Starting...')
    main()
    print('Done...')
