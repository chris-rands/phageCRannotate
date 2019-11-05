import sys
import os

MIN_NUM_ORFS = 2
MIN_ANNOTATED_ORFS = 2

# Only looking for tailed prophages
TERMS_EITHER_OR = [{'head ', 'capsid', 'clp protease', 'prohead protease'}, # HEAD
                   {'head-tail', 'head tail', 'head/tail'}, # HEAD-TAIL
                   {'portal'},  # PORTAL
                   {' tail', '\ttail', 'tail fibre', 'tail tube', 'tail sheath',  # TAIL
                     'tail tape measure', 'base plate', 'base-plate', 'baseplate'},  # TAIL
                   {'terminase'},  # PACKAGING
                   {'helicase', 'primase'},  # REPLICATION
                   {'lysis', 'lysin', 'holin'}, # Excludes likely bacterial prots, 'lysm', 'lysr'},  # LYSIS
                   {'repressor', 'antirepressor', 'anti-repressor', 'integrase'}]  # LYSOGENY

TERMS_TO_INCLUDE = sorted(set.union(*map(set, TERMS_EITHER_OR)))  # used for filtering
COLS = ('blue', 'purple', 'maroon1', 'light blue', 'pink', 'orange', 'coral', 'grey')
assert len(TERMS_EITHER_OR) == len(COLS)

TERM2COL = {}
for terms, col in zip(TERMS_EITHER_OR, COLS):
    for term in terms:
        TERM2COL[term] = col

TERM2COL['trna'] = 'yellow'

#Bacterial gene annots- These have been moved to a seperate script
for term in ('toxin', 'virulence',  # TOXIN/VIRULENCE
             'beta-lactamase', 'beta lactamase',  # ANTIBIOTIC RESITANCE (likely highly incomplete)
             'aminoglycoside', 'chloramphenicol',
             'glycopeptide', 'quinolone',
             'tetracycline', 'macrolide',
             'ansamycin', 'streptogramin',
             'lipopeptide', 'vanomycin'
             'efflux-pump', 'efflux pump'):
    TERM2COL[term] = 'red'

# Other bacterial genes
for term in ('16S', '30S', '23S', '50S'):
    TERM2COL[term] = 'red'

for term in ('recombinase', 'transposase'):
    TERM2COL[term] = 'green'

# Preferentially annotated ARGs/toxin/virulence prots
TERM2COL_lst = sorted(TERM2COL, key= lambda x: TERM2COL[x] != 'red')
# Make holin first, since toxin genes that are in fact holins are NOT bacterial genes
TERM2COL_lst.insert(0, TERM2COL_lst.pop(TERM2COL_lst.index('holin')))
# Put helicase/primase before VapE proteins
TERM2COL_lst.insert(0, TERM2COL_lst.pop(TERM2COL_lst.index('helicase')))
TERM2COL_lst.insert(0, TERM2COL_lst.pop(TERM2COL_lst.index('primase')))

# Move tail generic terms to end
TERM2COL_lst.append(TERM2COL_lst.pop(TERM2COL_lst.index(' tail')))
TERM2COL_lst.append(TERM2COL_lst.pop(TERM2COL_lst.index('\ttail')))

# Should I do the same for head terms?

print(TERM2COL_lst)

'''
Input format:
gene_id scaffold_info   NR_hit_info VOG_hit VOG_desc
GCA_000024945.1|Veillonella_1_#_1_#_330_#_-1_#_ID=1_1   GCA_000024945.1|Veillonella_1   MULTISPECIES: thioredoxin [Veillonella] VOG5240_0.04    hypothetical protein [n=27]
GCA_000024945.1|Veillonella_2_#_788_#_1891_#_-1_#_ID=1_2    GCA_000024945.1|Veillonella_2   Nif3-like dinuclear metal center hexameric protein [Veillonella parvula]    NA_NA   NA
GCA_000024945.1|Veillonella_3_#_1912_#_2634_#_-1_#_ID=1_3   GCA_000024945.1|Veillonella_3   SAM-dependent methyltransferase [Veillonella parvula]   NA_NA   NA
GCA_000024945.1|Veillonella_4_#_2967_#_3530_#_-1_#_ID=1_4   GCA_000024945.1|Veillonella_4   GTP cyclohydrolase I FolE [Veillonella parvula] VOG3222_2.7e-55 GTP cyclohydrolase I [n=12]
GCA_000024945.1|Veillonella_5_#_3617_#_4360_#_-1_#_ID=1_5   GCA_000024945.1|Veillonella_5   4Fe-4S cluster-binding domain-containing protein [Veillonella parvula]  VOG0998_4.7e-33 queuosine biosynthesis QueE radical SAM [n=6]
GCA_000024945.1|Veillonella_6_#_4347_#_4685_#_-1_#_ID=1_6   GCA_000024945.1|Veillonella_6   6-carboxytetrahydropterin synthase QueD [Veillonella parvula]   VOG5156_2.3e-20 gp4 [n=5]
GCA_000024945.1|Veillonella_7_#_4891_#_6051_#_-1_#_ID=1_7   GCA_000024945.1|Veillonella_7   RNA polymerase, sigma 70 subunit, RpoD subfamily [Veillonella parvula DSM 2008] VOG0052_1.7e-56 RNA polymerase sigma factor [n=33]
GCA_000024945.1|Veillonella_8_#_6080_#_7852_#_-1_#_ID=1_8   GCA_000024945.1|Veillonella_8   DNA primase [Veillonella parvula]   VOG4551_7.6e-46 DNA primase [n=259]
GCA_000024945.1|Veillonella_9_#_8110_#_9654_#_-1_#_ID=1_9   GCA_000024945.1|Veillonella_9   ribonuclease Y [Veillonella montpellierensis]   VOG0299_0.0056  endopeptidase [n=16]
'''

def group_seqs(in_f, id_suffix):
    group, n = [], 0
    for line in in_f:
        id_ = line.split()[1]
        generic, num = id_.rsplit('_', 1)
        if num == '1' and group:
            yield 'id_#{}{}'.format(n, id_suffix), group
            n += 1
            group = [line]
        else:
            group.append(line)
    if group:
        yield 'id_#{}{}'.format(n, id_suffix), group


def get_groups(in_file, filter_phages, id_suffix):
    with open(in_file) as in_f:
        header = next(in_f)
        yield None, header
        for id_, group in group_seqs(in_f, id_suffix):
            len_ = len(group)
            if len_ < MIN_NUM_ORFS:
                print('Warning: group with id: {} :contains too few ({}) ORFs, skipping'.format(id_, len_))
                continue
            seen_terms = set()
            for line in group:
                line = line.lower()
                for term in TERMS_TO_INCLUDE:
                    if term in line:
                        seen_terms.add(term)
            print(seen_terms)
            if len(seen_terms) < MIN_ANNOTATED_ORFS:
                print('Warning: group with id: {} :contains too few annotated ({}) ORFs, skipping'.format(id_, len(seen_terms)))
                continue
            if filter_phages:
                if all(any(term in or_terms for term in seen_terms) 
                       for or_terms in TERMS_EITHER_OR):
                    print('Keeping group with id: {}'.format(id_))
                    for line in group:
                        yield id_, line
                else:
                    print('Warning group with id: {} does not contain neccesray annotated protins'.format(id_))
            else:
                print('Ouputing all phages, id: {}'.format(id_))
                for line in group:
                    yield id_, line
                    
def main(in_file, out_file, id_suffix='N', filter_phages=False):
    print('write out_files')
    out_file = out_file if filter_phages else out_file + '_nonfiltered' 
    out_dir_name = 'col_files' if filter_phages else 'col_files'+'_nonfiltered'
    print('Making: {}'.format(out_dir_name))
    os.makedirs(out_dir_name)
    with open(out_file,'w') as out_f:
        for id_, line in get_groups(in_file, filter_phages, id_suffix):
            out_f.write(line)
            if line.startswith('gene_id'):
                continue  # skip header
            try:
                gca, taxon = line.split('|')[0], line.split('|')[1].split('_')[0]
            except IndexError:  # not GCA format
                gca, taxon = line.split('_')[0], id_suffix

            # Coping with NC ids
            if gca.startswith('NC'):  # HACK TO COPE WITH BAD BEHAVING NC TREATMENT- this is Refseq records
                gca = '_'.join(line.split('_')[:2])
                assert len(gca.split()) == 1  # should be no whitepsace in the id
                taxon = id_suffix  # give up on taxon naming for RefSeq records

            #print('id_: {}, gca: {}, taxon: {}'.format(id_, gca, taxon))
            with open(os.path.join(out_dir_name, '{}.{}.{}.colors'.format(id_, gca, taxon)), 'a') as col_out_f:
                #  c("Start", "End", "Strand", "Colour", "Annot")
                start, end, strand = line.split()[0].split('_#_')[1:4]
                strand = {'1': '+', '-1': '-'}[strand]
                tmp_line = line.lower()
                annot = next((term for term in TERM2COL_lst if term in tmp_line), None)
                if annot is None:
                    continue  # skipping uninformative prots
                col = TERM2COL[annot]
                col_out_f.write('|'.join([start, end, strand, col, annot]) + '\n') 
                    

if __name__ == '__main__':
    print('starting...')
    main(*sys.argv[1:])
    print('done!')
