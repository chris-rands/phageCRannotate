# phageCRannotate

**Description**
Take as input NCBI nucleotide id(s) or FASTA sequences of phages and outputs phage 'genome map' annotations.

**Dependancies**
- python3 with Biopython
- R with viz
- diamond
- hmmer3
- prodigal (or metagenemark)

**Setup**
Edit the paths at the top of `python/faOrNCBINucIds_to_phageProtAnnots.py` to specify the correct installations. Create a an input file, either 1) suffixed `.ids` containing a list of NCBI nucleotide ids of phages one per line, or 2) suffixed `.fa` containing sequences in FASTA format, one phage per sequence. 

**Example test**
Example test of annotating phage lambda. See `test_outputs/` for outputs

1st step: download sequence and annotate via pVOG HMMs and NR protein sequences
(by far the slowest step due to diamond, if using the full NR database: uses ~12G and runs for several hours;
note it scales better you might expect though due to the diamond implementation)
```python3 python/faOrNCBINucIds_to_phageProtAnnots.py --input test_input_lambda_phage.ids --evalue 0.1 --threads 8 --genepred prodigal```

2nd step: collate the annotation information
```python3 python/phageAnnotsFilter.py test_input_lambda_phage.merged_final.out test_input_lambda_phage.merged_final.phageAnnotsFilter.out > test_input_lambda_phage.merged_final.phageAnnotsFilter.log```

3rd step: plot the genome maps
```Rscript R/gvizPlotPhages.R col_files_nonfiltered/*```

See: `test_outputs/` for ouputs from each of these steps

**Example plot genome map with own annotation**

You can use your own annotations and then just use the Rscript to plot the genome map viz. In this case create a directory of files (one per phage) with the following format
`start|end|strand|colour|annotation` for each file, e.g. see `annots_example/`:

```Rscript R/gvizPlotPhages.R annots_example/*```

See: annots_example_GENOME_MAP.pdf 

**Disclaimer**
This in largely untested software licensed under the MIT License Copyright (c) 2019, Chris Rands.
