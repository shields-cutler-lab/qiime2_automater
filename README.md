# QIIME 2 Automater
Robin Shields-Cutler

Inspired by the frustrations of continuously going through the tedious commandline steps required to go from tab-delimited OTU table to a suite of alpha & beta diversity results, rank-collapsed taxonomy tables, taxa bar plots, etc., and then again to convert those diversity metrics, taxonomy tables, etc., back into a text format that other software (e.g. R) can easily read.

### Version note:
Occassionally QIIME 2 command arguments change. This script version is tested with QIIME2 v2018.11, after some significant changes to common commands from a previous version. Will try to test on newer versions periodically.

## The inputs:
To use: You must have:
* an OTU table in tab-delimited format (or .qza artifact)
* a metadata table in tab-delimited format (where sample names match the OTU table exactly)
* a headerless taxonomy map (tab-delimited _or_ .qza artifact) for your database
* if phylogenetic, a QIIME-friendly Newick format tree (or .qza rooted tree artifact) that matches your OTU references
* your rarefaction level, in counts

## Run the script:
```shell
$ python qiime2_automater_2018.11_automater.py --help`
python qiime2_automater.py otu_table_fp metadata_fp taxonomy_fp rarefaction_level_counts phylo_tree_OR_NONE PROK_OR_NONE'`
# Positional arguments - make sure they are in the right order. (What. I got lazy.)
# If not running phylogenetic, enter NONE instead of the tree file path
# If using the PROK database from Knights Lab, enter PROK here (issue with tree tip reformatting), otherwise NONE.
```

## The outputs:
In your current working directory, will generate and save:
* OTU table, biom format
* Taxa tables collapsed to levels 2 through 7
* core div alpha diversity metrics, merged into one text file with a column for each metric
* QIIME 2 visualizations for core div beta diversity metrics (Emperor plots), and the interactive Taxa Barplot
* All QIIME 2 intermediate artifacts, for future use
