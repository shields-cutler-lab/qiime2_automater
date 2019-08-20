#!/usr/bin/env Python

import os
import sys
import pandas as pd

usage = 'python qiime2_automater.py otu_table_fp metadata_fp taxonomy_fp rarefaction_level_counts phylo_tree_OR_NONE PROK_OR_NONE'

## Refactored for QIIME 2018.11 changes ##
## Not compatible with previous versions ##

if sys.argv[1] == "--help":
    print(usage)
    sys.exit()

otusfp = sys.argv[1]
metadatafp = sys.argv[2]
taxonomyfp = sys.argv[3]
rarefy = int(sys.argv[4])

with open('automater_cmd_run.txt', 'w') as cmdsave:
    cmdsave.write(' '.join(sys.argv))
os.mkdir('qiime2_artifacts')
os.mkdir('qiime2_viz')
os.mkdir('alpha_diversity')
os.mkdir('beta_diversity')

if not taxonomyfp.endswith('.qza'):
    os.system(' '.join(['qiime tools import --input-path', taxonomyfp, "--type 'FeatureData[Taxonomy]' --output-path qiime2_artifacts/taxonomy_map.qza --input-format 'HeaderlessTSVTaxonomyFormat'"]))
    taxonomyfp = 'qiime2_artifacts/taxonomy_map.qza'

phylo_go = False
phylo = sys.argv[5]
if not str(phylo) == 'NONE':
    if phylo.endswith('.tre') or phylo.endswith('.tree') or phylo.endswith('.nwk'):
        phylo_go = True
        os.system(' '.join(['qiime tools import --input-path', phylo, "--output-path qiime2_artifacts/rooted_tree.qza --input-format 'NewickFormat' --type 'Phylogeny[Rooted]'"]))
    elif phylo.endswith('.qza'):
        phylo_go = True
        os.system(' '.join(['cp',phylo,'qiime2_artifacts/rooted_tree.qza']))
    else:
        print('Error, you must provide a *.tree, *.tre, or *.nwk file for argument 5, or enter "NONE" if non-phylogenetic analysis')
        sys.exit()

if not otusfp.endswith('.qza'):
    # Convert the table to hdf5 BIOM
    os.system(' '.join(['biom convert -i', otusfp, '-o otu_table.biom --table-type="OTU table" --to-hdf5']))
    # Convert to qiime2 artifact
    os.system(' '.join(["qiime tools import --input-path otu_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path qiime2_artifacts/otu_table.qza"]))
else:
    os.system(' '.join(['cp',otusfp,'qiime2_artifacts/otu_table.qza']))
os.system(' '.join(["qiime taxa barplot --i-table qiime2_artifacts/otu_table.qza --i-taxonomy",taxonomyfp,"--m-metadata-file", metadatafp, "--o-visualization qiime2_viz/barplot.qzv"]))

# Generate and export tables collapsed to phylum, class, order, family, genus, and species
for L in range(2,8,1):
    table_nameq = ''.join(['otu_table_L',str(L),'.qza'])
    outpathq = os.path.join('qiime2_artifacts',table_nameq)
    table_name = ''.join(['otu_table_L',str(L),'.txt'])
    os.system(' '.join(["qiime taxa collapse --i-table qiime2_artifacts/otu_table.qza --i-taxonomy", taxonomyfp, "--p-level", str(L), "--o-collapsed-table", outpathq]))
    os.system(' '.join(["qiime tools export --input-path", outpathq, "--output-path ."]))
    os.system(' '.join(["biom convert -i feature-table.biom -o", table_name, "--to-tsv"]))
    os.system("rm feature-table.biom")

if phylo_go:  # Do the phylogenetic core metrics
    # If using PROK database, you need to replace the underscores in the OTU names to match the tree once imported
    if sys.argv[6] == 'PROK':
        outotu = 'otu_table_phylo.txt'
        with open(otusfp, 'r') as inf, open(outotu, 'w') as outf:
            df = pd.read_table(inf, header=0, index_col=0, engine='c')
            otu_labels = list(df.index)
            new_labels = []
            for i in otu_labels:
                new_label = i.replace('_',' ')  # Replacing underscores with spaces
                new_labels.append(new_label)
            df.index = new_labels
            df.index.name = '#OTU ID'
            df.to_csv(outf, sep='\t')
        # Convert the table to hdf5 BIOM
        os.system(' '.join(['biom convert -i', outotu, '-o otu_table_phylo.biom --table-type="OTU table" --to-hdf5']))
        # Convert to qiime2 artifact
        os.system(' '.join(["qiime tools import --input-path otu_table_phylo.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path qiime2_artifacts/otu_table_phylo.qza"]))
        phylotu = 'qiime2_artifacts/otu_table_phylo.qza'
    else:
        phylotu = otusfp
    os.system(' '.join(["qiime diversity core-metrics-phylogenetic --i-table qiime2_artifacts/otu_table.qza --i-phylogeny qiime2_artifacts/rooted_tree.qza --m-metadata-file", metadatafp, "--output-dir qiime2_analysis/core_metrics_phylo --p-sampling-depth", str(rarefy)]))
    dists = ['bray_curtis_distance_matrix.qza','jaccard_distance_matrix.qza','unweighted_unifrac_distance_matrix.qza','weighted_unifrac_distance_matrix.qza']
    alphas = ['evenness_vector.qza','faith_pd_vector.qza','observed_otus_vector.qza','shannon_vector.qza']
    vizzs = ['jaccard_emperor.qzv','unweighted_unifrac_emperor.qzv','weighted_unifrac_emperor.qzv','bray_curtis_emperor.qzv']
    for d in dists:
        namer = os.path.join('beta_diversity','.'.join([d.split('.')[0], 'txt']))
        inpath = os.path.join('qiime2_analysis/core_metrics_phylo',d)
        os.system(' '.join(["qiime tools export --input-path", inpath, "--output-path beta_diversity"]))
        os.rename('beta_diversity/distance-matrix.tsv', namer)
    for a in alphas:
        namer = os.path.join('alpha_diversity', '.'.join([a.split('.')[0], 'txt']))
        inpath = os.path.join('qiime2_analysis/core_metrics_phylo',a)
        os.system(' '.join(["qiime tools export --input-path", inpath, "--output-path alpha_diversity"]))
        os.rename('alpha_diversity/alpha-diversity.tsv', namer)
    combined_alpha = os.path.join('alpha_diversity','alphadiv.txt')
    newalphas = ['evenness_vector.txt','faith_pd_vector.txt','observed_otus_vector.txt','shannon_vector.txt']
    with open('alpha_diversity/evenness_vector.txt', 'r') as ev_in, open('alpha_diversity/faith_pd_vector.txt', 'r') as pd_in, open('alpha_diversity/observed_otus_vector.txt', 'r') as ob_in, open('alpha_diversity/shannon_vector.txt', 'r') as sh_in:
        alphains = [ev_in, pd_in, ob_in, sh_in]
        alphasers = [pd.read_csv(ainf, header=0, index_col=0, delimiter='\t', squeeze=False) for ainf in alphains]
    alphadfs = [pd.DataFrame(ser) for ser in alphasers]
    alphadf = alphasers[0].join(alphadfs[1:])
    alphadf.to_csv(combined_alpha, index_label='sampleID', sep='\t')
    for v in vizzs:
        inpath = os.path.join('qiime2_analysis/core_metrics_phylo',v)
        os.system(' '.join(["mv", inpath, 'qiime2_viz/.']))
    print('To visualize, run: qiime tools view [filepath]')
else:  # Do the non-phylogenetic metrics
    os.system(' '.join(["qiime diversity core-metrics --i-table qiime2_artifacts/otu_table.qza --m-metadata-file", metadatafp, "--output-dir qiime2_analysis/core_metrics --p-sampling-depth", str(rarefy)]))
    dists = ['bray_curtis_distance_matrix.qza','jaccard_distance_matrix.qza']
    alphas = ['evenness_vector.qza','observed_otus_vector.qza','shannon_vector.qza']
    vizzs = ['jaccard_emperor.qzv','bray_curtis_emperor.qzv']
    for d in dists:
        namer = os.path.join('beta_diversity','.'.join([d.split('.')[0], 'txt']))
        inpath = os.path.join('qiime2_analysis/core_metrics',d)
        os.system(' '.join(["qiime tools export --input-path", inpath, "--output-path beta_diversity"]))
        os.rename('beta_diversity/distance-matrix.tsv', namer)
    for a in alphas:
        namer = os.path.join('alpha_diversity','.'.join([a.split('.')[0], 'txt']))
        inpath = os.path.join('qiime2_analysis/core_metrics',a)
        os.system(' '.join(["qiime tools export --input-path", inpath, "--output-path alpha_diversity"]))
        os.rename('alpha_diversity/alpha-diversity.tsv', namer)
    combined_alpha = os.path.join('alpha_diversity','alphadiv.txt')
    newalphas = ['evenness_vector.txt','faith_pd_vector.txt','observed_otus_vector.txt','shannon_vector.txt']
    with open('alpha_diversity/evenness_vector.txt', 'r') as ev_in, open('alpha_diversity/observed_otus_vector.txt', 'r') as ob_in, open('alpha_diversity/shannon_vector.txt', 'r') as sh_in:
        alphains = [ev_in, ob_in, sh_in]
        alphasers = [pd.read_csv(ainf, header=0, index_col=0, delimiter='\t', squeeze=False) for ainf in alphains]
    alphadfs = [pd.DataFrame(ser) for ser in alphasers]
    alphadf = alphasers[0].join(alphadfs[1:])
    alphadf.to_csv(combined_alpha, index_label='sampleID', sep='\t')
    for v in vizzs:
        inpath = os.path.join('qiime2_analysis/core_metrics',v)
        os.system(' '.join(["mv", inpath, 'qiime2_viz/.']))
    print('To visualize, run: qiime tools view [filepath]')



