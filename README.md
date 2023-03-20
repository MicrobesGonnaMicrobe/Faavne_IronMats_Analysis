# Microbial iron mats at Fåvne hydrothermal vent

Tools and commands used for analyses of microbial iron mats at Fåvne hydrothermal vent field. Article in preparation.

## Genome-resolved metagenomics

## MAGs
* `GTDB-Tk v1.3.0`: https://github.com/Ecogenomics/GTDBTk
* `ncbi-genome-download v0.3.1`: https://github.com/kblin/ncbi-genome-download
* `dRep v3.2.2`: https://github.com/MrOlm/drep

### Manual refinement of MAGs
* `anvi’o v7.1`: https://anvio.org

### Download reference genomes 
### Get genome metadata

### Calculate ANI
### Calculate AAI

## Phylogenomics

* `anvi’o v7.1`: https://anvio.org
* `IQ-TREE v2.1.4`: http://www.iqtree.org

### Make list of genes
List available HMM sources in the contigs database
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --list-hmm-sources

Alternative HMM sources (add your own external database)
- Get HMMs from GTDB single copy marker genes database repository, r202: https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_all/bac120_msa_marker_genes_all_r202.tar.gz
- Import to anvi'o in a specific format: https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Bacteria_71
- To get descriptions of marker genes from GTDB, you can find the markers and their descriptions in gtdbtk output under "align/intermediate_results"

Build HMMs
* `hmmbuild 3.2.1`
```bash
#!/bin/bash
for bins in *.faa
do
	basename=${bins%.faa}
	id=${basename}
   hmmbuild --amino ${id}.hmm "${bins}"
done
```

Run GTDB marker genes HMMs
```bash
for i in *.db; do anvi-run-hmms -c $i -H /export/work_cgb/Shared/Anvio_GTDB_markers_r202/GTDB_bac120_r202 --num-threads 6; done
```

### Get amino acid sequences
```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --hmm-source GTDB_bac120_r202 --list-available-gene-names
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --hmm-source Bacteria_71 --list-available-gene-names
```

Decide which set of markers to use
- only select the single copy marker genes (check with a matrix if they are in duplicates)
- select the ones that are mostly present in the genomes
- more good quality marker genes, the better

To check presence of marker genes in all genomes (matrix):
```bash
anvi-script-gen-hmm-hits-matrix-across-genomes --external-genomes external_selected_genomes.txt --hmm-source GTDB_bac120_r202 -o Zeta_GTDB_bac120_markers_matrix.txt
```

Choose a subset of those genes and get them in one fasta file (without alignment and concatenation)
```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt -o Zetaproteobacteria_selectedmarkers_Bacteria71.fa --hmm-source Bacteria_71 --gene-names Zetaproteobacteria_selectedmarkers_Bacteria71.txt --return-best-hit --get-aa-sequences

anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt -o Zetaproteobacteria_selectedmarkers_GTDB_bac120.fa --hmm-source GTDB_bac120_r202 --gene-names Zetaproteobacteria_selectedmarkers_GTDB_bac120.txt --return-best-hit --get-aa-sequences
```

### Trim with trimal

### Align with mafft

### Concatenate

### Build tree

## 16S rRNA gene phylogeny

## Phylogeny of NiFe uptake hydrogenase and cyc2

## Annotation: Metabolism and other genes

### Iron oxidation

#### FeGenie: iron genes and metabolism
- Tutorial: https://github.com/Arkadiy-Garber/FeGenie/wiki/Tutorial
```bash
    FeGenie.py -bin_dir /export/work_cgb/Petra/Zetaproteobacteria/fasta_files_bins/ -bin_ext fa -out fegenie_output 
```

Sample command (if providing gene-calls):

```bash
    FeGenie.py -bin_dir Anvio_genecalls_aa -bin_ext fa -out Anvio_genecalls_aa_fegenie3 --orfs
```

### BacMet: Heavy metal resistance genes and biocides
Experimentally Verified Resistance Gene Information:
http://bacmet.biomedicine.gu.se/database_browse.pl

Predicted:

```bash
    for i in *.faa; do blastp -query $i -db BacMet2_predicted_database.fasta -out BacMet_Results/${i%.faa}_BacMet_Pred -outfmt 6 -evalue 1e-6 -num_threads 20; done
    
    cd BacMet_Results
    
    for i in *_BacMet_Pred; do perl best_blast.pl $i ${i%}_best; done
    
    cat *Pred_best > Faavne_BS4_allZetaproteobacteria_BacMet_Pred.tsv
```

What database and threshold criteria to use?
"... for screening of metagenomes and applying a uniform, more relaxed cut-off criteria, we recommend to use predicted database with a criterion of about 85-90% sequence identity with the full length coverage of the short reads (75-300+ bps) against resistance genes. Similarly, the predicted database is useful for investigating the presence of resistance genes in genomes of species that are not closely related to the ones our experimental knowledge about the genes resistance function are derived from." - http://bacmet.biomedicine.gu.se/FAQS.html#experimetal

## Other

### Predicting optimal growth temperatures (OGT)
https://github.com/DavidBSauer/OGT_prediction

1. Place each genome in a folder in the prediction directory called "genomes/XXX/" where XXX is the name of the species.  Genomes should be ziped in .gz format.
2. Create a tab separated file of the genomes and species pairs. Provide this file in place of genomes_retrieved.txt. 
3. Create a file for the taxonomic classification of each species. 
4. Run the script: Run the prediction script to predict the OGT of each species. If you want to account for absence of 16S rRNA gene and genome incompleteness (usual with MAGs), use specified regression models modified to exclude 16S rRNA gene and genome size: https://github.com/DavidBSauer/OGT_prediction/tree/master/data/calculations/prediction/regression_models

    ```
    python3 ~/OGT_prediction/OGT_prediction/prediction/prediction_pipeline.py regression_models_Bacteria_uncultivated/ Faavne_MAGs.txt Faavne_MAGs_taxonomic3.txt
    ```

5. Results in newly_predicted_OGTs.txt

