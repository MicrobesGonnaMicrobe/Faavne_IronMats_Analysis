# Microbial iron mats at Fåvne hydrothermal vent

Tools and commands used for analyses of microbial iron mats at Fåvne hydrothermal vent field. Article in preparation.

## Genome-resolved metagenomics

## MAGs
* `GTDB-Tk v1.3.0`: https://github.com/Ecogenomics/GTDBTk
* `ncbi-genome-download v0.3.1`: https://github.com/kblin/ncbi-genome-download
* `dRep v3.2.2`: https://github.com/MrOlm/drep

### Manual refinement of MAGs
* `anvi’o v7.0`: https://anvio.org

### Download reference genomes 
### Get genome metadata

### Calculate ANI
### Calculate AAI

## Phylogenomics

* `anvi’o v7.0`: https://anvio.org
* `IQ-TREE v2.1.4`: http://www.iqtree.org

### Make list of genes
### Get amino acid sequences
```bash
anvi-get-sequences-for-hmm-hits 
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

