# Microbial iron mats at Fåvne hydrothermal vent

Tools and commands used for analyses of microbial iron mats at Fåvne hydrothermal vent field. Article in preparation.

## Genome-resolved metagenomics

## MAGs

### Manual refinement of MAGs
* `anvi’o v7.1`: https://anvio.org

### Dereplication
* `dRep v3.2.2`: https://github.com/MrOlm/drep
Chosing parameters: https://drep.readthedocs.io/en/latest/choosing_parameters.html

At 95 ANI clustering:
```bash
dRep dereplicate -p 8 -g genomes_50_10_paths.txt -comp 50 -con 10 -sa 0.95 genomes_50_10_dRep_95
```

### Download reference genomes with metadata
* `ncbi-genome-download v0.3.1`: https://github.com/kblin/ncbi-genome-download

Make a taxid list for ncbi-genome-download
- "On first use, a small sqlite database will be created in your home directory by default (change the location with the --database flag). You can update this database by using the --update flag. What it's doing it's downloading the new version of this: wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
- This script lets you find out what TaxIDs to pass to ngd:
```bash
python gimme_taxa.py -o zetataxafile_2021_10_19.txt 580370 --update
```

Get unique taxids from the second column (descendent_taxid), extract second column and turn to comma separated list:
```bash
awk '{print $2}' file.txt | paste -s -d, -
```

Download taxid genomes with metadata
```bash
ncbi-genome-download --section genbank --format fasta,assembly-report,assembly-stats -t 314345,1188231,904974,1921010,1921086,1921087,2528639,933853,999537,1353261,1871319,1896269,1896270,1896271,1896272,1896273,1896274,1896276,1896277,1896278,1896279,1896280,1896281,1896282,1896283,1896284,1896297,1896298,1896300,1896301,1896303,1896304,1946828,1946829,1946830,1946831,1946832,2608715,2608716,2650971,2306050,2614254,587656,933866,933867,933868,1871100,1131282,1131283,1131284,1131285,1131286,1131287,1131288,1131289,1131290,1313310,1313311,1314827,1314828,1314829,1314830,1314831,1314832,1314833,1314834,1740636,1740639,1805426,1805427,1805428,1805429,1805430,1849588,1850256,1856032,1974110,1974111,1974112,1974113,1974114,1974115,1974116,1974117,1974118,1974119,1974120,1974121,1974122,1974123,1974124,1974125,1974126,1974127,1974128,1974129,1974130,2026807,1485545,2045305 --flat-output -v -o Zetaproteobacteria_NCBI_2021_10_19 -m Zetaproteobacteria_NCBI_2021_10_19_metadata.txt all
```

### Calculate ANI
### Calculate AAI

### Taxonomy
* `GTDB-Tk v1.3.0`: https://github.com/Ecogenomics/GTDBTk
```bash
gtdbtk classify_wf --genome_dir MAGs --out_dir MAGs_gtdbtk -x fa --cpus 20
```

### Genome quality and statistics
* `CheckM`: [https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/CheckM)
```bash
checkm lineage_wf allbins_only/ allbins_only/checkm -x fa -t 6 --pplacer_threads 6 --tmpdir /export/work_cgb/Petra/tmp -f BS4_checkM_113012021.tsv --tab_table
```

If you want to have a thorough bin statistics .tsv metadata file can be generated using the qa option in CheckM
```bash
checkm qa lineage.ms . --tmpdir /export/work_cgb/Petra/tmp -o 2 -f Phylogenomics_references_alltaxa_GB_checkm_stats.tsv --tab_table
```

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

Make separate trees for all separate proteins (to check if monophyletic/good marker genes)
```bash
anvi-gen-phylogenomic-tree -f Zetas_ribosomal_L15.fa -o Zetas_ribosomal_L15_tree.txt

or

FastTree protein_alignment > tree_file
```

Split multifasta markers from anvio into single marker files
- The file with all markers is divided based on the header before the "___" symbol into one file per marker:

```bash
seqkit split -i --id-regexp "^(\\S+)\___\s?" Zeta_ribosomal_markers20_proteins.fa
```

### Align individual sequences with mafft
* `MAFFT L-INS-i v7.397`
```bash
mkdir individual_mafft
for i in *.fa; do mafft-linsi $i | awk 'BEGIN{FS=":|[|]"}{if(/^>/){print ">"$2}else{print $0}}' > individual_mafft/${i%.fa}_mafft.fa; done
```

- Inspect the alignment manually
* `AliView v1.26`

### Trimming
* `trimAl v1.4.rev15`

Before using trimal, remove spaces in >fasta headers, so that trimal does not remove the taxonomic classification reported in the header
```bash
sed -i 's/ /_/' *mafft.fa
sed -i 's/:/_/' *mafft.fa

mkdir trimal
for i in *mafft.fa; do trimal -in $i -gt 0.5 -cons 60 |cut -f 1 -d ' ' > trimal/${i%.fa}_trimal.fa; done

-gt 0.5 -cons 60
#### Removes all positions in the alignment with gaps in 50% or more of the sequences, unless this leaves less than 60%. In such case, print the 60% best (with less gaps) positions.
```

### Concatenate
* `catfasta2.phyml v07.04.20`: (https://github.com/nylander/catfasta2phyml)
```bash
catfasta2phyml -v -c -f *mafft_trimal.fa > Zetaproteobacteria_concat_mafft_trimal.fa
```

### Build tree
* `IQ-TREE v2.0.3`: https://github.com/Cibiv/IQ-TREE
Choose the best substitution model (Best-fit model)
- "By default, substitution models are not included in these tests. If we want to test them we have to add them. Generally, it is recommended to include them in the test and the following selection would be quite comprehensive for testing models."

```bash
iqtree -s Zetaproteobacteria_concat_mafft_trimal.fa -m MFP -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60,LG+C10+R+F,LG+C20+R+F,LG+C30+R+F,LG+C40+R+F,LG+C50+R+F,LG+C60+R+F -v -nt 4
```

#### The standard non-parametric bootstraping with 1000 replicates
```bash
iqtree -s Zetaproteobacteria_concat_mafft_trimal.fa -m LG+F+R7 -b 1000 -nt 4 -pre Zetaproteobacteria_b1000_LG_F_R7
```

Visualise the tree by importing to iTOL and use templates for tree annotation: https://itol.embl.de/help.cgi#annot

## Phylogeny of NiFe uptake hydrogenase and cyc2

Dereplicate identical sequences
* `CD-HIT v4.8.1`
Clustering at 100% identity using CD-HIT

```bash
cd-hit -i Cyc2.fa -o Cyc2_drep100.fa -c 1.00
```
## Annotation: Metabolism and other genes

- Based on this annotation workflow: https://ndombrowski.github.io/Annotation_workflow/

Gene calling and functional annotation of MAGs was performed with an automated pipeline conducting separate searches against:
* `Prokka v1.14`
* NCBI COG (downloaded from NCBI webserver in February 2021)
* arCOG (version from 2018)
* KEGG (downloaded in February 2021)
* Pfam (release 33.0)
* TIGRFAM (release 15.0)
* CAZy (dbCAN v9)
* Transporter Classification Database (downloaded from TCDB webserver in February 2021)
* HydDB (downloaded from HydDB webserver in February 2021)
* NCBI_nr (downloaded from NCBI webserver in February 2021)

### Iron oxidation

#### FeGenie: iron genes and metabolism
* `FeGenie`
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

### Predicting gene expression using codon bias
* `coRdon v1.8.0`: R package https://github.com/BioinfoHR/coRdon
- GUI: `INCA`

## Other

### Zetaproteobacteria taxonomy based on OTUs
* `ZetaHunter v1.0.11`: https://github.com/mooreryan/ZetaHunter/

1. Extract 16S rRNA genes from MAGs
* `barrnap`: https://github.com/tseemann/barrnap
* `anvi’o v7.1`: https://anvio.org

```bash
for i in *.fa; do barrnap --kingdom bac --lencutoff 0.2 --reject 0.3 --evalue 1e-05 --threads 20 --outseq barrnap/${i%.fa}_barrnap.fa $i; done
```

or

```bash
anvi-get-sequences-for-hmm-hits --external-genomes external_selected_genomes.txt --hmm-source Ribosomal_RNA_16S --list-available-gene-names
```

2. Check for chimeras
3. Make a SINA alignment with them: https://www.arb-silva.de/aligner/
3. Run ZetaHunter (with a --no-check-chimeras tag, otherwise it doesn't want to run)

    ```
run_zeta_hunter --inaln Zetas_SINAalign.fasta --outdir Faavne_biofilms_Zetahunter --no-check-chimeras
    ```
    
### Predicting optimal growth temperatures (OGT)
* https://github.com/DavidBSauer/OGT_prediction

1. Place each genome in a folder in the prediction directory called "genomes/XXX/" where XXX is the name of the species.  Genomes should be ziped in .gz format.
2. Create a tab separated file of the genomes and species pairs. Provide this file in place of genomes_retrieved.txt. 
3. Create a file for the taxonomic classification of each species. 
4. Run the script: Run the prediction script to predict the OGT of each species. If you want to account for absence of 16S rRNA gene and genome incompleteness (usual with MAGs), use specified regression models modified to exclude 16S rRNA gene and genome size: https://github.com/DavidBSauer/OGT_prediction/tree/master/data/calculations/prediction/regression_models

    ```
    python3 ~/OGT_prediction/OGT_prediction/prediction/prediction_pipeline.py regression_models_Bacteria_uncultivated/ Faavne_MAGs.txt Faavne_MAGs_taxonomic3.txt
    ```

5. Results in newly_predicted_OGTs.txt

### Viruses
* `VIBRANT v.1.2.1`: https://github.com/AnantharamanLab/VIBRANT

    ```
VIBRANT_run.py -i unbinned_contigs_refined.fa
    ```
    
VIBRANT uses:
* KEGG (March release): https://www.genome.jp/kegg/ (FTP: ftp://ftp.genome.jp/pub/db/kofam/archives/2019-03-20/)
* Pfam (v32): https://pfam.xfam.org (FTP: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/)
* VOG (release 94): http://vogdb.org/ (FTP: http://fileshare.csb.univie.ac.at/vog/vog94/)

* `CheckV v.0.8.1`: https://bitbucket.org/berkeleylab/checkv/src/master/

    ```
checkv end_to_end input_file.fna output_directory -t 4
    ```
    
* `VirMatcher`: https://bitbucket.org/MAVERICLab/virmatcher/src/master/

- Tools that make VirMatcher possible: Minced, tRNAscan-SE, WIsH, BLAST
- NOTE! All fasta files need to end in .fasta!
- "--virus-fp": The putative/suspected/known viruses in a single FASTA-formatted file.
- Taxonomy files: A tab-delimited archaeal genome taxonomy file. Should be in the format of "host name \<TAB> taxon". The taxon can be any level.

    ```
VirMatcher --virus-fp Faavne_Nanopore_viruses.fasta --bacteria-host-dir Faavne_potential_hosts --bacteria-taxonomy Faavne_potential_hosts_taxonomy.txt --threads 8 -o Faavne_Nanopore_VirMatcher --python-aggregator
    ```


