# Microbial iron mats at Fåvne hydrothermal vent

Tools and commands used for analyses of microbial iron mats at Fåvne hydrothermal vent field. Article in preparation.

## Genome-resolved metagenomics

## Phylogenomics

## Phylogeny of NiFe uptake hydrogenase and cyc2

## Metabolism and other genes

## Iron oxidation

## Other

### Predicting optimal growth temperatures (OGT)
https://github.com/DavidBSauer/OGT_prediction

1. Place each genome in a folder in the prediction directory called "genomes/XXX/" where XXX is the name of the species.  Genomes should be ziped in .gz format.
2. Create a tab separated file of the genomes and species pairs. Provide this file in place of genomes_retrieved.txt. 
3. Create a file for the taxonomic classification of each species. 
4. Run the script: Run the prediction script to predict the OGT of each species. If you want to account for absence of 16S rRNA gene and genome incompleteness (usual with MAGs), use specified regression models modified to exclude 16S rRNA gene and genome size: https://github.com/DavidBSauer/OGT_prediction/tree/master/data/calculations/prediction/regression_models

    ```
    python3 ~/OGT_prediction/OGT_prediction/prediction/prediction_pipeline.py regression_models_Bacteria_uncultivated/     Faavne_MAGs.txt Faavne_MAGs_taxonomic3.txt
    ```

5. Results in newly_predicted_OGTs.txt

