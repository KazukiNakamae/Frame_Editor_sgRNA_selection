# Frame_Editor_sgRNA_selection

We selected the gRNA dataset using the below script. 

We referred to the following reports to develop it.
1. Shalem, Ophir, et al. Genome-scale CRISPR-Cas9 knockout screening in human cells. Science 343, 84-87 (2014).
2. Sondka, Zbyslaw, et al. The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nature Reviews Cancer 18 696-705 (2018).
3. Shen, Max W., et al. Predictable and precise template-free CRISPR editing of pathogenic variants. Nature 563, 646-651 (2018).

---

## Script

### Get target sequences (Upstream 13bp + protospacer + PAM + Downstream 24bp) using Blat

```bash
# create enviroment
conda create --name ExtendSeq_env;
conda activate ExtendSeq_env;
conda install -c bioconda blat bedtools ucsc-psltobed;
conda install -c anaconda wget python==3.8.5;
pip install pandas biopython;

# download reference genome (hg38)
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz --directory-prefix ./hg38;
cd ./hg38;
gunzip -c hg38.analysisSet.chroms.tar.gz | tar xopf -;
cd ..;
for i in `seq 1 22` X Y;do cat ./hg38/hg38.analysisSet.chroms/chr$i.fa >> ./hg38/hg38.fa;done;

# search target sequence using Blat
python ExtendSeq.py -i sgRNALibrary.csv -o sgRNALibrary_extended.csv -s "sgRNA sequence" -u "Number" -l 13 -r 27;

conda deactivate;
```

### Get the sgRNA set

```bash
# create enviroment
conda create --name ClassifyMutPatt_env;
conda activate ClassifyMutPatt_env;
conda install -c anaconda python scikit-learn==0.20.0;
pip install pandas biopython numexpr;
conda install -c pytorch pytorch torchvision;
git clone https://github.com/maxwshen/inDelphi-model.git;
# rename inDelphi-model -> inDelphi-HEK293-model
git clone https://github.com/maxwshen/inDelphi-model.git;
# rename inDelphi-model -> inDelphi-HCT116-model
git clone https://github.com/maxwshen/inDelphi-model.git;
# rename inDelphi-model -> inDelphi-K562-model
git clone https://github.com/maxwshen/inDelphi-model.git;
# rename inDelphi-model -> inDelphi-U2OS-model

# download data of COSMIC Cancer Gene Census (CMC; release v92, 27th August 2020) from https://cancer.sanger.ac.uk/census

# get guide RNA selection to target non-essential and non-cancer-causing genes and estimate the indel pattern using inDelphi 
python ClassifyMutPatt.py -i sgRNALibrary_extended.csv -o ClassifiedGenes -s "extended_sequence" -g "Gene name" -m "Number of 1bp mismatches" -m "Number of 2bp mismatches" -e A375_HUES62_essential_genes.txt -c "CMC_export.v92/cmc_export.tsv";
conda deactivate;
# NOTE: This program does not export data dominant in G insertion, NHEJ-mediated +2 frameshift, and NHEJ-mediated frame-in deletion because there are no applicants in this search

# check specificity using GGGenome(https://gggenome.dbcls.jp/ja/)

```