# Frame_Editor_sgRNA_selection

## Selection process

We compiled our test set of several sgRNAs (Table 1) from 64,751 SpCas9 guide set to reasonably assess the feasibility of Frame Proof in various sequence contexts. The set consisted of SpCas9 guide groups (A, B). Our paper showed the result of the partial sgRNA set because it is difficult to perform PCR amplification in some sgRNAs.

Group A. We selected 14 sgRNAs used in ref1 (Table_S1) that identified genes essentiality for cell viability in cancer and pluripotent stem cells using the genome-scale CRISPR-Cas9 knockout (GeCKO) library. The sgRNAs targeting a gene ranked within the top-1,000 essential genes in any cell types (Table 2) were discarded from the list of candidates to avoid bias resulting from the lethality of targeting gene mutations. Moreover, we excluded sgRNA targeting a gene known to have the genetic variation driving cancer according to the COSMIC Cancer Gene Census(ref2) (CMC; release v92, 27th August 2020).
 We applied the machine learning model for the indel prediction(ref3) (inDelphi, https://github.com/maxwshen/inDelphi-model, scikit-learn v0.20.0 models, for HEK293) to select the only guides that would show definite distributions on indels. The sgRNA targets whose “Precision” value was higher than 0.35 were selected as the final candidates. Lastly, we chose two sgRNAs from seven groups (a-g) that were predicted to show a specific indel distribution in inDelphi.

The targeting genes are here:
	A-insertion dominant (ZNF498 and C1orf51).  The predicted “1-bp ins frequency” was higher than the third quartile of the predicted “1-bp ins frequency” in the final candidates (Q31-bp ins frequency). The most frequent one-base insertion was estimated to be A-insertion.
	T-insertion dominant (ODZ3 and FAM159B).  The predicted “1-bp ins frequency” was higher than the Q31-bp ins frequency. The most frequent one-base insertion was estimated to be T-insertion.
	C-insertion dominant (A2LD1 and TMEM57).  The predicted “1-bp ins frequency” was higher than the Q31-bp ins frequency. The most frequent one-base insertion was estimated to be C-insertion.
	MMEJ-assisted frame-in deletion (C16orf93 and C1orf204). The predicted “1-bp ins frequency” was not higher than the Q31-bp ins frequency. The predicted “MH del frequency” is higher than the “MHless del frequency”. The total predicted frequency of indels resulting in +0 frame change (Freq+0) was the most frequent in the Freq+0, Freq+1, and Freq+2.
	MMEJ-assisted +1 frameshift deletion (C18orf1 and HBXIP). The predicted “1-bp ins frequency” was not higher than the Q31-bp ins frequency. The predicted “MH del frequency” is higher than the “MHless del frequency”. The Freq+1 was the most frequent in the Freq+0, Freq+1, and Freq+2.
	MMEJ-assisted +2 frameshift deletion (RP4-697K14.7 and BRP44). The predicted “1-bp ins frequency” was not higher than the Q31-bp ins frequency. The predicted “MH del frequency” is higher than the “MHless del frequency”. The Freq+2 was the most frequent in the Freq+0, Freq+1, and Freq+2.
	NHEJ-assisted +1 frameshift deletion (FAM26D and ATP5G1). The predicted “1-bp ins frequency” was not higher than the Q31-bp ins frequency. The predicted “MH del frequency” is not higher than the “MHless del frequency”. The Freq+1 was the most frequent in the Freq+0, Freq+1, and Freq+2.

We also evaluated the predicted indel distribution using the other human cell models (U2OS, HCT116, and K562). To calculate the divergence of the indel distributions, we used the Kullback–Leibler (KL) divergence. KL divergence is calculated as

D_KL(P||Q)=∑_i(P_i * log⁡(P_i/Q_i ))

where i indexes the different predicted indels, and Pi and Qi are their normalized rate of total predicted mutations in the compared cells. We estimate the weighted average of KL divergences between HEK293 and each cell line.
 The target sequences for ZNF498, ODZ3, A2LD1, C16orf93, C18orf1, RP4-697K14.7, and FAM26D have a high value of the average divergence, which means that the indel pattern could be affected by the cell-specific factors. The other target sequences (C1orf51, FAM159B, TMEM57, C1orf204, HBXIP, BRP44, and ATP5G1) have the low value. These indel distributions are predicted to be stable in the cell environment.

Group B. We chose six sgRNAs from the final candidate's list of Group A. These sgRNAs were predicted to present a low frequency of 1 bp insertion (<7%) and a high frequency of MMEJ-assisted deletion (>87%).

The Python script used for the selection process was the following section.

### References

1. Shalem, Ophir, et al. Genome-scale CRISPR-Cas9 knockout screening in human cells. Science 343, 84-87 (2014).
2. Sondka, Zbyslaw, et al. The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nature Reviews Cancer 18 696-705 (2018).
3. Shen, Max W., et al. Predictable and precise template-free CRISPR editing of pathogenic variants. Nature 563, 646-651 (2018).

---

## Script used for selection process

### Get target sequence (Upstream 13bp + protospacer + PAM + Downstream 24bp) using Blat

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

### Get a guide RNA selection

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