# Sc-candida paper

This is the code and data associated with our study published in eLife in 2023.

Please cite: 
Dumeaux V, Massahi S, Bettauer V, Khurdia S, Costa ACBP, Omran RP, Simpson S, Xie JL, Whiteway M, Berman J, Hallett MT. Candida albicans exhibits adaptive cytoprotective responses to anti-fungal compounds. Elife 2023;12:e81406 [https://dx.doi.org/10.7554/eLife.81406](https://dx.doi.org/10.7554/eLife.81406)

We use OSF as a [special remote](https://git-annex.branchable.com/special_remotes/) to store data in the annex of this repository : we published the repository to GitHub and have the data published to the [OSF](https://osf.io/5tpk3/) (via a publication dependency). 

## To clone repository and get the data

You need to install datalad 

```bash
conda create -n datalad-osf -c conda-forge datalad python=3.6

conda activate datalad-osf
```

You then can clone the repository and get the data from osf

```bash
git clone https://github.com/vdumeaux/sc-candida_paper.git
cd sc-candida_paper
datalad get .
```


## To run the code 

### Software

You need to have the following software installed and change the paths in `src/init.R` : 

* sra-toolkit
* Salmon 1.6.0
* alevin-fry 0.4.3 (bioconda)
* STAR 2.7.9a
* dropEst 0.8.6
* TrimGalore 0.6.6
* cutadapt 1.18 (bioconda)


### scvi environment 

To run scvi notebook, create a new scvi-env as follow:

```
conda create --name scvi-env scvi-tools rpy2=3.4.2 scanpy anndata2ri anndata bioconductor-singlecellexperiment leidenalg phate gprofiler-official scprep graphtools meld -c bioconda -c conda-forge

pip install --user magic-impute
```

OR

```
conda env create -f scvi-env.yml
```

