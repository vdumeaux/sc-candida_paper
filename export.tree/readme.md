



## Software requirements

### Software

You need to have the following software installed and change the paths in `src/init.R` : 

sra-toolkit

Salmon 1.6.0
alevin-fry 0.4.3 (bioconda)
STAR 2.7.9a
dropEst 0.8.6
TrimGalore 0.6.6
cutadapt 1.18 (bioconda)


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

