# AffinitySPI: IC50 cutoffs
Process public MHC-I immunopeptidome data to infer the distribution of binding affinites.

## Data sources
IEDB - download and unzip the folder
`https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip`

SysteMHC Atlas
`https://systemhcatlas.org/`

Human HLA Atlas
```
https://github.com/CaronLab/MHCIatlas/blob/master/inst/extdata/HumanHLAatlas_reprocessed1.csv
https://github.com/CaronLab/MHCIatlas/blob/master/inst/extdata/HumanHLAatlas_reprocessed2.csv
https://github.com/CaronLab/MHCIatlas/blob/master/inst/extdata/HumanHLAatlas_reprocessed3.csv
https://github.com/CaronLab/MHCIatlas/blob/master/inst/extdata/HumanHLAatlas_reprocessed4.csv
https://github.com/CaronLab/MHCIatlas/blob/master/inst/extdata/HumanHLAatlas_reprocessed5.csv
```

## Options
`netMHCpan` - local path to netMHCpan-4.1
`alleles` - if `TRUE`, use all aleleles, otherwise - a character vector with allele names
`Use_strictly_monoallelic` - if `TRUE`, use only IEDB monoallelic data, otherwise all data will be used
`min_peptides_per_allele.length` - minimal number of peptide sequences per allele and peptide length to estimate the cutoffs
`n_cores` - number of CPUs to use

