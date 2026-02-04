---
title: "ID Mapping - Worked Tutorial"
output: html_document
---

# Chapter 4: biomaRt step-by-step

``` r
library(biomaRt) # Load the biomaRt package in
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:biomaRt':
## 
##     select
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
listMarts() # Available biomaRts
```

```
##                biomart                version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 115
## 2   ENSEMBL_MART_MOUSE      Mouse strains 115
## 3     ENSEMBL_MART_SNP  Ensembl Variation 115
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 115
```

``` r
ensembl <- useMart("ensembl") # Not pinning a specific version

datasets <- listDatasets(ensembl) # Datasets part of the ensembl biomaRt

human <- datasets[grep(datasets$dataset, pattern = "hsapiens"), ]
human # List all datasets where the sample was collected from humans
```

```
##                  dataset              description    version
## 80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14
```

``` r
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl) # Use the dataset name that was listed in the previous dataset

all_filters <- listFilters(ensembl) # Filters are input ID type

all_filters[grep(all_filters$name, pattern = "ensembl_gene"), ] # Find the filters that are related to ensembl gene IDs
```

```
##                       name
## 54         ensembl_gene_id
## 55 ensembl_gene_id_version
##                                                 description
## 54                 Gene stable ID(s) [e.g. ENSG00000000003]
## 55 Gene stable ID(s) with version [e.g. ENSG00000000003.17]
```

``` r
all_attr <- listAttributes(ensembl) # Attributes are the output columns that we want to get

searchAttributes(ensembl, "hgnc") # Find the available outputs that are related to HGNC
```

```
##               name        description         page
## 63     hgnc_symbol        HGNC symbol feature_page
## 64         hgnc_id            HGNC ID feature_page
## 95 hgnc_trans_name Transcript name ID feature_page
```

``` r
strip_ensembl_version <- function(x) sub("\\..*$", "", x)
ids <- c("ENSG00000141510.17", "ENSG00000157764.2")
ids_clean <- strip_ensembl_version(ids)
ids_clean # Clean the IDs (remove version suffixes) so they can be queried
```

```
## [1] "ENSG00000141510" "ENSG00000157764"
```

``` r
cache_file <- "id_conversion.rds"

if(file.exists(cache_file)) {
  map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ids_clean,
    mart = ensembl
  )
  saveRDS(map, cache_file)
} # Cache file

df <- tibble(ensembl_gene_id = ids_clean, values = c(1, 2))
df_annot <- left_join(df, map, by = "ensembl_gene_id")
df_annot # Merge all data together! Now we have ensembl_gene_id and the hgnc_symbol.
```

```
## # A tibble: 2 × 3
##   ensembl_gene_id values hgnc_symbol
##   <chr>            <dbl> <chr>      
## 1 ENSG00000141510      1 TP53       
## 2 ENSG00000157764      2 BRAF
```

# Chapter 6.1 (Exercise 1)

``` r
df <- readr::read_table("data/GSE119732/GSE119732_count_table_RNA_seq.txt.gz")
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   gene_id = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

``` r
filtered_df <- head(df, 20) %>%
  dplyr::select(gene_id) # Extract the first 20 IDs

version_counts <- sum(grepl("\\.", filtered_df$gene_id))
version_counts # Number of IDs that have a version code
```

```
## [1] 20
```

``` r
filtered_df <- filtered_df %>%
  mutate(stripped = sub("\\..*$", "", gene_id)) # Create column with stripped IDs

ids <- filtered_df$stripped[1:20] # Put stripped IDs into a vector

# Load in ensembl again (from above)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Map to HGNC symbols
map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ids,
  mart = ensembl
)
map
```

```
##    ensembl_gene_id hgnc_symbol
## 1  ENSG00000186092       OR4F5
## 2  ENSG00000222623  RNU6-1100P
## 3  ENSG00000223972     DDX11L1
## 4  ENSG00000227232      WASH7P
## 5  ENSG00000233750      CICP27
## 6  ENSG00000237613     FAM138A
## 7  ENSG00000239906            
## 8  ENSG00000239945            
## 9  ENSG00000240361     OR4G11P
## 10 ENSG00000241599            
## 11 ENSG00000241860            
## 12 ENSG00000243485 MIR1302-2HG
## 13 ENSG00000268020      OR4G4P
## 14 ENSG00000268903            
## 15 ENSG00000269981            
## 16 ENSG00000273874   MIR6859-2
## 17 ENSG00000278267   MIR6859-1
## 18 ENSG00000279457      WASH9P
## 19 ENSG00000279928    DDX11L17
```

# Chapter 6.2 (Exercise 2)

``` r
# Get files from GEO
source("./fetch_geo_supp.R")
fetch_geo_supp(gse = "GSE122380")
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Using locally cached version of supplementary file(s) GSE122380 found here:
## data/GSE122380/GSE122380_Supplementary_Data_Table_S1.xlsx
```

```
## Using locally cached version of supplementary file(s) GSE122380 found here:
## data/GSE122380/GSE122380_raw_counts.txt.gz
```

``` r
# Repeating code from above
df <- readr::read_table("data/GSE122380/GSE122380_raw_counts.txt.gz")
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   Gene_id = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

``` r
filtered_df <- head(df, 20) %>%
  dplyr::select(Gene_id) # Extract the first 20 IDs

version_counts <- sum(grepl("\\.", filtered_df$Gene_id))
version_counts # Number of IDs that have a version code
```

```
## [1] 0
```

``` r
filtered_df <- filtered_df %>%
  mutate(stripped = sub("\\..*$", "", Gene_id)) # Create column with stripped IDs

ids <- filtered_df$stripped[1:20] # Put stripped IDs into a vector

# Load in ensembl again (from above)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Map to HGNC symbols
map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ids,
  mart = ensembl
)
map
```

```
##    ensembl_gene_id hgnc_symbol
## 1  ENSG00000000419        DPM1
## 2  ENSG00000000457       SCYL3
## 3  ENSG00000000460       FIRRM
## 4  ENSG00000000938         FGR
## 5  ENSG00000000971         CFH
## 6  ENSG00000001036       FUCA2
## 7  ENSG00000001084        GCLC
## 8  ENSG00000001167        NFYA
## 9  ENSG00000001460       STPG1
## 10 ENSG00000001461      NIPAL3
## 11 ENSG00000001561       ENPP4
## 12 ENSG00000001617      SEMA3F
## 13 ENSG00000001626        CFTR
## 14 ENSG00000001629      ANKIB1
## 15 ENSG00000001630     CYP51A1
## 16 ENSG00000001631       KRIT1
## 17 ENSG00000002016       RAD52
## 18 ENSG00000002330         BAD
## 19 ENSG00000002549        LAP3
## 20 ENSG00000002587      HS3ST1
```
