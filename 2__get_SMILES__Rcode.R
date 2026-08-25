library(tidyverse)
library(rcdk)
library(itertools)
## Input
struct_path <- "..."
# Read the data
# Table
data_raw <- read_tsv(str_glue("{struct_path}/gathered__hazard_labels_PC_CIDP_NLM_ghs.tsv"))

## Process
# get SMILES using CDK for R, see the https://cran.r-project.org/web/packages/rcdk/vignettes/rcdk.html
i <- 1
iter <- iload.molecules(str_glue("{struct_path}/gathered__hazard_labels_PC_CIDP_NLM_ghs.SDF"), type="sdf")
while(hasNext(iter)) {
 mol <- nextElem(iter)
 data_raw[i, 7] <- get.property(mol, "CID")
 data_raw[i, 8] <- get.smiles(mol, flavor = smiles.flavors("Canonical"))
 i <- i + 1
}
# check IDs
n_prblms <- data_raw |> filter(cid != `...7`) |> nrow()
# should be OK

## Process
data <- data_raw |> rename(smiles = `...8`) |>
					select(smiles, cid, cid_all, code, label)

## Export the results
# to CSV / TAB
write_tsv(data, str_glue("{struct_path}/SMILES_gathered__hazard_labels_PC_CIDP_NLM_ghs.tsv"), quote = "all")