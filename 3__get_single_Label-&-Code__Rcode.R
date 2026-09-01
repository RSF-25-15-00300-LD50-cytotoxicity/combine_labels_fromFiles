library(tidyverse)
## Input
struct_path <- "..."
# Read the data
# Table
data_raw <- read_tsv(str_glue("{struct_path}/SMILES_gathered__hazard_labels_PC_CIDP_NLM_ghs.tsv"))

## Process, divide into the distinct records according to the associated codes (one code per record) maintaining the structure to code relations
#  The same thing for the text labels
#  Text labels and codes are processed separately, since one of the files processed using the first script did not have codes,
#  Thus, storing codes and text labels is not an easy option
data_codes  <- data_raw |> separate_longer_delim(code, delim = "\r\n") |>
						select(smiles, cid, cid_all, code) |>
						distinct()
data_labels <- data_raw |> separate_longer_delim(label, delim = "\r\n") |>
						select(smiles, cid, cid_all, label) |>
						distinct()


## Export the results
# to CSV / TAB
write_tsv(data_codes, str_glue("{struct_path}/SMILES-&-codes_divided__hazard_labels_PC_CIDP_NLM_ghs.tsv"), quote = "all")
write_tsv(data_labels, str_glue("{struct_path}/SMILES-&-labels_divided__hazard_labels_PC_CIDP_NLM_ghs.tsv"), quote = "all")