library(tidyverse)
## Input & Process
path <- "..."
data <- read_file(str_glue("{path}/SDF_gathered__hazard_labels_PC_CIDP_NLM_ghs_Training_for_Supplementary.sdf")) |> str_replace_all("\r\n", "\n") |>
							  str_replace_all("\n", "\r\n") |>
							  as_tibble() |>
							  separate_longer_delim(value, delim = "$$$$") |>
							  mutate(value = str_trim(value)) |>
							  filter(value != "") |>
							  separate_wider_delim(value, delim = "\r\n>  <Molecule Name>\r\n", names=c("mol", "data"), too_few = "align_start") |>
							  separate_wider_delim(data, delim = "\r\n>  <CID>\r\n", names=c("name", "data"), too_few = "align_start") |>
							  separate_wider_delim(data, delim = "\r\n>  <CID_all>\r\n", names=c("cid", "data"), too_few = "align_start") |>
							  separate_wider_delim(data, delim = "\r\n>  <GHS code>\r\n", names=c("cid_all", "data"), too_few = "align_start") |>
							  separate_wider_delim(data, delim = "\r\n>  <GHS label>\r\n", names=c("code", "data"), too_few = "align_start") |>
							  separate_wider_delim(data, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names=c("label", "mna"), too_few = "align_start") |>
							  mutate_all(str_trim) |>
							  select(mol, cid_all, cid, name, code, label) |>
							  group_by(code) |>
							  slice_sample(prop = .25) |>
							  ungroup()
## Export the results to SDF
sdf_prep <- data |> mutate(id_rec = "\r\n\r\n>  <CID>\r\n", allid_rec = "\r\n>  <CID_all>\r\n",
							code_rec = "\r\n\r\n>  <GHS code>\r\n", label_rec = "\r\n\r\n>  <GHS label>\r\n", end_rec = "\r\n\r\n$$$$") |>
					select(mol, allid_rec, cid_all, id_rec, cid, code_rec, code, label_rec, label, end_rec) |>
					unite("record", mol:end_rec, sep = "")
write_lines(str_c("", sdf_prep[[1]]), str_glue("{path}/sample__SDF_gathered__hazard_labels_PC_CIDP_NLM_ghs_Training_for_Supplementary.sdf"))