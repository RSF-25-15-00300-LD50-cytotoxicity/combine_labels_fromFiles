library(tidyverse)
## Input
struct_path <- "..."
# List files from AAL
raw_fn <- list.files(struct_path, full.names = TRUE, pattern = "*.SDF")

# Calculate descriptors to be able to merge the data according to the PASS internal representation
for (i in seq(1:length(raw_fn))) {
	cmd_line <- str_glue("{struct_path}/makemna.exe {raw_fn[i]}")
	system(cmd_line)
}
# List SDFs with descriptors
mna_fn <- list.files(struct_path, full.names = TRUE, pattern = "*_SD.SDF")
# Read the files with descriptors and parse them to get access to structures, descriptors, codes and labels
data_raw <- tibble(mol = "example", cid = "example", cid_all = "example", code = "example", label = "example", mna = "example")
for (i in seq(1:length(mna_fn))) {
	data_i <- read_file(mna_fn[i]) |>
				str_replace_all(">  <hazard label>", ">  <GHS label>") |>
				str_trim()
	if (!str_detect(data_i, "\r\n>  <GHS code>\r\n")) {
		data_i <- str_replace_all(data_i, "\r\n>  <GHS label>\r\n", "\r\n>  <GHS code>\r\nunknown\r\n\r\n>  <GHS label>\r\n")
	}
	data_i	<-	as_tibble(data_i)		|>
				separate_longer_delim(value, delim = "$$$$") |>
				mutate(value = str_trim(value)) |>
				filter(value != "") |>
				separate_wider_delim(value, delim = "\r\n>  <CID>\r\n", names=c("mol", "data"), too_few = "align_start") |>
				separate_wider_delim(data, delim = "\r\n>  <CID_all>\r\n", names=c("cid", "data"), too_few = "align_start") |>
				separate_wider_delim(data, delim = "\r\n>  <GHS code>\r\n", names=c("cid_all", "data"), too_few = "align_start") |>
				separate_wider_delim(data, delim = "\r\n>  <GHS label>\r\n", names=c("code", "data"), too_few = "align_start") |>
				separate_wider_delim(data, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names=c("label", "mna"), too_few = "align_start") |>
				mutate_all(str_trim)
	# Slow solution, but OK considering the moderate data volume
	data_raw <- bind_rows(data_raw, data_i)
}

# Process the tibble
data <- data_raw |> filter(mol != "example") |>
					filter(!is.na(mna)) |>
					group_by(mna) |>
					mutate(cid_all = str_c(cid_all, collapse = ", "),
						   code = str_c(code, collapse = "\r\n"),
						   label = str_c(label, collapse = "\r\n")) |>
					slice_head(n = 1) |>
					ungroup() |>
					rowwise() |>
					mutate(cid_all = cid_all |> str_split(", ", simplify = TRUE) |> as.integer() |> sort() |> unique() |> str_c(collapse = ", "),
						   code    = code |> str_split("\r\n",  simplify = TRUE) |> sort() |> unique() |> str_c(collapse = "\r\n") |> str_replace("unknown", "") |> str_trim(),
						   label   = label |> str_split("\r\n", simplify = TRUE) |> sort() |> unique() |> str_c(collapse = "\r\n") |> str_trim()) |>
					ungroup()
# CID 222 -> no MNA

# Export the results to SDF
sdf_prep <- data |> mutate(id_rec = "\r\n>  <CID>\r\n", allid_rec = "\r\n\r\n>  <CID_all>\r\n",
							code_rec = "\r\n\r\n>  <GHS code>\r\n", label_rec = "\r\n\r\n>  <GHS label>\r\n", end_rec = "\r\n\r\n$$$$") |>
					select(mol, id_rec, cid, allid_rec, cid_all, code_rec, code, label_rec, label, end_rec) |>
					unite("record", mol:end_rec, sep = "")
write_lines(str_c("", sdf_prep[[1]]), str_glue("{struct_path}/gathered__hazard_labels_PC_CIDP_NLM_ghs.SDF"))
# Export the results to CSV / TAB
write_tsv(data, str_glue("{struct_path}/gathered__hazard_labels_PC_CIDP_NLM_ghs.tsv"), quote = "all")