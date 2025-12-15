
library(dplyr)
library(tidyverse)
load(file.path("compileduids.RData"))

moore_list <- read.csv("fromSteveMoore_metabolite_list.csv")


#########################################
#     MetLinkR Harm of COMETS + Moore   #
#########################################
# 
# install.packages("devtools")
# library(devtools)
# install_github("ncats/RAMP-DB")
# 
# library(RaMP)
# rampDB <- RaMP()
# rampDB <- RaMP(version = "2.5.4")
# install.packages("remotes")
# remotes::install_github("ncats/MetLinkR")

write.csv(mastermetid, "COMETS_metabolite_df.csv")
COMETS_file_name <- "COMETS_metabolite_df.csv"
moore_file_name <- "fromSteveMoore_metabolite_list.csv"
metlinkr_input_file <- data.frame("FileNames" = c(COMETS_file_name, moore_file_name),
                                  "ShortFileName" = c("COMETS", "Moore"),
                                  "HMDB" = c("hmdb_id", "input_hmdb_id"),
                                  "Metabolite_Name" = c("biochemical", "input_metabolite_name"),
                                  "PubChem_CID" = c(NA, "input_pubchem"),
                                  "KEGG" = c(NA, "input_kegg"),
                                  "LIPIDMAPS" = c(NA, NA),
                                  "chebi" = c(NA, "input_chebi"))

# library(devtools)
# load_all("/Users/heysiri/MetLinkR")
# library(dplyr)
# library(readr)
# write.csv(metlinkr_input_file, "metlinkr_input_file_for_merging.csv")
# metLinkR::harmonizeInputSheets("metlinkr_input_file_for_merging.csv")

moore_comets <- readxl::read_excel("mapping_library_COMETS_Moore.xlsx")
moore_comets_overlapping <- moore_comets %>% filter(`Input name (COMETS)` != "-" & `Input name (Moore)` != "-") # 1688 rows


###################################################
#     Merged Dataset Creation (COMETS + MOORE)    #
###################################################

mastermetid_clean <- mastermetid %>%
  select(-uidsource, -main_class, -chemical_id, -comp_id)

moore_list_clean <- moore_list %>%
  select(-dataset_id, -study_name, -platform, -merge_variable,
         -input_chemical_id, -input_comp_id, -input_cas) %>%
  mutate(
    biochemical_lower = tolower(input_metabolite_name),
    input_hmdb_id = trimws(input_hmdb_id),
    metabolite_id = trimws(metabolite_id)
  )

mastermetid_clean <- mastermetid_clean %>%
  mutate(biochemical_lower = tolower(biochemical))

collapse_vals <- function(x) {
  x <- x[!is.na(x) & x != ""]        # remove NAs and blanks
  if (length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = " | ")
}

# --- Collapse moore list by identifiers ---
moore_metid <- moore_list_clean %>%
  filter(!is.na(metabolite_id) & metabolite_id != "") %>%  # skip blank IDs
  group_by(metabolite_id) %>%
  summarise(across(everything(), collapse_vals), .groups = "drop")

moore_hmdb <- moore_list_clean %>%
  filter(!is.na(input_hmdb_id) & input_hmdb_id != "") %>%  # skip blank IDs
  group_by(input_hmdb_id) %>%
  summarise(across(everything(), collapse_vals), .groups = "drop")

moore_biochem <- moore_list_clean %>%
  filter(!is.na(biochemical_lower) & biochemical_lower != "") %>%
  group_by(biochemical_lower) %>%
  summarise(across(everything(), collapse_vals), .groups = "drop")

mastermetid_clean <- mastermetid_clean %>% mutate(source_comets = TRUE)
moore_metid <- moore_metid %>% mutate(source_moore_metid = TRUE)
moore_hmdb <- moore_hmdb %>% mutate(source_moore_hmdb = TRUE)
moore_biochem <- moore_biochem %>% mutate(source_moore_biochem = TRUE)

# --- Join sequentially ---
merged <- mastermetid_clean %>%
  full_join(moore_metid, by = c("uid_01" = "metabolite_id")) %>%
  full_join(moore_hmdb, by = c("hmdb_id" = "input_hmdb_id")) %>%
  full_join(moore_biochem, by = c("biochemical_lower" = "biochemical_lower")) %>%
  mutate(both_moore_comets = ifelse(
    !is.na(uid_01) & (!is.na(input_metabolite_name.x) | !is.na(input_metabolite_name.y) | !is.na(input_metabolite_name)), 
    TRUE, FALSE
  ))

merged_clean <- merged %>%
  mutate(
    biochemical_lower = coalesce(biochemical_lower.x, biochemical_lower.y, biochemical_lower),
    input_metabolite_name = coalesce(input_metabolite_name.x, input_metabolite_name.y, input_metabolite_name),
    input_metabid = coalesce(input_metabid.x, input_metabid.y, input_metabid),
    input_hmdb_id = coalesce(input_hmdb_id.x, input_hmdb_id.y),
    input_pubchem = coalesce(input_pubchem.x, input_pubchem.y, input_pubchem),
    input_chemspider = coalesce(input_chemspider.x, input_chemspider.y, input_chemspider),
    input_kegg = coalesce(input_kegg.x, input_kegg.y, input_kegg),
    input_inchikey = coalesce(input_inchikey.x, input_inchikey.y, input_inchikey),
    input_chebi = coalesce(input_chebi.x, input_chebi.y, input_chebi),
    metabolite_id = coalesce(metabolite_id.x, metabolite_id.y)
  ) %>%
  select(
    -ends_with(".x"),
    -ends_with(".y")
  )


# Step 0: Ensure all dfs have the same columns
force_character <- function(df) {
  df[] <- lapply(df, as.character)
  return(df)
}

merged_collapsed <- force_character(merged_clean)
mastermetid_clean <- force_character(mastermetid_clean)
moore_list_clean <- force_character(moore_list_clean)

all_cols <- union(
  names(merged_collapsed),
  union(names(mastermetid_clean), names(moore_list_clean))
)

add_missing_cols <- function(df, all_cols) {
  missing_cols <- setdiff(all_cols, names(df))
  if (length(missing_cols) > 0) {
    df <- df %>% add_column(!!!setNames(rep(list(NA), length(missing_cols)), missing_cols))
  }
  return(df)
}

merged_collapsed <- add_missing_cols(merged_collapsed, all_cols)
mastermetid_clean <- add_missing_cols(mastermetid_clean, all_cols)
moore_list_clean <- add_missing_cols(moore_list_clean, all_cols)

# Step 1: Set flag for rows not in merged_collapsed
# For mastermetid_clean
master_only <- mastermetid_clean %>%
  anti_join(merged_collapsed, by = c("metid", "uid_01", "hmdb_id")) %>%
  mutate(both_moore_comets = "FALSE")

# For moore_list_clean
moore_only <- moore_list_clean %>%
  anti_join(
    merged_collapsed,
    by = c(
      "metabolite_id" = "uid_01",
      "input_hmdb_id" = "hmdb_id",
      "biochemical_lower" = "biochemical_lower"
    )
  ) %>%
  mutate(
    source_comets = NA,
    both_moore_comets = FALSE
  )

# Ensure all have same columns
all_cols <- Reduce(union, list(
  names(merged_collapsed),
  names(master_only),
  names(moore_only)
))

# Add missing columns filled with NA_character_
standardize_df <- function(df, all_cols) {
  df %>%
    mutate(across(everything(), as.character)) %>%  # ensure consistent types
    mutate(across(all_of(setdiff(all_cols, names(df))), ~NA_character_)) %>%
    select(all_of(all_cols))
}

merged_collapsed_std <- standardize_df(merged_collapsed, all_cols)
master_only_std <- standardize_df(master_only, all_cols)
moore_only_std <- standardize_df(moore_only, all_cols)

# Bind everything together
final_merged <- bind_rows(
  merged_collapsed_std,
  master_only_std,
  moore_only_std
)

nrow(final_merged)

##########################################################
#     Merged Dataset Creation (Adding MetLinkR Names)    #
##########################################################
moore_comets <- readxl::read_excel("mapping_library_COMETS_Moore.xlsx")
moore_comets_overlapping <- moore_comets %>% filter(`Input name (COMETS)` != "-" & `Input name (Moore)` != "-") # 1688 rows

# Make sure identifiers are character
final_merged_metlinkr <- final_merged %>%
  mutate(across(c(uid_01, hmdb_id, metid, input_metabid), as.character))

# Prepare metlinkr wide: gather all identifiers into one column for joining
metlinkr_long <- moore_comets %>%
  pivot_longer(
    cols = starts_with("Input name"),
    names_to = "origin_type",
    values_to = "input_id"
  ) %>%
  separate_rows(input_id, sep = ";") %>%   # split multiple ids in one cell
  filter(input_id != "-" & input_id != "")

# Join to final_merged on multiple possible identifiers
final_merged_metlinkr <- final_merged_metlinkr %>%
  left_join(
    metlinkr_long,
    by = c("uid_01" = "input_id"), relationship="many-to-many"
  ) %>%
  left_join(
    metlinkr_long %>% rename(input_hmdb_id_mlink = input_id),
    by = c("hmdb_id" = "input_hmdb_id_mlink"), relationship="many-to-many"
  ) %>%
  left_join(
    metlinkr_long %>% rename(metid_mlink = input_id),
    by = c("metid" = "metid_mlink"), relationship="many-to-many"
  ) %>%
  mutate(
    Harmonized_name = coalesce(`Harmonized name.x`, `Harmonized name.y`, `Harmonized name`),
    in_metlinkr = !is.na(Harmonized_name)
  ) %>%
  select(-`Harmonized name.x`, -`Harmonized name.y`, -`Harmonized name`)

final_merged_metlinkr_filt <- final_merged_metlinkr %>% select(
  -`Origin (COMETS).x`, -`Origin (Moore).x`, -`origin_type.x`,
  -`Origin (COMETS).y`, -`Origin (Moore).y`, -`origin_type.y`, -`Origin (COMETS)`,
  -`Origin (Moore)`, -`origin_type`
) %>% distinct()

write.csv(final_merged_metlinkr_filt, "final_merged_metlinkr_draft.csv")
save(final_merged_metlinkr_filt, file="final_merged_metlinkr_filt.Rdata")

#########################################
#     Analysis of merged dataset        #
#########################################

final_merged_metlinkr_filt.all_3 <- final_merged_metlinkr_filt %>% filter(
  in_metlinkr == "TRUE" &
    both_moore_comets == "TRUE"
)


cols_to_check <- c("hmdb_id", "biochemical_lower", "Harmonized_name")

final_merged_metlinkr_filt.all_3.unique <- final_merged_metlinkr_filt.all_3 %>%
  distinct(across(any_of(cols_to_check)))

final_merged_metlinkr_filt.all_3.unique_all_col <- final_merged_metlinkr_filt.all_3 %>%
  distinct(across(any_of(cols_to_check)), .keep_all = TRUE)

################################################
#     Sample harmonization (ALSPAC): Draft     #
################################################

alspac.orig <-readxl::read_xlsx("fromEGG_alspac_metabolite_list.xlsx")
alspac_biochem_name <- "CHEMICAL_NAME"
alspac_hmdb_name <- "HMDB"
alspac_kegg_name <- "KEGG"
alspac_pubchem_name <- "PUBCHEM"
alpsac_inchikey_name <- "INCHIKEY"


expand_dataframes <- function(input_df=final_merged_metlinkr_filt, input_column, punctuation=c("|", ";")) {
  escape_regex <- function(x) {
    if (x %in% c(".", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "?")) {
      paste0("\\", x)
    } else {
      x
    }
  }
  
  regex_sep <- str_c(sapply(punctuation, escape_regex), collapse = "|")
  
  merged_df_expanded <- input_df %>%
    separate_rows({{input_column}}, sep = regex_sep) %>%
    mutate({{ input_column }} := str_trim({{ input_column }}))
  
  return (merged_df_expanded)
  
}


standardize_hmdb <- function(input_df=final_merged_metlinkr_filt, input_hmdb_col){
  input_df.hmdb_fixed <- input_df %>%
    mutate(
      !!sym(input_hmdb_col) := ifelse(
        !is.na(.data[[input_hmdb_col]]) &
          .data[[input_hmdb_col]] != "" &
          nchar(.data[[input_hmdb_col]]) != 11,
        str_replace(.data[[input_hmdb_col]], "HMDB", "HMDB00"),
        .data[[input_hmdb_col]]
      )
    )
  
  return(input_df.hmdb_fixed)
}


metlinkr_merge_hmdb_expanded <- expand_dataframes(input_column = hmdb_id, punctuation=c("#", ",", "_", ";"))
metlinkr_merge_comets_hmdb_std <- standardize_hmdb(input_df=metlinkr_merge_hmdb_expanded, input_hmdb_col = "hmdb_id")
metlinkr_merge_hmdb_expanded.2 <- expand_dataframes(input_df=metlinkr_merge_comets_hmdb_std, input_column = input_hmdb_id, punctuation=c("#", ",", "_", ";"))
metlinkr_merge_hmdb_std <- standardize_hmdb(input_df = metlinkr_merge_hmdb_expanded.2, input_hmdb_col = "input_hmdb_id")


alspac_hmdb_expanded <- expand_dataframes(input_df = alspac.orig, input_column = HMDB, punctuation = c(","))
alspac_hmdb_expanded <- standardize_hmdb(input_df = alspac_hmdb_expanded, input_hmdb_col = "HMDB")

alspac_harm.hmdb.1 <- metlinkr_merge_hmdb_std %>%
  filter(!is.na(hmdb_id), hmdb_id != "") %>%
  inner_join(
    alspac_hmdb_expanded %>%
      filter(!is.na(HMDB), HMDB != "", HMDB != "NA"),
    by = c("hmdb_id" = "HMDB"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_remain.hmdb.1 <- alspac_hmdb_expanded %>% filter(!(HMDB %in% alspac_harm.hmdb.1$HMDB))

alspac_harm.hmdb.2 <- metlinkr_merge_hmdb_std %>%
  filter(!is.na(input_hmdb_id), hmdb_id != "") %>%
  inner_join(
    alspac_remain.hmdb.1 %>%
      filter(!is.na(HMDB), HMDB != "", HMDB != "NA"),
    by = c("input_hmdb_id" = "HMDB"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_harm.hmdb <- bind_rows(alspac_harm.hmdb.1, alspac_harm.hmdb.2)

alspac_remain.hmdb <- alspac_hmdb_expanded %>% filter(!(HMDB %in% alspac_harm.hmdb$HMDB))

metlinkr_merge_kegg_expanded <- expand_dataframes(input_df =metlinkr_merge_hmdb_std, input_column = input_kegg, punctuation = c("|"))
alspac_kegg_expanded <- expand_dataframes(input_df=alspac_remain.hmdb, input_column=KEGG, punctuation=c(","))


alspac_harm.kegg <- metlinkr_merge_kegg_expanded %>%
  filter(!is.na(input_kegg), input_kegg != "") %>%
  inner_join(
    alspac_kegg_expanded %>%
      filter(!is.na(KEGG), KEGG != "", KEGG != "NA"),
    by = c("input_kegg" = "KEGG"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_harm <- bind_rows(alspac_harm.hmdb, alspac_harm.kegg)
alspac_remain.kegg <- alspac_kegg_expanded %>% filter(!(KEGG %in% alspac_harm$KEGG))

metlinkr_merge_pubchem_expanded <- expand_dataframes(input_df=metlinkr_merge_kegg_expanded, input_column = input_pubchem, punctuation = c("|"))
alspac_pubchem_expanded <- expand_dataframes(input_df=alspac_remain.kegg, input_column=PUBCHEM, punctuation=c(","))

alspac_harm.pubchem <- metlinkr_merge_pubchem_expanded %>%
  filter(!is.na(input_pubchem), input_pubchem != "") %>%
  inner_join(
    alspac_pubchem_expanded %>%
      filter(!is.na(PUBCHEM), PUBCHEM != "", PUBCHEM != "NA"),
    by = c("input_pubchem" = "PUBCHEM"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_harm <- bind_rows(alspac_harm, alspac_harm.pubchem)
alspac_remain.pubchem <- alspac_pubchem_expanded %>% filter(!(PUBCHEM %in% alspac_harm$PUBCHEM))

metlinkr_merge_inchikey_expanded <- expand_dataframes(input_df=metlinkr_merge_pubchem_expanded, input_column = input_inchikey, punctuation = c("|"))
alspac_inchikey_expanded <- expand_dataframes(input_df=alspac_remain.pubchem, input_column=INCHIKEY, punctuation=c(","))

alspac_harm.inchikey <- metlinkr_merge_inchikey_expanded %>%
  filter(!is.na(input_inchikey), input_inchikey != "") %>%
  inner_join(
    alspac_inchikey_expanded %>%
      filter(!is.na(INCHIKEY), INCHIKEY != "", INCHIKEY != "NA"),
    by = c("input_inchikey" = "INCHIKEY"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_harm <- bind_rows(alspac_harm, alspac_harm.inchikey)
alspac_remain.inchikey <- alspac_pubchem_expanded %>% filter(!(INCHIKEY %in% alspac_harm$INCHIKEY))

metlinkr_merge_chemspider_expanded <- expand_dataframes(input_df=metlinkr_merge_inchikey_expanded, input_column = input_chemspider, punctuation = c("|"))
alspac_chemspider_expanded <- expand_dataframes(input_df=alspac_remain.inchikey, input_column=CHEMSPIDER, punctuation=c(","))

alspac_harm.chemspider <- metlinkr_merge_chemspider_expanded %>%
  filter(!is.na(input_chemspider), input_chemspider != "") %>%
  inner_join(
    alspac_chemspider_expanded %>%
      filter(!is.na(CHEMSPIDER), CHEMSPIDER != "", CHEMSPIDER != "NA"),
    by = c("input_chemspider" = "CHEMSPIDER"),
    suffix = c("", "_alspac"),
    relationship="many-to-many",
    keep=TRUE
  )

alspac_harm <- bind_rows(alspac_harm, alspac_harm.chemspider)
alspac_remain.chemspider <- alspac_chemspider_expanded %>% filter(!(CHEMSPIDER %in% alspac_harm$CHEMSPIDER))


metlinkr_merge_biochem_name_expanded <- expand_dataframes(input_df=metlinkr_merge_chemspider_expanded, input_column = input_metabolite_name, punctuation = c("|", "or"))

alspac_harm.biochem.moore <- metlinkr_merge_biochem_name_expanded %>%
  filter(!is.na(input_metabolite_name), input_metabolite_name != "") %>%  
  inner_join(
    alspac_remain.chemspider %>%
      filter(!is.na(CHEMICAL_NAME), CHEMICAL_NAME != "", CHEMICAL_NAME != "NA"),
    by = c("input_metabolite_name" = "CHEMICAL_NAME"),
    suffix = c("", "_alspac"),
    relationship = "many-to-many",
    keep = TRUE
  )


alspac.biochem.1 <- alspac_remain.chemspider %>% filter(
  !(CHEMICAL_NAME %in% alspac_harm.biochem.moore$CHEMICAL_NAME)
)

alspac_harm <- bind_rows(alspac_harm, alspac_harm.biochem.moore)


alspac_harm.biochem.comets <- metlinkr_merge_biochem_name_expanded %>% 
  filter(!is.na(biochemical_lower), biochemical_lower != "") %>%  
  inner_join(
    alspac.biochem.1 %>%
      filter(!is.na(CHEMICAL_NAME), CHEMICAL_NAME != "", CHEMICAL_NAME != "NA"),
    by = c("biochemical_lower" = "CHEMICAL_NAME"),
    suffix = c("", "_alspac"),
    relationship = "many-to-many",
    keep = TRUE
  )


alspac.biochem.2 <- alspac.biochem.1 %>% filter(
  !(CHEMICAL_NAME %in% alspac_harm.biochem.moore$CHEMICAL_NAME)
)

alspac_harm <- bind_rows(alspac_harm, alspac_harm.biochem.comets)


alspac_harm.biochem.metlinkr <- metlinkr_merge_biochem_name_expanded %>%
  filter(!is.na(Harmonized_name), Harmonized_name != "") %>%
  mutate(Harmonized_name_lower = tolower(Harmonized_name)) %>%
  inner_join(
    alspac.biochem.2 %>%
      filter(!is.na(CHEMICAL_NAME), CHEMICAL_NAME != "", CHEMICAL_NAME != "NA"),
    by = c("Harmonized_name_lower" = "CHEMICAL_NAME"),
    suffix = c("", "_alspac"),
    relationship = "many-to-many",
    keep = TRUE
  ) %>%
    select(-Harmonized_name_lower)

  
alspac_harm <- bind_rows(alspac_harm, alspac_harm.biochem.metlinkr)
  
alspac_remain <- alspac.biochem.2 %>% filter(
  !(CHEMICAL_NAME %in% alspac_harm.biochem.moore$CHEMICAL_NAME)
)




