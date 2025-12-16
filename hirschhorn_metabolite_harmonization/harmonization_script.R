################################################
#     Sample Harmonization (Generic): Draft     #
################################################

library(tidyverse)
library(readxl)

#TODO: Abstraction of harmonization
#TODO: Comments added soon

# -------------------------
# Input configuration
# -------------------------
load("final_merged_metlinkr_filt.RData")
target_df <- readxl::read_xlsx("fromEGG_alspac_metabolite_list.xlsx")

target_cols <- list(
  biochem = "CHEMICAL_NAME",
  hmdb = "HMDB",
  kegg = "KEGG",
  pubchem = "PUBCHEM",
  inchikey = "INCHIKEY",
  chemspider = "CHEMSPIDER")

# -------------------------
# Helper functions
# -------------------------
expand_dataframes <- function(input_df, input_column, punctuation = c("|", ";")) {
  escape_regex <- function(x) {
    if (x %in% c(".", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "?")) paste0("\\", x) else x
  }
  regex_sep <- str_c(sapply(punctuation, escape_regex), collapse = "|")
  
  input_df %>%
    separate_rows({{input_column}}, sep = regex_sep) %>%
    mutate({{input_column}} := str_trim({{input_column}}))
}

standardize_hmdb <- function(input_df, input_hmdb_col) {
  input_df %>%
    mutate(
      !!sym(input_hmdb_col) := ifelse(
        !is.na(.data[[input_hmdb_col]]) &
          .data[[input_hmdb_col]] != "" &
          nchar(.data[[input_hmdb_col]]) != 11,
        str_replace(.data[[input_hmdb_col]], "HMDB", "HMDB00"),
        .data[[input_hmdb_col]]
      )
    )
}

# -------------------------
# HMDB harmonization
# -------------------------
metlinkr_merge_hmdb_expanded <- expand_dataframes(input_df = final_merged_metlinkr_filt, input_column = hmdb_id, punctuation = c("#", ",", "_", ";"))
metlinkr_merge_comets_hmdb_std <- standardize_hmdb(input_df = metlinkr_merge_hmdb_expanded, input_hmdb_col = "hmdb_id")
metlinkr_merge_hmdb_expanded.2 <- expand_dataframes(input_df = metlinkr_merge_comets_hmdb_std, input_column = input_hmdb_id, punctuation = c("#", ",", "_", ";"))
metlinkr_merge_hmdb_std <- standardize_hmdb(input_df = metlinkr_merge_hmdb_expanded.2, input_hmdb_col = "input_hmdb_id")

target_hmdb_expanded <- expand_dataframes(input_df = target_df, input_column = !!sym(target_cols$hmdb), punctuation = c(","))
target_hmdb_expanded <- standardize_hmdb(input_df = target_hmdb_expanded, input_hmdb_col = target_cols$hmdb)

harm.hmdb.1 <- metlinkr_merge_hmdb_std %>%
  filter(!is.na(hmdb_id), hmdb_id != "") %>%
  inner_join(
    target_hmdb_expanded %>%
      filter(!is.na(.data[[target_cols$hmdb]]), .data[[target_cols$hmdb]] != "", .data[[target_cols$hmdb]] != "NA"),
    by = setNames(target_cols$hmdb, "hmdb_id"),
    suffix = c("", "_target"),
    relationship = "many-to-many",
    keep = TRUE
  )

remain.hmdb.1 <- target_hmdb_expanded %>% filter(!(!!sym(target_cols$hmdb) %in% harm.hmdb.1[[target_cols$hmdb]]))

harm.hmdb.2 <- metlinkr_merge_hmdb_std %>%
  filter(!is.na(input_hmdb_id), input_hmdb_id != "") %>%
  inner_join(
    remain.hmdb.1 %>%
      filter(!is.na(.data[[target_cols$hmdb]]), .data[[target_cols$hmdb]] != "", .data[[target_cols$hmdb]] != "NA"),
    by = setNames(target_cols$hmdb, "input_hmdb_id"),
    suffix = c("", "_target"),
    relationship = "many-to-many",
    keep = TRUE
  )

harm.hmdb <- bind_rows(harm.hmdb.1, harm.hmdb.2)
remain.hmdb <- target_hmdb_expanded %>% filter(!(!!sym(target_cols$hmdb) %in% harm.hmdb[[target_cols$hmdb]]))

# -------------------------
# KEGG harmonization
# -------------------------

if (!is.na(target_cols$kegg)) {
  metlinkr_merge_kegg_expanded <- expand_dataframes(input_df = metlinkr_merge_hmdb_std, input_column = input_kegg, punctuation = c("|"))
  target_kegg_expanded <- expand_dataframes(input_df = remain.hmdb, input_column = !!sym(target_cols$kegg), punctuation = c(","))
  
  harm.kegg <- metlinkr_merge_kegg_expanded %>%
    filter(!is.na(input_kegg), input_kegg != "") %>%
    inner_join(
      target_kegg_expanded %>%
        filter(!is.na(.data[[target_cols$kegg]]), .data[[target_cols$kegg]] != "", .data[[target_cols$kegg]] != "NA"),
      by = setNames(target_cols$kegg, "input_kegg"),
      suffix = c("", "_target"),
      relationship = "many-to-many",
      keep = TRUE
    )
  
  harm.all <- bind_rows(harm.hmdb, harm.kegg)
  remain.kegg <- target_kegg_expanded %>% filter(!(!!sym(target_cols$kegg) %in% harm.all[[target_cols$kegg]]))
  
} else {
  harm.all <- harm.hmdb
  remain.kegg <- remain.hmdb
}


# -------------------------
# PUBCHEM harmonization
# -------------------------
if (!is.na(target_cols$pubchem)){
  metlinkr_merge_pubchem_expanded <- expand_dataframes(input_df = metlinkr_merge_kegg_expanded, input_column = input_pubchem, punctuation = c("|"))
  target_pubchem_expanded <- expand_dataframes(input_df = remain.kegg, input_column = !!sym(target_cols$pubchem), punctuation = c(","))
  
  harm.pubchem <- metlinkr_merge_pubchem_expanded %>%
    filter(!is.na(input_pubchem), input_pubchem != "") %>%
    inner_join(
      target_pubchem_expanded %>%
        filter(!is.na(.data[[target_cols$pubchem]]), .data[[target_cols$pubchem]] != "", .data[[target_cols$pubchem]] != "NA"),
      by = setNames(target_cols$pubchem, "input_pubchem"),
      suffix = c("", "_target"),
      relationship = "many-to-many",
      keep = TRUE
    )
  
  harm.all <- bind_rows(harm.all, harm.pubchem)
  remain.pubchem <- target_pubchem_expanded %>% filter(!(!!sym(target_cols$pubchem) %in% harm.all[[target_cols$pubchem]]))
  
} else {
  remain.pubchem <- remain.kegg
}

# -------------------------
# INCHIKEY harmonization
# -------------------------
if (!is.na(target_cols$inchikey)){
  metlinkr_merge_inchikey_expanded <- expand_dataframes(input_df = metlinkr_merge_pubchem_expanded, input_column = input_inchikey, punctuation = c("|"))
  target_inchikey_expanded <- expand_dataframes(input_df = remain.pubchem, input_column = !!sym(target_cols$inchikey), punctuation = c(","))
  
  harm.inchikey <- metlinkr_merge_inchikey_expanded %>%
    filter(!is.na(input_inchikey), input_inchikey != "") %>%
    inner_join(
      target_inchikey_expanded %>%
        filter(!is.na(.data[[target_cols$inchikey]]), .data[[target_cols$inchikey]] != "", .data[[target_cols$inchikey]] != "NA"),
      by = setNames(target_cols$inchikey, "input_inchikey"),
      suffix = c("", "_target"),
      relationship = "many-to-many",
      keep = TRUE
    )
  
  harm.all <- bind_rows(harm.all, harm.inchikey)
  remain.inchikey <- target_inchikey_expanded %>% filter(!(!!sym(target_cols$inchikey) %in% harm.all[[target_cols$inchikey]]))
  
}else{
  remain.inchikey <- remain.pubchem
}


# -------------------------
# CHEMSPIDER harmonization
# -------------------------
if (!is.na(target_cols$chemspider)) {
  metlinkr_merge_chemspider_expanded <- expand_dataframes(input_df = metlinkr_merge_inchikey_expanded, input_column = input_chemspider, punctuation = c("|"))
  target_chemspider_expanded <- expand_dataframes(input_df = remain.inchikey, input_column = !!sym(target_cols$chemspider), punctuation = c(","))
  
  harm.chemspider <- metlinkr_merge_chemspider_expanded %>%
    filter(!is.na(input_chemspider), input_chemspider != "") %>%
    inner_join(
      target_chemspider_expanded %>%
        filter(!is.na(.data[[target_cols$chemspider]]), .data[[target_cols$chemspider]] != "", .data[[target_cols$chemspider]] != "NA"),
      by = setNames(target_cols$chemspider, "input_chemspider"),
      suffix = c("", "_target"),
      relationship = "many-to-many",
      keep = TRUE
    )
  
  harm.all <- bind_rows(harm.all, harm.chemspider)
  remain.chemspider <- target_chemspider_expanded %>% filter(!(!!sym(target_cols$chemspider) %in% harm.all[[target_cols$chemspider]]))
}else{
  remain.chemspider <- remain.inchikey
}



# -------------------------
# Biochemical name harmonization
# -------------------------
metlinkr_merge_biochem_name_expanded <- expand_dataframes(input_df = metlinkr_merge_chemspider_expanded, input_column = input_metabolite_name, punctuation = c("|", "or"))

harm.biochem.moore <- metlinkr_merge_biochem_name_expanded %>%
  filter(!is.na(input_metabolite_name), input_metabolite_name != "") %>%  
  inner_join(
    remain.chemspider %>%
      filter(!is.na(.data[[target_cols$biochem]]), .data[[target_cols$biochem]] != "", .data[[target_cols$biochem]] != "NA"),
    by = setNames(target_cols$biochem, "input_metabolite_name"),
    suffix = c("", "_target"),
    relationship = "many-to-many",
    keep = TRUE
  )

biochem_1 <- remain.chemspider %>% filter(!(!!sym(target_cols$biochem) %in% harm.biochem.moore[[target_cols$biochem]]))
harm.all <- bind_rows(harm.all, harm.biochem.moore)

harm.biochem.comets <- metlinkr_merge_biochem_name_expanded %>% 
  filter(!is.na(biochemical_lower), biochemical_lower != "") %>%  
  inner_join(
    biochem_1 %>%
      filter(!is.na(.data[[target_cols$biochem]]), .data[[target_cols$biochem]] != "", .data[[target_cols$biochem]] != "NA"),
    by = setNames(target_cols$biochem, "biochemical_lower"),
    suffix = c("", "_target"),
    relationship = "many-to-many",
    keep = TRUE
  )

biochem_2 <- biochem_1 %>% filter(!(!!sym(target_cols$biochem) %in% harm.biochem.comets[[target_cols$biochem]]))
harm.all <- bind_rows(harm.all, harm.biochem.comets)

harm.biochem.metlinkr <- metlinkr_merge_biochem_name_expanded %>%
  filter(!is.na(Harmonized_name), Harmonized_name != "") %>%
  mutate(Harmonized_name_lower = tolower(Harmonized_name)) %>%
  inner_join(
    biochem_2 %>%
      filter(!is.na(.data[[target_cols$biochem]]), .data[[target_cols$biochem]] != "", .data[[target_cols$biochem]] != "NA"),
    by = setNames(target_cols$biochem, "Harmonized_name_lower"),
    suffix = c("", "_target"),
    relationship = "many-to-many",
    keep = TRUE
  ) %>%
  select(-Harmonized_name_lower)

harm.all <- bind_rows(harm.all, harm.biochem.metlinkr)

remain.final <- biochem_2 %>%
  filter(!(!!sym(target_cols$biochem) %in% harm.all[[target_cols$biochem]]))

# -------------------------
# Output
# -------------------------
output <- list(
  harmonized = harm.all,
  remaining = remain.final
)

output

save(output, file="harmonized_output_alspac.Rdata")
