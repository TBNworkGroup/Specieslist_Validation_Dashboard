
library(dplyr)

# 1. 設定預期父分類階層
rank_rules <- tibble::tribble(
  ~taxonRank,      ~expected_parentRank,
  "phylum",       "kingdom",
  "class",        "phylum",
  "order",        "class",
  "family",       "order",
  "genus",        "family",
  "species",      "genus",
  "infraspecies", "species"
)

# 2. 整理分類群資料
df_TT_taxon <- df_TTsplist %>%
  dplyr::select(
    taxonUUID, taxonRank, kingdom,
    simplifiedScientificName, parentUUID
  ) %>%
  mutate(
    taxonRank = tolower(trimws(as.character(taxonRank)))
  )

# 3. 建立父分類群資料，沿用你重新命名的欄位
df_TT_parent <- df_TT_taxon %>%
  transmute(
    parentUUID = taxonUUID,
    parent_taxonRank = taxonRank,
    parent_name = simplifiedScientificName
  )

# 4. 篩出父分類階層不符合規則的紀錄
df_TT_parent_rank_mismatch <- df_TT_taxon %>%
  inner_join(
    df_TT_parent,
    by = "parentUUID",
    relationship = "many-to-one"
  ) %>%
  inner_join(
    rank_rules,
    by = "taxonRank"
  ) %>%
  filter(parent_taxonRank != expected_parentRank)




fwrite(df_TT_parent_rank_mismatch, "../../data/output/TT_infrataxoncheck.csv")
fwrite(df_TT_parent_rank_mismatch, "../../data/output/TT_to_TC_beta/TT_infrataxoncheck.csv")
