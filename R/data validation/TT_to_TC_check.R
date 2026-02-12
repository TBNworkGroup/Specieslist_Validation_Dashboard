rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table","progressr")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)

# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20260128"  # 舉例

# (2) 讀取檔案 & 篩選欄位
df_TTsplist <- fread(sprintf("../../data/input/TT/TTsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE)

df_TCsplist <- fread(sprintf("../../data/input/TC/TCsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE) %>% 
  filter(is_in_taiwan%in%"TRUE")



df_TT_select <- df_TTsplist %>%
  select(taxonUUID, taxonRank, kingdom, taiCOLNameCode, simplifiedScientificName) %>% 
  setnames(., c("taxonUUID", "taxonRank", "kingdom", "taiCOLNameCode", "simplifiedScientificName"), c("TT_taxonUUID", "TT_taxonRank", "TT_kingdom", "TT_taiCOLNameCode", "TT_simplifiedScientificName"))

df_TC_select <- df_TCsplist %>% 
  select(taxon_id, rank, kingdom, simple_name)%>%
  setnames(., c("taxon_id", "simple_name", "rank", "kingdom"), c("TC_taxon_id", "TC_simple_name", "TC_rank", "TC_kingdom")) %>% 
  .[, TC_rank := tolower(TC_rank)]


# ------------------------------------------------------------------
# Part A: TT與TC，taxon_id、simplifiedScientificName、rank、kindgdom完全一致的分類群
# ------------------------------------------------------------------
TT_TC_all_same <- df_TT_select %>%
  inner_join(
    df_TC_select,
    by = c(
      "TT_simplifiedScientificName" = "TC_simple_name",
      "TT_taxonRank"               = "TC_rank",
      "TT_kingdom"                 = "TC_kingdom",
      "TT_taiCOLNameCode"          = "TC_taxon_id"
    )
  )

fwrite(TT_TC_all_same, "../../data/output/TT_to_TC/TT_nochange.csv")

# ------------------------------------------------------------------
# Part B: TT與TC，simplifiedScientificName、rank、kindgdom完全一致的分類群
# ------------------------------------------------------------------

TT_checkrank_taxon_id <- df_TT_select[!(TT_taxonUUID %in% as.vector(TT_TC_all_same$TT_taxonUUID))]%>%
  inner_join(
    df_TC_select,
    by = c(
      "TT_simplifiedScientificName" = "TC_simple_name",
      "TT_taxonRank"               = "TC_rank",
      "TT_kingdom"                 = "TC_kingdom"
    )
  )

fwrite(TT_checkrank_taxon_id, "../../data/output/TT_to_TC/TT_checkrank_taxon_id.csv")
# ------------------------------------------------------------------
# Part C: TT與TC，simplifiedScientificName、rank、kindgdom完全一致的分類群
# ------------------------------------------------------------------

TT_checkrank_scientificName <- df_TT_select[!(TT_taxonUUID %in% as.vector(TT_TC_all_same$TT_taxonUUID))]%>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_taxon_id$TT_taxonUUID))] %>% 
  inner_join(
    df_TC_select,
    by = c(
      "TT_taxonRank"               = "TC_rank",
      "TT_kingdom"                 = "TC_kingdom",
      "TT_taiCOLNameCode"          = "TC_taxon_id"
    )
  )

fwrite(TT_checkrank_scientificName, "../../data/output/TT_to_TC/TT_checkrank_scientificName.csv")
# ------------------------------------------------------------------
# Part D: TT已經是非接受名的分類群
# ------------------------------------------------------------------

TT_add_taxonid <- df_TT_select[!(TT_taxonUUID %in% as.vector(TT_TC_all_same$TT_taxonUUID))]%>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_taxon_id$TT_taxonUUID))] %>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_scientificName$TT_taxonUUID))] %>% 
  filter(TT_taiCOLNameCode%in%"")

handlers(global = TRUE)
handlers("txtprogressbar")  # console 版進度條

query_taicol_one <- function(sciname){
  
  url <- sprintf(
    "https://api.taicol.tw/v2/name?scientific_name=%s",
    URLencode(sciname, reserved = TRUE)
  )
  
  res <- tryCatch(
    fromJSON(url),
    error = function(e) NULL
  )
  
  if (is.null(res) || res$info$total == 0) {
    return(
      tibble(
        taicol_taxon_id = NA_character_,
        taicol_status   = NA_character_
      )
    )
  }
  
  res$data %>%
    as_tibble() %>%
    separate(
      taxon,
      into = c("taicol_taxon_id", "taicol_status", "taicol_flag"),
      sep = ",\\s*",
      fill = "right"
    ) %>%
    select(taicol_taxon_id, taicol_status)
}

with_progress({
  p <- progressor(along = TT_add_taxonid$TT_simplifiedScientificName)
  
  taicol_lookup <- TT_add_taxonid %>%
    distinct(TT_simplifiedScientificName) %>%
    mutate(
      taicol_result = map(
        TT_simplifiedScientificName,
        ~{
          p(sprintf("Querying: %s", .x))
          Sys.sleep(0.1)  # API 節流
          query_taicol_one(.x)
        }
      )
    ) %>%
    unnest(taicol_result)
})




TT_add_taxonid_checked <- TT_add_taxonid %>%
  left_join(
    taicol_lookup,
    by = "TT_simplifiedScientificName"
  ) %>% .[!(is.na(.$taicol_taxon_id))]
fwrite(TT_add_taxonid_checked, "../../data/output/TT_to_TC/TT_add_taxonid.csv")

TT_unknow <- df_TT_select[!(TT_taxonUUID %in% as.vector(TT_TC_all_same$TT_taxonUUID))]%>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_taxon_id$TT_taxonUUID))] %>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_scientificName$TT_taxonUUID))] %>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_add_taxonid_checked$TT_taxonUUID))]
fwrite(TT_unknow, "../../data/output/TT_to_TC/TT_unknow.csv")
