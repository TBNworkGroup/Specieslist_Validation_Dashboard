rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table","progressr")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)

# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20260225"  # 舉例

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


query_taicol_raw <- function(sciname){
  url <- sprintf(
    "https://api.taicol.tw/v2/name?scientific_name=%s",
    URLencode(sciname, reserved = TRUE)
  )
  
  res <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  
  if (is.null(res) || is.null(res$info$total) || res$info$total == 0) {
    return(tibble::tibble(
      hit_index = NA_integer_,
      taxon_raw = list(NULL)
    ))
  }
  
  dat <- tibble::as_tibble(res$data)
  
  dat %>%
    dplyr::mutate(
      hit_index = dplyr::row_number(),
      taxon_raw = purrr::map(taxon, ~ .x)
    ) %>%
    dplyr::select(hit_index, taxon_raw)
}



handlers(global = TRUE)
handlers("txtprogressbar")   # console 進度條

with_progress({
  
  p <- progressor(along = unique(TT_add_taxonid$TT_simplifiedScientificName))
  
  taicol_raw <- TT_add_taxonid %>%
    distinct(TT_simplifiedScientificName) %>%
    mutate(
      raw = map(
        TT_simplifiedScientificName,
        function(x){
          
          p(sprintf("Querying: %s", x))   # ← 進度更新
          Sys.sleep(0.1)                  # API 節流
          
          query_taicol_raw(x)
        }
      )
    ) %>%
    unnest(raw)
  
})



taicol_shape <- taicol_raw %>%
  mutate(
    taxon_class = map_chr(taxon_raw, ~ paste(class(.x), collapse="|")),
    taxon_type  = map_chr(taxon_raw, ~ typeof(.x)),
    taxon_len   = map_int(taxon_raw, ~ if (is.null(.x)) 0L else length(.x)),
    taxon_names = map_chr(taxon_raw, ~ {
      if (is.null(.x)) return("")
      n <- names(.x)
      if (is.null(n)) "" else paste(n, collapse="|")
    }),
    # 判斷是否含向量（例如 taxon_id 是 c("t1","t2")）
    has_vector_inside = map_lgl(taxon_raw, ~{
      if (!is.list(.x)) return(FALSE)
      any(map_lgl(.x, ~ is.atomic(.x) && length(.x) > 1))
    })
  ) %>%
  count(taxon_class, taxon_type, taxon_names, taxon_len, has_vector_inside, sort = TRUE)



normalize_taxon <- function(x){
  # x 可能是 NULL / character / list
  if (is.null(x) || length(x) == 0) {
    return(tibble(taicol_taxon_id = NA_character_, taicol_status = NA_character_))
  }
  
  # Case 1: 字串 "t0000001, accepted, ..."
  if (is.character(x) && length(x) >= 1) {
    # 若是一個向量字串，就逐一處理
    out <- map_dfr(x, ~{
      parts <- strsplit(.x, ",\\s*")[[1]]
      tibble(
        taicol_taxon_id = parts[1] %||% NA_character_,
        taicol_status   = parts[2] %||% NA_character_
      )
    })
    return(out)
  }
  
  # Case 2: list，有 named elements
  if (is.list(x)) {
    # 常見欄位名：taxon_id / taicol_name_status（依你截圖）
    taxon_id <- x$taxon_id %||% x$taicol_taxon_id %||% NA_character_
    status   <- x$taicol_name_status %||% x$taicol_status %||% x$status %||% NA_character_
    
    # 這裡處理「內含向量」：把 taxon_id / status 做成對齊的長度
    # 若 status 長度=1 但 taxon_id 多個，status 會 recycle
    n <- max(length(taxon_id), length(status), 1)
    taxon_id <- rep(taxon_id, length.out = n)
    status   <- rep(status,   length.out = n)
    
    return(tibble(
      taicol_taxon_id = as.character(taxon_id),
      taicol_status   = as.character(status)
    ))
  }
  
  # 其他型態：先保底
  tibble(taicol_taxon_id = NA_character_, taicol_status = NA_character_)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b



taicol_lookup_long <- taicol_raw %>%
  mutate(norm = map(taxon_raw, normalize_taxon)) %>%
  select(TT_simplifiedScientificName, hit_index, norm) %>%
  unnest(norm)

TT_add_taxonid_checked <- TT_add_taxonid %>%
  left_join(taicol_lookup_long, by = "TT_simplifiedScientificName")

TT_add_taxonid_hit <- TT_add_taxonid_checked %>%
  filter(!is.na(taicol_taxon_id) & taicol_taxon_id != "")

TT_add_taxonid_miss <- TT_add_taxonid_checked %>%
  filter(is.na(taicol_taxon_id) | taicol_taxon_id == "")

data.table::fwrite(taicol_shape, "../../data/output/TT_to_TC/taicol_taxon_shape_summary.csv")
data.table::fwrite(taicol_lookup_long, "../../data/output/TT_to_TC/taicol_lookup_long.csv")
data.table::fwrite(TT_add_taxonid_checked, "../../data/output/TT_to_TC/TT_add_taxonid_checked_long.csv")

TT_unknow <- df_TT_select[!(TT_taxonUUID %in% as.vector(TT_TC_all_same$TT_taxonUUID))]%>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_taxon_id$TT_taxonUUID))] %>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_checkrank_scientificName$TT_taxonUUID))] %>%
  .[!(.$TT_taxonUUID %in%  as.vector(TT_add_taxonid_checked$TT_taxonUUID))]
fwrite(TT_unknow, "../../data/output/TT_to_TC/TT_unknow.csv")
