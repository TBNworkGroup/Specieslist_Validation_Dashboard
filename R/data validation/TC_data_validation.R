rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)



# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20250528"  # 舉例


# 先抓第一頁
splist <- fromJSON(sprintf("https://api.taicol.tw/v2/taxon?&limit=300"))
df_TCsplist <- splist$data
total <- splist$info$total

# 計算頁碼
pg <- floor(total / 300)
if (total %% 300 == 0) {
  pg <- pg - 1
}
sequence <- seq(300, pg * 300, by = 300)

# loop 從第二頁開始
for (i in sequence) {
  splist <- fromJSON(sprintf("https://api.taicol.tw/v2/taxon?&limit=300&offset=%d", i))
  df_TCsplist <- rbind(df_TCsplist, splist$data)
  print(i)
}


# highertaxon <- fromJSON(sprintf("https://api.taicol.tw/v2/higherTaxa?taxon_id=%s", df_TCsplist$taxon_id[33609]))
# df_highertaxon <- highertaxon$data
# filter(
#   taxonRank %in% c("species", "infraspecies"),
#   alien_type == ""
# )


fwrite(df_TCsplist, sprintf("../../data/input/TCsplist_%s.csv", modified_date))






# (2) 讀取檔案 & 篩選欄位
df_TCsplist <- fread(sprintf("../../data/input/TCsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE) %>% 
  filter(taxon_status%in%"accepted")



# ------------------------------------------------------------------
# Part A: 重複學名比對（重複的 simplifiedScientificName）
# ------------------------------------------------------------------
# 這部份只要偵測學名重複的分類群
# 第一階段：學名重複
# 第二階段：同一界（Kingdom）+ 學名重複
# 第三階段：同一界（Kingdom）+ 命名者 + 學名重複
# 最後輸出一張表 TT_repeated
# Step 1: 篩選欄位
df_TCrepeated <- df_TCsplist %>%
  select(
    taxon_id, rank,
    simple_name, name_author
  )

# Step 2: 判斷重複樣態，並標記為 reason
df_taxa_duplicates <- df_TCrepeated %>%
  group_by(simple_name) %>%
  mutate(is_dup_simple_name = n() > 1) %>%
  ungroup() %>%
  
  group_by(simple_name, name_author) %>%
  mutate(is_dup_name_author = n() > 1) %>%
  ungroup()
  


# Step 3: 建立 reason 欄位（依照條件逐一指定）
# 先定義三種條件的篩選與理由
dup_name <- df_taxa_duplicates %>%
  filter(is_dup_simple_name) %>%
  mutate(reason = "學名重複")


dup_name_author <- df_taxa_duplicates %>%
  filter(is_dup_name_author) %>%
  mutate(reason = "學名加命名者重複")

# 再把它們合併起來
df_duplicates_reasoned <- bind_rows(dup_name, dup_name_author) %>%
  distinct()  # 去掉完全重複列（避免重複列入）



# Step 4: 加入連結欄位
df_duplicates_reasoned$TC_URL <- sprintf("https://taicol.tw/zh-hant/taxon/%s", df_duplicates_reasoned$taxon_id)

# Step 5: 輸出結果
fwrite(df_duplicates_reasoned, "../../data/output/TC_duplicates_result.csv")




# ------------------------------------------------------------------
# Part B: 學名欄位錯誤樣態檢核
# ------------------------------------------------------------------
# 先定義兩個檢查函式 check_string_lower()、check_string_vers()。
df_TCsubset <- df_TCsplist %>%
  select(
    taxon_id, rank, simple_name
  )

check_string_vers_TC <- function(string) {
  reasons <- character(0)
  
  if (str_detect(string, "&") || str_detect(string, "_") || str_detect(string, "\\.")) {
    reasons <- c(reasons, "特殊符號錯誤")
  }
  if (str_detect(string, " \\)") || str_detect(string, "\\( ")) {
    reasons <- c(reasons, "括號前後空白")
  }
  if (str_starts(string, " ")) {
    reasons <- c(reasons, "文字前空格")
  }
  if (str_ends(string, " ")) {
    reasons <- c(reasons, "文字後空格")
  }
  if (str_detect(string, "  ")) {
    reasons <- c(reasons, "連續空格")
  }
  
  return(reasons)
}

check_higher_rank_word_error <- function(string, rank) {
  if (!rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) return(FALSE)
  if (str_detect(str_to_lower(string), "incertae sedis")) return(FALSE)
  return(str_count(string, "\\S+") > 1)
}

check_higher_rank_case_error <- function(string, rank) {
  if (!rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) return(FALSE)
  if (str_detect(str_to_lower(string), "incertae sedis")) return(FALSE)
  
  first <- substr(string, 1, 1)
  rest  <- substr(string, 2, nchar(string))
  return(!(str_to_upper(first) == first && str_to_lower(rest) == rest))
}


errors_list_TC <- list()
target_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")

for (i in seq_len(nrow(df_TCsubset))) {
  row_data <- df_TCsubset[i, ]
  val <- row_data[["simple_name"]]
  rank <- row_data[["rank"]]
  
  if (rank %in% target_ranks) {
    # 高階層多詞
    if (check_higher_rank_word_error(val, rank)) {
      errors_list_TC <- append(errors_list_TC, list(
        data.frame(
          taxon_id = row_data[["taxon_id"]],
          rank = rank,
          simple_name = val,
          reason = "高階層多詞格式錯誤",
          stringsAsFactors = FALSE
        )
      ))
    }
    # 高階層大小寫
    if (check_higher_rank_case_error(val, rank)) {
      errors_list_TC <- append(errors_list_TC, list(
        data.frame(
          taxon_id = row_data[["taxon_id"]],
          rank = rank,
          simple_name = val,
          reason = "高階層大小寫格式錯誤",
          stringsAsFactors = FALSE
        )
      ))
    }
    
    # 一般空白與符號錯誤
    error_reasons <- check_string_vers_TC(val)
    if (length(error_reasons) > 0) {
      for (err in error_reasons) {
        errors_list_TC <- append(errors_list_TC, list(
          data.frame(
            taxon_id = row_data[["taxon_id"]],
            rank = rank,
            simple_name = val,
            reason = err,
            stringsAsFactors = FALSE
          )
        ))
      }
    }
  }
}


df_TC_errors <- do.call(rbind, errors_list_TC)
df_TC_errors$TC_URL <- sprintf("https://taicol.tw/zh-hant/taxon/%s", df_TC_errors$taxon_id)


# Output
fwrite(df_TC_errors, "../../data/output/TC_errortypes_result.csv")


# ------------------------------------------------------------------
# Part C: 原生性與敏感狀態比對
# ------------------------------------------------------------------
# 這部份只要偵測原生性與敏感狀態有問題的分類群
# 第一階段：檢查「種」與「種下」階層是否有原生性（TC）是空白的分類群
# 第二階段：挑出敏感狀態 = 無的保育類
# 第三階段：挑出敏感狀態 = 無的國內紅皮書等級高於「VU（含）」的物種
# 第四階段：挑出敏感狀態 = 無的國際紅皮書等級高於「VU（含）」的物種
# 第五階段：檢查外來種的敏感狀態（TT專屬）：挑出敏感狀態 != 無的外來種
# 最後輸出一張表df_TT_attribute_error

df_TC_attribute <- df_TCsplist %>%
  select(
    taxon_id, rank, simple_name, 
    alien_type, is_in_taiwan,               
    protected, redlist, iucn, cites, sensitive
  )

df_TC_undertaxon <- df_TC_attribute %>%
  filter(
    rank %in% c("Species", "Subspecies"),
    alien_type == "",
    is_in_taiwan == "TRUE",
  )
df_TC_undertaxon$reason <- "種與種下原生性空白"

df_TC_protected <- df_TC_attribute %>%
  filter(
    !(alien_type %in% c("cultured", "invasive", "naturalized")),
    is_in_taiwan == "TRUE",
    sensitive =="",
    protected != ""
  )

df_TC_protected$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"

df_TC_redlist <- df_TC_attribute %>%
  filter(
    !(alien_type %in% c("cultured", "invasive", "naturalized")),
    is_in_taiwan == "TRUE",
    sensitive == "",
    redlist %in% c(
      "VU",
      "RE",
      "EW",
      "CR",
      "EX",
      "EN"
    )
  )


df_TC_redlist$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"

df_TC_IUCN <- df_TC_attribute %>%
  filter(
    !(alien_type %in% c("cultured", "invasive", "naturalized")),
    is_in_taiwan == "TRUE",
    sensitive == "",
    iucn %in% c(
      "VU",
      "EX",
      "CR",
      "EN"
    )
  )

df_TC_IUCN$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"



df_TC_invasive <- df_TC_attribute %>%
  filter(
    sensitive != "",
    alien_type %in% c("cultured", "invasive", "naturalized"),
    is_in_taiwan == "TRUE",
  )

df_TC_invasive$reason <- "敏感狀態不等於無的外來種"

df_TC_attribute_error <- rbind(df_TC_undertaxon, df_TC_redlist, df_TC_protected, df_TC_IUCN, df_TC_invasive)



df_TC_attribute_error$TC_URL <- sprintf("https://taicol.tw/zh-hant/taxon/%s", df_TC_attribute_error$taxon_id)
fwrite(df_TC_attribute_error, "../../data/output/TC_attributeerror_result.csv")





