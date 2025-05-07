rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)



# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20250507"  # 舉例


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
#   nativeness == ""
# )


fwrite(df_TCsplist, sprintf("../../data/input/TCsplist_%s.csv", modified_date))






# (2) 讀取檔案 & 篩選欄位
df_TCsplist <- fread(sprintf("../../data/input/TCsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE)



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
  mutate(is_dup_global = n() > 1) %>%
  ungroup() %>%
  
  group_by(rank, simple_name) %>%
  mutate(is_dup_kingdom = n() > 1) %>%
  ungroup() %>%
  
  group_by(rank, simple_name, name_author) %>%
  mutate(is_dup_kingdom_author = n() > 1) %>%
  ungroup()


df_duplicates_reasoned <- df_taxa_duplicates %>%
  mutate(reason = case_when(
    is_dup_kingdom_author ~ "分類階層學名加命名者重複",
    is_dup_kingdom ~ "分類階層學名重複",
    is_dup_global ~ "學名重複",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(reason))



# Step 4: 加入連結欄位
df_duplicates_reasoned$TC_URL <- sprintf("https://taicol.tw/zh-hant/taxon/%s", df_duplicates_reasoned$taxon_id)

# Step 5: 輸出結果
fwrite(df_duplicates_reasoned, "../../data/output/TT_duplicates_result.csv")
