rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)

# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20260702"  # 舉例

# (2) 讀取檔案 & 篩選欄位
df_TTsplist <- fread(sprintf("../../data/input/TT/TTsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE)%>%
  filter(
    kingdom %in% "Animalia"
  )






# ------------------------------------------------------------------
# Part A: 重複學名比對（重複的 simplifiedScientificName）
# ------------------------------------------------------------------
# 這部份只要偵測學名重複的分類群
# 第一階段：學名重複
# 第二階段：同一界（Kingdom）+ 學名重複
# 第三階段：同一界（Kingdom）+ 命名者 + 學名重複
# 最後輸出一張表 TT_repeated
# Step 1: 篩選欄位
df_TTrepeated <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, kingdom, taiCOLNameCode, tfNameCode, 
    simplifiedScientificName, scientificName, 
  )

TT_authorship <- df_TTsplist %>%
  select(
    taxonUUID, authorship, subspeciesAuthorship, varietyAuthorship, formAuthorship
  )

# 假設 TT_authorship 已經讀入並為 data.table
cols_to_join <- c("authorship", "subspeciesAuthorship", "varietyAuthorship", "formAuthorship")

TT_authorship[, authorshipCheck := apply(.SD, 1, function(x) {
  # 過濾空字串
  non_empty <- x[nzchar(x)]
  # 如果全是空值回傳空字串，否則用"-"串接
  if (length(non_empty) == 0) "" else paste(non_empty, collapse = "-")
}), .SDcols = cols_to_join]

df_TTrepeated <- df_TTrepeated %>%
  left_join(
    TT_authorship[, .(taxonUUID, authorshipCheck)], 
    by = "taxonUUID"
  )

# Step 2: 判斷重複樣態，並標記為 reason
df_taxa_duplicates <- df_TTrepeated %>%
  group_by(simplifiedScientificName) %>%
  mutate(is_dup_global = n() > 1) %>%
  ungroup() %>%
  
  group_by(kingdom, simplifiedScientificName) %>%
  mutate(is_dup_kingdom = n() > 1) %>%
  ungroup() %>%
  
  group_by(kingdom, simplifiedScientificName, authorshipCheck) %>%
  mutate(is_dup_kingdom_author = n() > 1) %>%
  ungroup() %>%
  
  group_by(kingdom, simplifiedScientificName, authorshipCheck, taxonRank) %>%
  mutate(is_dup_kingdom_author_rank = n() > 1) %>%
  ungroup()




# Step 3: 建立 reason 欄位（依照條件逐一指定）
# 先定義三種條件的篩選與理由
dup_global <- df_taxa_duplicates %>%
  filter(is_dup_global) %>%
  mutate(reason = "學名重複")

dup_kingdom <- df_taxa_duplicates %>%
  filter(is_dup_kingdom) %>%
  mutate(reason = "相同Kingdom學名重複")

dup_kingdom_author <- df_taxa_duplicates %>%
  filter(is_dup_kingdom_author) %>%
  mutate(reason = "相同Kingdom學名加命名者重複")

dup_kingdom_author_rank <- df_taxa_duplicates %>%
  filter(is_dup_kingdom_author_rank) %>%
  mutate(reason = "相同Kingdom學名加命名者階層重複")

# ===== 新增：重複識別碼欄位檢查 =====
# 補充欄位（必要時先加入）
df_id_check <- df_TTsplist %>%
  select(taxonUUID, taiCOLNameCode, tfNameCode, taxonRank, kingdom, simplifiedScientificName, scientificName) %>%
  mutate(across(everything(), as.character))  # 保險起見全轉成字串

# 檢查各欄位是否重複
dup_taxonUUID <- df_id_check %>%
  group_by(taxonUUID) %>%
  filter(n() > 1) %>%
  mutate(reason = "taxonUUID重複")

dup_taiCOLNameCode <- df_id_check %>%
  filter(!is.na(taiCOLNameCode) & taiCOLNameCode != "") %>%
  group_by(taiCOLNameCode) %>%
  filter(n() > 1) %>%
  mutate(reason = "taiCOLNameCode重複")

dup_tfNameCode <- df_id_check %>%
  filter(!is.na(tfNameCode) & tfNameCode != "") %>%
  group_by(tfNameCode) %>%
  filter(n() > 1) %>%
  mutate(reason = "tfNameCode重複")

# 合併為一份資料框
df_id_duplicates <- bind_rows(
  dup_taxonUUID,
  dup_taiCOLNameCode,
  dup_tfNameCode) %>%
  distinct() %>%
  mutate(TT_URL = sprintf("https://taxatree.tbn.org.tw/taxa/%s", taxonUUID)) %>%
  select(taxonUUID, taxonRank, kingdom, taiCOLNameCode, tfNameCode, simplifiedScientificName, scientificName, reason)





df_duplicates_reasoned <-bind_rows(
  df_id_duplicates,
  dup_global,
  dup_kingdom,
  dup_kingdom_author,
  dup_kingdom_author_rank) %>%
  select(taxonUUID, taxonRank, kingdom, taiCOLNameCode, tfNameCode, simplifiedScientificName, scientificName, reason)

# 
# 
# for (i in 1:nrow(df_duplicates_reasoned)) {
#   
#   while(TRUE) {
#     tryCatch({
#       TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_duplicates_reasoned$taxonUUID[i]))
#       break
#     }, error = function(e) {
#       message("Error occurred: ", e)
#       message("Retrying after 5 seconds")
#       Sys.sleep(5)
#     })
#   }
#   
#   if (TBN_result$meta$status == "SUCCESS") {
#     df_duplicates_reasoned$number_of_occurrence[i] <- TBN_result$meta$total
#   } else if (TBN_result$meta$status == "NOT FOUND") {
#     df_duplicates_reasoned$number_of_occurrence[i] <- 0
#   } else {
#     df_duplicates_reasoned$number_of_occurrence[i] <- TBN_result$meta$status
#   }
#   
#   print(paste("finish i=", i, " download"))
# }

# Step 4: 加入連結欄位
df_duplicates_reasoned$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_duplicates_reasoned$taxonUUID)

# Step 5: 輸出結果
fwrite(df_duplicates_reasoned, "../../data/output/TT_duplicates_result_versAn.csv")




df_TT_attribute <- df_TTsplist %>%
  select(
    taxonUUID, taiCOLNameCode, taxonRank, kingdom, simplifiedScientificName, 
    endemism, nativeness,               
    protectedStatusTW, categoryRedlistTW, categoryIUCN, sensitiveCategory
  )

df_TT_withouttcnamecode <- df_TT_attribute %>%
  filter(
    taiCOLNameCode == ""
  )
df_TT_withouttcnamecode$reason <- "taiCOLNameCode空白"

df_TT_undertaxon <- df_TT_attribute %>%
  filter(
    taxonRank %in% c("species", "infraspecies"),
    nativeness == ""
  )
df_TT_undertaxon$reason <- "種與種下原生性空白"

df_TT_protected <- df_TT_attribute %>%
  filter(
    !(nativeness %like% "外"),
    sensitiveCategory =="",
    protectedStatusTW != ""
  )

df_TT_protected$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"

df_TT_redlist <- df_TT_attribute %>%
  filter(
    !(nativeness %like% "外"),
    sensitiveCategory == "",
    categoryRedlistTW %in% c(
      "易危（VU, Vulnerable）",
      "極危（CR, Critically Endangered）",
      "瀕危（EN, Endangered）"
    )
  )


df_TT_redlist$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"

df_TT_IUCN <- df_TT_attribute %>%
  filter(
    !(nativeness %like% "外"),
    sensitiveCategory == "",
    categoryIUCN %in% c(
      "易危（VU, Vulnerable）",
      "極危（CR, Critically Endangered）",
      "瀕危（EN, Endangered）"
    )
  )


df_TT_IUCN$reason <- "敏感狀態=無的保育類or國內紅皮書VU以上or國際IUCN VU以上的原生種"

df_TT_invasive <- df_TT_attribute %>%
  filter(
    sensitiveCategory != "",
    nativeness %like% "外"
  )

df_TT_invasive$reason <- "敏感狀態不等於無的外來種"

df_TT_animalia <- df_TT_attribute %>%
  filter(
    kingdom == "Animalia",
    nativeness %in% c("外來栽培 Cultivated (non-native)", "外來歸化 Naturalized (non-native)")
  )

df_TT_animalia$reason <- "原生性等於外來栽培or外來歸化的動物"


df_TT_attribute_error <- rbind(df_TT_withouttcnamecode, df_TT_undertaxon, df_TT_redlist, df_TT_protected, df_TT_IUCN, df_TT_invasive, df_TT_animalia)
# ✅ 在這裡加判斷：完全包住 for 迴圈


# if (nrow(df_TT_attribute_error) > 0) {
#   
#   for (i in 1:nrow(df_TT_attribute_error)) {
#     
#     while(TRUE) {
#       tryCatch({
#         TBN_result <- fromJSON(sprintf(
#           "https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20",
#           df_TT_attribute_error$taxonUUID[i]
#         ))
#         break
#       }, error = function(e) {
#         message("Error occurred: ", e)
#         message("Retrying after 5 seconds")
#         Sys.sleep(5)
#       })
#     }
#     
#     if (TBN_result$meta$status == "SUCCESS") {
#       df_TT_attribute_error$number_of_occurrence[i] <- TBN_result$meta$total
#     } else if (TBN_result$meta$status == "NOT FOUND") {
#       df_TT_attribute_error$number_of_occurrence[i] <- 0
#     } else {
#       df_TT_attribute_error$number_of_occurrence[i] <- TBN_result$meta$status
#     }
#     
#     print(paste("finish i=", i, " download"))
#   }
#   
# } else {
#   message("🛑 df_TT_attribute_error 沒有資料，跳過整個查詢迴圈。")
# }

df_TT_attribute_error$number_of_occurrence <- NA


df_TT_attribute_error$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_TT_attribute_error$taxonUUID)
fwrite(df_TT_attribute_error, "../../data/output/TT_attributeerror_result__versAn.csv")
