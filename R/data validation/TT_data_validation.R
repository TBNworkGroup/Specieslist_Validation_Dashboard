rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)



# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20260121"  # 舉例

# (2) 讀取檔案 & 篩選欄位
df_TTsplist <- fread(sprintf("../../data/input/TT/TTsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE)



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



for (i in 1:nrow(df_duplicates_reasoned)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_duplicates_reasoned$taxonUUID[i]))
      break
    }, error = function(e) {
      message("Error occurred: ", e)
      message("Retrying after 5 seconds")
      Sys.sleep(5)
    })
  }
  
  if (TBN_result$meta$status == "SUCCESS") {
    df_duplicates_reasoned$number_of_occurrence[i] <- TBN_result$meta$total
  } else if (TBN_result$meta$status == "NOT FOUND") {
    df_duplicates_reasoned$number_of_occurrence[i] <- 0
  } else {
    df_duplicates_reasoned$number_of_occurrence[i] <- TBN_result$meta$status
  }
  
  print(paste("finish i=", i, " download"))
}

# Step 4: 加入連結欄位
df_duplicates_reasoned$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_duplicates_reasoned$taxonUUID)

# Step 5: 輸出結果
fwrite(df_duplicates_reasoned, "../../data/output/TT_duplicates_result_vers2.csv")




# ------------------------------------------------------------------
# Part B: 學名欄位錯誤樣態檢核
# ------------------------------------------------------------------
# 先定義兩個檢查函式 check_string_lower()、check_string_vers()。
# ------------------------
# 輔助檢查函式（回傳錯誤標籤）
# ------------------------

check_string_lower <- function(string, columnname) {
  # 如果是 under species 欄位，必須全部小寫或以 '×' 開頭
  if (columnname %in% c("specificEpithet", "subspecies", "variety", "form", "cultigen", "cultivar")) {
    if (str_to_lower(string) == string || str_starts(string, "×")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (!columnname %in% c("taxonUUID", "taxonRank", "simplifiedScientificName")) {
      first_letter <- substr(string, 1, 1)
      remainder    <- substr(string, 2, nchar(string))
      if ((str_to_upper(first_letter) == first_letter &&
           str_to_lower(remainder) == remainder) ||
          str_starts(string, "×")) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    return(TRUE)
  }
}



check_all_errors <- function(string, columnname, taxonRank) {
  reasons <- character(0)
  
  if (str_detect(string, " \\)") || str_detect(string, "\\( ") ||
      str_detect(string, "&") || str_detect(string, "_") || str_detect(string, "\\.")) {
    reasons <- c(reasons, "錯誤符號與括號前後空格")
  }
  if (str_starts(string, " ")) reasons <- c(reasons, "文字前空格")
  if (str_ends(string, " "))   reasons <- c(reasons, "文字後空格")
  if (str_detect(string, "  ")) reasons <- c(reasons, "連續空格")
  
  # 大小寫錯誤（僅適用於非 authorship 類）
  if (!(columnname %in% c("authorship", "subspeciesAuthorship", "varietyAuthorship", "formAuthorship"))) {
    if (!check_string_lower(string, columnname)) {
      reasons <- c(reasons, "大小寫錯誤")
    }
  }
  
  # 高階層多詞檢查（僅 simplifiedScientificName + 高階層 taxonRank）
  if (columnname == "simplifiedScientificName" &&
      taxonRank %in% c("kingdom", "phylum", "class", "order", "family", "genus") &&
      !str_detect(str_to_lower(string), "incertae sedis") &&
      str_count(string, "\\S+") > 1) {
    reasons <- c(reasons, "高階層欄位出現多詞格式")
  }
  
  return(unique(reasons))
}

# ------------------------
# 主檢查邏輯
# ------------------------
# ============================
# Initialize 儲存錯誤記錄的清單
# ============================
error_records <- list()

for (i in seq_len(nrow(df_TTsplist))) {
  row <- df_TTsplist[i, ]
  uuid <- row$taxonUUID
  rank <- row$taxonRank
  kingdom <- row$kingdom
  
  # 基本結構（其他欄位初始化為 NA）
  base <- list(
    taxonUUID = uuid,
    taxonRank = rank,
    kingdom = kingdom,
    simplifiedScientificName = NA_character_,
    specificEpithet = NA_character_,
    subspecies = NA_character_,
    variety = NA_character_,
    form = NA_character_,
    cultigen = NA_character_,
    authorship = NA_character_,
    subspeciesAuthorship = NA_character_,
    varietyAuthorship = NA_character_,
    formAuthorship = NA_character_,
    errortypes = NA_character_,
    TT_URL = sprintf("https://taxatree.tbn.org.tw/taxa/%s", uuid)
  )
  
  # ===== 類型 1：高階層 simplifiedScientificName =====
  if (rank %in% c("kingdom", "phylum", "class", "order", "family", "genus")) {
    val <- row$simplifiedScientificName
    if (!is.na(val) && val != "") {
      reasons <- check_all_errors(val, "simplifiedScientificName", rank)
      for (reason in reasons) {
        record <- base
        record$simplifiedScientificName <- val
        record$errortypes <- reason
        error_records[[length(error_records) + 1]] <- as.data.frame(record, stringsAsFactors = FALSE)
      }
    }
  }
  
  # ===== 類型 2：下階層字串欄位 =====
  for (col in c("specificEpithet", "subspecies", "variety", "form", "cultigen")) {
    val <- row[[col]]
    if (!is.na(val) && val != "") {
      reasons <- check_all_errors(val, col, rank)
      for (reason in reasons) {
        record <- base
        record[[col]] <- val
        record$errortypes <- reason
        error_records[[length(error_records) + 1]] <- as.data.frame(record, stringsAsFactors = FALSE)
      }
    }
  }
  
  # ===== 類型 3：命名者欄位空白檢查 =====
  authorship_cols <- c()
  if (rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species ")) {
    authorship_cols <- c("authorship")
  } else if (rank == "subspecies") {
    authorship_cols <- c("subspeciesAuthorship", "varietyAuthorship", "formAuthorship")
  }
  
  for (col in authorship_cols) {
    val <- row[[col]]
    if (!is.na(val) && val != "") {
      # 只檢查空白錯誤
      reasons <- character(0)
      if (str_starts(val, " ")) reasons <- c(reasons, "文字前空格")
      if (str_ends(val, " "))   reasons <- c(reasons, "文字後空格")
      if (str_detect(val, "  ")) reasons <- c(reasons, "連續空格")
      
      for (reason in reasons) {
        record <- base
        record[[col]] <- val
        record$errortypes <- reason
        error_records[[length(error_records) + 1]] <- as.data.frame(record, stringsAsFactors = FALSE)
      }
    }
  }
}

# 合併結果為資料框
df_errors <- do.call(rbind, error_records)


for (i in 1:nrow(df_errors)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_errors$taxonUUID[i]))
      break
    }, error = function(e) {
      message("Error occurred: ", e)
      message("Retrying after 5 seconds")
      Sys.sleep(5)
    })
  }
  
  if (TBN_result$meta$status == "SUCCESS") {
    df_errors$number_of_occurrence[i] <- TBN_result$meta$total
  } else if (TBN_result$meta$status == "NOT FOUND") {
    df_errors$number_of_occurrence[i] <- 0
  } else {
    df_errors$number_of_occurrence[i] <- TBN_result$meta$status
  }
  
  print(paste("finish i=", i, " download"))
}




df_errors$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_errors$taxonUUID)


# (B3) 輸出到csv
fwrite(df_errors, "../../data/output/TT_errortypes_result.csv")

# ------------------------------------------------------------------
# Part C: 原生性與敏感狀態比對
# ------------------------------------------------------------------
# 這部份只要偵測原生性與敏感狀態有問題的分類群
# 第一階段：檢查「種」與「種下」階層是否有原生性（TT）是空白的分類群
# 第二階段：挑出敏感狀態 = 無的保育類
# 第三階段：挑出敏感狀態 = 無的國內紅皮書等級高於「VU（含）」的物種
# 第四階段：挑出敏感狀態 = 無的國際紅皮書等級高於「VU（含）」的物種
# 第五階段：檢查外來種的敏感狀態（TT專屬）：挑出敏感狀態 != 無的外來種
# 第六階段：挑出原生性=「外來栽培」「外來歸化」的動物
# 最後輸出一張表df_TT_attribute_error

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


if (nrow(df_TT_attribute_error) > 0) {
  
  for (i in 1:nrow(df_TT_attribute_error)) {
    
    while(TRUE) {
      tryCatch({
        TBN_result <- fromJSON(sprintf(
          "https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20",
          df_TT_attribute_error$taxonUUID[i]
        ))
        break
      }, error = function(e) {
        message("Error occurred: ", e)
        message("Retrying after 5 seconds")
        Sys.sleep(5)
      })
    }
    
    if (TBN_result$meta$status == "SUCCESS") {
      df_TT_attribute_error$number_of_occurrence[i] <- TBN_result$meta$total
    } else if (TBN_result$meta$status == "NOT FOUND") {
      df_TT_attribute_error$number_of_occurrence[i] <- 0
    } else {
      df_TT_attribute_error$number_of_occurrence[i] <- TBN_result$meta$status
    }
    
    print(paste("finish i=", i, " download"))
  }
  
} else {
  message("🛑 df_TT_attribute_error 沒有資料，跳過整個查詢迴圈。")
}




df_TT_attribute_error$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_TT_attribute_error$taxonUUID)
fwrite(df_TT_attribute_error, "../../data/output/TT_attributeerror_result.csv")


# ------------------------------------------------------------------
# Part D: 沒有種以下階層的分類觀點
# ------------------------------------------------------------------
# 這部份只要偵測分類群沒有子階層的「屬」以上階層的分類群
# 第一階段：檢查沒有「種」與「種下」階層的分類群
# 最後輸出一張表df_TT_without_species

df_TT_taxon <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, parentUUID, simplifiedScientificName, kingdom, class
  )


df_TT_without_species <- df_TT_taxon %>%
  # 1. 留下 taxonUUID 不在 parentUUID 集合裡
  filter(! taxonUUID %in% df_TT_taxon$parentUUID) %>%
  # 2. 先排除 rank 是 species、subspecies
  filter(! taxonRank %in% c("species", "infraspecies")) 

for (i in 1:nrow(df_TT_without_species)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_TT_without_species$taxonUUID[i]))
      break
    }, error = function(e) {
      message("Error occurred: ", e)
      message("Retrying after 5 seconds")
      Sys.sleep(5)
    })
  }
  
  if (TBN_result$meta$status == "SUCCESS") {
    df_TT_without_species$number_of_occurrence[i] <- TBN_result$meta$total
  } else if (TBN_result$meta$status == "NOT FOUND") {
    df_TT_without_species$number_of_occurrence[i] <- 0
  } else {
    df_TT_without_species$number_of_occurrence[i] <- TBN_result$meta$status
  }
  
  print(paste("finish i=", i, " download"))
}



df_TT_without_species$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_TT_without_species$taxonUUID)


fwrite(df_TT_without_species, "../../data/output/TT_without_species.csv")


# ------------------------------------------------------------------
# Part E: 種與種下階層的命名法規
# ------------------------------------------------------------------
# 這部份只要偵測種與種下階層的命名法規
# 第一階段：檢查植物「種」與「種下」階層的命名法規
# 第二階段：檢查動物「種」與「種下」階層的命名法規
# 第三階段：檢查「種」與「種下」階層的命名法規有不同的
# 最後輸出一張表df_TT_nomenclaturalCode

df_TT_species_nomenclaturalCode <- df_TTsplist %>%
  filter( taxonRank %in% c("species", "infraspecies")) %>% 
  select(
    taxonUUID, taxonRank, parentUUID, kingdom, simplifiedScientificName, 
    variety, form, cultigen, nomenclaturalCode 
  )

df_TT_plant_nomenclaturalCode <- df_TT_species_nomenclaturalCode %>%
  # 1. 留下 taxonUUID 不在 parentUUID 集合裡
  filter( kingdom %in% "Plantae") %>%
  # 2. 先排除 rank 是 species、subspecies
  filter(! nomenclaturalCode %in% c("ICN", "ICNCP"))
df_TT_plant_nomenclaturalCode$reason <- "命名法規錯誤的植物"

df_TT_animal_nomenclaturalCode <- df_TT_species_nomenclaturalCode %>%
  # 1. 留下 taxonUUID 不在 parentUUID 集合裡
  filter( kingdom %in% "Animalia") %>%
  # 2. 先排除 rank 是 species、subspecies
  filter(! nomenclaturalCode %in% c("三名法、二名法"))
df_TT_animal_nomenclaturalCode$reason <- "命名法規錯誤的動物"

# 新增 groupID 欄位：species 用 taxonUUID，自身是 infraspecies 的用 parentUUID
df_species_grouped <- df_TT_species_nomenclaturalCode %>%
  filter(taxonUUID %in% parentUUID | taxonRank == "infraspecies") %>%
  mutate(
    groupID = case_when(
      taxonRank == "species" ~ taxonUUID,
      taxonRank == "infraspecies" ~ parentUUID
    )
  ) %>%
  split(., .$groupID)

# 定義一個函式：檢查 nomenclaturalCode 是否一致
has_nomenclaturalCode_difference <- function(df) {
  length(unique(df$nomenclaturalCode)) > 1
}

# 找出有差異的 group，並加入 reason
conflict_groups <- Filter(has_nomenclaturalCode_difference, df_species_grouped)
df_species_nomenclaturalCode_mismatch <- bind_rows(conflict_groups, .id = "groupID") %>%
  mutate(reason = "種與種下命名法規不同")


df_TT_animaliaerror <- df_TT_species_nomenclaturalCode %>%
  filter(
    kingdom == "Animalia" &
      (
        !is.na(variety) & variety != "" |
          !is.na(form) & form != "" |
          !is.na(cultigen) & cultigen != ""
      )
  ) %>% 
  mutate(reason = "出現變種名、型名、栽培類名的動物")




df_TT_nomenclaturalCode <- bind_rows(df_TT_animal_nomenclaturalCode, df_TT_plant_nomenclaturalCode, df_species_nomenclaturalCode_mismatch, df_TT_animaliaerror)

for (i in 1:nrow(df_TT_nomenclaturalCode)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_TT_nomenclaturalCode$taxonUUID[i]))
      break
    }, error = function(e) {
      message("Error occurred: ", e)
      message("Retrying after 5 seconds")
      Sys.sleep(5)
    })
  }
  
  if (TBN_result$meta$status == "SUCCESS") {
    df_TT_nomenclaturalCode$number_of_occurrence[i] <- TBN_result$meta$total
  } else if (TBN_result$meta$status == "NOT FOUND") {
    df_TT_nomenclaturalCode$number_of_occurrence[i] <- 0
  } else {
    df_TT_nomenclaturalCode$number_of_occurrence[i] <- TBN_result$meta$status
  }
  
  print(paste("finish i=", i, " download"))
}




df_TT_nomenclaturalCode$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_TT_nomenclaturalCode$taxonUUID)

fwrite(df_TT_nomenclaturalCode, "../../data/output/TT_nomenclaturalCode.csv")

# ------------------------------------------------------------------
# Part F: 種與種下階層的屬性資料是否一致
# ------------------------------------------------------------------
# 這部份只要偵測種與種下階層的命名法規
# 第一階段：檢查「種」與「種下」階層的命名法規有不同的
# 第二階段：檢查「種」與「種下」階層的敏感狀態有不同的
# 第三階段：檢查「種」與「種下」階層的保育等級有不同的
# 第四階段：檢查「種」與「種下」階層的國際紅皮書有不同的
# 第五階段：檢查「種」與「種下」階層的國內紅皮書有不同的
# 第六階段：檢查「種」與「種下」階層的原生性有不同的
# 最後輸出一張表df_TT_nomenclaturalCode

df_TT_speciesinfraspecies_attribute <- df_TTsplist %>%
  filter( taxonRank %in% c("species", "infraspecies")) %>% 
  select(
    taxonUUID, taxonRank, parentUUID, kingdom, simplifiedScientificName, 
    nativeness, protectedStatusTW, categoryRedlistTW, categoryIUCN
  )
  
df_species_list <- df_TT_speciesinfraspecies_attribute %>%
  # 1. 留下 taxonUUID 不在 parentUUID 集合裡
  filter(taxonUUID %in% df_TT_speciesinfraspecies_attribute$parentUUID|taxonRank %in% "infraspecies") %>%
  # 新增 groupID
  mutate(
    groupID = case_when(
      taxonRank == "species" ~ taxonUUID,
      taxonRank == "infraspecies" ~ parentUUID
    )
  ) %>%
  # 依 groupID 分組並 split
  split(., .$groupID)

# 🔁 每個 group 做檢查：每個欄位的差異產生一筆紀錄
records <- list()

# 要比對的欄位與對應原因
check_columns <- list(
  protectedStatusTW = "保育等級不同",
  categoryRedlistTW = "國內紅皮書不同",
  categoryIUCN = "IUCN紅皮書不同"
)

# 定義紅皮書「不比較」的值
ignore_vals <- c("不適用", "暫無危機（LC, Least Concern）")

# 遍歷每一組 group
for (group_id in names(df_species_list)) {
  group_df <- df_species_list[[group_id]]
  
  for (col in names(check_columns)) {
    distinct_vals <- unique(group_df[[col]])
    #distinct_vals <- distinct_vals[distinct_vals != ""]  # 移除空白
    
    # 對於紅皮書類別，要先排除 ignore 值
    if (col %in% c("categoryRedlistTW", "categoryIUCN")) {
      filtered_vals <- setdiff(distinct_vals, ignore_vals)
    } else {
      filtered_vals <- distinct_vals
    }
    
    # 如果排除後還有超過1個不同值，才視為差異
    if (length(filtered_vals) > 1) {
      group_df$reason <- check_columns[[col]]
      group_df$check_column <- col
      group_df$filtered_values <- paste(filtered_vals, collapse = ";")
      records[[length(records) + 1]] <- group_df
    }
  }
}

# 將所有有問題的 group 綁在一起
df_speciesinfraspecies_attribute_mismatch <- bind_rows(records)



for (i in 1:nrow(df_speciesinfraspecies_attribute_mismatch)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_speciesinfraspecies_attribute_mismatch$taxonUUID[i]))
      break
    }, error = function(e) {
      message("Error occurred: ", e)
      message("Retrying after 5 seconds")
      Sys.sleep(5)
    })
  }
  
  if (TBN_result$meta$status == "SUCCESS") {
    df_speciesinfraspecies_attribute_mismatch$number_of_occurrence[i] <- TBN_result$meta$total
  } else if (TBN_result$meta$status == "NOT FOUND") {
    df_speciesinfraspecies_attribute_mismatch$number_of_occurrence[i] <- 0
  } else {
    df_speciesinfraspecies_attribute_mismatch$number_of_occurrence[i] <- TBN_result$meta$status
  }
  
  print(paste("finish i=", i, " download"))
}





df_speciesinfraspecies_attribute_mismatch$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_speciesinfraspecies_attribute_mismatch$taxonUUID)
fwrite(df_speciesinfraspecies_attribute_mismatch, "../../data/output/TT_speciesinfraspecies_attribute_mismatch.csv")

 

