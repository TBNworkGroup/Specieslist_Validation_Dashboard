rm(list=ls())

usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)



# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20250521"  # 舉例

# (2) 讀取檔案 & 篩選欄位
df_TTsplist <- fread(sprintf("../../data/input/TTsplist_%s.csv", modified_date), sep = ",", fill=TRUE, encoding = "UTF-8", colClasses="character", header=TRUE)



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
    taxonUUID, taxonRank, kingdom,
    simplifiedScientificName, scientificName
  )

# Step 2: 判斷重複樣態，並標記為 reason
df_taxa_duplicates <- df_TTrepeated %>%
  group_by(simplifiedScientificName) %>%
  mutate(is_dup_global = n() > 1) %>%
  ungroup() %>%
  
  group_by(kingdom, simplifiedScientificName) %>%
  mutate(is_dup_kingdom = n() > 1) %>%
  ungroup() %>%
  
  group_by(kingdom, scientificName) %>%
  mutate(is_dup_kingdom_author = n() > 1) %>%
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

# 再把它們合併起來
df_duplicates_reasoned <- bind_rows(dup_global, dup_kingdom, dup_kingdom_author) %>%
  distinct()  # 去掉完全重複列（避免重複列入）



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
fwrite(df_duplicates_reasoned, "../../data/output/TT_duplicates_result.csv")




# ------------------------------------------------------------------
# Part B: 學名欄位錯誤樣態檢核
# ------------------------------------------------------------------
# 先定義兩個檢查函式 check_string_lower()、check_string_vers()。

check_string_lower <- function(string, columnname) {
  # 如果是 under species 欄位，必須全部小寫或以 '×' 開頭
  if (columnname %in% c("specificEpithet", "subspecies", "variety", "form", "cultigen", "cultivar")) {
    # 注意：原 Python 似乎是 "cultivar"；如有 "cultigen"，自行確認。
    # 判定方式：全部小寫或以 '×' 開頭
    if (str_to_lower(string) == string || str_starts(string, "×")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    # 其他欄位(非 taxonUUID, taxonRank, simplifiedScientificName)，
    # 要求：首字大寫 + 後面小寫；或以 '×' 開頭
    if (!columnname %in% c("taxonUUID", "taxonRank", "simplifiedScientificName")) {
      # 先擷取字串首字母、其餘部分
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
    # 若是 taxonUUID, taxonRank, simplifiedScientificName，就不檢查大小寫
    return(TRUE)
  }
}

check_string_vers <- function(string, columnname) {
  # 預設沒有問題 (FALSE 表示"沒有問題")
  check <- FALSE
  
  # 1. 括號前後空白 & 特定符號
  #    - ' )', '( ', '&', '_', '.' 均視為錯誤
  if (str_detect(string, " \\)")  ||
      str_detect(string, "\\( ")  ||
      str_detect(string, "&")     ||
      str_detect(string, "_")     ||
      str_detect(string, "\\.")) {
    check <- TRUE
  }
  
  # 2. 檢查頭尾空白、連續空白
  if (str_starts(string, " ")) {
    check <- TRUE
  }
  if (str_ends(string, " ")) {
    check <- TRUE
  }
  if (str_detect(string, "  ")) {
    check <- TRUE
  }
  
  # 3. 大小寫檢查
  if (!check_string_lower(string, columnname)) {
    check <- TRUE
  }
  
  return(check)
}

check_higher_rank_format <- function(string, columnname) {
  # 只對 kingdom 到 genus 檢查
  if (!columnname %in% c("kingdom", "phylum", "class", "order", "family", "genus")) {
    return(FALSE)  # 不檢查
  }
  
  # 包含 incertae sedis 則視為合法
  if (str_detect(str_to_lower(string), "incertae sedis")) {
    return(FALSE)
  }
  
  # 字串內若出現多個詞（空白分隔 > 1），則不合法
  word_count <- str_count(string, "\\S+")
  if (word_count > 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# (B1) 對非 simplifiedScientificName 欄位執行檢查
# 模仿 Python 的 iterrows()，逐列 + 逐欄檢查
errors_list <- list()  # 收集所有錯誤紀錄


df_TTsubset <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, kingdom, phylum, class, order, family,
    genus, specificEpithet, subspecies, variety, form, cultigen,
    simplifiedScientificName
  )


for (i in seq_len(nrow(df_TTsubset))) {
  row_data <- df_TTsubset[i, ]
  
  # 準備一個空的紀錄（同 df_TTsubset 的所有欄位，先設 NA）
  error_record <- as.list(rep(NA, ncol(df_TTsubset)))
  names(error_record) <- colnames(df_TTsubset)
  
  # 依序檢查每個欄位
  for (colname in colnames(df_TTsubset)) {
    val <- row_data[[colname]]
    
    if (is.na(val)) {
      next
    }
    
    # simplifiedScientificName 不檢查
    if (colname == "simplifiedScientificName") {
      next
    }
    
    # 執行檢查
    check_result <- FALSE
    
    check_result <- tryCatch({
      # 加入新的檢查邏輯
      if (check_higher_rank_format(val, colname)) {
        error_record[["taxonUUID"]] <- row_data[["taxonUUID"]]
        error_record[["taxonRank"]] <- row_data[["taxonRank"]]
        error_record[[colname]]     <- val
        errors_list <- append(errors_list, list(error_record))
        break
      }
      check_string_vers(val, colname)
    }, error = function(e) {
      # 如果函式本身執行時出錯，就回傳 NA 讓後面好判斷
      NA
    })
    
    # 若出現 NA，代表檢查函式壞掉，這邊就直接跳出該列
    if (is.na(check_result)) {
      break
    }
    if (check_result) {
      # 一旦發現錯誤，記錄 taxonUUID, taxonRank, 問題欄位的值
      error_record[["taxonUUID"]] <- row_data[["taxonUUID"]]
      error_record[["taxonRank"]] <- row_data[["taxonRank"]]
      error_record[[colname]]     <- val
      
      # 放入錯誤清單
      errors_list <- append(errors_list, list(error_record))
      # 因為該列已確定有問題，不檢查其他欄位，直接跳出
      break
    }
  }
}

# (B2) 將錯誤清單轉成 DataFrame
df_errors <- do.call(rbind, lapply(errors_list, as.data.frame))
df_errors <- as.data.frame(df_errors, stringsAsFactors = FALSE)

check_string_vers_detail_extended <- function(string, columnname) {
  reasons <- character(0)
  
  # 舊有邏輯
  if (str_detect(string, " \\)") || str_detect(string, "\\( ") ||
      str_detect(string, "&") || str_detect(string, "_") || str_detect(string, "\\.")) {
    reasons <- c(reasons, "錯誤符號與括號前後空格")
  }
  if (str_starts(string, " ")) reasons <- c(reasons, "文字前空格")
  if (str_ends(string, " ")) reasons <- c(reasons, "文字後空格")
  if (str_detect(string, "  ")) reasons <- c(reasons, "連續空格")
  if (!check_string_lower(string, columnname)) reasons <- c(reasons, "大小寫錯誤")
  
  # 新增邏輯
  if (check_higher_rank_format(string, columnname)) {
    reasons <- c(reasons, "高階層欄位出現多詞格式")
  }
  
  if (length(reasons) == 0) {
    return("")
  } else {
    return(paste(unique(reasons), collapse = ";"))
  }
}


# 先新增欄位 errortypes
df_errors$errortypes <- NA_character_

# 遍歷 df_errors 的每一列, 找出哪個欄位是出錯欄位(即有值), 
# 然後用 check_string_vers_detail() 取得錯誤種類
for (i in seq_len(nrow(df_errors))) {
  # 這行數據
  row_data <- df_errors[i, ]
  
  # 假設只有1個欄位(除了 taxonUUID, taxonRank, TT_URL... ) 會存到值
  # 先找出 "非 NA" 的欄位
  non_na_cols <- colnames(row_data)[which(!is.na(row_data) & row_data != "")]
  
  # 排除不需要檢查的欄位 (taxonUUID, taxonRank, TT_URL, simplifiedScientificName等)
  # 你可自行增減排除清單
  exclude_cols <- c("taxonUUID","taxonRank","TT_URL","simplifiedScientificName")
  flagged_cols <- setdiff(non_na_cols, exclude_cols)
  
  if (length(flagged_cols) == 1) {
    # 就用這個欄位為 "出錯欄位"
    colname <- flagged_cols[1]
    val <- row_data[[colname]]
    
    # 執行加強版檢查 -> 回傳一串錯誤描述
    error_str <- check_string_vers_detail_extended(val, colname)
    df_errors$errortypes[i] <- error_str
  } else if (length(flagged_cols) > 1) {
    # 若不只1個欄位(理論上不該發生, 因為你 break 了),
    # 這裡可以自行決定怎麼處理, 例如只檢查第一個
    colname <- flagged_cols[1]
    val <- row_data[[colname]]
    error_str <- error_str <- check_string_vers_detail_extended(val, colname)
    df_errors$errortypes[i] <- error_str
  } else {
    # flagged_cols 長度是 0 => 找不到出錯欄位, 可能都 NA => 不做事
  }
}
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
# 最後輸出一張表df_TT_attribute_error

df_TT_attribute <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, simplifiedScientificName, 
    endemism, nativeness,               
    protectedStatusTW, categoryRedlistTW, categoryIUCN, sensitiveCategory
    )

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
      "區域滅絕（RE, Regionally Extinct）",
      "野外滅絕（EW, Extinct in the Wild）",
      "極危（CR, Critically Endangered）",
      "滅絕（EX, Extinct）",
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
      "野外滅絕（EW, Extinct in the Wild）",
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

df_TT_attribute_error <- rbind(df_TT_undertaxon, df_TT_redlist, df_TT_protected, df_TT_IUCN, df_TT_invasive)

for (i in 1:nrow(df_TT_attribute_error)) {
  
  while(TRUE) {
    tryCatch({
      TBN_result <- fromJSON(sprintf("https://www.tbn.org.tw/api/v25/occurrence?taxonUUID=%s&limit=20", df_TT_attribute_error$taxonUUID[i]))
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
    nomenclaturalCode 
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






df_TT_nomenclaturalCode <- bind_rows(df_TT_animal_nomenclaturalCode, df_TT_plant_nomenclaturalCode, df_species_nomenclaturalCode_mismatch)

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

# 遍歷每一組 group
for (group_id in names(df_species_list)) {
  group_df <- df_species_list[[group_id]]
  
  for (col in names(check_columns)) {
    distinct_vals <- unique(group_df[[col]])
    # 保留 NA 與空字串以便檢查完整差異
    if (length(distinct_vals) > 1) {
      group_df$reason <- check_columns[[col]]
      group_df$check_column <- col  # 可選擇是否保留這個輔助欄
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


#### === 確認結果 ===
# A: 重複學名 -> df_duplicates_result
# B: 欄位錯誤 -> df_errors
# C: 屬性資料錯誤 <- df_TT_attribute_error
# D: 沒有種階層分類群 <- df_TT_without_species
# E: 種與種下階層的命名法規與保育等級 <- df_TT_species_attribute
# F: 種與種下階層的屬性資料是否一致 <- df_speciesinfraspecies_attribute_mismatch
print(df_duplicates_reasoned)
print(df_errors)
print(df_TT_attribute_error)
print(df_TT_without_species)
print(df_TT_species_attribute)
print(df_speciesinfraspecies_attribute_mismatch)
