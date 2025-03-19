
usepackage <- c("jsonlite", "tidyverse", "data.table")
install.packages(usepackage[!(usepackage %in% installed.packages()[,1])])
sapply(usepackage, library, character.only = TRUE)



# (1) 假設你有一個 modified_date 變數；如果沒有，就直接指定檔名。
modified_date <- "20250319"  # 舉例

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

df_TTrepeated <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, kingdom,
    simplifiedScientificName, scientificName
  )

df_taxa_duplicates <- df_TTrepeated %>%
  # (A) simplifiedScientificName (不分 kingdom)
  group_by(simplifiedScientificName) %>%
  mutate(dup_simplifiedGlobal = n() > 1) %>%
  ungroup() %>%
  
  # (B) 在相同 kingdom 下，simplifiedScientificName 重複
  group_by(kingdom, simplifiedScientificName) %>%
  mutate(dup_simplifiedKingdom = n() > 1) %>%
  ungroup() %>%
  
  # (C) 在相同 kingdom 下，scientificName 重複
  group_by(kingdom, scientificName) %>%
  mutate(dup_scientificKingdom = n() > 1) %>%
  ungroup()


df_duplicates_result <- df_taxa_duplicates %>%
  filter(
    dup_simplifiedGlobal |
      dup_simplifiedKingdom |
      dup_scientificKingdom
  )

df_duplicates_result$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_duplicates_result$taxonUUID)

fwrite(df_duplicates_result, "../../data/output/TT_duplicates_result.csv")


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

check_string_vers_detail <- function(string, columnname) {
  # 回傳違反規則的字串 (若沒錯誤 => "")
  reasons <- character(0)
  
  # 1. 括號前後空白 & 特定符號
  if (str_detect(string, " \\)") ||
      str_detect(string, "\\( ") ||
      str_detect(string, "&")    ||
      str_detect(string, "_")    ||
      str_detect(string, "\\.")) {
    reasons <- c(reasons, "錯誤符號與括號前後空格")
  }
  
  # 2. 頭尾或連續空白
  if (str_starts(string, " ")) reasons <- c(reasons, "文字前空格")
  if (str_ends(string, " "))   reasons <- c(reasons, "文字後空格")
  if (str_detect(string, "  ")) reasons <- c(reasons, "連續空格")
  
  # 3. 大小寫檢查
  if (!check_string_lower(string, columnname)) {
    reasons <- c(reasons, "大小寫錯誤")
  }
  
  if (length(reasons) == 0) {
    return("")
  } else {
    return(paste(unique(reasons), collapse=";"))
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
    error_str <- check_string_vers_detail(val, colname)
    df_errors$errortypes[i] <- error_str
  } else if (length(flagged_cols) > 1) {
    # 若不只1個欄位(理論上不該發生, 因為你 break 了),
    # 這裡可以自行決定怎麼處理, 例如只檢查第一個
    colname <- flagged_cols[1]
    val <- row_data[[colname]]
    error_str <- check_string_vers_detail(val, colname)
    df_errors$errortypes[i] <- error_str
  } else {
    # flagged_cols 長度是 0 => 找不到出錯欄位, 可能都 NA => 不做事
  }
}
df_errors$TT_URL <- sprintf("https://taxatree.tbn.org.tw/taxa/%s", df_errors$taxonUUID)


# (B3) 輸出到csv
fwrite(df_errors, "../../data/output/TT_errortypes_result.csv")

# ------------------------------------------------------------------
# Part C: 原生性與敏感狀態比對（重複的 simplifiedScientificName）
# ------------------------------------------------------------------
# 這部份只要偵測原生性與敏感狀態有問題的分類群
# 第一階段：檢查「種」與「種下」階層是否有原生性（TT）是空白的分類群
# 第二階段：挑出敏感狀態 = 無的保育類
# 第三階段：挑出敏感狀態 = 無的國內紅皮書等級高於「VU（含）」的物種
# 第四階段：挑出敏感狀態 = 無的國際紅皮書等級高於「VU（含）」的物種
# 第五階段：檢查外來種的敏感狀態（TT專屬）：挑出敏感狀態 /= 無的外來種
# 最後輸出一張表 TT_repeated

df_TTrepeated <- df_TTsplist %>%
  select(
    taxonUUID, taxonRank, kingdom,
    simplifiedScientificName, scientificName
  )



# === 確認結果 ===
# A: 重複學名 -> df_duplicates_result
# B: 欄位錯誤 -> df_errors
print(df_duplicates_result)
print(df_errors)
