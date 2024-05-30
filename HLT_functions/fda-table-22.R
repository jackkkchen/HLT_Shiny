make_table_22 <- function(df,
                          alt_counts_df = NULL,
                          show_colcounts = TRUE,
                          arm_var = "ARM",
                          saffl_var = "SAFFL",
                          vars = c("SEX", "AGEGR1", "RACE", "ETHNIC"),
                          denom = c("N_s", "N_col", "n"),
                          lbl_overall = NULL,
                          lbl_vars = formatters::var_labels(df, fill = TRUE)[vars],
                          prune_0 = FALSE,
                          annotations = NULL) {
  checkmate::assert_subset(c(vars, arm_var, saffl_var), names(df))
  
  df <- df %>%
    filter(.data[[saffl_var]] == "Y") %>%
    df_explicit_na()
  
  # For percentages calculations in case of N_s, add the overall observations
  denom <- match.arg(denom)
  if (!is.null(lbl_overall) && denom == "N_s") {
    df_ovrl <- df
    df_ovrl[[arm_var]] <- lbl_overall
    df <- rbind(df, df_ovrl)
    
    if (!is.null(alt_counts_df)) {
      alt_df_ovrl <- alt_counts_df
      alt_df_ovrl[[arm_var]] <- lbl_overall
      alt_counts_df <- rbind(alt_counts_df, alt_df_ovrl)
    }
  }
  
  alt_counts_df <- alt_counts_df_preproc(alt_counts_df, arm_var, saffl_var)
  
  lyt <- basic_table_annot(show_colcounts, annotations)
  
  lyt <- if (!is.null(lbl_overall) && denom != "N_s") {
    lyt %>% split_cols_by_arm(arm_var, lbl_overall)
  } else {
    lyt %>% split_cols_by_arm(arm_var)
  }
  
  lyt <- lyt %>%
    count_patients_with_event(
      "USUBJID",
      filters = c("TRTEMFL" = "Y"),
      .labels = c(count_fraction = "Any AE, n (%)")
    ) %>%
    analyze(
      vars = vars,
      var_labels = paste0(lbl_vars, ", n (%)"),
      afun = a_count_occurrences_trtem_ae,
      extra_args = list(
        denom = denom,
        arm_var = arm_var,
        df_denom = if (!is.null(alt_counts_df)) alt_counts_df else df
      ),
      show_labels = "visible"
    ) %>%
    append_topleft(c("", "Characteristic"))
  
  tbl <- build_table(lyt, df = df, alt_counts_df = alt_counts_df)
  if (prune_0) tbl <- prune_table(tbl)
  
  tbl
}

#' Analysis Function to Calculate Count/Fraction of Any Adverse Event Occurrences
#'
#' @inheritParams tern::s_count_occurrences
#' @inheritParams argument_convention
#' @param df_denom (`data.frame`)\cr Full data frame used to calculate denominator subgroup counts
#'   when `denom = "N_s"`.
#' @param denom (`character`)\cr Denominator to use to calculate fractions. Can be `"N_s"` (total `df_denom`
#'   subgroup/row counts), `"N_col"` (total `df` column counts), or `"n"` (total `df` overall patient count).
#'   Note that `df` is filtered to only include treatment-emergent adverse events (`TRTEMFL == "Y"`).
#'
#' @keywords internal
a_count_occurrences_trtem_ae <- function(df,
                                         .var,
                                         .N_col, # nolint
                                         df_denom = NULL,
                                         denom = c("N_s", "N_col", "n"),
                                         id_var = "USUBJID",
                                         arm_var = "ARM") {
  df <- df %>% filter(TRTEMFL == "Y")
  occurrences <- df[[.var]]
  ids <- factor(df[[id_var]])
  has_occurrence_per_id <- table(occurrences, ids) > 0
  n_ids_per_occurrence <- as.list(rowSums(has_occurrence_per_id))
  lvls <- names(n_ids_per_occurrence)
  
  denom <- match.arg(denom)
  if (denom == "N_s" && is.null(df_denom)) {
    stop("If using subgroup population counts, `df_denom` must be specified.") # nocov
  }
  dn <- switch(denom,
               N_s = lapply(
                 lvls,
                 function(x) {
                   df_denom %>%
                     filter(.data[[.var]] == x, .data[[arm_var]] == df[[arm_var]][1]) %>%
                     select(USUBJID) %>%
                     distinct() %>%
                     nrow()
                 }
               ),
               n = nlevels(ids),
               N_col = .N_col
  )
  if (denom == "N_s") names(dn) <- lvls
  
  x_stats <- lapply(
    lvls,
    function(x) {
      i <- n_ids_per_occurrence[[x]]
      denom <- if (denom == "N_s") dn[[x]] else dn
      if (i == 0 && denom == 0) c(0, 0) else c(i, i / denom)
    }
  )
  names(x_stats) <- names(n_ids_per_occurrence)
  
  in_rows(
    .list = x_stats,
    .formats = tern::format_count_fraction
  )
}

load_libraries_and_data <- function() {
  library(dplyr)
  library(falcon)
  library(random.cdisc.data)
  library(formatters)
  library(sassy)
  library(stringr)
  
  adsl <- random.cdisc.data::cadsl %>%
    mutate(AGEGR1 = as.factor(case_when(
      AGE >= 17 & AGE < 65 ~ ">=17 to <65",
      AGE >= 65 ~ ">=65",
      AGE >= 65 & AGE < 75 ~ ">=65 to <75",
      AGE >= 75 ~ ">=75"
    )) %>% formatters::with_label("Age Group, years")) %>%
    formatters::var_relabel(AGE = "Age, years")
  
  adae <- random.cdisc.data::cadae
  
  df <- left_join(adsl, adae, by = intersect(names(adsl), names(adae)))
  df
  list(df = df, adsl = adsl)
}

create_table_22 <- function(df, adsl) {
  tbl <- make_table_22(df = df, alt_counts_df = adsl, denom = "N_s")
  ft <- tt_to_flextable(tbl)
  df_tbl <- ft$body$dataset
  
  df_tbl <- df_tbl %>%
    mutate(V1 = ifelse(V1 %in% c("F", "M", ">=17 to <65", ">=65", "ASIAN", "BLACK OR AFRICAN AMERICAN", "WHITE", "AMERICAN INDIAN OR ALASKA NATIVE", "MULTIPLE", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "OTHER", "UNKNOWN", "HISPANIC OR LATINO", "NOT HISPANIC OR LATINO", "NOT REPORTED", "UNKNOWN"), paste0("  ", V1), V1))
  
  df_tbl  # Return df_tbl
}

create_sassy_report <- function(df_tbl) {
  library(sassy)
  
  tbl <- create_table(df_tbl, first_row_blank = TRUE, borders = "bottom") %>%
    column_defaults(from = "V1", to = "V4", align = "left", width = 1.5) %>%
    spanning_header(2, 4, label = "基线情况") %>% 
    define(var = "V1", label = "访视 \n 治疗组", width = 2.5) %>%  # Set the width of the first column to 2.5
    define("V2", label = "正常 \nN = 134", align = "right") %>%
    define("V3", label = "异常无临床意义 \nN = 134", align = "right") %>%
    define("V4", label = "异常有临床意义 \nN = 132", align = "right") %>%
    titles("表格14.3.4.2 用药前后体格检查参数临床评估的交叉表 — 安全性分析集", borders = "bottom", bold = TRUE) %>%
    footnotes("数据来源：列表16.2.x", "注：百分比计算的分母基于安全性分析集的受试者人数。", "[1] 依从性 =实际用药量/计划用药量x100%。")
  
  rpt <- create_report("fda-table-22.rtf", 
                       paper_size = "RD4",
                       output_type = "RTF", 
                       font = "Times", font_size = 10) %>%
    page_header(width = 2,
                right = c("[Mock-up TLF shells (CN)", "Statistical Analysis Tables and Figures, List (Chinese Version)]", "[TP-HLT-BS-004,V1.0,15Mar2024]"), 
                blank_row = "below") %>%
    set_margins(top = 1.26, bottom = 1.26, left = 1, right = 1) %>%
    add_content(tbl, align = "centre") %>%
    page_footer(left = c("Associated Process:", "Associated Process: SOP-HLT-BS-001", "Confidentiality保密"), right = "Page [pg] of [tpg]")
  
  res <- write_report(rpt)
  
  res$modified_path
}

library(dplyr)
library(falcon)
library(random.cdisc.data)
library(formatters)
library(sassy)
library(stringr)

# Load the libraries and data
data <- load_libraries_and_data()
print(data) 
# Create the table
df_tbl <- create_table_22(data$df, data$adsl)
class(df_tbl)
print(df_tbl)  


# Create the report
report_path <- create_sassy_report(df_tbl)  # Pass df_tbl as an argument

# Uncomment to view the report
file.show(report_path)