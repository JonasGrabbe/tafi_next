#' Data loader / preloader for Archology taphonomy datasets
#'
#' Reads a site-by-individual (or site aggregate) table and validates it against
#' a schema described in `meta.xlsx` (field names, types, and constraints).
#'
#' The loader is designed to be *warning-first*: it will try to coerce columns
#' into the expected types and report problems as warnings rather than stopping.
#'
#' @param data Either a data.frame / tibble OR a path to a data file
#'   (csv/tsv/xlsx).
#' @param meta_path Path to the schema Excel file (default: `inst/extdata/meta.xlsx`).
#' @param meta_sheet Excel sheet name that contains the schema table.
#' @param coerce Logical; if TRUE (default), try to coerce columns into the expected types.
#' @param na_to_zero_skeletal Logical; if TRUE (default), replace NA in skeletal count
#'   columns with 0 and warn if any NA were present.
#' @param check_enums Logical; if TRUE (default), validate categorical values when the
#'   schema provides a parseable enumeration.
#' @param check_cross_fields Logical; if TRUE (default), run cross-field consistency
#'   checks (e.g., start/end dates, karstic vs open-air).
#' @param treat_unknown_numeric_as_skeletal Logical; if TRUE (default), columns that
#'   are not in the schema and are numeric are validated as skeletal counts (>= 0 integers).
#'
#' @return A data.frame with cleaned column names and coerced columns.
#'   A validation report (list of issues) is attached as attribute `validation_report`.
#'
#' @export
archo_preload <- function(
  data,
  meta_path = system.file("extdata", "meta.xlsx", package = "archology"),
  meta_sheet = "Site Name",
  coerce = TRUE,
  na_to_zero_skeletal = TRUE,
  check_enums = TRUE,
  check_cross_fields = TRUE,
  treat_unknown_numeric_as_skeletal = TRUE
) {
  schema <- archo_read_schema(meta_path = meta_path, meta_sheet = meta_sheet)

  df <- if (is.character(data) && length(data) == 1L) {
    archo_read_data_file(data)
  } else {
    as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Clean names to be robust against trailing tabs/spaces in excel exports
  names(df) <- archo_clean_names(names(df))
  schema$field_name <- archo_clean_names(schema$field_name)

  # Detect duplicate columns after cleaning
  if (any(duplicated(names(df)))) {
    dups <- unique(names(df)[duplicated(names(df))])
    warning(
      paste0(
        "Duplicate columns detected after name cleaning: ",
        paste(dups, collapse = ", "),
        ". Consider renaming upstream to unique names."
      ),
      call. = FALSE
    )
  }

  out <- archo_validate_and_coerce(
    df = df,
    schema = schema,
    coerce = coerce,
    na_to_zero_skeletal = na_to_zero_skeletal,
    check_enums = check_enums,
    check_cross_fields = check_cross_fields,
    treat_unknown_numeric_as_skeletal = treat_unknown_numeric_as_skeletal
  )

  out$data
}

# --------- Schema reading ---------

#' @noRd
archo_read_schema <- function(meta_path, meta_sheet = "Site Name") {
  if (!file.exists(meta_path)) {
    stop("meta_path does not exist: ", meta_path, call. = FALSE)
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read meta.xlsx. Please install it.", call. = FALSE)
  }

  sch <- readxl::read_excel(meta_path, sheet = meta_sheet)
  sch <- as.data.frame(sch, stringsAsFactors = FALSE)

  # expected columns in the schema spreadsheet
  # Field Name | Description | Data Type | Allowed Values/Constraints | Example | Explanation
  names(sch) <- archo_clean_headers(names(sch))

  # Normalize
  if (!"field_name" %in% names(sch)) {
    stop("Schema sheet must contain a 'Field Name' column.", call. = FALSE)
  }

  # Drop empty rows
  sch$field_name <- archo_clean_names(sch$field_name)
  sch <- sch[!is.na(sch$field_name) & nzchar(sch$field_name), , drop = FALSE]

  # Normalize Data Type
  if ("data_type" %in% names(sch)) {
    sch$data_type <- archo_clean_dtype(sch$data_type)
  } else {
    sch$data_type <- NA_character_
  }

  # Helper: classify skeletal fields (in your schema these rows have the NA=0 explanation)
  sch$is_skeletal <- FALSE
  if ("explanation" %in% names(sch)) {
    sch$is_skeletal <- !is.na(sch$explanation) & grepl("NA\\s*=\\s*0", sch$explanation, ignore.case = TRUE)
  }

  sch
}

# --------- Data reading ---------

#' @noRd
archo_read_data_file <- function(path) {
  if (!file.exists(path)) {
    stop("Data file does not exist: ", path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel input data.", call. = FALSE)
    }
    df <- readxl::read_excel(path)
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    return(df)
  }

  if (ext %in% c("csv")) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    return(df)
  }

  if (ext %in% c("tsv", "txt")) {
    df <- utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    return(df)
  }

  stop("Unsupported file extension: .", ext, " (supported: csv/tsv/txt/xlsx).", call. = FALSE)
}

# --------- Validation ---------

#' @noRd
archo_validate_and_coerce <- function(
  df,
  schema,
  coerce = TRUE,
  na_to_zero_skeletal = TRUE,
  check_enums = TRUE,
  check_cross_fields = TRUE,
  treat_unknown_numeric_as_skeletal = TRUE
) {
  issues <- list()
  add_issue <- function(type, msg) {
    issues[[length(issues) + 1L]] <<- paste0("[", type, "] ", msg)
    invisible(NULL)
  }

  expected <- unique(schema$field_name)
  present <- names(df)

  missing_expected <- setdiff(expected, present)
  if (length(missing_expected) > 0) {
    add_issue(
      "missing_columns",
      paste0("Missing schema columns: ", paste(missing_expected, collapse = ", "))
    )
  }

  extra_cols <- setdiff(present, expected)

  # identify schema-defined skeletal fields (first few were uploaded, but constraints are the same for all)
  skeletal_fields_schema <- schema$field_name[schema$is_skeletal]
  meta_fields_schema <- setdiff(expected, skeletal_fields_schema)

  # Validate schema-defined fields
  for (i in seq_len(nrow(schema))) {
    field <- schema$field_name[i]
    dtype <- schema$data_type[i]
    constraints <- if ("allowed_values_constraints" %in% names(schema)) schema$allowed_values_constraints[i] else NA_character_
    is_skeletal <- isTRUE(schema$is_skeletal[i])

    if (!field %in% names(df)) next

    x <- df[[field]]

    # Coercion + type checks
    if (dtype == "integer") {
      conv <- archo_as_integerish(x, field = field, coerce = coerce)
      df[[field]] <- conv$x
      if (length(conv$warn) > 0) add_issue("type", conv$warn)

      # basic non-negative constraint
      neg <- which(!is.na(df[[field]]) & df[[field]] < 0L)
      if (length(neg) > 0) {
        add_issue("constraint", paste0(field, ": has ", length(neg), " negative values (must be >= 0)."))
      }

      if (is_skeletal && na_to_zero_skeletal) {
        na_n <- sum(is.na(df[[field]]))
        if (na_n > 0) {
          df[[field]][is.na(df[[field]])] <- 0L
          add_issue("skeletal_na", paste0(field, ": replaced ", na_n, " NA with 0 (NA treated as not identified)."))
        }
      }
    } else if (dtype == "boolean") {
      conv <- archo_as_logicalish(x, field = field, coerce = coerce)
      df[[field]] <- conv$x
      if (length(conv$warn) > 0) add_issue("type", conv$warn)
    } else if (dtype == "string") {
      if (!is.character(x) && coerce) {
        df[[field]] <- as.character(x)
        add_issue("type", paste0(field, ": coerced to character."))
      }
    } else if (dtype == "categorical") {
      if (!is.character(x) && !is.factor(x) && coerce) {
        df[[field]] <- as.character(x)
        add_issue("type", paste0(field, ": coerced to character (categorical)."))
      }

      if (check_enums) {
        allowed <- archo_parse_enum(constraints)
        if (!is.null(allowed) && length(allowed) > 0) {
          bad <- which(!is.na(df[[field]]) & !df[[field]] %in% allowed)
          if (length(bad) > 0) {
            # show a few distinct bad values
            bad_vals <- unique(df[[field]][bad])
            bad_vals <- bad_vals[seq_len(min(length(bad_vals), 8))]
            add_issue(
              "constraint",
              paste0(field, ": has values outside allowed set: ", paste(bad_vals, collapse = ", "))
            )
          }
        }
      }
    } else {
      # Unknown dtype in schema; no validation performed
      if (!is.na(dtype) && nzchar(dtype)) {
        add_issue("schema", paste0(field, ": unsupported data type in schema: ", dtype))
      }
    }
  }

  # Validate extra numeric columns as skeletal counts (useful because only the first few skeletal fields were included in meta.xlsx)
  if (treat_unknown_numeric_as_skeletal && length(extra_cols) > 0) {
    for (field in extra_cols) {
      x <- df[[field]]
      if (is.numeric(x) || is.integer(x)) {
        conv <- archo_as_integerish(x, field = field, coerce = coerce)
        df[[field]] <- conv$x
        if (length(conv$warn) > 0) add_issue("type", paste0(field, ": ", conv$warn, " (treated as skeletal)."))

        neg <- which(!is.na(df[[field]]) & df[[field]] < 0L)
        if (length(neg) > 0) {
          add_issue("constraint", paste0(field, ": has ", length(neg), " negative values (must be >= 0)."))
        }

        if (na_to_zero_skeletal) {
          na_n <- sum(is.na(df[[field]]))
          if (na_n > 0) {
            df[[field]][is.na(df[[field]])] <- 0L
            add_issue("skeletal_na", paste0(field, ": replaced ", na_n, " NA with 0 (treated as skeletal count)."))
          }
        }
      }
    }
  }

  # Cross-field checks that matter for your comparisons
  if (check_cross_fields) {
    # start/end BP: start >= end (older BP number should be greater)
    if ("start_date_cal_bp" %in% names(df) && "end_date_cal_bp" %in% names(df)) {
      s <- df[["start_date_cal_bp"]]
      e <- df[["end_date_cal_bp"]]
      bad <- which(!is.na(s) & !is.na(e) & s < e)
      if (length(bad) > 0) {
        add_issue("cross_field", paste0("start_date_cal_bp < end_date_cal_bp in ", length(bad), " rows (expected start >= end)."))
      }
    }

    # karstic vs open-air mutual exclusivity (NA where not applicable)
    if ("karstic_system" %in% names(df) && "open_air_site" %in% names(df)) {
      k <- df[["karstic_system"]]
      o <- df[["open_air_site"]]
      bad <- which(!is.na(k) & nzchar(k) & !is.na(o) & nzchar(o))
      if (length(bad) > 0) {
        add_issue("cross_field", paste0("Both karstic_system and open_air_site are filled in ", length(bad), " rows; expected one to be NA."))
      }
    }
  }

  # Emit one aggregated warning (much friendlier than 200 warnings)
  if (length(issues) > 0) {
    warning(
      paste0(
        "Archology data validation produced ", length(issues), " issue(s):\n",
        paste0(" - ", unique(issues), collapse = "\n")
      ),
      call. = FALSE
    )
  }

  attr(df, "validation_report") <- unique(issues)
  list(data = df, issues = unique(issues))
}

# --------- Helpers ---------

#' @noRd
# Clean field names (data column names / schema field names) while preserving bone names with spaces
archo_clean_names <- function(x) {
  x <- as.character(x)
  x <- gsub("[\\t\\r\\n]+", " ", x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x
}

#' @noRd
# Clean *header* names (schema column headers like "Field Name") to lower_snake_case
archo_clean_headers <- function(x) {
  x <- as.character(x)
  x <- gsub("[\\t\\r\\n]+", " ", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

#' @noRd
archo_clean_dtype <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("int", "integer", "whole number")] <- "integer"
  x[x %in% c("logical", "bool", "boolean")] <- "boolean"
  x[x %in% c("cat", "categorical", "factor")] <- "categorical"
  x[x %in% c("str", "string", "character", "text")] <- "string"
  x
}

#' @noRd
archo_as_integerish <- function(x, field, coerce = TRUE) {
  warn <- character(0)

  # allow blanks
  if (all(is.na(x))) {
    return(list(x = as.integer(x), warn = warn))
  }

  if (is.integer(x)) {
    return(list(x = x, warn = warn))
  }

  if (is.logical(x)) {
    # excel sometimes stores 0/1 as logical TRUE/FALSE
    if (coerce) {
      warn <- paste0(field, ": logical values coerced to integer (TRUE->1, FALSE->0).")
      return(list(x = as.integer(x), warn = warn))
    }
    warn <- paste0(field, ": expected integer but found logical.")
    return(list(x = x, warn = warn))
  }

  if (is.character(x)) {
    # try numeric parse
    if (!coerce) {
      warn <- paste0(field, ": expected integer but found character.")
      return(list(x = x, warn = warn))
    }
    suppressWarnings(num <- as.numeric(x))
    if (all(is.na(num) & !is.na(x))) {
      warn <- paste0(field, ": could not parse character values as numeric.")
      return(list(x = x, warn = warn))
    }
    x <- num
  }

  if (!is.numeric(x)) {
    warn <- paste0(field, ": expected integer but found type ", class(x)[1], ".")
    return(list(x = x, warn = warn))
  }

  # numeric -> integer check
  frac <- x - floor(x)
  bad_frac <- which(!is.na(frac) & abs(frac) > 1e-8)
  if (length(bad_frac) > 0) {
    if (coerce) {
      warn <- paste0(field, ": has non-integer numeric values; rounded to nearest integer.")
      x <- round(x)
    } else {
      warn <- paste0(field, ": has non-integer numeric values.")
      return(list(x = x, warn = warn))
    }
  }

  list(x = as.integer(x), warn = warn)
}

#' @noRd
archo_as_logicalish <- function(x, field, coerce = TRUE) {
  warn <- character(0)

  if (is.logical(x)) {
    return(list(x = x, warn = warn))
  }

  if (!coerce) {
    warn <- paste0(field, ": expected boolean but found type ", class(x)[1], ".")
    return(list(x = x, warn = warn))
  }

  if (is.numeric(x) || is.integer(x)) {
    bad <- which(!is.na(x) & !x %in% c(0, 1))
    if (length(bad) > 0) {
      warn <- paste0(field, ": numeric boolean has values other than 0/1; coerced non-zero to TRUE.")
    }
    return(list(x = as.logical(x), warn = warn))
  }

  if (is.character(x)) {
    z <- tolower(trimws(x))
    out <- rep(NA, length(z))
    out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[z %in% c("false", "f", "0", "no", "n")] <- FALSE
    if (any(is.na(out) & !is.na(z) & nzchar(z))) {
      warn <- paste0(field, ": some values could not be interpreted as TRUE/FALSE.")
    }
    return(list(x = out, warn = warn))
  }

  warn <- paste0(field, ": could not coerce to logical from type ", class(x)[1], ".")
  list(x = x, warn = warn)
}

#' @noRd
archo_parse_enum <- function(constraints) {
  if (is.null(constraints) || is.na(constraints) || !nzchar(constraints)) return(NULL)
  s <- as.character(constraints)

  # If numeric constraints are present, do not parse as enum
  if (grepl("[≥<=]", s) || (grepl("\\bNA\\b", s, ignore.case = TRUE) && grepl("[0-9]", s))) {
    return(NULL)
  }

  # remove parenthetical notes like "NA if unknown"
  s <- gsub("\\(.*?\\)", "", s)
  s <- gsub("[\"“”]", "", s)
  s <- gsub("…", "", s)
  parts <- unlist(strsplit(s, ","))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts <- parts[!grepl("^NA\\b", parts, ignore.case = TRUE)]
  # If only one part, it's probably not a real enum
  if (length(parts) <= 1) return(NULL)
  parts
}
