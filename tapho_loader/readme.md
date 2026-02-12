
# tapho_loader_new.R

1. summarize_bone_sections_zero_na
2. recal_mne_six
3. align_counts_to_skeleton
4. tapho.load
5. tapho.load.db
6. read_excel_sheets



# `tapho_loader_new.R`

## `summarize_bone_sections_zero_na(input_df, reference_df = elements_data)`

**Inputs**

- `input_df`: data.frame/tibble with columns **`element`** and **`count`** (skeletal rows only)
- `reference_df`: defaults to `elements_data`; must support `add_bone_info()` (i.e., provides element → `max_Level` and `Skeleton_Section` mapping)

**Output**

- **unnamed** numeric vector length **11**, ordered as:
  `(cranium, teeth, hyoid, vertebrae, thorax, shoulder, arms, hands, pelvic, legs, foot)`

**Dependencies**

- `add_bone_info()` (from `data_loader_new.R`)
- `elements_data` (from `meta_new.R`)

**What it does**

- Drops rows where `count` is `NA` or `0` (so “no observed bone” doesn’t drive levels).
- Adds bone metadata (`max_Level`, `Skeleton_Section`) using `add_bone_info()`.
- For each skeleton section, takes the **minimum `max_Level`** observed.
- Initializes the result with a hard-coded “most detailed” fallback vector:
  `c(5, 5, 1, 3, 4, 3, 3, 5, 4, 3, 5)`
  and only replaces sections that are present in the input.
- Returns the vector **without names** (`unname(...)`).

------

## `recal_mne_six(current_mne, target_levels, skeleton_hierarchy, bone_names, sections_order)`

**Inputs**

- `current_mne`: numeric vector of counts (same order as `bone_names`)
- `target_levels`: numeric vector length **11** (target levels in the fixed section order)
- `skeleton_hierarchy`: *passed in but not used in the current implementation*
- `bone_names`: character vector of element names (same length/order as `current_mne`)
- `sections_order`: character vector of section names (length 11), matching `target_levels`

**Output**

- tibble with columns:
  - `element`: standardized element name (**spaces replaced with dots**)
  - `count`: recalculated count at the target level

**Dependencies**

- `add_bone_info()` (adds `Skeleton_Section`)

- `category_hierarchy` (must exist in scope; typically from `meta_new.R`)

  

**What it does**

- Treats the input counts as being at fixed “current levels”:
  `c(5, 5, 1, 3, 4, 3, 3, 5, 4, 3, 5)` (aligned with `sections_order`).
- For each skeleton section:
  - pulls only rows in that `Skeleton_Section`
  - iteratively “rolls up” counts from `current_level` down to `target_level` using mappings in `category_hierarchy`
  - mapping name pattern: `"{part}_{from}_to_{to}"` (e.g., `arms_3_to_2`)
  - elements that map are aggregated to their parent (sum)
  - non-mappable elements are carried forward
  - if a mapping is missing, warns and stops rolling-up that section
- Returns a tibble of the rolled-up element counts (element names are the standardized `std_element`, not prefixed by section).

------

## `align_counts_to_skeleton(one_skeleton_element, recalc_element, recalc_count, fill_value = NA, aggregate = c("sum", "first"))`

**Inputs**

- `one_skeleton_element`: character vector of the **reference** element order (the “one skeleton” inventory)
- `recalc_element`: character vector of element names from recalculated data
- `recalc_count`: numeric vector of counts aligned to `recalc_element`
- `fill_value`: value used when an element in `one_skeleton_element` is missing from `recalc_element` (default `NA`)
- `aggregate`: how to handle duplicates in `recalc_element`
  - `"sum"`: sum duplicate element counts
  - `"first"`: keep the first occurrence only

**Output**

- data.frame with columns:
  - `element` (from `one_skeleton_element`)
  - `count` (aligned and filled)

**What it does**

- Builds a lookup table from `recalc_element`/`recalc_count`, resolving duplicates by `aggregate`.

- Matches counts into the **exact order** of `one_skeleton_element`.

- Fills missing matches with `fill_value`.

- Returns a 2-column aligned data.frame ready for `loaded_tapho()`.

  

------

## `tapho.load(new_site_data, strict = FALSE)`

**Inputs**

- `new_site_data`: site table data.frame in the Archo format:
  - column 1 = `element`, column 2 = counts (column name becomes `ind_name`)
  - rows **1:31** meta rows; rows **32+** skeletal rows
- `strict`: logical; forwarded to validation and final packaging

**Output**

- a single-row tibble produced by `loaded_tapho()` with:
  - parsed meta fields as columns
  - `sum_count`
  - list-cols: `bone.level`, `element`, `count`

**Dependencies**

- `archo_preload_site_table()` (from `archo_preload_2.R`)
- `summarize_bone_sections_zero_na()`
- `recal_mne_six()`
- `get_bones_dataframe()` (from `data_loader_new.R`)
- `loaded_tapho()` (from `archo_preload_2.R`)
- Requires in-scope objects:
  - `skeleton_hierarchy`
  - `sections_order`

**What it does**

1. Validates/coerces `new_site_data` with `archo_preload_site_table()` (NA/invalid skeletal counts become 0, etc.).
2. Extracts `ind_name` from the second column name, then renames columns to `element/count`.
3. Splits into:
   - `new_site_meta` = rows 1:31
   - `new_site_raw` = rows 32+
4. Computes `target_levels` from observed skeletal data (`summarize_bone_sections_zero_na()`).
5. Recalculates skeletal counts down to `target_levels` with `recal_mne_six()`.
6. Builds the “one skeleton at those levels” element list (`get_bones_dataframe(target_levels)`).
7. Aligns recalculated counts to that skeleton inventory (`align_counts_to_skeleton()`).
8. Packages everything into the standard dataset row via `loaded_tapho(...)`.

------

## `tapho.load.db(data)`

**Inputs**

- `data`: data.frame/tibble where:
  - column 1 = element names
  - column 2 = first individual/site counts
  - columns 3+ (optional) = additional individuals’ counts

**Output**

- tibble with **one or more rows** (each row = result of `tapho.load()` on one count column)

**Dependencies**

- `tapho.load()`

**What it does**

- Loads the first two columns using `tapho.load(data[, 1:2])`.
- If the returned record has `count_type == FALSE`, it assumes columns 3+ are additional individuals and:
  - iterates `i = 3:ncol(data)`
  - loads `data[, c(1, i)]` with `tapho.load()`
  - appends each result as a new row via `add_row(...)`
- Returns the combined tibble.

------

## `read_excel_sheets(file_path)`

**Inputs**

- `file_path`: path to an Excel file (`.xlsx`, etc.) with **one or more sheets** in the same layout expected by `tapho.load.db()`

**Output**

- tibble combining all processed sheets (row-bind of each sheet’s `tapho.load.db()` output)

**Dependencies**

- `readxl::excel_sheets()`
- `readxl::read_excel()`
- `tapho.load.db()`
- `dplyr::bind_rows()`

**What it does**

- Lists all sheet names and prints them.
- For each sheet:
  - reads it into a data.frame
  - processes it with `tapho.load.db()`
- Binds all processed sheets together with `bind_rows()` and returns the combined dataset.

# select_data_new.R

1. process_data
2. select_data
3. get_min
4. recalculate_counts
5. group_and_summarize
6. sumby_site



## `process_data(data)`

**Inputs**

- `data`: dataset tibble with `sum` column (0=individual, 1=site-aggregated) plus the standard list-columns.

**Output**

- tibble where:
  - all `sum==0` rows are first harmonized, then collapsed to site-level (`sum=1`)
  - merged with the original `sum==1` rows

**Dependencies**

- `recal_mne_selected_data()` (from `recalMNE.R`)

**What it does**

- If there are individual rows (`sum==0`):

  1. harmonize them to a common bone.level
  2. group by `site_index`
  3. sum `number_of_individuals`
  4. mark metadata columns as `"mixed"` when not consistent
  5. sum the `count` vectors **by position** (pads shorter vectors with zeros)
  6. strips the `"Section_"` prefix from element names (everything before `_`)

- Then binds with any pre-aggregated rows.

  

------

## `select_data(data, ...filters..., sum_filter=1, hominin_species_filter=NULL)`

**Inputs**

- `data`: dataset tibble
- filters: vectors to keep matching rows (time_range, context, age, gender, site_name, etc.)
- `bone_parts_filter`: keeps rows where any `element` matches the provided list
- `sum_filter`: if 0 → keep only individual rows; else → runs `process_data()`

**Output**

- filtered dataset tibble

**Dependencies**

-  `process_data()`

  

**What it does**

- Applies metadata filters.

  

------

## `get_min(selected_data)`

**Inputs**

- `selected_data`: tibble with list-column `bone.level`

**Output**

- unnamed numeric vector: minimum bone level per bone section (across all rows)
- 

**What it does**

- Expands `bone.level` into columns, takes per-column minima, returns numeric vector.

------

## `recalculate_counts(data_filtered)`

**Inputs**

- `data_filtered`: dataset tibble in internal format

**Output**

- tibble where list-columns `element/count` are converted into **wide columns**, one per element, with 0 fill.

**Dependencies**

- `get_min()`

- `recal_mne()` (from `recalMNE_new.R`)

- `skeleton_hierarchy`, `sections_order` 

  

**What it does**

1. finds minimum bone level across the dataset
2. recalculates each row to that minimum using `recal_mne()`
3. builds a consistent set of all elements across rows
4. expands into a matrix-like wide dataframe (elements become columns)

------

## `group_and_summarize(data, group_cols)`

**Inputs**

- `data`: dataset tibble
- `group_cols`: character vector of column names to group by

**Output**

- grouped/summarized tibble



**What it does**

- Groups rows and merges metadata (marks “mixed” when inconsistent).

- Sums `number_of_individuals`.

- Keeps first `bone.level`, first `element`.

- Sums `count` list-columns element-wise

  

----

## `sumby_site(data)`

**Inputs**

- `data`: dataset tibble with (at least):
  - `site_id`, `ind_name`, `individual_id`, `site_name`, `layer_name`
  - demographic columns: `mni`, `age_category`, `gender`, `n_infant`, `n_juvenil`, `n_subadult`, `n_adult`
  - `observation`, `reference`, `sum_count`
  - list-cols: `bone.level`, `element`, `count`

**Output**

- tibble with **one row per `site_id`**, aggregating individuals into a single “site-level” record:
  - `ind_name` and `individual_id` concatenated across individuals
  - `layer_name` kept if uniform; otherwise a comma-separated list of unique layers
  - `mni`, `sum_count`, and the `n_*` demographics summed
  - `age_category` / `gender` merged with special logic (see below)
  - `observation` concatenated into a single string
  - `count` summed into one vector
  - `bone.level` / `element` taken from the first record (after harmonization)

**Dependencies**

- `recal_mne_selected_data()` (from `recalMNE_new.R`)

**What it does**

1. Aggregates MNE data into one aggregated max shared level 



----

---

---



recalcMNE_new.R

1) recal_mne
2) recal_mne_selected_data



------

# `recalMNE_new.R` (harmonize/aggregate counts to a target level)

## `recal_mne(current_mne, current_levels, target_levels, skeleton_hierarchy, bone_names, sections_order)`

**Inputs**

- `current_mne`: numeric vector of counts
- `current_levels`: numeric vector length 11 (current bone.level)
- `target_levels`: numeric vector length 11 (desired bone.level)
- `skeleton_hierarchy`
- `bone_names`: character vector aligned to `current_mne`
- `sections_order`: character vector of section names in the same order as levels

**Output**

- tibble with columns:
  - `element`: names like `"Skull_frontal.bone"` or `"Skull_neurocranium_bone"`
  - `count`: aggregated counts

**Dependencies**

- `add_bone_info()` (from `data_loader_new.R`)
- `category_hierarchy` (from `meta_new.R`)

**What it does**

- Adds section/level metadata to each element using `add_bone_info()`.
- For each skeletal section:
  - if current level == target: copies counts through
  - else: looks up a mapping like `"Skull_5_to_3"` in `category_hierarchy`,
    maps fine elements → coarse categories, then **sums counts into target categories**

------

## `recal_mne_selected_data(selected_data)`

**Inputs**

- `selected_data`: tibble with list-columns `count`, `bone.level`, `element`

**Output**

- same tibble but with:
  - `bone.level` replaced by the minimum levels across all rows
  - `element` / `count` replaced by recalculated lists at that min level

**Dependencies**

- `get_min()` from `select_data_new.R`
- `recal_mne()`
- Needs `sections_order`, `skeleton_hierarchy` 

**What it does**

- Computes the dataset-wide minimum bone.level vector.

- Re-aggregates every row to match that minimum.

  

----

---



# data_loader_new.R

1. add_bone_info 

2. summarize_bone_sections

3. get_bones_dataframe

4. compareBoneLevels

   

# `data_loader_new.R`

## `add_bone_info(input_df, reference_df = elements_data)`

**Inputs**

- `input_df`: dataframe with columns **`element`** and **`count`**
- `reference_df`: defaults to `elements_data` and must have `ElementsOne`, `max_Level`, `Skeleton_Section`

**Output**

- dataframe with columns: `element`, `count`, `max_Level`, `Skeleton_Section`

**Dependencies**

- `elements_data` from `meta_new.R`

**What it does**

- Adds two columns  `max_Level` of bone element and  `Skeleton_Section` (cranium, teeth, hyoid, vert...)

------

## `summarize_bone_sections(input_df, reference_df = elements_data)`

**Inputs**

- `input_df`: dataframe with `element`, `count`
- `reference_df`: `elements_data`

**Output**

- **unnamed** numeric vector length 10: (Skull, Teeth, Hyoid, Vertebrae, Thorax, Shoulder, Arms, Hands, Pelvic, Legs, Foot)

**Dependencies**

- `add_bone_info()`

**What it does**

- Adds bone info then takes the **minimum `max_Level` per `Skeleton_Section`**.
- Missing sections default to **1**.
- Returns vector of those minima in a fixed order.

------

## `get_bones_dataframe(levels_vector)`

**Inputs**

- `levels_vector`: numeric vector length **11**, ordered to match sections (Skull..Foot)

**Output**

- dataframe with columns: `element`, `count` (expected per one skeleton at those levels)

**Dependencies**

- `skeleton_hierarchy` from `meta_new.R`

**What it does**

- For each section, looks up `skeleton_hierarchy[[section]][[paste0("Level ", level)]]`.
- Unlists nested lists and returns the flattened element inventory.

------

## `compareBoneLevels(...)`  *(duplicate name, different expectation!)*

**Inputs**

- `...`: 2+ lists, each expected to have **`$bone.level`** (not `$bone_level`)

**Output**

- numeric vector: element-wise minimum.

**Dependencies**

- Base R `Reduce(pmin, ...)`

**What it does**

-  Calculates shared maximum level.





# archo_preload_2.R

1. archo_preload_site_table
2. archo_preload
3. tapho.load_2
4. loaded_tapho



---





# `archo_preload_2.R`

## `archo_preload_site_table(x, strict = FALSE)`

**Inputs**

- `x`: data.frame with **at least 2 columns** and **at least 32 rows**, expected layout:
  - Column 1 = element names (character)
  - Column 2 = one site/individual column (counts; column name is treated as `ind_name`)
  - Rows **1:31** = meta rows (key/value stored as `element` / `ind_name`)
  - Rows **32+** = skeletal element rows (counts)
- `strict`: logical
  - `FALSE` (default) → collect issues and emit a single warning report
  - `TRUE` → stop immediately on the first validation issue

**Output**

- data.frame with **2 columns**:

  

**What it does**

- Enforces the “31 meta rows + skeletal rows” table shape and keeps **only the first two columns**.
- Normalizes:
  - trims whitespace in `element`
  - preserves/repairs the 2nd column name as the individual name
- Validates meta rows contain required keys:
  - `site_id`, `individual_id`, `site_name`, `count_type`, `geological_period`, `cultural_period`, `reference`
- Validates / coerces skeletal counts (rows 32+):
  - non-numeric values → coerced to `NA`, then set to **0**
  - `NA` → replaced with **0** 
  - non-integers → **rounded** to nearest integer
  - negatives → set to **0**

------

## `archo_preload(x, strict = FALSE)`

**Inputs**

- `x`: expected to be a two-column (or more) site table data.frame in the format required by `archo_preload_site_table()`
- `strict`: logical, forwarded to `archo_preload_site_table()`

**Output**

- Validated/preloaded site table returned by `archo_preload_site_table()` (2-column data.frame with meta rows + skeletal rows)

**Dependencies**

- `archo_preload_site_table()`

**What it does**

- Possible wrapper for `archo_preload_site_table()`

------

## `tapho.load_2(new_site_data, strict = FALSE)`

**Inputs**

- `new_site_data`: a “site table” data.frame in the same expected format as `archo_preload_site_table()` (meta rows 1:31 + skeletal rows)

**Output**

- A tibble (from `loaded_tapho()`) with:
  - `ind_name`
  - all parsed meta fields (`site_id` … `reference`)
  - `sum_count`
  - list-cols: `bone.level`, `element`, `count`

**Dependencies**

- Internal to this file:
  - `archo_preload()`
  - `loaded_tapho()`
- Must already exist in the calling environment (not defined in this file):
  - `summarize_bone_sections_zero_na()`
  - `recal_mne_six()`
  - `skeleton_hierarchy`
  - `sections_order`
  - `get_bones_dataframe()`
  - `align_loaded_data()`

**What it does**

- Validates `new_site_data` as a properly-formed 2-column site table via `archo_preload()`.
- Extracts:
  - `ind_name` from the original second column name
  - `new_site_meta` = rows 1:31
  - `new_site_raw` = rows 32+ (skeletal elements)
- Computes target bone preservation/levels from the raw skeletal section data:
  - `target_levels <- summarize_bone_sections_zero_na(new_site_raw)`
- Ensures skeletal counts are usable (NA → 0) and recalculates/derives counts using:
  - `recal_mne_six(counts, target_levels, skeleton_hierarchy, elements, sections_order)`
- Builds the “expected one skeleton” element inventory at those levels (`get_bones_dataframe(target_levels)`)
- Aligns recalculated counts to the expected element list (`align_loaded_data(...)`)
- Packages everything into a validated output record using `loaded_tapho(...)`



------

## `loaded_tapho(ind_name, meta, level, col_element, col_count, strict = FALSE)`

**Inputs**

- `ind_name`: character scalar (identifier for the individual / column name)
- `meta`: data.frame with columns **`element`** and **`count`**
  - expected to represent the **31 meta rows** as key/value pairs
- `level`: numeric vector, expected length **11** (bone section levels)
- `col_element`: character vector of skeletal element names
- `col_count`: numeric/integer vector of skeletal counts (same length as `col_element`)

**Output**

- A `tibble::tibble` with columns:

  - `ind_name`
  - meta fields (explicit columns):
    - `site_id`, `individual_id`, `site_name`, `layer_name`, `count_type`,
      `geological_period`, `cultural_period`, `culture_period_2`,
      `start_date_cal_bp`, `end_date_cal_bp`, `mis_stage`,
      `karstic_system`, `karstic_system_2`, `open_air_site`, `open_air_site_2`,
      `region`, `country`,
      `mni`, `age_category`, `n_infant`, `n_juvenil`, `n_subadult`, `n_adult`,
      `gender`,
      `taphonomic_context`, `funerary_context`, `discovery_context`, `degree_of_excavation`,
      `hominin_species`, `observation`, `reference`
  - `sum_count` (sum of skeletal counts)
  - list-cols:
    - `bone.level` (list containing the `level` vector)
    - `element` (list containing `col_element`)
    - `count` (list containing validated integer counts)

  

**What it does**

- Validates / standardizes `level`:
  - expects length 11
  - in non-strict mode: pads/truncates and records a warning
  - in strict mode: stops on mismatch
- Validates skeletal vectors:
  - `col_element` and `col_count` must be same length (otherwise stop)
  - trims whitespace in `col_element`
  - coerces `col_count` to numeric → `NA` becomes 0
  - rounds non-integers, clamps negatives to 0, returns **integer** counts
- Parses the 31-row meta block into named fields using hard-coded coercers:
  - `site_id`, `start_date_cal_bp`, etc. coerced to integer (with warnings on bad values)
  - text fields trimmed / blank → `NA`
  - `count_type` parsed as boolean-ish (`true/false`, `1/0`, `yes/no`, etc.)
- Applies “soft” meta rules (warn/stop depending on `strict`):
  - BP date consistency: `start_date_cal_bp` should be **>=** `end_date_cal_bp`
  - mutual exclusivity: `karstic_system` and `open_air_site` shouldn’t both be filled
  - allowed-value checks for:
    - `geological_period`
    - `karstic_system`
    - `open_air_site`
    - `gender`
    - `age_category`
  - `individual_id` pattern check: expects `"1-0"` or `"1-0-3"` style
  - if both present: checks `site_id` matches the prefix of `individual_id`
- Computes `sum_count` from the validated skeletal counts.
- Returns a single tibble row with meta as columns and skeletal data as list-cols, plus a consolidated validation warning report in non-strict mode.



















