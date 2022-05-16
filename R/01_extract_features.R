##-----------------------------------------------------##
##          Part 1: Extracting features                ##
##-----------------------------------------------------##

# load library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, readxl, openxlsx, reshape2, glue, stargazer, corrr, 
               corrplot, zoo, lubridate, tidyverse)

# set the current working directory the same as shell
wd_path = Sys.getenv("PWD") 
if(getwd() != wd_path) setwd(wd_path)

# set reference path 
if(!file.exists(".here")) here::set_here(path= ".")

# load functions
source(here::here("R/utils.R"))

# set system time zone
Sys.setenv(TZ='Europe/Zurich')

# files
file3 <- here::here("data/Gesundheitsbericht_2021_BProt_EBJ.xlsx")
file4_head <- here::here("data/Job354_Behandlungen_20-4 bis 21-4_20220121.csv")
file4_3 <- here::here("data/job_noheader_cleaned_quote_bis_213.csv")
file4_4 <- here::here("data/job_noheader_cleaned_quote_bis_214.csv")
file5 <- here::here("data/Qualitas_SUS-675_LSG_EXPORT_Q3_2021_LDat 20-4 bis 21-3.csv")
file6 <- here::here("data/Qualitas_SUS-675_TBi_EXPORT_Q3_2021_TBI 20-4 bis 21-3.csv")
file7 <- here::here("data/Impfstoffe Schwein_MHA_25.03.22.xls")

# Headers: Gesundheitsbericht column names is declared in 2 rows. We cleaned it first on excel file.
gbr_raw <- read_excel(here::here("data/variable_names.xlsx"), 
                      sheet = "Gesundheitsbericht", 
                      n_max = 280)

##-----------------------------------------------------##
## 1) Job354_Behandlungen  : vaccine & pathogens count ##
##-----------------------------------------------------##

# read the head
bhd_head <- readr::read_delim(file4_head,
                              ";", escape_double = FALSE, locale = locale(encoding = "ISO-8859-1"), 
                              col_names = FALSE, n_max = 1, trim_ws = TRUE, show_col_types = FALSE)

# clean up column names
bhd_head_clean <- data.frame(
  ori = cleanup_colname(as.character(bhd_head), minlen=16)[[1]],
  new = cleanup_colname(as.character(bhd_head), minlen=16)[[2]]
)

# update betriebs id col name
bhd_head_clean[which(bhd_head_clean$new == "betriebstvdnummr"), "new"] <- "tvdnummer"

# read col spec
bhd_spec <- readr::spec_delim(file4_head, ";")

# update the spec
names(bhd_spec$cols) <- bhd_head_clean[[2]]

# select column that contains date
selected <- grep("dtm$", names(bhd_spec$cols), value = TRUE)
bhd_spec$cols[selected] <- list(col_date(format = "%Y-%m-%d"))

# id column as character
selected_id <- grep("id$|nummer$|nr$", names(bhd_spec$cols), value = TRUE)
bhd_spec$cols[selected_id] <- list(col_character())

# index as integer
selected_idx <- grep("^idx", names(bhd_spec$cols), value = TRUE)
bhd_spec$cols[selected_idx] <- list(col_integer())

# preview
bhd_spec

# select which column to ingest
f_delim <- function(x, pos) {x}

# ingest file
bhd <- readr::read_delim_chunked(
  file = file4_4, DataFrameCallback$new(f_delim), 
  delim = ";", 
  escape_double = FALSE, 
  locale = locale(encoding = "ISO-8859-1"), 
  col_names = bhd_head_clean$new, 
  col_types = bhd_spec, 
  trim_ws = TRUE,  
  chunk_size = 100) 

# replace character
ReplaceChar <- function(x) { stringr::str_replace_all(x, c(
  "Ã¼" = "ü", 
  "Ã¶" = "ö", 
  "Ã¤" = "ä", 
  "Â®" = "®",
  "Ã©" = "é")) }

# convert date, and years > 2021 to be recorded as 2020
bhd <- bhd %>%
  # add rowid
  mutate(id = row_number()) %>%
  # date format
  mutate(
    erstsbhndlngsdtm = as.POSIXct(erstsbhndlngsdtm, format = '%Y-%m-%d'),
    ltztsbhndlngsdtm = as.POSIXct(ltztsbhndlngsdtm, format = '%Y-%m-%d'),
    fleischfreigbdtm = as.POSIXct(fleischfreigbdtm, format = '%Y-%m-%d'),
    organfreigabedtm = as.POSIXct(organfreigabedtm, format = '%Y-%m-%d'),
    injstellfregbdtm = as.POSIXct(injstellfregbdtm, format = '%Y-%m-%d')
  ) %>%
  mutate(
    jahr_ersdtm = lubridate::year(erstsbhndlngsdtm),
    jahr_ltzdtm = lubridate::year(ltztsbhndlngsdtm),
    jahr_flfdtm = lubridate::year(fleischfreigbdtm),
    jahr_orfdtm = lubridate::year(organfreigabedtm),
    jahr_injdtm = lubridate::year(injstellfregbdtm)
  ) %>%
  mutate(
    jahr_ltzdtm = ifelse(jahr_ltzdtm > 2021, 2020, jahr_ltzdtm),
    jahr_flfdtm = ifelse(jahr_ltzdtm > 2021, 2020, jahr_flfdtm),
    jahr_orfdtm = ifelse(jahr_ltzdtm > 2021, 2020, jahr_orfdtm),
    jahr_injdtm = ifelse(jahr_ltzdtm > 2021, 2020, jahr_injdtm),
  ) %>%
  mutate_at(vars(tiergruppe, identifikation, bhndlngsgrndkrnk, herkunftarznmttl), ReplaceChar)

# all medicine names
medicine_name <- bhd %>% 
  select(handlsnmarznmttl) %>%
  drop_na() %>%
  unique() %>%
  pull


#### Com & Sp Impfstoff: Using reference to calculate vaccine ######
# use table ref as of email 25.3.2022
# 3. sheet: merge of commercial and specific vaccines.
# With the 3rd sheet you can create one variable for the total number of pathogens covered, 
# as we discussed it by telephone. As you can see we might lose some farms (those who use the vaccines 
# in grey) and we lose some detail information (e.g. different types of E.coli or C. perfringens needed 
# to be summarized due to missing information).
# And let me mention once again – because we did not discuss it during our telephone call – 
# that we don’t need the vaccination information by animal category but in general. (I expect that 
# this will facilitate the interpretation.😉)

# vaccine
impfstoff_ref <- readxl::read_excel(file7, sheet = 'Comm. & spec. vaccines')
impfstoff_refsp <- readxl::read_excel(file7, sheet = 'Specific vaccines')

# cleaup col names
names(impfstoff_ref) <- cleanup_colname(names(impfstoff_ref), minlen = 12)[[2]]
names(impfstoff_refsp) <- cleanup_colname(names(impfstoff_refsp), minlen = 12)[[2]]

# marker for specific impfstoff
marker_spvac <- impfstoff_refsp  %>% 
  select(vaccinename) %>%
  mutate( vactype = "spec")

ReplaceX <- function(x) { ifelse(x == "x", 1, NA)}

# table 
impfstoff_tab <- impfstoff_ref %>%
  # change "x" to 1
  mutate(across(c(eclf18dmorf4:rabies), ReplaceX)) %>%
  mutate(nmbrofpthgns = ifelse(nmbrofpthgns == "NA", NA, nmbrofpthgns)) %>%
  # remove duplicate
  distinct() %>%
  # name for pattern search
  mutate(vacname_ptrn = gsub( "\\+", "\\\\\\+", vaccinename )) %>%
  left_join(marker_spvac, by = "vaccinename") %>%
  mutate(vactype = replace_na(vactype, "com"))

# Prep in 2 steps: if a column has the terms, write it down or else leave it blank. do left join
vac_ptrn <- paste(impfstoff_tab$vacname_ptrn, collapse = "|")

# to cal pathogen column names
pathogens <- names(impfstoff_ref)[-1:-2]

# pathogen covered by vaccine
# To use the number of covered pathogens I consider as best practice for the commercial vaccines. 
# Since the “spezifische” vaccines do also cover different pathogens (which we unfortunately 
# don’t know for all of them), I would suggest to rather code this variable by 1/0  (does use or 
# doesn’t use specific vaccine) rather than generating an artificial number of pathogens (1).

# step 1: match vaccine name from ref table, with medicine name in the handlsnmarznmttl
vac_coverage <- bhd %>%
  
  # subset 1 year
  dplyr::filter(erstsbhndlngsdtm >= as.Date('2021-01-01') &
                  erstsbhndlngsdtm <= as.Date('2021-12-31')) %>%
  
  # match vaccine name 
  dplyr::mutate(has_pos1 = regexpr(vac_ptrn, handlsnmarznmttl, ignore.case = TRUE)) %>%
  dplyr::mutate(has_pos2 = has_pos1 + attr(has_pos1, 'match.length') -1) %>%
  dplyr::mutate(has_vaccine = substr(handlsnmarznmttl, has_pos1, has_pos2)) %>%
  dplyr::left_join(impfstoff_tab, by=c("has_vaccine" = "vaccinename")) %>%
  mutate(vactype = replace_na(vactype, "notvaccinated"))

# step 2: count pathogens covered by each vaccine
ConvertToOne <- function(x) { ifelse(x > 0, 1, 0)}

vac_npathogens <- vac_coverage %>%
  # count pathogen
  group_by(betriebsid) %>%
  summarise(across(c(eclf18dmorf4:rabies), sum, na.rm = TRUE )) %>%
  ungroup() %>%
  # total pathogen as sum
  mutate(pathogen_ttl = base::rowSums( select(., eclf18dmorf4:rabies), na.rm = TRUE)) %>%
  # count pathogen coverage
  mutate(across(c(eclf18dmorf4:rabies), ConvertToOne, .names = "one_{col}" )) %>%
  mutate(pathogen_count = base::rowSums( select(., starts_with("one_")), na.rm = TRUE)) %>%
  select(-starts_with("one_")) %>%
  # drop pathogens 
  select(-(eclf18dmorf4:rabies))



##-----------------------------------------------------##
## 4) Gesundheitsbericht  : symptoms count             ##
##-----------------------------------------------------##

# Column names of the original file is stored in 2 rows. 
# We've preprocessed it in excel to make it more readable.
gbr_row1 <- cleanup_colname(gbr_raw$row1_name, minlen=12)[[2]]
gbr_row2 <- cleanup_colname(gbr_raw$row2_name, minlen=12)[[2]]

# combine
gbr_name <- paste(gbr_row1, gbr_row2, sep="_")

# mark duplicates
#length(col_name[duplicated(col_name)]) 
gbr_name[duplicated(gbr_name)] <- paste0(gbr_name[duplicated(gbr_name)], "_dup")

# check if there's no more duplicates
all.equal(sum(duplicated(gbr_name)), 0)

# update betriebsid col name
gbr_name[which(gbr_name == "betriebsid_info")] <- "betriebsid"
gbr_name[which(gbr_name == "tvdnr_info")] <- "tvdnummer"

# save headers
gbr_head_clean <- data.frame(
  ori = paste(gbr_raw$row1_name, gbr_raw$row2_name),
  new = gbr_name
)

# load data
gbr <- read_excel(file3, 
                  sheet = "BP Alle", col_names = gbr_name, 
                  skip = 3)

# symptoms
symt_names <- c("naehrzustand", "kuemmerer", "fruchtbarket", "milchfieber", "fieber", "durchfall", 
                "znsstoerungn", "mortalitaet", "kannibalisms", "niesen", "husten", "lahmheiten", 
                "hautvrndrngn", "juckreiz")
symt_ptrn <- paste(symt_names, collapse = '|')

## define farm groups
group_züchter <- c("Züchter", "Züchter-Mäster", "Ehem. Züchter")
group_mäster <- c("Mäster", "Ehem. Mäster")
group_afpring <- c("AFP Ring", "Betriebsgemeinschaft", "Narkosegerätegemeinschaft")

# subset column for symptoms & mark hybrid farms and afrpring
symptoms_cnt <- gbr %>% 
  select(betriebsid, produktnstyp_info, ringart_info, besuchsdatum_beschsprtkll, 
         sauen_info, absetzjager_info, remonten_info, eber_info, mastplaetze_info,
         naehrzustand_absetzjager:juckreiz_endmast, naehrzustand_saeugendsaun:juckreiz_mastjager) %>%
  # Betriebs tpye: Züchter or Mäster
  mutate(type2_info = ifelse(produktnstyp_info %in% group_züchter, "Züchter", 
                             ifelse(produktnstyp_info %in% group_mäster, "Mäster", "Hybrid"))) %>%
  # Ring type: drop and keep
  mutate(afpring_info = ifelse(ringart_info %in% group_afpring, "drop", "keep")) %>%
  
  # total animal
  rowwise() %>%
  mutate(ttlpig_info = sum(across(sauen_info:mastplaetze_info)), na.rm = TRUE) %>%
  ungroup %>%
  rename_with( ~gsub('$_dup', '', .x)) %>%
  
  # filter: ringart == keep & prodtype = züchter or mäster
  filter(type2_info %in% c("Züchter", "Mäster") & afpring_info == "keep")

# check data : Have I already filter the data correctly?
check_filtered_gbr <- symptoms_cnt %>% select(ends_with("_info")) %>% 
  gather("pigtype", "pign",  sauen_info, absetzjager_info, remonten_info, eber_info, mastplaetze_info) %>%
  group_by(produktnstyp_info, ringart_info, afpring_info, type2_info, pigtype) %>%
  summarise(nrows=n(), sum_animals=sum(pign), .groups = "keep")

#readr::write_csv(check_filtered_gbr, here::here("analysis/overview_filteredrows_gbr.csv"))
View(check_filtered_gbr)

# count symptoms per pig category
symt_saeugendsaun <- count_symptoms2(symptoms_cnt, piggrp = '_saeugendsaun')
symt_saugferkel <- count_symptoms2(symptoms_cnt, piggrp = '_saugferkel')
symt_galtsauen <- count_symptoms2(symptoms_cnt, piggrp = '_galtsauen')
symt_absetzjager <- count_symptoms2(symptoms_cnt, piggrp = '_absetzjager')
symt_remonten <- count_symptoms2(symptoms_cnt, piggrp = '_remonten')
symt_eber <- count_symptoms2(symptoms_cnt, piggrp = '_eber')
symt_mastjager <- count_symptoms2(symptoms_cnt, piggrp = '_mastjager')
symt_vormast <- count_symptoms2(symptoms_cnt, piggrp = '_vormast')
symt_endmast <- count_symptoms2(symptoms_cnt, piggrp = '_endmast')

# column to includes
gbr_cols <- c("betriebsid", "type2_info", "afpring_info", "besuchsdatum_beschsprtkll", "ttlpig_info",
               "sauen_info", "absetzjager_info", "remonten_info", "eber_info", "mastplaetze_info")

# join into 1 df
gbr_symptoms <- dplyr::bind_cols(
  symptoms_cnt[, gbr_cols], 
  symt_saeugendsaun[,4:5],
  symt_saugferkel[,4:5],
  symt_galtsauen[,4:5],
  symt_absetzjager[,4:5],
  symt_remonten[,4:5],
  symt_eber[,4:5],
  symt_mastjager[,4:5],
  symt_vormast[,4:5],
  symt_endmast[,4:5]
)

# collect pig group
pig_group <- c("saeugendsaun", "saugferkel", "galtsauen", "absetzjager", "remonten", "eber", "mastjager", "vormast", "endmast")
pig_smscore_ptrn <- paste0('smscore_', pig_group, collapse = '|')
pig_smpct_ptrn <- paste0('smpct_', pig_group, collapse = '|')

# add the sum across all pig group
gbr_symptoms <- gbr_symptoms %>%
  rowwise() %>%
  mutate(
    smscore_all = ifelse(
      all(is.na(c_across(matches( pig_smscore_ptrn )))), NA,
      mean(c_across(matches( pig_smscore_ptrn )), na.rm = TRUE)
    ),
    smpct_all = ifelse(
      all(is.na(c_across(matches( pig_smpct_ptrn )))), NA,
      mean(c_across(matches( pig_smpct_ptrn )), na.rm = TRUE)
    )
  ) %>%
  ungroup() 

##-----------------------------------------------------##
## 5) TBI                                              ##
##-----------------------------------------------------##

# read head
tbi_head <- readr::read_delim(file6,
                              ";", escape_double = FALSE, locale = locale(encoding = "ISO-8859-1"), 
                              col_names = FALSE, n_max = 1, trim_ws = TRUE, show_col_types = FALSE)

# clean up column names
tbi_head_clean <- data.frame(
  ori = cleanup_colname(as.character(tbi_head), minlen=16)[[1]],
  new = cleanup_colname(as.character(tbi_head), minlen=16)[[2]]
)

# update betriebs id col name
tbi_head_clean[which(tbi_head_clean$new == "tvdnr"), "new"] <- "tvdnummer"
tbi_head_clean[which(tbi_head_clean$new == "tbi"), "new"] <- "tbi_val"

# read col spec
tbi_spec <- readr::spec_delim(file6, ";")

# update the spec
names(tbi_spec$cols) <- tbi_head_clean[[2]]

# id column as character
selected_id <- c("tvdnummer")
tbi_spec$cols[selected_id] <- list(col_character())

# integer as integer
selected_db <- grep("nr$|jahr|quartal", names(tbi_spec$cols), value = TRUE)
tbi_spec$cols[selected_db] <- list(col_integer())

# preview
tbi_spec

# read the file
tbi <- readr::read_delim(file6,
                         ";", escape_double = FALSE, locale = locale(encoding = "ISO-8859-1"), 
                         col_names = tbi_head_clean$new, col_types = tbi_spec, skip = 1,
                         trim_ws = TRUE) 


# count by tvdnummer
qtrjahr_keep <- c("2020_4", "2021_1", "2021_2", "2021_3")

prep_tbi <- tbi %>%
  # drop unused
  select(-maxbetriebsnr, -anzahlbetriebsnr) %>%
  
  # simplify Tierkategorie naming before column spread
  mutate(tierkategorie = cleanup_colname(tierkategorie, minlen = 12)[[2]] ) %>%
  # use all qtr
  unite("tierkqtr", c(jahr, quartal, tierkategorie), remove = TRUE ) %>%
  spread(tierkqtr, tbi_val) %>%
  # add prefix to colname
  rename_with( ~ paste0('tbi_', .x), where(is.numeric) )


prep_tbi <- tbi %>%
  # drop unused
  select(-maxbetriebsnr, -anzahlbetriebsnr) %>%
  
  # simplify Tierkategorie naming before column spread
  mutate(tierkategorie = cleanup_colname(tierkategorie, minlen = 12)[[2]] ) %>%
  # use all qtr
  unite("tierkqtr", c(jahr, quartal, tierkategorie), remove = TRUE ) %>%
  spread(tierkqtr, tbi_val) %>%
  # add prefix to colname
  rename_with( ~ paste0('tbi_', .x), where(is.numeric) )


# for züchter : use summary for 1 year
prep_tbi_yearly <- tbi %>%
  unite("qtrjahr", jahr:quartal, remove = FALSE ) %>%
  dplyr::filter(qtrjahr %in% qtrjahr_keep) %>%
  group_by(tvdnummer, tierkategorie) %>%
  summarise(tbiavg = mean( tbi_val, na.rm = TRUE )) %>%
  spread(tierkategorie, tbiavg) %>%
  rename_with( ~ paste0('tbiavg_', .x), where(is.numeric) )




##-----------------------------------------------------##
## 6) LSG                                              ##
##-----------------------------------------------------##

# read head
lsg_head <- readr::read_delim(file5,
                              ";", escape_double = FALSE, locale = locale(encoding = "ISO-8859-1"), 
                              col_names = FALSE, n_max = 1, trim_ws = TRUE, show_col_types = FALSE)

# clean up column names
lsg_head_clean <- data.frame(
  ori = cleanup_colname(as.character(lsg_head), minlen=16)[[1]],
  new = cleanup_colname(as.character(lsg_head), minlen=16)[[2]]
)

# update betriebs id col name
lsg_head_clean[which(lsg_head_clean$new == "tvdnr"), "new"] <- "tvdnummer"

# read col spec
lsg_spec <- readr::spec_delim(file5, ";")

# update the spec
names(lsg_spec$cols) <- lsg_head_clean[[2]]

# id column as character
selected_id <- c("tvdnummer")
lsg_spec$cols[selected_id] <- list(col_character())

# double as double
selected_db <- c("abgmsundeber")
lsg_spec$cols[selected_db] <- list(col_double())

# integer as integer
selected_db <- grep("nr$|jahr|quartal", names(lsg_spec$cols), value = TRUE)
lsg_spec$cols[selected_db] <- list(col_integer())

# preview
lsg_spec

# read the file
lsg <- readr::read_delim(file5,
                         ";", escape_double = FALSE, locale = locale(encoding = "ISO-8859-1"), 
                         col_names = lsg_head_clean$new, col_types = lsg_spec, skip = 1,
                         trim_ws = TRUE) 


# include all 4 quarters
# each measurement as column
prep_lsg_qtr <- lsg %>%
  # drop unused
  select(-maxbetriebsnr, -anzahlbetriebsnr) %>%
  
  # use all quarters
  unite("qtrjahr", jahr:quartal, remove = TRUE ) %>%
  
  # unite measurement and quarters
  gather("measurement", "val", lebendgebferkel:abgmsundeber) %>%
  unite("qtrlsg", c(qtrjahr, measurement), remove = TRUE) %>%
  
  # spread quarterly measurement vs value
  spread(qtrlsg, val ) %>%
  
  # add prefix to colname
  rename_with( ~ paste0('lsg_', .x), where(is.numeric) )

# subset to 4 quarters
qtrjahr_keep <- c("2020_4", "2021_1", "2021_2", "2021_3")

# each measurement in a column
prep_lsg_tall <- lsg %>%
  unite("qtrjahr", jahr:quartal, remove = FALSE ) %>%
  
  # consider 4 quarters
  dplyr::filter(qtrjahr %in% qtrjahr_keep) %>%

  # each quartal as column
  mutate(quartal = paste0("q", quartal)) %>%
  gather("measurement", "val", lebendgebferkel:abgmsundeber) %>%
  spread(quartal, val ) %>%
  
  # do not include zero in calculating the average
  #mutate(across(q1:q2, ~ifelse(.x == 0, NA, .x), .names = "{col}_zero")) %>%
  
  # calculate each measurement average over all quartal
  rowwise() %>%
  mutate(avglsg = mean(c(q1, q2, q4), na.rm = TRUE)) %>%
  mutate(avglsg = ifelse(is.nan(avglsg), NA, avglsg)) %>%
  ungroup() 

# each measurement as column
prep_lsg_wide <- lsg %>%
  # drop unused
  select(-maxbetriebsnr, -anzahlbetriebsnr) %>%
  
  # consider 4 quarters
  unite("qtrjahr", jahr:quartal, remove = FALSE ) %>%
  dplyr::filter(qtrjahr %in% qtrjahr_keep) %>%
  
  # do not include zero in calculating the average
  # mutate(across(lebendgebferkel:abgmsundeber, ~ifelse(.x == 0, NA, .x))) %>%
  
  # calculate average over all measurement
  group_by(tvdnummer) %>%
  summarise(across(lebendgebferkel:abgmsundeber, ~mean(.x, na.rm = TRUE) )) %>%
  
  # add prefix to colname
  rename_with( ~ paste0('lsg_', .x), where(is.numeric) )


##-----------------------------------------------------##
## 7) save data                                        ##
##-----------------------------------------------------##

save(vac_coverage, vac_npathogens, gbr_symptoms, 
     prep_tbi, prep_tbi_yearly, prep_lsg_qtr, prep_lsg_wide,
     gbr_head_clean, bhd_head_clean, lsg_head_clean, tbi_head_clean,
     file = here::here('analysis/feature_extracted_1.Rdata'))
