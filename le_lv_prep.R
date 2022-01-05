###################################################
#   Author: Isabel De Ramos                       #
#   Date Created: 4 January 2022                  #
#   SAMPLE CODE - DATA MANAGEMENT                 #
#   Function: Life Expectancy Data Preparation    #
###################################################

## this code along with its files and generated .rdata can be accessed via my GitHub repository
# https://github.com/isabelderamos/thesis-sample  

# loading in libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)


# NOTE: 1989-2019 raw mortality text files were obtained from NCHS 
# 1999-2019 were cleaned and exported as .csv's beforehand 
# these mortality .csv files can be accessed via Drexel BILAL_DP5 server and downloaded onto local computer 
# as they are too big to upload to GitHub


#################### PREP WORK #################### 
# (1) get mortality .csv files
# (2) create function that formats fips 
# (3) prepare population denominators 


#~~~~ (1) get mortality .csv files ~~~~#
### Get mortality file names
mort_file_list <- list.files(pattern = "mort[1|2]", full.names = TRUE)
df_mort_files = tibble(file = list.files(pattern = "mort[1|2]", full.names = TRUE)) %>% 
  mutate(year = str_sub(file,7,10) %>% parse_integer() ) %>% 
  arrange(year) 


#~~~~ (2) format fips ~~~~#
# create function that pads all fips (adds back leading zeroes)
format_5dig_fips <- function(fp) {
  ifelse(str_length(fp)==4,
         str_pad(fp, width=5, side="left", pad="0"),
         fp)
}
format_5digfips_v <- Vectorize(format_5dig_fips)


#~~~~ (3) preparing population denominators ~~~~#
# following population files can be accessed via Drexel BILAL_DP5 server

# pop_county_asrh
# years 1990-2019 available 
# split by age group, gender, and race
# race is 4 categories = H, NHB, NHW, NHO

# bridged_pop_county_asrh
# years 1990-2019 available 
# split by age group (0,5,10,15...), gender, and race
# race is 5 categories = H, NHB, NHW, NHAIAN, NHAPI

# bridged_pop_county_asrhLT 
# years 2000-2019 available
# split by age group (0,1,5,10...), gender, and race
# race is 5 categories = H, NHB, NHW, NHAIAN, NHAPI

# using bridged_pop_county_asrhLT as it's more appropriate to build life tables with age groups 0,1,5,10
load('../Data/bridged_pop_county_asrhLT.rdata')
pop_county <- pop_county %>% as_tibble() %>% 
  mutate(year=as.character(year),
         fips=as.character(fips),
         age_5yr_group=as.character(age_5yr_group),
         fips=format_5digfips_v(fips)) %>% # apply 5 dig fips function to pop_county_race
  rename(`sex`=`male`,
         `race`=`hispanic_bridged`)


#################### CREATING FUNCTIONS TO CLEAN AND PREPARE NCHS MORTALITY AND POPULATION DATA #################### 
# (1) clean_nchs 
# (2) clean_nchs1
# (3) clean_popdenoms
# (4) clean_popdenoms1


#~~~~ (1) clean_nchs ~~~~#
# clean_nchs extracts variables of interest from mortality data - year, fips, age, sex, and race
clean_nchs <- function(file) {
  #file <- df_mort_files %>% filter(year==2000) %>% pull(file) 
  
  # isolate all fips / all groups for given year 
  nchs_dta_tmp <- read_csv(file, n_max = Inf) %>% as_tibble() 
  
  # find counts per age group by race
  nchs_dta_tmp <- nchs_dta_tmp %>% 
    transmute(year=as.character(death_year),
              fips=res_fips_effective18,
              age=age_red,
              age_5yr_group=case_when(
                age == 0 ~ '0',
                age < 5 ~ '1',
                between(age, 5,84)~as.character(floor(age/5)*5),
                TRUE~"85"),
              sex=male,
              race=hispanic_bridged) %>%
    group_by(fips, year, age_5yr_group, sex, race) %>% 
    count(age_5yr_group) %>%
    rename(count=n) %>% 
    arrange(fips, year, age_5yr_group, sex, race) %>%
    ungroup()
  print(nchs_dta_tmp) %>% 
    return()
}

#~~~~ (2) clean_nchs1 ~~~~#
# clean_nchs1 fixes problematic FIPS  
## (a) filters out non-continental FIPS 
## (b) fixes FIPS that were renamed over time 
## (c) removes foreign fips overall (foreign residents)
clean_nchs1 <- function(nchs_dta_tmp) {
  
  # (a) filters out non-continental FIPS 
  nonUS_FIPS=c("^(02)","^(60)","^(66)","^(69)","^(72)","^(78)") 
  # Alaska, American Samoa, Guam, Northern Mariana Islands, Puerto Rico, Virgin Islands
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!grepl(paste(nonUS_FIPS, collapse="|"), fips))
  
  # (b) fixes FIPS that were renamed over time 
  nchs_dta_tmp <- nchs_dta_tmp %>% mutate(fips=ifelse(fips=="12025", "12086", fips)) %>% # fixing 12025
    mutate(fips=ifelse(fips=="51560", "51005", fips)) %>% # fixing 51560
    # special case: 30113 
    # pop_denoms %>% filter(fips%in%"30031") %>% count(fips, wt=pop_denom) >> population of 1863635
    # pop_denoms %>% filter(fips%in%"30067") %>% count(fips, wt=pop_denom) >> population of 332532
    # since FIPS 30031 has largest population, 30113, 30067 and 30031 will all merge into 30031 
    mutate(fips=ifelse(fips%in%c("30113", "30067"), "30031", fips))
  
  # (b) removes foreign fips overall (foreign residents, fips=00000) and missing FIPS
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!fips%in%c("00000", "00999")) %>% 
    filter(!is.na(fips)==TRUE) %>% 
    arrange(fips, year, age_5yr_group, sex, race) %>% 
    group_by(fips, year, age_5yr_group, sex, race) %>%
    summarise(count=sum(count)) %>% 
    ungroup()
  print(nchs_dta_tmp) %>% 
    return()
}


#~~~~ (3) clean_popdenoms ~~~~#
# clean_popdenoms extracts variables of interest - year, fips, age_5yr_group, sex, race, and pop_county
clean_popdenoms <- function(file) {
  #file <- df_mort_files %>% filter(year==2000) %>% pull(file) 
  
  # isolate all fips / all groups for given year 
  pop_county_tmp <- pop_county %>% filter(year%in%substr(file, 7,10)) %>%
    group_by(fips, year, age_5yr_group, sex, race) %>% 
    summarise(pop_denom=sum(pop_county)) %>%
    ungroup()
  print(pop_county_tmp) %>%
    return()
}


#~~~~ (4) clean_popdenoms1 ~~~~#
## clean_popdenoms1  fixes problematic FIPS  
## (a) fixes FIPS that were renamed over time 
clean_popdenoms1 <- function(pop_county_tmp) {

  # (a) fixes FIPS that were renamed over time 
  pop_county_tmp <- pop_county_tmp %>% mutate(fips=ifelse(fips=="46113", "46102", fips)) %>% # fixing 46102
    mutate(fips=ifelse(fips=="51560", "51005", fips)) %>%  # fixing 51560
    mutate(fips=ifelse(fips%in%c("30113", "30067"), "30031", fips)) %>% # fixing 30113 (see explanation above in clean_nchs1)
    arrange(year, age_5yr_group, race) %>% 
    group_by(fips, year, age_5yr_group, sex, race) %>%
    summarise(pop_denom=sum(pop_denom)) %>% 
    ungroup()
  print(pop_county_tmp) %>%
    return()
}




#################### CREATING FUNCTIONS FOR STRATIFYING FUNCTIONS #################### 
# (1) crosswalk_metro_regions


#~~~~ (1) crosswalk_metro_regions ~~~~#
# crosswalk_metro_regions joins NCHS mortality data to population denoms, 
# adds 2013 NCHS Urban-Rural Classification Scheme to NCHS mortality data (6 levels)
# adds binary urbanicity variable = metropolitan (1) or nonmetropolitan (0) (2 levels)
# adds U.S. Census Regions & Divisions (4 levels)

crosswalk_metro_regions <- function(nchs_dta_tmp, pop_county_tmp) {
  # join mortality (NCHS) to pops (POP DENOMS)
  final_dta_tmp <- nchs_dta_tmp %>% full_join(pop_county_tmp, by=c("fips", "year", "age_5yr_group", "sex", "race")) %>%
    arrange(fips, year, age_5yr_group, sex, race)
  # creating template to see full combination of fips / age groups / sexes/ races
  template <- expand_grid(fips=unique(final_dta_tmp$fips), 
                          age_5yr_group=unique(final_dta_tmp$age_5yr_group),
                          sex=unique(final_dta_tmp$sex),
                          race=unique(final_dta_tmp$race)) %>% 
    arrange(fips, age_5yr_group, sex, race)
  # replacing NA counts (deaths) with 0 
  final_dta_tmp <- template %>% full_join(final_dta_tmp, by=c("fips", "age_5yr_group", "sex", "race")) %>%
    fill(year) %>%
    mutate(count=replace_na(count,0)) %>% 
    arrange(fips, year, age_5yr_group, sex, race)
  #~~~~~ ADDING URBANICITY CLASSIFICATION ~~~~~#
  # reading in 2013 NCHS urban-rural crosswalk 
  xwalk_fips_urbanrural <- read.csv('../Crosswalks/US/Clean/fips_urbanrural_xwalk.csv',  
                                    sep = ',', header =  T, colClasses = "character") %>% as_tibble() %>%
    transmute(fips=format_5digfips_v(fips),
              `metro6`=`X2013_code`)
  # merging urban-rural codes to data
  final_dta_tmp <- final_dta_tmp %>% full_join(xwalk_fips_urbanrural, by="fips") %>% 
    # creating `metro_class` column to delineate metropolitan or nonmetropolitan
    mutate(metro2=ifelse(metro6%in%c("1", "2", "3", "4"), "1", "0"))
  #~~~~~ ADDING U.S. CENSUS REGIONS & DIVISIONS CLASSIFICATION ~~~~~#
  # reading in U.S. Census Regions dataset 
  xwalk_census_regions <- read.csv('../Crosswalks/US/Clean/fips_censusregion_xwalk.csv',
                                   sep = ',', header =  T, colClasses = "character") %>% as_tibble() 
  final_dta_tmp <- final_dta_tmp %>% mutate(state_code=substr(fips, 1, 2)) %>% 
    full_join(xwalk_census_regions, by="state_code") %>% 
    rename(`census_region`=`region`) %>% 
    select(-state_code, -division)
  # excluding NA data (caused by faulty fips)
  final_dta_tmp <- final_dta_tmp %>% filter(!is.na(year)==TRUE)
  print(final_dta_tmp) %>%
    return()
}




#################### EXECUTING FUNCTIONS #################### 


### For my MS thesis project, only examined years 2000-2019
mort_files_tmp = df_mort_files %>% filter(between(year, 2000, 2019)) %>% pull(file) 

# clean mortality data
nchs <- map_dfr(mort_files_tmp, ~clean_nchs(.x))
nchs <- clean_nchs1(nchs)

# prepare population denominators
pop_denoms <- map_dfr(mort_files_tmp, ~clean_popdenoms(.x))
pop_denoms <- clean_popdenoms1(pop_denoms)

# create master file by joining mortality and population denominators, and adding urbanicity and geography categorizations
master_dta <- crosswalk_metro_regions(nchs, pop_denoms)
save(master_dta, file='../Clean/00_19_nchs_mortality.rdata')


