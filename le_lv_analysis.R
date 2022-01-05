###################################################
#   Author: Isabel De Ramos                       #
#   Date Created: 4 January 2022                  #
#   SAMPLE CODE - DATA ANALYSIS                   #
#   Function: Life Expectancy Calculations        #
###################################################

## code, files, and final .rdata can be accessed via my GitHub repository
# https://github.com/isabelderamos/thesis-sample 

### loading libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)
library(scales)
library(multcomp)
library(RColorBrewer)
library(scales)
library(ggpubr)
#install.packages("devtools")
#library(devtools)
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("remotes")
#remotes::install_github("timriffe/DemoTools")
#install.packages("DemoTools")
library(DemoTools)


#################### PREP WORK #################### 
# (1) ggplot helpers
# (2) Camarda's functions for lifetables, life expectancy 
# (3) loading in .rdata
# (4) create 5-year pooled periods and factors
# (5) create le_lv function

#~~~~ (1) ggplot helpers ~~~~#
### prep for generating plots ###
select<-dplyr::select
fontsize<-16
isabel_theme <-theme_bw()+
  theme(axis.text=element_text(color="black", size=fontsize),
        axis.title=element_text(color="black", size=fontsize, face="bold"),
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=fontsize, face="bold"),
        legend.text=element_text(color="black", size=fontsize),
        legend.title=element_text(color="black", size=fontsize, face="bold"))
cols<-c(brewer_pal(type="qual", palette=2)(8),"blue", brewer_pal(type="qual", palette=2)(8), "blue", "red")
shapes<-c(rep(21, times=9), rep(22, times=9), 22)


#~~~~ (2) loading in Camarda's functions for lifetables, life expectancy ~~~~#
# see https://sites.google.com/site/carlogiovannicamarda/r-stuff/life-expectancy-confidence-interval 
## lifetable() calculates life expectancy
## CIex() calculates 95% CI 
source('Data/LifeTableFUN.R')


#~~~~ (3) loading in master datafile (see LE_data_prep.R) ~~~~#
load("Clean/00_19_nchs_mortality.rdata")
dta <- master_dta %>% select(-fips)


#~~~~ (4) create year5 variable that pools 5-year periods from 2000-2019 ~~~~#
dta <- master_dta %>% mutate(year5=case_when(
  year%in%c(2015:2019) ~ as.character('2015-2019'),
  year%in%c(2010:2014) ~ as.character('2010-2014'),
  year%in%c(2005:2009) ~ as.character('2005-2009'),
  year%in%c(2000:2004) ~ as.character('2000-2004'))) %>% 
  # factoring metro2, metro6, sex, census_region to help in ggplot later
  mutate(metro2=factor(metro2, levels=c(1,0),
                       labels=c("Metropolitan", "Non-metro")),
         metro6=factor(metro6, levels=c(1:6),
                       labels=c("Large central", "Large fringe", "Medium", "Small", "Micro", "Noncore")),
         sex=factor(sex, level=c(1,0),
                    labels=c("Male", "Female")),
         census_region=factor(census_region, levels=c(1:4),
                              labels=c("Northeast", "Midwest", "South", "West")))

#~~~~ (5) create le_lv function that calculates life expectancy and lifespan variation  ~~~~#
le_lv<-function(.x, .y){
  ## .x -> the data with one race and metro
  # example: .x<-data_baseline %>% filter(metro2=="Metropolitan", race=="H")
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=0, ns=1000, level=0.95)
  # adding variables needed to calculate lifespan variation 
  lt <- lt %>% mutate(xbar=NA_real_,
                      noname=NA_real_,
                      v=NA_real_,
                      sd=NA_real_)
  lt <- lt %>% mutate(xbar=x+ax,
                      noname=case_when(
                        x==0 ~ dx/lx*(xbar-ex)^2,
                        x!=0 ~ {
                          lx_forlv <- lt %>% filter(x==0) %>% pull(lx)
                          ex_forlv <- lt %>% filter(x==0) %>% pull(ex)
                          dx/lx_forlv*(xbar-ex_forlv)^2}),
                      v=rev(cumsum(rev(noname))),
                      sd=sqrt(v))
  # extracting LE
  le <- lt %>% filter(x==0) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting LV 
  lv <- lt %>% filter(x==0) %>% pull(sd)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv=lv)
}



############ BASELINE ANALYSIS: 2015-2019, GENDER POOLED, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# baseline analysis examines disparities in life expectancy and lifespan variation by race/ethnicity 
# and by urbanicity (2 categories) within the 5-year pooled period of 2015-2019
results_baseline <- dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, race, age_5yr_group) %>% 
  group_by(race, metro2) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(.)) 


####### table 1: BA_lelv_table ######## 
# le/lv and absolute/relative differences
# spreading LE and 95% CIs by metro2 
BA_lelv_table <- results_baseline %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_baseline %>% select(race, metro2, le) %>% 
      spread(metro2, le) %>% 
      rename(metro_le=`Metropolitan`,
             nonmetro_le=`Non-metro`) %>% 
      full_join(
        results_baseline %>% select(race, metro2, lv) %>% 
          spread(metro2, lv) %>% 
          rename(metro_lv=`Metropolitan`,
                 nonmetro_lv=`Non-metro`)) %>% 
      mutate(abs_le=metro_le-nonmetro_le,
             rel_le=metro_le/nonmetro_le,
             abs_lv=metro_lv-nonmetro_lv,
             rel_lv=metro_lv/nonmetro_lv)) %>% 
  ungroup() %>% 
  transmute(race, 
         metro_le_ci, 
         nonmetro_le_ci,
         abs_le=as.numeric(format(abs_le, digits=2, nsmall=2)), 
         rel_le=as.numeric(format(rel_le, digits=2, nsmall=2)), 
         metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1)), 
         nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1)), 
         abs_lv=as.numeric(format(abs_lv, digits=2, nsmall=2)), 
         rel_lv=as.numeric(format(rel_lv, digits=2, nsmall=2)))





############ ANALYSIS 1 LONGITUDINAL: 2000-2019, GENDER POOLED, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# analysis 1 examines disparities in life expectancy and lifespan variation by race/ethnicity 
# and by urbanicity (2 categories) longitudinally across 5-year pooled period from 2000-2019
results_longitudinal <- dta %>% group_by(year5, age_5yr_group, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(year5, metro2, race, age_5yr_group) %>% 
  group_by(year5, race, metro2) %>% 
  group_modify(~le_lv(.)) 



####### table 2: A1_le_table #### 
# LE longitudinal by race/urbanicity
df <- results_longitudinal %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(year5, race, metro2, le_ci) %>% 
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) 
# spreading metro LE, nonmetro LE, then combining
A1_le_table <-
  df %>% select(year5, race, metro_le_ci) %>% 
  spread(year5, metro_le_ci) %>% 
  rename(metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(year5, race, nonmetro_le_ci) %>% 
      spread(year5, nonmetro_le_ci) %>% 
      rename(nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`))

####### table 3: A1_lv_table #### 
# LV longitudinal by race/urbanicity
df <- results_longitudinal %>%
  transmute(year5, race, metro2, lv=as.numeric(format(lv, digits=1, nsmall=1))) %>% 
  spread(metro2, lv) %>% 
  rename(metro_lv=`Metropolitan`,
         nonmetro_lv=`Non-metro`) 
# spreading metro LV, nonmetro LV, then combining
A1_lv_table <-
  df %>% select(year5, race, metro_lv) %>% 
  spread(year5, metro_lv) %>% 
  rename(metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(year5, race, nonmetro_lv) %>% 
      spread(year5, nonmetro_lv) %>% 
      rename(nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`))



####### figure 1: A1_lelv_longtrends ####### 
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4))
df <- df %>% mutate(year5=factor(year5, 
                                 levels=c(2000, 2005, 2010, 2015),
                                 labels=c("2000-2004", 
                                          "2005-2009",
                                          "2010-2014",
                                          "2015-2019")))
# le longitudinal figure 
le_longtrend <- ggplot(df, aes(year5, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(color=race))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Year",
       y="Life Expectancy",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# lv longitudinal figure 
lv_longtrend <- ggplot(df, aes(year5, lv, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(color=race))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Year",
       y="Lifespan Variation",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# arranging figures above into one figure with two rows
A1_lelv_longtrends<- ggarrange(le_longtrend, lv_longtrend, nrow=2, 
                               common.legend = TRUE, 
                               legend="bottom")




####### figure 2: A1_le_vs_lv ####### 
# association b/w LE and LV by race/ethnicity, 2000-2019
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  select(year5, race, metro2, le, lv) %>% 
  arrange(race, metro2, year5)

A1_le_vs_lv <- ggplot(df, aes(lv, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_path(aes(color=race))+
  geom_text(aes(label=year5), vjust=-1) +
  facet_wrap(~metro2)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Lifespan Variation",
       y="Life Expectancy",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())




############ ANALYSIS 2 GENDER: 2015-2019, MALE V FEMALE, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# analysis 2 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (2 categories), and by gender (males v females) within the 5-year pooled period of 2015-2019
results_gender <- dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race, sex) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, race, sex, age_5yr_group) %>% 
  group_by(race, metro2, sex) %>% 
  group_modify(~le_lv(.)) 




######## figure 3: A2_le_bars ######## 
# le by race/urbanicity and sex
# making dataset appropriate for ggplot
df <- results_baseline %>% ungroup() %>% 
  select(race, metro2, le) %>% 
  full_join(results_gender, by=c("race", "metro2")) %>% 
  transmute(race, metro2, total_sex=`le.x`, sex, le.y) %>% 
  spread(sex, `le.y`) %>% 
  rename(male_le=`Male`,
         female_le=`Female`) %>% 
  gather("sex", "le", 3:5) %>% 
  mutate(sex=case_when(
    sex=="total_sex" ~ 0,
    sex=="male_le" ~ 1,
    sex=="female_le" ~ 2,
  )) %>% 
  mutate(sex=factor(sex, levels=c(0,1,2), labels=c("Total", "Male", "Female"))) %>% 
  arrange(race)

A2_le_bars <- ggplot(df, aes(x=race, y=le, fill=as.factor(sex))) + 
  geom_col(position=position_dodge(width=0.7), width=0.7) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-1) +
  scale_y_continuous(limits=c(df %>% ungroup() %>% select(le) %>% min,
                              df %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Race/Ethnicity",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  isabel_theme




######## table 4: A2_lelv_table ######## 
# le/lv by sex
df <- results_gender %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",                      
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_gender %>% select(race, metro2, sex, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`)) 

A2_lelv_table <- df %>% select(race, sex, metro_le_ci) %>% 
  spread(sex, metro_le_ci) %>% 
  rename(metro_male_le=`Male`,
         metro_female_le=`Female`) %>% 
  full_join(
    df %>% select(race, sex, nonmetro_le_ci) %>% 
      spread(sex, nonmetro_le_ci) %>% 
      rename(nonmetro_male_le=`Male`,
             nonmetro_female_le=`Female`)) %>% 
  full_join(
    df %>% transmute(race, sex, metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1))) %>% 
      spread(sex, metro_lv) %>% 
      rename(metro_male_lv=`Male`,
             metro_female_lv=`Female`)) %>% 
  full_join(
    df %>% transmute(race, sex, nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1))) %>% 
      spread(sex, nonmetro_lv) %>% 
      rename(nonmetro_male_lv=`Male`,
             nonmetro_female_lv=`Female`))


######## figure 4: A2_MvF_disparities ######## 
# MvF scatterplot (LE and LV)  

## first, creating LE scatterplot
# spreading LE and 95% CIs by sex 
df <- results_gender %>% full_join(
  results_gender %>% select(race, metro2, sex, le) %>% 
    spread(sex,le) %>% 
    rename(male_le=`Male`,
           female_le=`Female`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, sex, lci) %>% 
      spread(sex,lci) %>% 
      rename(male_lci=`Male`,
             female_lci=`Female`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, sex, uci) %>% 
      spread(sex,uci) %>% 
      rename(male_uci=`Male`,
             female_uci=`Female`)) %>% 
  select(race, metro2, sex, male_le, male_lci, male_uci, female_le, female_lci, female_uci)
# generating scatterplot
A2_le_scatter <- ggplot(df, aes(x=male_le, y=female_le)) +
  #geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(male_le=seq(70, 90, by=1),
                                   female_le=seq(70, 90, by=1)))+
  geom_hline(yintercept = 70, lty=2)+
  geom_vline(xintercept = 70, lty=2)+
  geom_linerange(aes(ymin=female_lci, ymax=female_uci))+
  geom_linerange(aes(xmin=male_lci, xmax=male_uci))+
  geom_point(aes(color=race, shape=metro2), size=3) +
  #scale_y_continuous(limits=c(df %>% ungroup() %>% select(female_lci) %>% min,
  #                            df %>% ungroup() %>% select(female_uci) %>% max))+
  #scale_x_continuous(limits=c(df %>% ungroup() %>% select(male_lci) %>% min,
  #                            df %>% ungroup() %>% select(male_uci) %>% max))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Male Life Expectancy (95% CI)",
       y="Female Life Expectancy (95% CI)",
       color="", fill="", shape="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# next, creating LV scatterplot
# spreading LV by sex 
df <- results_gender %>% select(race, metro2, sex, lv) %>% 
  spread(sex,lv) %>% 
  rename(male_lv=`Male`,
         female_lv=`Female`)
# generating scatterplot
A2_lv_scatter <- 
  ggplot(df, aes(x=male_lv, y=female_lv)) +
  # geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(male_lv=seq(14, 22, by=1),
                                   female_lv=seq(14, 22, by=1)))+
  geom_hline(yintercept = 14, lty=2)+
  geom_vline(xintercept = 14, lty=2)+
  geom_point(aes(color=race, shape=metro2), size=3) +
  # scale_y_continuous(limits=c(df %>% ungroup() %>% select(female_lv) %>% min,
  #                            df %>% ungroup() %>% select(female_lv) %>% max))+
  # scale_x_continuous(limits=c(df %>% ungroup() %>% select(male_lv) %>% min,
  #                            df %>% ungroup() %>% select(male_lv) %>% max))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Male Lifespan Variation",
       y="Female Lifespan Variation",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# arranging figures above into one figure with two panels
A2_MvF_disparities <- ggarrange(A2_le_scatter, A2_lv_scatter, ncol=2, 
                                common.legend = TRUE, 
                                legend="bottom")





############ ANALYSIS 3 URBANICITY: 2015-2019, GENDER POOLED, GEOGRAPHY POOLED, 6 METRO CATEGORIES ############ 
# analysis 3 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (6 categories) within the 5-year pooled period of 2015-2019
results_urbanicity <- dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro6, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro6, race, age_5yr_group) %>% 
  group_by(race, metro6) %>%
  group_modify(~le_lv(.)) 



######## table 5: A3_le_table ######## 
# le by race/urbanicity (6 metro categories)
A3_le_table <- results_urbanicity %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro6, le_ci) %>% 
  spread(metro6, le_ci)


######## table 6: A3_lv_table ######## 
# lv by race/urbanicity (6 metro categories)
A3_lv_table <- results_urbanicity %>% transmute(race, metro6, lv=as.numeric(format(lv, digits=1, nsmall=1))) %>% 
  spread(metro6, lv)


######## figure 5:  A3_le_bars ######## 
# one panel per race, one bar per metro category
A3_le_bars <- 
  ggplot(results_urbanicity, aes(x=metro6, y=le)) + 
  geom_col(aes(fill=metro6), color="black", position=position_dodge(width=0.7), width=0.7) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-0.5) +
  scale_y_continuous(limits=c(results_urbanicity %>% ungroup() %>% select(le) %>% min,
                              results_urbanicity %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_x_discrete(labels=function(results_urbanicity) str_wrap(results_urbanicity, width = 10))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Levels of Urbanization",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~race, nrow=1)+
  guides(fill=F)+
  isabel_theme +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))




############ ANALYSIS 4 GEOGRAPHY: 2015-2019, GENDER POOLED, 4 CENSUS REGIONS, 2 METRO CATEGORIES  ######
# analysis 4 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (2 categories), and by census region (NE, MW, S, W) within the 5-year pooled period of 2015-2019
results_geography <- dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race, census_region) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, race, census_region, age_5yr_group) %>% 
  group_by(race, metro2, census_region) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(.)) 


######## table 7:  A4_le_table ######## 
# le by race/urbanicity and census region 
df <- results_geography %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_geography %>% select(race, metro2, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`)) 

A4_le_table <- df %>% select(census_region, race, metro_le_ci) %>% 
  spread(census_region, metro_le_ci) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% select(census_region, race, nonmetro_le_ci) %>% 
      spread(census_region, nonmetro_le_ci) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro)


######## table 8:  A4_lv_table ######## 
# lv by race/urbanicity and census region 
A4_lv_table <- df %>% transmute(census_region, race, metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1))) %>% 
  spread(census_region, metro_lv) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% transmute(census_region, race, nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1))) %>% 
      spread(census_region, nonmetro_lv) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro)


######## figure 6:  A4_le_bars1 ######## 
A4_le_bars1 <- ggplot(results_geography, aes(x=census_region, y=le, fill=as.factor(race))) + 
  geom_col(position=position_dodge(width=0.8), width=0.8) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Census Region",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  isabel_theme


######## figure 7:  A4_le_bars2 ######## 
A4_le_bars2 <- ggplot(results_geography, aes(x=race, y=le, fill=as.factor(metro2))) + 
  geom_col(position=position_dodge(width=0.7), width=0.7) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-0.6) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Race/Ethnicity",
       y="Life Expectancy",
       color="", fill="") +
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~census_region, scales="free")+
  isabel_theme


######## figure 8:  A4_le_bars3 ######## 
A4_le_bars3 <- ggplot(results_geography, aes(x=metro2, y=le, fill=as.factor(race))) + 
  geom_col(position=position_dodge(width=0.7), width=0.7) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-0.6) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Metropolitan/Non-metropolitan",
       y="Life Expectancy",
       color="", fill="") +
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~census_region, scales="free")+
  isabel_theme


######## figure 9:  A4_le_bars4 ######## 
A4_le_bars4 <- ggplot(results_geography, aes(x=race, y=le, fill=as.factor(census_region))) + 
  geom_col(position=position_dodge(width=0.8), width=0.8) +
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(x="Census Region",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  isabel_theme