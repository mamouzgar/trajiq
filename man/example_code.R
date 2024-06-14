
library(tidyverse)
library(trajiq)

## load example data
glm_input=trajiq::glm_input

## features of interest
my_features = colnames(glm_input)[-c(1:3)]
covariates_in_model = c('mouse') ## covariates of interest
contrast_variables = c('tissue') ## variable you want to compare

###############################################################################################
## (A) example for running analysis on glm_input features that look like celltype abundance fractions
###############################################################################################
## note: be sure to make covariates as factor variables when it's not a continuous variable. Eg, mouse ids as a covariate for matched analysis
glm_input = glm_input%>%
     mutate(mouse = factor(mouse))
dif_res =trajiq::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE)
head(dif_res)

###############################################################################################
## (B) to run interactions, you would run something like this:
###############################################################################################
contrast_variables = c('tissue','ed')
dif_res =trajiq::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE)

###############################################################################################
## (C) example for running analysis on glm_input features that look like marker expression
###############################################################################################
glm_input =trajiq::ct %>% mutate(mouse = factor(mouse)) %>%
     bind_rows(trajiq::ct %>% mutate(mouse = factor(mouse)))
contrast_variables = c('tissue')
my_features = c('CD4','CD40','FOXP3') ## make sure these are in the object (no spelling errors)
# my_features %in% colnames(glm_input)
dif_res = trajiq::differential_analysis_program(glm_input = glm_input,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE)

###############################################################################################
# (D) example code to split across multiple subsets of interest
###############################################################################################
dif_res_by_ct = trajiq::differential_analysis_program(glm_input = x,outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = 'patient_id',intercept = TRUE,
                                                      SPLIT_BY_NAMES = c('annotation_lin_sub','annotation_subset'))












