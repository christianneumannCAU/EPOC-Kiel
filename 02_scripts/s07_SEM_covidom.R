library(lavaan)

# read covidom data 
merged_clean <- read.csv("~/GitHub/EPOC/02_data/01_prep/covidom_data/covidom_subset.csv")

merged_scaled <- merged_clean

merged_scaled$tmt_a_time <- as.numeric(gsub(",", ".", merged_scaled$tmt_a_time))
merged_scaled$tmt_b_time <- as.numeric(gsub(",", ".", merged_scaled$tmt_b_time))
merged_scaled$pvt_reaction_time <- as.numeric(gsub(",", ".", merged_scaled$pvt_reaction_time))

merged_scaled[, c("moca", "tmt_a_time", "tmt_b_time", "pvt_reaction_time", "facit_f_FS", "age")] <- 
  scale(merged_scaled[, c("moca", "tmt_a_time", "tmt_b_time", "pvt_reaction_time", "facit_f_FS", "age")])


model <- "
  # fatigue through sex, age, education
  facit_f_FS ~ s1*sex + age + education_over_12

  # group through fatigue
  SCD ~ a*facit_f_FS + age + education_over_12

  # outcomes directly through fatigue 
  moca              ~ f1*facit_f_FS + age + education_over_12
  tmt_a_time        ~ f2*facit_f_FS + age + education_over_12
  tmt_b_time        ~ f3*facit_f_FS + age + education_over_12
  pvt_reaction_time ~ f4*facit_f_FS + age + education_over_12

  # covariance TMT A ~~ B
  tmt_a_time ~~ tmt_b_time

  # indirect effects sex → fatigue → outcome
  ind1 := s1*f1
  ind2 := s1*f2
  ind3 := s1*f3
  ind4 := s1*f4
"

fit_sem <- sem(model, data = merged_scaled, se = "bootstrap", bootstrap = 10000)
fitMeasures(fit_sem, c("cfi", "rmsea", "srmr"))
summary(fit_sem, standardized = TRUE, ci = TRUE)

