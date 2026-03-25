library(lavaan)

merged_clean <- read.csv("~/GitHub/EPOC/02_data/01_prep/behavioral/behavioral_clean.csv")

merged_scaled <- merged_clean
merged_scaled[, c("moca", "tmt_a_time", "tmt_b_time", "pvt_reaction_time", "nback_miss_1", "facit_f_FS", "age", "years_of_education")] <- 
  scale(merged_scaled[, c("moca", "tmt_a_time", "tmt_b_time", "pvt_reaction_time", "nback_miss_1", "facit_f_FS", "age", "years_of_education")])


model <- "
  # Fatigue durch sex, age, education
  facit_f_FS ~ s1*sex + age + years_of_education

  # Gruppe durch fatigue
  group ~ a*facit_f_FS + age + years_of_education

  # Outcomes direkt durch fatigue (nicht durch group!)
  moca                ~ f1*facit_f_FS + age + years_of_education
  tmt_a_time        ~ f2*facit_f_FS + age + years_of_education
  tmt_b_time        ~ f3*facit_f_FS + age + years_of_education
  pvt_reaction_time ~ f4*facit_f_FS + age + years_of_education
  nback_miss_1      ~ f5*facit_f_FS + age + years_of_education

  # Kovarianzen   
  tmt_a_time ~~ tmt_b_time

  # Indirekte Effekte sex → fatigue → outcome
  ind1 := s1*f1
  ind2 := s1*f2
  ind3 := s1*f3
  ind4 := s1*f4
  ind5 := s1*f5
"

fit_sem <- sem(model, data = merged_scaled, se = "bootstrap", bootstrap = 10000)
fitMeasures(fit_sem, c("cfi", "rmsea", "srmr"))
summary(fit_sem, standardized = TRUE, ci = TRUE)

