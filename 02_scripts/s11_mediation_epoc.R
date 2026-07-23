library(readr)
library(mediation)
library(car)

# read data
merged_clean <- readr::read_tsv(
  "C:/Users/chris/Desktop/EPOC/EPOC-Kiel/01_data/00_bids/participants.tsv"
)

merged <- merged_clean

# numerical variables 
vars_num <- c(
  "facit_f_FS",
  "hads_a_total_score",
  "hads_d_total_score",
  "psqi_total_score",
  "pvt_reaction_time",
  "age",
  "years_of_education"
)

merged[vars_num] <- lapply(
  merged[vars_num],
  as.numeric
)

# mediator
model.m <- lm(
  facit_f_FS ~ group +
    age +
    sex +
    years_of_education +
    hads_a_total_score +
    hads_d_total_score +
    psqi_total_score,
  data = merged
)

# outcome
model.y <- lm(
  pvt_reaction_time ~
    group +
    facit_f_FS +
    age +
    sex +
    years_of_education +
    hads_a_total_score +
    hads_d_total_score +
    psqi_total_score,
  data = merged
)

# mediation
set.seed(500)

med.out <- mediate(
  model.m,
  model.y,
  treat = "group",
  mediator = "facit_f_FS",
  boot = TRUE,
  sims = 10000
)

summary(med.out)
summary(model.m)
summary(model.y)

vif(model.m)
vif(model.y)
