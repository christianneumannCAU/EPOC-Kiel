library(readr)
library(ggplot2)

# read participants 
merged_clean <- readr::read_tsv("C:/Users/chris/Desktop/EPOC/01_data/00_bids/participants.tsv")

# number of NA per column
colSums(is.na(merged_clean))

# pick variables
corr_vars <- c(
  "age", "years_of_education",
  "facit_f_FS","hads_d_total_score","hads_a_total_score","psqi_total_score",
  "pvt_reaction_time","nback_miss_1","tmt_a_time","tmt_b_time","moca"
)

merged_clean[corr_vars] <- lapply(merged_clean[corr_vars], function(x) as.numeric(as.character(x)))

# calculate spearman correlations
ct2 <- psych::corr.test(merged_clean[, corr_vars],
                        method = "spearman", adjust = "none")

# give out correlations
R2 <- ct2$r
P2 <- ct2$p

ct2$r["facit_f_FS", "pvt_reaction_time"]    # Spearman rho
ct2$p["facit_f_FS", "pvt_reaction_time"]    # p value
ct2$r["facit_f_FS", "hads_d_total_score"]   # Spearman rho
ct2$p["facit_f_FS", "hads_d_total_score"]   # p value
ct2$r["facit_f_FS", "hads_a_total_score"]   # Spearman rho
ct2$p["facit_f_FS", "hads_a_total_score"]   # p value
ct2$r["facit_f_FS", "psqi_total_score"]     # Spearman rho
ct2$p["facit_f_FS", "psqi_total_score"]     # p value
ct2$r["hads_d_total_score", "hads_a_total_score"]   # Spearman rho
ct2$p["hads_d_total_score", "hads_a_total_score"]   # p value
ct2$r["hads_d_total_score", "psqi_total_score"]   # Spearman rho
ct2$p["hads_d_total_score", "psqi_total_score"]   # p value
ct2$r["hads_a_total_score", "psqi_total_score"]   # Spearman rho
ct2$p["hads_a_total_score", "psqi_total_score"]   # p value




# Long format
dat2 <- merge(
  reshape2::melt(R2, varnames = c("Var1","Var2"), value.name = "rho"),
  reshape2::melt(P2, varnames = c("Var1","Var2"), value.name = "p"),
  by = c("Var1","Var2")
)

# only get one half without the diagonal
dat2 <- dat2[as.numeric(factor(dat2$Var1, levels = corr_vars)) >
               as.numeric(factor(dat2$Var2, levels = corr_vars)), ]

# highlight significances 
dat2$label <- paste0(sprintf("%.2f", dat2$rho), ifelse(dat2$p < (.05/55), "*", ""))

# rename for plot
var_labels <- c(
  "age" = "Age",
  "years_of_education" = "Years of Education",
  "facit_f_FS" = "FACIT-Fatigue Subscale",
  "hads_d_total_score" = "HADS-D Total Score",
  "hads_a_total_score" = "HADS-A Total Score",
  "psqi_total_score" = "PSQI Total Score",
  "pvt_reaction_time" = "PVT Reaction Time",
  "nback_miss_1" = "1-Back Misses",
  "tmt_a_time" = "TMT A Time",
  "tmt_b_time" = "TMT B Time",
  "moca" = "MoCA Total score"
)


# Plot
ggplot(dat2, aes(Var2, Var1, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_x_discrete(labels = var_labels) +
  scale_y_discrete(labels = var_labels) +
  scale_fill_gradient2(
    low  = RColorBrewer::brewer.pal(11, "RdBu")[11],
    mid  = "white",
    high = RColorBrewer::brewer.pal(11, "RdBu")[1],
    midpoint = 0, limits = c(-1, 1), name = expression(rho)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)

ggsave("C:/Users/chris/Desktop/EPOC/03_figures/06_heatmaps/behavioral.png", width = 7, height = 5, units = "in", dpi = 300)
