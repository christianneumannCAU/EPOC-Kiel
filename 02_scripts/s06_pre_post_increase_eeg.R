library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(psych)
library(tidyr)
library(ggtext)

merged_clean <- readr::read_tsv("C:/Users/chris/Desktop/EPOC/EPOC-Kiel/01_data/00_bids/participants.tsv")

eeg_vars  <- c("theta_peak_power_Pre","alpha_peak_power_Pre",
               "theta_peak_power_Post","alpha_peak_power_Post")
test_vars <- c("facit_f_FS","hads_d_total_score","hads_a_total_score","psqi_total_score",
               "pvt_reaction_time","nback_miss_1","tmt_a_time","tmt_b_time","moca")

merged_clean[eeg_vars] <- lapply(merged_clean[eeg_vars], function(x) as.numeric(as.character(x)))
merged_clean[test_vars] <- lapply(merged_clean[test_vars], function(x) as.numeric(as.character(x)))

# remove participants with missing EEG data
merged_clean <- merged_clean[!is.na(merged_clean$alpha_peak_power_Pre),]

########################## outlier analyses ###################################

# Scatterplots
plot(merged_clean$alpha_peak_power_Pre, 
     main = "Scatterplot of Alpha Peak Power (Pre)", 
     xlab = "Index", 
     ylab = "Alpha Peak Power (Pre)", 
     pch = 16,  
     col = "blue")
plot(merged_clean$alpha_peak_power_Post, 
     main = "Scatterplot of Alpha Peak Power (Post)", 
     xlab = "Index", 
     ylab = "Alpha Peak Power (Post)", 
     pch = 16,  
     col = "blue")
plot(merged_clean$theta_peak_power_Pre, 
     main = "Scatterplot of Theta Peak Power (Pre)", 
     xlab = "Index", 
     ylab = "Theta Peak Power (Pre)", 
     pch = 16,  
     col = "blue")
plot(merged_clean$theta_peak_power_Post, 
     main = "Scatterplot of Theta Peak Power (Post)", 
     xlab = "Index", 
     ylab = "Theta Peak Power (Post)", 
     pch = 16,  
     col = "blue")



pre_lim_a <- mean(merged_clean$alpha_peak_power_Pre) + 5 * sd(merged_clean$alpha_peak_power_Pre)
post_lim_a <- mean(merged_clean$alpha_peak_power_Post) + 5 * sd(merged_clean$alpha_peak_power_Post)
pre_lim_t <- mean(merged_clean$theta_peak_power_Pre) + 5 * sd(merged_clean$theta_peak_power_Pre)
post_lim_t <- mean(merged_clean$theta_peak_power_Post) + 5 * sd(merged_clean$theta_peak_power_Post)

original <- merged_clean

to_remove <- which(
  original$alpha_peak_power_Pre  >= pre_lim_a |
    original$alpha_peak_power_Post >= post_lim_a |
    original$theta_peak_power_Pre  >= pre_lim_t |
    original$theta_peak_power_Post >= post_lim_t
)

to_remove

merged_clean <- merged_clean[merged_clean$alpha_peak_power_Pre  < pre_lim_a, ]
merged_clean <- merged_clean[merged_clean$alpha_peak_power_Post < post_lim_a, ]
merged_clean <- merged_clean[merged_clean$theta_peak_power_Pre  < pre_lim_t, ]
merged_clean <- merged_clean[merged_clean$theta_peak_power_Post < post_lim_t, ]

####################################################################

# colors und labels
merged_clean <- merged_clean %>%
  filter(group_Fatigue != "n/a")

group_colors <- c(
  "clinically significant fatigue" = "#CC79A7",
  "no clinically significant fatigue" = "#82325E"
)

###############################################################################
# THETA

theta_long <- merged_clean %>%
  select(group_Fatigue, theta_peak_power_Pre, theta_peak_power_Post) %>%
  pivot_longer(
    cols = starts_with("theta_peak_power"),
    names_to = "Time",
    values_to = "ThetaPower"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c("theta_peak_power_Pre", "theta_peak_power_Post"),
      labels = c("Pre-Stimulus", "Post-Stimulus")
    ),
    group = factor(
      group_Fatigue,
      levels = c(
        "clinically relevant fatigue",
        "no clinically relevant fatigue"
      ),
      labels = c(
        "clinically significant fatigue",
        "no clinically significant fatigue"
      )
    )
  )

theta_summary <- theta_long %>%
  group_by(group, Time) %>%
  summarise(
    mean = mean(ThetaPower, na.rm = TRUE),
    se = sd(ThetaPower, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_theta <- ggplot(theta_summary,
                  aes(Time, mean,
                      group = group,
                      color = group,
                      linetype = group)) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "B",
    hjust = -0.5,
    vjust = 1.5,
    size = 6,
    fontface = "bold"
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = .1, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(x = NULL,
       y = "<span style='color:#E69F00;'>Theta Peak Power</span> (µV²)",
       color = NULL,
       linetype = NULL) +
  scale_x_discrete(labels = c(
    "Pre-Stimulus" = "<span style='color:#006047;'>Pre-Stimulus</span>",
    "Post-Stimulus" = "<span style='color:#009E73;'>Post-Stimulus</span>"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(),
    axis.title.y = element_markdown()
  ) +
  coord_cartesian(ylim = c(0, 0.06))

###############################################################################
# ALPHA

alpha_long <- merged_clean %>%
  select(group_Fatigue, alpha_peak_power_Pre, alpha_peak_power_Post) %>%
  pivot_longer(
    cols = starts_with("alpha_peak_power"),
    names_to = "Time",
    values_to = "AlphaPower"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c("alpha_peak_power_Pre", "alpha_peak_power_Post"),
      labels = c("Pre-Stimulus", "Post-Stimulus")
    ),
    group = factor(
      group_Fatigue,
      levels = c(
        "clinically relevant fatigue",
        "no clinically relevant fatigue"
      ),
      labels = c(
        "clinically significant fatigue",
        "no clinically significant fatigue"
      )
    )
  )

alpha_summary <- alpha_long %>%
  group_by(group, Time) %>%
  summarise(
    mean = mean(AlphaPower, na.rm = TRUE),
    se = sd(AlphaPower, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_alpha <- ggplot(alpha_summary,
                  aes(Time, mean,
                      group = group,
                      color = group,
                      linetype = group)) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "A",
    hjust = -0.5,
    vjust = 1.5,
    size = 6,
    fontface = "bold"
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = .1, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(x = NULL,
       y = "<span style='color:#C08500;'>Alpha Peak Power</span> (µV²)",
       color = NULL,
       linetype = NULL) +
  scale_x_discrete(labels = c(
    "Pre-Stimulus" = "<span style='color:#006047;'>Pre-Stimulus</span>",
    "Post-Stimulus" = "<span style='color:#009E73;'>Post-Stimulus</span>"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(),
    axis.title.y = element_markdown()
  ) +
  coord_cartesian(ylim = c(0, 0.035))

combined_plot <- (p_alpha + p_theta) + plot_layout(guides = "collect") &
  theme(legend.position = "top")


combined_plot

ggsave(
  "C:/Users/chris/Desktop/EPOC/EPOC-Kiel/03_figures/06_eeg_peaks/pre_post_alpha_theta.pdf",
  plot = combined_plot,
  width = 210,
  height = 180,
  units = "mm",
  dpi = 600
)
