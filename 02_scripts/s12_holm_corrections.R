# Applying Holm-correction to all relevant p-values

# correlation analysis for association between trait fatigue and task-related cognitive performance
p <- c(0.008, 0.257, 0.442, 0.525)
p.adjust(p, method = "holm")

# Mann-whitney U tests for neurophysiological correlats of trait fatigue
p <- c(0.001, 0.004)
p.adjust(p, method = "holm")

# Mediation analysis for differences in task-related performance between participants with and without subjective cognitive symptoms
p <- c(0.0188, 0.2056)
p.adjust(p, method = "holm")
