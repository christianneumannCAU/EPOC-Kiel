# Applying Holm-correction to all relevant p-values


# Mann-whitney U tests for neurophysiological correlats of trait fatigue
p <- c(0.002, 0.038)
p.adjust(p, method = "holm")

