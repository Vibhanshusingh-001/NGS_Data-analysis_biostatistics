# load package 
library(ggplot2)

variant_data <- read.csv("variant_metrics.csv")
# Calculate overall Ti/Tv ratio
ti_count <- sum(variant_data$Ti.Tv == "Transition")
tv_count <- sum(variant_data$Ti.Tv == "Transversion")
ti_tv_ratio <- ti_count / tv_count
cat("Transition/Transversion Ratio:", ti_tv_ratio, "\n")

# Plot VAF distribution
ggplot(variant_data, aes(x = AF)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Variant Allele Frequency Distribution", x = "VAF", y = "Count")
