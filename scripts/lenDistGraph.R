# Load libraries
library(ggplot2)
library(reshape2)

# Get command line argument
args <- commandArgs(TRUE)

# Load Data
lenData <- read.csv(args[1], header = T, sep = "\t")

# Convert data for ggplot input
bit <- melt(lenData, id.vars = "X")

# Determine number of samples per role in output (based on total # of samples)
wrap = round(length(lenData) * .33)

# Write output to lenDistHistogram.png
png("lenDistHistogram.png")
ggplot(bit, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat="identity") +
  facet_wrap(~variable, ncol = wrap) +
  scale_fill_discrete(guide=FALSE) +
  labs(list(title="Length Distribution", x="Read length", y="Ratio of total reads")) 
dev.off()
