# load your normalized RRM file for your samples
data = as.data.frame(read.csv("RPMM_miRs_only.tsv", sep="\t", row.names = 1))
dim(data)
head(data)

# Re-name samples
names(data)
new_names = c()
for (x in names(data)) {
  first = strsplit(x, 'mouse.')[[1]][2]
  second = strsplit(first, '..TAB')[[1]][1]
  new_names <- append(new_names, second)
}
names(data) <- new_names
names(data)

# Get correlation values
data <- data.matrix(data)
head(data)
# get the correlation values for the data
data_corr <- cor(as.matrix(data))

# write correlation values to file
write.table(data_corr, file = "sample_correlation_values.tsv", sep = "\t", quote = F)

# make heat map using the correlation values
library(RColorBrewer)
library(gplots)
png("Sample_correlation_heatmap.png")
heatmap.2(data_corr, key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          main="Correlation Table",
          cexRow=1,
          cexCol=1,
          margins=c(12,12),
          srtCol=45)
dev.off()
