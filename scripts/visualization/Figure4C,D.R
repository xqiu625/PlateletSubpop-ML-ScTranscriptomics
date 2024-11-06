library(ggplot2)
library(dplyr)

df <- read.csv("Supp_Table2.csv")
df$Clusters <-  factor(df$Clusters,levels= c("C0", "C1", "C2", "C3", "C4", "C5", "C6", 
                                             "C7","C8", "C9", "C10", "C11", "C12"))


dpi = 300
png(file = "Fig4D_Severity_cluster_density.png", width = 6, height = 4, units = 'in',res = 300)
ggplot(df, aes(fill=Severity, y=Frequency, x=Clusters)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = 
                      c('#8dd3c7','#ffffb3','#bebada','#fdb462','#80b1d3','#fb8072','#b3de69')
                    ) + 
  theme(axis.text=element_text(size=16, face="bold")) + xlab("Cluster") + ylab("Percentage") + 
  guides(fill=guide_legend(title="Severity")) + 
  theme_minimal()
invisible(dev.off())



df <- read.csv("Supp_Table3.csv")
df$Outcome <- factor(df$Outcome,
                            levels = c("Unknown", "HC", "S", "FT"),
                            labels = c("Unknown", "Healthy Control", "Survivor", "Fatal"))
df$Clusters <-  factor(df$Clusters,levels= c("C0", "C1", "C2", "C3", "C4", "C5", "C6", 
                                                           "C7","C8", "C9", "C10", "C11", "C12"))

dpi = 300
png(file = "Fig4D_Outcome_cluster_density.png", width = 6, height = 4, units = 'in',res = 300)
ggplot(df, aes(fill=Outcome, y=Frequency, x=Clusters)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c('#FFD16F','#A1CEED', '#456681', '#ad3c53')) + 
  theme(axis.text=element_text(size=16, face="bold")) + xlab("Cluster") + ylab("Percentage") + 
  guides(fill=guide_legend(title="Outcome")) + 
  theme_minimal()
invisible(dev.off())
