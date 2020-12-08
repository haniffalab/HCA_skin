library(ggplot2)
a <- read.csv("../pct_inf_df_proportion.csv", header=T, stringsAsFactors = F)

a$Celltype <- factor(a$Celltype, levels = c(
  "post_inf",
  "post",
  "pro",
  "pre"
  )) # specify order

a$Status_site <- factor(a$Status_site, levels = c(
  "HS_NL",
  "E_NL",
  "E_L",
  "P_NL",
  "P_L"
  ))

colour_celltypes <- c(
    '#E87D72',
    '#0e6c8b',
    '#87AC34',
    '#b8bae5'
    )

p <- ggplot(data=a, aes(x=Status_site, y=Proportion, fill=Celltype)) +
  geom_bar(stat="identity", position="fill",color="black") +
  theme(axis.line =element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),strip.text.x = element_text(size = 12)) +
  theme_bw() +
  scale_fill_manual(values=colour_celltypes) +
  theme(aspect.ratio = 2/1) +
  theme(axis.text.x = element_text(angle = 90)) +
  # theme(axis.text.x=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank()
  )

pdf("KC_proportions.pdf", width = 3, height = 3)
p
dev.off()
