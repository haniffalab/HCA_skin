library(ggplot2)
a <- read.csv("../pct_exp_only_df.csv", header=T, stringsAsFactors = F)
a$Celltype <- factor(a$Celltype, levels = c(
  "Tc",
  "Tc IL13 IL22",
  "Tc17_Th17",
  "Th",
  "Treg"
  )) # specify order
a$Status_site <- factor(a$Status_site, levels = c(
  "Healthy_non_lesion",
  "Eczema_non_lesion",
  "Eczema_lesion",
  "Psoriasis_non_lesion",
  "Psoriasis_lesion"
  ))

colour_celltypes <- c(
    '#a3dc82',  # Tc
    'red',
    'black',  # Tc17_Th17
    '#a02b9d',
    '#B482D6'
    )

ggplot(data=a, aes(fill=condition, y=Proportion, x=specie)) +
    geom_bar(position="dodge", stat="identity")

p <- ggplot(data=a, aes(x=Status_site, y=Proportion, fill=Celltype)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  theme(axis.line=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 12)) +
  theme_bw() +
  scale_fill_manual(values=colour_celltypes) +
  # theme(aspect.ratio = 2/1) +
  theme(axis.text.x = element_text(angle = 90)) +
  # theme(axis.text.x=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.y=element_blank(),
        # axis.text.y=element_blank()
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

p <- p + ylab("Expanded clonotype [%]")
p <- p + ylim(0, 100)

pdf("lymphoid_expanded_proportions_side.pdf", width = 8, height = 6)
p
dev.off()
