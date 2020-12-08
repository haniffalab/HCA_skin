library(ggplot2)
shi_df <- read.csv("../shi_df_all_t_cells.csv", row.names = 1, stringsAsFactors = F)

shi_df$Status_Site <- factor(shi_df$Status_Site, levels = c(
  "Healthy_non_lesion",
  "Eczema_non_lesion",
  "Eczema_lesion",
  "Psoriasis_non_lesion",
  "Psoriasis_lesion"
  ))


plot_list <- list()
for(i_cell in 1:length(unique(shi_df[["celltype"]]))) {
  celltype <- unique(shi_df[["celltype"]])[i_cell]

  shi_df_sub <- shi_df
  shi_df_sub <- shi_df_sub[shi_df_sub["celltype"] == celltype, ]

  col_plot <- as.character(shi_df_sub[["colour"]][1])

  p <- ggplot(
    shi_df_sub,
    aes(x=Status_Site, y=shi)) +

    geom_dotplot(
      binaxis='y',
      stackdir='center',
      colour=col_plot,
      fill=col_plot,
      dotsize=0.7
    ) +
    # scale_color_manual(values = col_plot)
    stat_summary(fun.y=mean, geom="point", shape="-",
                   size=10, color=shi_df_sub[["colour"]][1]) +
    theme_bw() +
    # theme(axis.title.x=element_blank(),
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank(),
      # axis.title.y=element_blank(),
      # axis.text.y=element_blank(),
      # axis.ticks.y=element_blank()) +

    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +

    ylim(3.5, 6.5)

  # p + stat_summary(
  #       fun.data=mean_sdl,
  #       fun.args = list(mult=1),
  #       geom="pointrange",
  #       color=col_plot)

  p <- p + stat_summary(
        fun.data="mean_sdl",
        fun.args = list(mult=1),
        geom="crossbar",
        width=0.75,
        color=col_plot)

  plot(p)
  ggsave(filename = paste("Figure4G_dotplot_", celltype, ".png", sep = ""), width = 1.5, height = 3, units = "in")

  # ggsave(filename = paste("Figure4G_dotplot_", celltype, ".ps", sep = ""), width = 1.5, height = 2, units = "in")

  ggsave(filename = paste("Figure4G_dotplot_narrow_", celltype, ".png", sep = ""), width = 1, height = 3, units = "in")

  # ggsave(filename = paste("Figure4G_dotplot_narrow_", celltype, ".ps", sep = ""), width = 1, height = 2, units = "in")

}
