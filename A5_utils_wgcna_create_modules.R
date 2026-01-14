library(ggplot2)
# input: sample x gene matrix

createModules <- function(dat, wgcnaRes) {
  genes <- colnames(dat)
  colors <- wgcnaRes$colors
  mods <- list()
  for (i in colors) mods[[i]] <- genes[colors == i]
  stopifnot(sort(table(colors)) == sort(unlist(lapply(mods, length))))

  df <- data.frame(mod = names(mods), size = as.numeric(lapply(mods, length)), stringsAsFactors = FALSE) %>%
    dplyr::mutate(prop = scales::percent(size / sum(size))) %>%
    dplyr::mutate(label = paste(mod, " (", prop, ")", sep = "")) %>%
    arrange(desc(size))

  color.values <- setNames(df$mod, df$label)
  p <- ggplot(df, aes(x = "", y = size, fill = label)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = color.values, limits = df$label) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    )
  return(list(p = p, mods = mods))
}
