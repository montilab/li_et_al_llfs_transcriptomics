library(Biobase)
library(WGCNA)
library(flashClust)
library(dynamicTreeCut)
library(doParallel)

# input: sample x gene matrix

sft.plot <- function(sft, powers) {
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = 1,
    col = "black"
  )
  abline(h = 0.85, col = "red")
}

sft.check <- function(sft) {
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- 6 # Default
    cat("Using the following power:", beta, "\n")
  } else {
    cat("Optimal power selected:", beta, "\n")
  }
  return(beta)
}

# Note both the adjacency matrix and topological overlap functions
# use `cor.type` but have their own sets of additional options.
# Please see their respective documentation if you need to use something
# beyond simply signed or unsigned.
wgcna.wrapper <- function(dat, # sample x gene matrix
                          min.size = 10,
                          min.sft = 0.85,
                          beta = NULL,
                          cores = 1,
                          cor.fn = c("cor", "bicor"),
                          powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                          merging = FALSE,
                          merging.cut = 0.2,
                          cor.type = c("unsigned", "signed"),
                          hclust.method = "average",
                          do.plot = TRUE) {
  # Handle arguments
  args <- as.list(environment())
  cor.fn <- match.arg(cor.fn)
  cor.type <- match.arg(cor.type)

  # Correlation options
  if (cor.fn == "cor") cor.options <- list(use = "p")
  if (cor.fn == "bicor") cor.options <- list(use = "p", pearsonFallback = "individual")

  # Set parallel computing environment
  doParallel::registerDoParallel(cores = cores)

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(
      data = dat,
      corFnc = cor.fn,
      networkType = cor.type,
      RsquaredCut = min.sft,
      powerVector = powers,
      corOptions = cor.options
    )

    if (do.plot) sft.plot(sft, powers)

    # Check selected power
    beta <- sft.check(sft)
  }

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(
    datExpr = dat,
    power = beta,
    corFnc = cor.fn,
    type = cor.type,
    corOptions = cor.options
  )

  # Topological overlap dissimilarity transformation
  dis <- WGCNA::TOMdist(adjMat = adj, TOMType = cor.type)

  # Fast hierarchical clustering of dissimilarity
  dendro <- flashClust::flashClust(d = as.dist(dis), method = hclust.method)

  # Module identification using dynamic tree cut algorithm
  modules <- dynamicTreeCut::cutreeDynamic(
    dendro = dendro,
    method = "hybrid",
    distM = dis,
    deepSplit = 4,
    pamRespectsDendro = FALSE,
    minClusterSize = min.size
  )

  # Assign module colours
  colors <- WGCNA::labels2colors(labels = modules, zeroIsGrey = TRUE)

  # Merging close modules
  if (merging) {
    merged <- WGCNA::mergeCloseModules(
      exprData = dat,
      colors = colors,
      corFnc = cor.fn,
      corOptions = cor.options,
      cutHeight = merging.cut
    )

    if (do.plot) {
      WGCNA::plotDendroAndColors(
        dendro = dendro,
        colors = cbind(colors, merged$colors),
        groupLabels = c("Original Modules", "Merged Modules"),
        dendroLabels = FALSE,
        addGuide = TRUE
      )
    }
    # Merged data
    colors <- merged$colors
    eigengenes <- merged$newMEs
  } else {
    eigengenes <- WGCNA::moduleEigengenes(expr = dat, colors = colors)$eigengenes

    if (do.plot) {
      WGCNA::plotDendroAndColors(
        dendro = dendro,
        colors = colors,
        groupLabels = "Modules",
        dendroLabels = FALSE,
        addGuide = TRUE
      )
    }
  }

  # Formatting
  colnames(eigengenes) <- substr(colnames(eigengenes), 3, 500)

  # Define modules
  genes <- colnames(dat)
  mods <- list()
  for (i in unique(colors)) {
    mods[[i]] <- genes[colors == i]
  }

  # Check
  stopifnot(sort(table(colors)) == sort(unlist(lapply(mods, length))))

  # For downstream analysis
  return(list(
    dat = dat,
    beta = beta,
    genes = genes,
    colors = colors,
    mods = mods,
    dendro = dendro,
    eigengenes = eigengenes,
    args = args
  ))
}
