#' Identifies MRs for given data using stouffer integration.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param weights If included, will compute MRs using a weighted stouffer integration.
#' @return Returns the stouffer integrated scores for each protien.
StoufferMRs <- function(dat.mat, weights) {
  # generate dummy weights if missing
  if (missing(weights)) {
    weights = rep(1, ncol(dat.mat))
  }
  # stouffer integrate and return
  sInt <- rowSums(t(t(dat.mat) * weights)) / sqrt(sum(weights ** 2))
  return(sInt)
}

#' Identifies MRs based on ANOVA analysis for a given clustering.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering vector
#' @return A named vector of p-values for each protein
AnovaMRs <- function(dat.mat, clustering) {
  pVals <- c()
  group.vec <- clustering[colnames(dat.mat)]
  # perform an anova for each protein, storing pValues in a vector
  for (i in 1:nrow(dat.mat)) {
    aov.df <- data.frame('weights' = dat.mat[i,], 'group' = group.vec)
    #print(aov.df)
    aov.test <- aov(weights ~ group, aov.df)
    pVal <- summary(aov.test)[[1]][1,5]
    pVals <- c(pVals, pVal)
  }
  # name and return the vector
  names(pVals) <- rownames(dat.mat)
  return(pVals)
}

#' Identifies MRs based on a bootstraped Ttest between clusters.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Vector of cluster labels.
#' @param bootstrapNum Number of bootstraps to use. Default of 10 
#' @return Returns a list of lists; each list is a vector of sorted p-values for all proteins in the matrix.
BTTestMRs <- function(dat.mat, clustering, bootstrapNum = 100) {
  # set initial variables
  clustering <- clustering
  k <- length(table(clustering))
  mrs <- list()
  # identify MRs for each cluster
  for (i in 1:k) {
    clust <- names(table(clustering))[i]
    clust.vect <- which(clustering == clust)
    mList <- bootstrapTtest(dat.mat[, clust.vect ], dat.mat[, -clust.vect ], per = bootstrapNum)
    mrs[[clust]] <- sort(rowMeans(mList), decreasing = TRUE) # return mean p-value for the bootstraps
  }
  # return 
  return(mrs)
}

#' Returns the master regulators for the given data.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param method 'Stouffer' or 'ANOVA'
#' @param clustering Optional argument for a vector of cluster labels.
#' @param numMRs Number of MRs to return per cluster. Default of 50.
#' @param bottom Switch to return downregulated proteins in MR list. Default FALSE>
#' @param weights Optional argument for weights, which can be used in the Stouffer method.
#' @return Returns a list of master regulators, or a list of lists if a clustring is specified.
GetMRs <- function(dat.mat, clustering, method, numMRs = 50, bottom = FALSE, weights, ...) {
  if (method == 'ANOVA') {
    mr.vals <- AnovaMRs(dat.mat, clustering)
  } else if (method == 'Stouffer') {
    # generate dummy weights if not specified
    if (missing(weights)) {
      weights <- rep(1, ncol(dat.mat))
      names(weights) <- colnames(dat.mat)
    }
    # recursive calls for each cluster
    if (missing(clustering)) { # no clustering specified
      mr.vals <- StoufferMRs(dat.mat, weights)
    } else {
      k <- length(table(clustering))
      mrs <- list()
      for (i in 1:k) {
        # get cluster specific matrix and weights
        clust.cells <- names(which(clustering == i))
        clust.mat <- dat.mat[, clust.cells]
        clust.weights <- weights[clust.cells]
        # find mrs and add to list
        clust.mrs <- GetMRs(clust.mat, method = method, weights = clust.weights, numMRs = numMRs, bottom = bottom)
        mrs[[paste('c', i, sep = '')]] <- clust.mrs
      }
      return(mrs)
    }
  } else {
    print('Invalid method: must be "Stouffer" or "ANOVA".')
  }
  # return appropriate portion of MR list
  mr.vals <- sort(mr.vals, decreasing = TRUE)
  if (bottom) {
    return(c(mr.vals[1:numMRs], tail(mr.vals, numMRs)))
  } else {
    return(mr.vals[1:numMRs])
  }
}

#' Identifies MRs on a cell-by-cell basis and returns a merged, unique list of all such MRs.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param numMRs Default number of MRs to identify in each cell. Default of 25.
#' @return Returns a list of master regulators, the unique, merged set from all cells.
CBCMRs <- function(dat.mat, numMRs = 25) {
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:numMRs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # return
  return(cbc.mrs)
}

#' Make Cluster Metacells for ARACNe. Will take a clustering and produce saved meta cell matrices.
#'
#' @param dat.mat Matrix of raw gene expression (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor calculation. We recommend using a viper similarity matrix.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 5.
#' @param clustering Vector of cluster labels. 
#' @param subSize Size to subset the data too. Since 200 cells is adequate for ARACNe runs, this allows for speedup. Default of 200.
#' @param out.dir Directory for sub matrices to be saved in.
#' @param out.name Optional argument for preface of file names. 
MakeCMfA <- function(dat.mat, dist.mat, numNeighbors = 5, clustering, subSize = 200, out.dir, out.name = '') {
  # generate cluster matrices
  clust.mats <- ClusterMatrices(dat.mat, clustering)
  # produce metaCell matrix and save for each cluster matrix
  k <- length(clust.mats)
  meta.mats <- list()
  for (i in 1:k) {
    mat <- clust.mats[[i]]
    meta.mat <- MetaCells(mat, dist.mat, numNeighbors, subSize)
    meta.mat <- CPMTransform(meta.mat)
    file.name <- paste(out.dir, out.name, '_clust-', i, '-metaCells.tsv', sep = '')
    ARACNeTable(meta.mat, file.name, subset = FALSE)
    meta.mats[[i]] <- meta.mat
  }
  return(meta.mats)
}

#' Generates a meta cell matrix for given data.
#' 
#' @param dat.mat Raw gene expression matrix (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor inference.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 10.
#' @param subSize If specified, number of metaCells to be subset from the final matrix. No subsetting occurs if not incldued.
#' @return A matrix of meta cells (genes X samples).
MetaCells <- function(dat.mat, dist.mat, numNeighbors = 10, subSize) {
  # prune distance matrix if necessary
  dist.mat <- as.matrix(dist.mat)
  dist.mat <- dist.mat[colnames(dat.mat), colnames(dat.mat)]
  dist.mat <- as.dist(dist.mat)
  # KNN function
  KNN <- function(dist.mat, k){
    dist.mat <- as.matrix(dist.mat)
    n <- nrow(dist.mat)
    neighbor.mat <- matrix(0L, nrow = n, ncol = k)
    for (i in 1:n) {
      neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
    }
    return(neighbor.mat)
  }
  knn.neighbors <- KNN(dist.mat, numNeighbors)
  # create imputed matrix
  imp.mat <- matrix(0, nrow = nrow(dat.mat), ncol = ncol(dat.mat))
  rownames(imp.mat) <- rownames(dat.mat); colnames(imp.mat) <- colnames(dat.mat)
  for (i in 1:ncol(dat.mat)) {
    neighbor.mat <- dat.mat[,c(i, knn.neighbors[i,])]
    imp.mat[,i] <- rowSums(neighbor.mat)
  }
  # subset if requested and return
  if (missing(subSize)) {
    return(imp.mat)
  } else if (subSize > ncol(imp.mat)) {
    return(imp.mat)
  } else {
    return(imp.mat[, sample(colnames(imp.mat), subSize) ])
  }
}

#' Merges two viper matrices, giving priority to one over the other.
#' 
#' @param p.mat Priority viper matrix (proteins X samples). Proteins here will override those in the other matrix.
#' @param q.mat Secondary viper matrix (proteins X samples). Proteins here will fin in for gaps in the priority matrix.
#' @return A merged viper matrix.
ViperMerge <- function(p.mat, q.mat) {
  fill.genes <- setdiff(rownames(q.mat), rownames(p.mat))
  merged.mat <- rbind(p.mat, q.mat[fill.genes,])
  return(merged.mat)
}

#' Processes ARACNe results into a regulon object compatible with VIPER.
#'
#' @param a.file ARACNe final network .tsv.
#' @param exp.mat Matrix of expression from which the network was generated (genes X samples).
#' @param out.dir Output directory for networks to be saved to.
#' @param out.name Optional argument for prefix of the file name.
RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  require(viper)
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}

#' Performs a stouffer integration of a list of viper matrices.
#' 
#' @param vip.mats List of viper matrices (proteins X samples).
#' @param weights Vector of weights for the stouffer integration. If not included, all matrices are weighted equally.
#' @return An integrated viper matrix (proteins X samples).
VIPIntegrate <- function(vip.mats, weights) {
  # set weights, if not specified
  if (missing(weights)) {
    weights <- rep(1, length(vip.mats))
  }
  # identify set of regulators and create final matrix
  regs <- unique(Reduce(union, lapply(vip.mats, rownames)))
  int.mat <- matrix(0L, nrow = length(regs), ncol = ncol(vip.mats[[1]]))
  rownames(int.mat) <- regs; colnames(int.mat) <- colnames(vip.mats[[1]])
  # integrate each reg
  for (i in 1:length(regs)) {
    # generate matrix of all instances of this regulon
    reg.mat <- as.data.frame( lapply(vip.mats, function(x) {
      if (regs[i] %in% rownames(x)) {
        return( x[ regs[i] ,] )
      } else {
        return( rep(NA, ncol(x)) )
      }
    } ) )
    incl.vec <- which(!is.na(reg.mat[1,]))
    # integrata
    c.weights <- weights[incl.vec]
    sInt <- rowSums( t(t(reg.mat[, incl.vec]) * c.weights) ) / sqrt(sum(c.weights ** 2))
    int.mat[regs[i] ,] <- sInt
  }
  # return
  return(int.mat)
}

#' Performs VIPER analysis with all the GTEx Networks, then identifies and uses the top set for a metaVIPER analysis.
#' 
#' @param dat.mat Gene expression signature in matrix format (genes X samples) with ENSG ids.
#' @param gtex.path Path to the directory containing the GTEx networks.
#' @param num.nets Number of top networks to use. Default of 3.
#' @return A metaVIPER integration of the VIPER results for the top num.nets of GTEx networks (proteins X samples).
GTExVIPER <- function(dat.mat, gtex.path, num.nets = 3) {
  ## convert to entrez
  convert.dict <- readRDS(paste(gtex.path, 'gene-convert-dict.rds', sep = ''))
  dat.mat <- dat.mat[ which(rownames(dat.mat) %in% convert.dict$Ensembl.Gene.ID) ,] # remove rows with no ENSG match
  rname.match <- match(rownames(dat.mat), convert.dict$Ensembl.Gene.ID) # match remaining ENSG names
  entrez.names <- convert.dict$Entrez.Gene.ID[rname.match] # get Entrez names
  na.inds <- which(is.na(entrez.names)) # remove NA entrez names from matrix
  dat.mat <- dat.mat[ -na.inds ,]; entrez.names <- entrez.names[ -na.inds]
  rownames(dat.mat) <- entrez.names
  print(dim(dat.mat))
  ## load in all gtex networks
  net.files <- dir(gtex.path, pattern = '*.rda')
  gtex.nets <- list()
  for (i in 1:length(net.files)) {
    gtex.nets[[i]] <- get(load( paste(gtex.path, net.files[i], sep = '') ))
  }
  ## compute VIPER for all the networks
  viper.mats <- list()
  for (i in 1:length(gtex.nets)) {
    viper.mats[[i]] <- viper(dat.mat, gtex.nets[[i]], method = 'none')
  }
  ## identify the most important networks (as defined by num.nets)
  net.counts <- rep(0, length(gtex.nets)); names(net.counts) <- 1:length(gtex.nets)
  shared.regs <- Reduce(intersect, lapply(viper.mats, rownames))
  for (i in 1:length(shared.regs)) {
    reg <- shared.regs[i]
    reg.mat <- as.data.frame(lapply(viper.mats, function(x) { x[reg,] }))
    max.inds <- apply(reg.mat, 1, function(x) { which.max(abs(x)) } )
    for (j in 1:length(max.inds)) {
      net.counts[ max.inds[j] ] <- net.counts[ max.inds[j] ] + 1
    }
  }
  net.counts <- sort(net.counts, decreasing = TRUE)
  top.nets <- as.numeric(names(net.counts)[1:num.nets])
  ## clean up for memory purposes
  rm(viper.mats)
  ## integrate the results of the three selected networks
  mVip.mat <- viper(dat.mat, regulon = gtex.nets[top.nets], method = 'none')
  ## convert to ensemble
  mVip.mat <- mVip.mat[ which(rownames(mVip.mat) %in% convert.dict$Entrez.Gene.ID) ,]
  rname.match <- match(rownames(mVip.mat), convert.dict$Entrez.Gene.ID)
  ensg.names <- convert.dict$Ensembl.Gene.ID[rname.match]
  na.inds <- which(ensg.names == '')
  mVip.mat <- mVip.mat[ -na.inds ,]; ensg.names <- ensg.names[ -na.inds ]
  rownames(mVip.mat) <- ensg.names
  ## return
  return(mVip.mat)
}