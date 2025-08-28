#!/usr/bin/env Rscript

library(ggplot2)
library(grid)
library(optparse)
library(parallel)

option_list=list(
  make_option(c("-p","--proteinortho"), action="store", default=NA, type="character"),
  make_option(c("-i", "--iterations"), action="store", default=10, type="integer"),
  make_option(c("-d","--draw_only"), action="store", default=FALSE, type="logical"),
  make_option(c("-t","--table"), action="store", default=NA, type="character"),
  make_option(c("-o","--output"), action="store", default="proteinortho_curves", type="character"),
  make_option(c("--groups"), action="store", default=NA, type="character", 
              help="Two-column file: Species<TAB>Group. Species must match proteinortho columns."),
  make_option(c("--plot_width"), action="store", default=10, type="numeric"),
  make_option(c("--plot_height"), action="store", default=7, type="numeric"),
  make_option(c("--pan_color"), action="store", default="#82c6b8", type="character"),
  make_option(c("--core_color"), action="store", default="#c35b8e", type="character"),
  make_option(c("--plot_title"), action="store", default="Proteinortho Pan- and Core-Genome Size", type="character")
)

opt_parser=OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# check non-optional input
if (is.na(opt$proteinortho) && opt$draw_only==FALSE){
  stop("Proteinortho output is needed when -d FALSE", call.=FALSE)
}
if (is.na(opt$table) && opt$draw_only==TRUE){
  stop("proteinortho_curves.r output is needed when -d TRUE", call.=FALSE)
}

if (opt$draw_only==FALSE){
  iterations=as.integer(opt$iterations)
  
  # read in proteinortho output
  proteinortho=read.table(opt$proteinortho, header=F, sep="\t",stringsAsFactors=FALSE)
  table_all=proteinortho[,4:ncol(proteinortho)]
  
  # presence/absence matrix
  pa_matrix <- (table_all != "*") * 1  
  
  # If you have a groups file, read it.
  if (!is.na(opt$groups)){
    groups_df <- read.table(opt$groups, header=F, sep="\t", stringsAsFactors=FALSE)
    colnames(groups_df) <- c("Species","Group")
    
    # get col name(Assuming the column names correspond to species names)
    colnames(pa_matrix) <- groups_df$Species[match(1:ncol(pa_matrix), 1:ncol(pa_matrix))]
    
    group_levels <- unique(groups_df$Group)
    if (length(group_levels)!=2){
      stop("Currently only supports exactly two groups", call.=FALSE)
    }
    
    group_cols <- lapply(group_levels, function(g){
      species <- groups_df$Species[groups_df$Group==g]
      species[species %in% colnames(pa_matrix)]
    })
    names(group_cols) <- group_levels
  } else {
    group_levels <- "All"
    group_cols <- list(All=colnames(pa_matrix))
  }
  
  # detect cores
  ncores <- max(1, detectCores() - 1)
  cat("Using", ncores, "cores\n")
  
  # parallel wrapper
  do_iteration <- function(iteration){
    res_all <- list()
    
    for (grp in names(group_cols)){
      submat <- pa_matrix[, group_cols[[grp]], drop=FALSE]
      shuffled <- submat[, sample(ncol(submat)), drop=FALSE]
      
      res_iter <- data.frame(
        nrGenomes = sprintf("%03d", 1:ncol(shuffled)),
        sizeCore  = NA,
        sizePan   = NA,
        Group     = grp
      )
      for (n in 1:ncol(shuffled)){
        sub <- shuffled[,1:n, drop=FALSE]
        res_iter$sizeCore[n] <- sum(rowSums(sub) == n)
        res_iter$sizePan[n]  <- sum(rowSums(sub) > 0)
      }
      res_all[[grp]] <- res_iter
    }
    do.call(rbind, res_all)
  }
  
  # run parallel
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(ncores)
    results <- parLapply(cl, 1:iterations, do_iteration)
    stopCluster(cl)
  } else {
    results <- mclapply(1:iterations, do_iteration, mc.cores=ncores)
  }
  
  out <- do.call(rbind, results)
  
  write.table(out, file=paste0(opt$output,".txt"), sep="\t", row.names=F, quote=F)
  
} else if (opt$draw_only==TRUE){
  out=read.table(opt$table, sep="\t", header=T, stringsAsFactors=FALSE)
}

# plot
gg=ggplot(out, aes(x=nrGenomes, y=sizePan, fill=Group)) +
  geom_boxplot(alpha=0.5, position=position_dodge(width=0.8)) +
  geom_boxplot(aes(y=sizeCore), alpha=0.5, position=position_dodge(width=0.8)) +
  ggtitle(opt$plot_title) +
  theme_bw()

ggsave(filename=paste0(opt$output,".pdf"), scale=1, width=opt$plot_width, height=opt$plot_height)

cat("\nDone!\n")

