#' @param samplefiles
#' @param postMax.out
#' @param contexts
#' @param input.dir
#' @param out.dir
#' @param magrittr
#' @importFrom data.table fread fwrite
#' @importFrom dplyr inner_join
#' @export
#'
split.groups <- function(samplefiles, postMax.out=FALSE, contexts=c("CG","CHG","CHH"), input.dir, out.dir){
  ft <- data.table::fread(samplefiles)
  ft$name <- paste0(ft$sample,"_", ft$replicate)
  gps <- ft$group[!ft$group %in% c('control')]
  gps <- unique(gps)
  for (m in seq_along(gps)){
    myvec <- c("control", gps[m])
    gp1 <- ft$name[which(ft$group==myvec[1])]
    gp2 <- ft$name[which(ft$group==myvec[2])]
    gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
    gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
    out.name <- paste0(gp1.sample, "_", gp2.sample)
    for (cn in seq_along(contexts)){
      fn1 <- paste0(input.dir, '/', contexts[cn], "_StateCalls.txt")
      if (file.exists(fn1)) {
        StateCall <- data.table::fread(fn1)
        df1 <- subset(StateCall,, which(colnames(StateCall) %in% c('seqnames', 'start', 'end', gp1, gp2)))
        data.table::fwrite(x=df1,
                           file=paste0(out.dir, "/", contexts[cn], "_", out.name, "_StateCalls.txt"),
                           quote=FALSE,
                           row.names=FALSE,
                           col.names=TRUE,
                           sep="\t")
      }

      fn2 <- paste0(input.dir, '/', contexts[cn], "_rcMethlvl.txt")
      if (file.exists(fn2)) {
        rcMethlvl <- data.table::fread(fn2)
        df2 <- subset(rcMethlvl,, which(colnames(rcMethlvl) %in% c('seqnames', 'start', 'end', gp1, gp2)))
        data.table::fwrite(x=df2,
                           file=paste0(out.dir, "/", contexts[cn],"_", out.name, "_rcMethlvl.txt"),
                           quote=FALSE,
                           row.names=FALSE,
                           col.names=TRUE,
                           sep="\t")
      }
      if (postMax.out==TRUE){
        fn3 <- paste0(input.dir, '/', contexts[cn], "_postMax.txt")
        if (file.exists(fn3)) {
          postMax <- data.table::fread(fn3)
          df3 <- subset(postMax,, which(colnames(postMax) %in% c('seqnames', 'start', 'end', gp1, gp2)))
          data.table::fwrite(x=df3,
                             file=paste0(out.dir, "/", contexts[cn], "_", out.name, "_postMax.txt"),
                             quote=FALSE,
                             row.names=FALSE,
                             col.names=TRUE,
                             sep="\t")
        }
      }
    }
  }
}
