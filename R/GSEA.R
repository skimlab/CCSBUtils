#' A function to read GMT file
#'
#' @param gmt.file a text file in 'gmt' format.
#' @param min.size The minimum gene set size to read.  Defaults to 0.
#' @param max.size The maximum gene set size to read.  Defaults to no limit.
#' @param simplify Defaults to TRUE.  See 'output' for detail.
#' @param sep1 a character string separating 'gene set name', 'url' and 'gene set'
#' @param sep2 a character string separating each gene in 'gene set'
#' @return a list of gene set.
#'        If TRUE, the return structure is simplified to a list of genesets
#'        each of which is a list of genes and the name of each geneset is
#'        the name of gene set from MSigDB file, the first element of each line.
#'        If FALSE, it returns a list of three element list each element of which
#'        is \code{name}, \code{url}, and \code{geneset}.
#' @export
read_gmt <- function(gmt.file, min.size = 1, max.size = -1, simplify = TRUE, sep1 = "\t", sep2 = sep1) {
  parse_one_line_gmt <-
    function(one.line,
             split1 = "\t",
             split2 = split1) {
      gset <- unlist(strsplit(as.character(one.line), split = split1))

      if (split1 != split2) {
        gs <- trimws(unlist(strsplit(gset[3], split = split2)))
      } else {
        gs <- gset[-c(1:2)]
      }

      return(list(
        name = gset[1],
        url = gset[2],
        gene.set = gs
      ))
    }

  gs.temp <- readLines(con.gmt <- file(gmt.file))
  close(con.gmt)
  gene.sets <- lapply(gs.temp, parse_one_line_gmt, split1 = sep1, split2 = sep2)

  list.name <- array("", length(gene.sets))

  # convert to named list -- this will make the access to the list much easier.
  for (ii in 1:length(gene.sets)) {
    list.name[ii] <- gene.sets[[ii]]$name
  }
  names(gene.sets) <- list.name

  # Prune out ones that are too long/short
  for (ii in length(gene.sets):1) {
    if (length(gene.sets[[ii]]$gene.set) < min.size ||
        (max.size > 0 && length(gene.sets[[ii]]$gene.set) > max.size))
      gene.sets[[ii]] <- NULL
  }

  if (simplify) {
    gene.sets <- sapply(gene.sets, "[", "gene.set")
    names(gene.sets) <- sub(".gene.set$", "", names(gene.sets))
  }

  # rid of empty list
  return(gene.sets[!is.na(gene.sets)])
}
