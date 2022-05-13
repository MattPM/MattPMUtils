



# write wrapper for chea3 
RunChea3 = function(gene.list){
#' @param gene.list - a named list 
#' @importFrom httr POST with_config config content
#' @importFrom jsonlite fromJSON 
#' @return 
#' @export
#'
#' @examples
#' '\dontrun{
#' res = list()
#   for (i in 1:length(li.joint)) {
#     gene.use = li.joint[i]
#     res[[i]] = RunChea3(gene.list = li.joint[i])
#   }
# }
 get.res.1 <- function() {
    print('if using this function cite Keenan et al Nucleic Acids Research (2019) doi.org/10.1093/nar/gkz446')
    payload = list(query_name = names(gene.list), gene_set = gene.list[[1]])
    response = POST(url = "https://maayanlab.cloud/chea3/api/enrich/",
                    body = payload, 
                    encode = "json")
    return(response)
  }
  get.response = function() {
    httr::with_config(
      config = httr::config(ssl_verifypeer = FALSE),
      get.res.1()
    )
  }
  response2 = get.response()
  json = content(response2, "text")
  results = fromJSON(json)
  return(results)
}


#' p.adjust.cormat adjust a matrix of p values from a correlation matrix 
#'
#' @param hmisc.cor - object returned by `Hmisc::rcorr()`
#' @param method - argument to p.adjust
#' @importFrom stats padjust
#'
#' @return
#' @export
#'
#' @examples
p.adjust.cormat = function(hmisc.cor, method = 'fdr'){ 
  stopifnot(isTRUE(isSymmetric(hmisc.cor$P)))
  p.adj =  p.adjust(hmisc.cor$P[lower.tri(hmisc.cor$P)], method = method)
  p.adj.mx <- matrix(rep(0,ncol(hmisc.cor$P)*ncol(hmisc.cor$P)), nrow = ncol(hmisc.cor$P))
  p.adj.mx[lower.tri(p.adj.mx)] <- p.adj
  p.adj.mx[upper.tri(p.adj.mx)] <- p.adj
  diag(p.adj.mx) = 1
  colnames(p.adj.mx) = rownames(p.adj.mx) = colnames(hmisc.cor$P)
  return(p.adj.mx)
}






# Visualizations

# add alpha transparency to color 
# taken from rethinking package
# https://github.com/rmcelreath/rethinking/
#' col.alpha - taken from Rethinking library - add alpha to color vector 
#'
#' @param col - color  
#' @param alpha - amount of transparency lower means more transparent 
#'
#' @return
#' @export
#'
#' @examples
col.alpha <- function( col , alpha=0.2 ) {
  col <- grDevices::col2rgb(col)
  col <- grDevices::rgb(col[1]/255,col[2]/255,col[3]/255,alpha)
  col
}








# blue red palette white in middle for low mid high diverging at mid point
blue.red = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3", 
             "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")



