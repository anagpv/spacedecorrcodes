

get_smoothed <- function(sample_coords,Y,distneighb){
  num_cells <- nrow(sample_coords)

  expr_values <- as.matrix(t(Y))
  #colnames(expr_values) <- 1:num_cells

  distmat <- as.matrix(dist(sample_coords, diag = T))
  distances <- reshape2::melt(distmat) %>% filter(value < distneighb) %>% dplyr::rename(from=Var1, to=Var2)
  sample_coords2 <- sample_coords %>% mutate(from = rownames(sample_coords),
                                             to = rownames(sample_coords))
  nb <- merge(distances, sample_coords2[,c("sdimx","sdimy","from")], by = "from")
  nb <- as.data.table(merge(nb,  sample_coords2[,c(1,2,4)], by = "to", suffixes=c("_from", "_to")))
  nb$to <- as.factor(nb$to)
  nb$from <- as.factor(nb$from)

  smooth <- do_spatial_knn_smoothing(spatial_network = nb,expr_values=expr_values)
  smooth
}



convert_to_full_spatial_network =  function(reduced_spatial_network_DT) {

  # data.table variables
  distance = rank_int = NULL

  # find location coordinates
  coordinates = grep('sdim', colnames(reduced_spatial_network_DT), value = T)

  # create normal source --> target
  part1 = data.table::copy(reduced_spatial_network_DT)
  part1 = part1[, c('from', 'to'), with = F]
  colnames(part1) = c('source', 'target')

  # revert order target (now source) --> source (now target)
  part2 = data.table::copy(reduced_spatial_network_DT[, c('to', 'from'), with = F])
  colnames(part2) = c('source', 'target')

  # combine and remove duplicates
  full_spatial_network_DT = rbind(part1, part2)
  full_spatial_network_DT = unique(full_spatial_network_DT)

  # create ranking of interactions
  data.table::setorder(full_spatial_network_DT, source)
  full_spatial_network_DT[, rank_int := 1:.N, by = 'source']

  # create unified column
  #full_spatial_network_DT = sort_combine_two_DT_columns(full_spatial_network_DT, 'source', 'target', 'rnk_src_trgt')

  return(full_spatial_network_DT)

}
dt_to_matrix <- function(x) {
  rownames = as.character(x[[1]])
  mat = methods::as(Matrix::as.matrix(x[,-1]), 'Matrix')
  rownames(mat) = rownames
  return(mat)
}
do_spatial_knn_smoothing = function(spatial_network,expr_values,
                                    subset_genes = NULL,
                                    spatial_network_name = 'Delaunay_network',
                                    b = NULL) {

  # checks
  if(!is.null(b)) {
    if(b > 1 | b < 0) {
      stop('b needs to be between 0 (no spatial contribution) and 1 (only spatial contribution)')
    }
  }

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }

  # data.table variables
  gene_ID = value = NULL

  # merge spatial network with expression data
  expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
  expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')

  ## test ##
  spatial_network = convert_to_full_spatial_network(spatial_network)
  ## stop test ##

  #print(spatial_network)

  spatial_network_ext = data.table::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)

  #print(spatial_network_ext)

  # calculate mean over all k-neighbours
  # exclude 0's?
  # trimmed mean?
  spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]

  # convert back to matrix
  spatial_smooth_dc = data.table::dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
  spatial_smooth_matrix = dt_to_matrix(spatial_smooth_dc)


  # if network was not fully connected, some cells might be missing and are not smoothed
  # add the original values for those cells back
  all_cells = colnames(expr_values)
  smoothed_cells = colnames(spatial_smooth_matrix)
  missing_cells = all_cells[!all_cells %in% smoothed_cells]
  if(length(missing_cells) > 0) {
    missing_matrix = expr_values[, missing_cells]
    spatial_smooth_matrix = cbind(spatial_smooth_matrix[rownames(expr_values),], missing_matrix)
  }

  spatial_smooth_matrix = spatial_smooth_matrix[rownames(expr_values), colnames(expr_values)]

  # combine original and smoothed values according to smoothening b
  # create best guess for b if not given
  if(is.null(b)) {
    k = stats::median(table(spatial_network$source))
    smooth_b = 1 - 1/k
  } else {
    smooth_b = b
  }

  expr_b = 1 - smooth_b
  spatsmooth_expr_values = ((smooth_b*spatial_smooth_matrix) + (expr_b*expr_values))

  return(t(as.matrix(spatsmooth_expr_values)))

}

binSpectSingle = function(gobject,
                          bin_method = c('kmeans', 'rank'),
                          expression_values = c('normalized', 'scaled', 'custom'),
                          subset_genes = NULL,
                          spatial_network_name = 'Delaunay_network',
                          reduce_network = FALSE,
                          kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                          nstart = 3,
                          iter_max = 10,
                          extreme_nr = 50,
                          sample_nr = 50,
                          percentage_rank = 30,
                          do_fisher_test = TRUE,
                          adjust_method = 'fdr',
                          calc_hub = FALSE,
                          hub_min_int = 3,
                          get_av_expr = TRUE,
                          get_high_expr = TRUE,
                          implementation = c('data.table', 'simple', 'matrix'),
                          group_size = 'automatic',
                          do_parallel = TRUE,
                          cores = NA,
                          verbose = T,
                          set.seed = NULL,
                          bin_matrix = NULL) {

  if(verbose == TRUE) cat('\n This is the single parameter version of binSpect')


  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # data.table: set global variable
  genes = p.value = estimate = score = NULL

  # set binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # kmeans algorithm
  kmeans_algo = match.arg(kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))

  # implementation
  implementation = match.arg(implementation, choices = c('data.table', 'simple', 'matrix'))


  # spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial_network_name: ', spatial_network_name, ' does not exist, create a spatial network first')
  }

  # convert to full network
  if(reduce_network == FALSE) {
    spatial_network = convert_to_full_spatial_network(spatial_network)
    data.table::setnames(spatial_network, old = c('source', 'target'), new = c('from', 'to'))
  }




  ## start binarization ##
  ## ------------------ ##

  if(!is.null(bin_matrix)) {
    bin_matrix = bin_matrix
  } else {
    if(bin_method == 'kmeans') {

      bin_matrix = kmeans_binarize_wrapper(gobject = gobject,
                                           expression_values = expression_values,
                                           subset_genes = subset_genes,
                                           kmeans_algo = kmeans_algo,
                                           nstart = nstart,
                                           iter_max = iter_max,
                                           extreme_nr = extreme_nr,
                                           sample_nr = sample_nr,
                                           set.seed = set.seed)

    } else if(bin_method == 'rank') {

      # expression
      values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
      expr_values = select_expression_values(gobject = gobject, values = values)

      if(!is.null(subset_genes)) {
        expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
      }

      max_rank = (ncol(expr_values)/100)*percentage_rank
      bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = rank_binarize, max_rank = max_rank))
    }
  }

  if(verbose == TRUE) cat('\n 1. matrix binarization complete \n')

  ## start with enrichment ##
  ## --------------------- ##

  if(implementation == 'simple') {
    if(do_parallel == TRUE) {
      warning('Parallel not yet implemented for simple. Enrichment will default to serial.')
    }

    if(calc_hub == TRUE) {
      warning('Hub calculation is not possible with the simple implementation, change to matrix if requird.')
    }


    result = calc_spatial_enrichment_minimum(spatial_network = spatial_network,
                                             bin_matrix = bin_matrix,
                                             adjust_method = adjust_method,
                                             do_fisher_test = do_fisher_test)


  } else if(implementation == 'matrix') {

    result = calc_spatial_enrichment_matrix(spatial_network = spatial_network,
                                            bin_matrix = bin_matrix,
                                            adjust_method = adjust_method,
                                            do_fisher_test = do_fisher_test,
                                            do_parallel = do_parallel,
                                            cores = cores,
                                            calc_hub = calc_hub,
                                            hub_min_int = hub_min_int)

  } else if(implementation == 'data.table') {

    result = calc_spatial_enrichment_DT(bin_matrix = bin_matrix,
                                        spatial_network = spatial_network,
                                        calc_hub = calc_hub,
                                        hub_min_int = hub_min_int,
                                        group_size = group_size,
                                        do_fisher_test = do_fisher_test,
                                        adjust_method = adjust_method,
                                        cores = cores)
  }

  if(verbose == TRUE) cat('\n 2. spatial enrichment test completed \n')





  ## start with average high expression ##
  ## ---------------------------------- ##

  if(get_av_expr == TRUE) {

    # expression
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)

    if(!is.null(subset_genes)) {
      expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
    }

    sel_expr_values = expr_values * bin_matrix
    av_expr = apply(sel_expr_values, MARGIN = 1, FUN = function(x) {
      mean(x[x > 0])
    })
    av_expr_DT = data.table::data.table(genes = names(av_expr), av_expr = av_expr)
    result = merge(result, av_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 3. (optional) average expression of high expressing cells calculated \n')
  }



  ## start with number of high expressing cells ##
  ## ------------------------------------------ ##

  if(get_high_expr == TRUE) {
    high_expr = rowSums(bin_matrix)
    high_expr_DT = data.table::data.table(genes = names(high_expr), high_expr = high_expr)
    result = merge(result, high_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 4. (optional) number of high expressing cells calculated \n')
  }


  # sort
  if(do_fisher_test == TRUE) {
    data.table::setorder(result, -score)
  } else {
    data.table::setorder(result, -estimate)
  }


  return(result)

}
