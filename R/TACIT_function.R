#' Plot Marker Based on Quantile
#'
#' This function creates a spatial plot of `x` coordinate and `y` coordinate with points colored based on
#' the `val` parameter. Points are colored according to a gradient scale from grey to dark blue,
#' with colors representing the value of `val` up to the quantile specified by `q`.
#'
#' @param x numeric vector of x coordinate values.
#' @param y numeric vector of y coordinate values.
#' @param val numeric vector of values used to color points.
#' @param q single numeric value between 0 and 1 indicating the quantile of `val` to use for coloring.
#'
#' @return A ggplot object representing the scatter plot with colored points.
#' @export
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' val <- rnorm(100)
#' q <- 0.95
#' plot_marker(x, y, val, q)
#'
#' @import ggplot2
#' @importFrom stats quantile
plot_marker <- function(x, y, val, q) {

  upper_value <- as.numeric(quantile(val, q))
  data <- data.frame(x, y, val)

  ### plot scatter plot
  p1 <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = val), size = 0.5) +
    theme_classic(base_size = 15) +
    scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, upper_value), na.value = "darkblue")
  p1
}




#' Find blocks from PROTEINxCELLTYPE matrix
#'
#' This function analyzes a PROTEINxCELLTYPE matrix to extract main cell types and their associated markers.
#' It also identifies subset cell types that share all markers of a main cell type and have additional unique markers.
#'
#' @param signature a matrix where rows correspond to cell types and columns to markers.
#'        The first column should be the name of the cell types.
#'
#' @return A list where each main cell type is the name of the list element.
#'         Each element contains a character vector of additional subset markers not present in the main cell type.
#'
#' @importFrom stats setdiff
#'
#' @export

find_main_cell_types_with_subsets_and_markers <- function(signature) {
  # Extract cell types and markers
  cell_types <- signature[,1]
  markers <- colnames(signature)[-1] # Exclude the first column which is Cell_type

  # Initialize a list to store the main cell types and their subset markers
  main_cell_types <- list()

  # Loop over each cell type
  for (i in 1:length(cell_types)) {
    main_type <- cell_types[i]
    main_markers <- colnames(signature)[signature[i, ] == 1]

    # Store additional markers for subsets
    additional_subset_markers <- character()

    # Check against all other cell types for a subset
    for (j in 1:length(cell_types)) {
      if (i != j) {
        subset_type <- cell_types[j]
        subset_markers <- colnames(signature)[signature[j, ] == 1]

        # Check if subset has all markers of main and additional ones
        if (all(main_markers %in% subset_markers) && length(subset_markers) > length(main_markers)) {
          additional_subset_markers <- unique(c(additional_subset_markers, setdiff(subset_markers, main_markers)))
        }
      }
    }

    # Add to main cell types with their additional subset markers
    if (length(additional_subset_markers) > 0) {
      main_cell_types[[main_type]] <- additional_subset_markers
    }
  }

  return(main_cell_types)
}



#' Determine Thresholds for Cell Types Score and Markers
#'
#' This function identifies optimal threshold values for cell type score and markers in a CELLxPROTEIN and CELLxCELLTYPE matrix based on cluster analysis.
#' It uses segmented regression to find breakpoints in the median intensity values of markers across clusters.
#' The function then assigns intensity groups based on these breakpoints and calculates a threshold to
#' binary marker expressions into high (1) and low (0) intensity groups.
#'
#' @param k the index of the marker column in the `data_anb` dataframe.
#' @param data_anb a dataframe containing marker intensity data.
#' @param clusters a vector or factor of cluster assignments for each cell in `data_anb`.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{Group_Threshold_new}{A binary matrix representing high (1) and low (0) expressions of markers.}
#'   \item{threshold_point_final}{The calculated threshold value for marker binarization.}
#' }
#'
#' @importFrom dplyr data_frame
#' @importFrom stats setdiff reorder lm
#' @importFrom segmented segmented
#' @examples
#' # Example usage:
#' # Assuming 'data_anb' is your dataframe and 'clusters' is a vector of cluster assignments.
#' k <- 2 # Index of the marker column in 'data_anb'
#' result <- threshold_function(k, data_anb, clusters)
#' print(result)
#'
#' @export

threshold_function=function(k,data_anb,clusters){

  ## Save the name of antibodies
  name = colnames(data_anb[k])

  ## add louvain cluster to data intensity
  data_intensity_cluster=data_frame(data_anb,lv_clusters=clusters)

  # Order louvain cluster by median
  group_ordered <- with(data_intensity_cluster,
                        reorder(lv_clusters,
                                get(name),median))

  # Create data with reordered group levels
  data_ordered <- data_intensity_cluster  # keep order of orginial data
  data_ordered$group <- factor(data_ordered$lv_clusters,
                               levels = levels(group_ordered))

  # Create data with median value in each cluster
  data_median=data.frame(median_value=attr(group_ordered, "scores"),lv_clusters=names(attr(group_ordered, "scores")))
  data_median=data_median[order(data_median$median_value,decreasing = F),]


  if(sum(data_median$median_value)==0){
    threshold_point_final=max(data_ordered[,name]) # if all median louvain = 0, we do not have positive cells (threshold at max)


  }else{
    x = 1:nrow(data_median)
    y = data_median$median_value
    fit <- lm(y ~ x)

    #fit piecewise regression model to original model, estimating a breakpoint
    AIC_value=NULL
    for (i in 1:3) {
      out <-tryCatch( { segmented.fit <- segmented(fit, seg.Z = ~x,npsi=i) }
                      , error = function(e) {e})
      if(any(class(out) == "error")==T){
        AIC_value[i]=100000
      }else{
        segmented.fit <- segmented(fit, seg.Z = ~x,npsi=i)
        AIC_value[i]=AIC(segmented.fit)
      }
    }

    ##Choose the smallest AIC for better fit
    segmented.fit <- segmented(fit, seg.Z = ~x,npsi=which.min(AIC_value))

    ##get the breakpoint value
    breakpoint_point=data.frame(x=segmented.fit$psi[,2],y=segmented.fit$fitted.values[floor(segmented.fit$psi[,2])])


    ##assign group base on multiple breakpoints
    Group_expression=threshold_groups(thresholds=breakpoint_point$x, input_vector=1:nrow(data_median))
    data_median_update=data.frame(data_median,Group_expression)
    new_data_anb=merge(data_intensity_cluster,data_median_update,"lv_clusters")






    ## assign lowest intensity group (0) and highest intensity group (4)
    if(which.min(AIC_value)==3){
      new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==3,4,new_data_anb$Group_expression)
      slope=as.numeric(segmented.fit$coefficients[c(2,3,4,5)])
      if(slope[2]<0.05){
        new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==1,0,new_data_anb$Group_expression)
      }else{
        new_data_anb$Group_expression=new_data_anb$Group_expression
      }
    }

    if(which.min(AIC_value)==2){
      new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==2,4,new_data_anb$Group_expression)
      slope=as.numeric(segmented.fit$coefficients[c(2,3,4)])
      if(slope[2]<0.05){
        new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==1,0,new_data_anb$Group_expression)
      }else{
        new_data_anb$Group_expression=new_data_anb$Group_expression
      }
    }

    if(which.min(AIC_value)==1){
      new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==1,4,new_data_anb$Group_expression)
      slope=as.numeric(segmented.fit$coefficients[c(2,3)])
      if(slope[2]<0.05){
        new_data_anb$Group_expression=ifelse(new_data_anb$Group_expression==1,0,new_data_anb$Group_expression)
      }else{
        new_data_anb$Group_expression=new_data_anb$Group_expression
      }
    }





    ## create a vector of threshold
    cut_seq=seq(as.numeric(quantile(new_data_anb[,k+1],0.05)),as.numeric(quantile(new_data_anb[,k+1],0.99)),length.out=1000)
    cut_seq=cut_seq[-1] #sometime it 0, I would remove
    threshold_point_final=NULL
    low_high = high_low = high_high = low_low = NULL
    for (g in 1:length(cut_seq)) {

      high_data = new_data_anb[which(new_data_anb$Group_expression==4),k+1]
      low_data = new_data_anb[which(new_data_anb$Group_expression==0),k+1]

      # low_high: low expr in high group
      # high_low: high expr in low group

      low_high[g] = length(high_data[which(high_data < cut_seq[g] | high_data == cut_seq[g])])/length(high_data)
      high_low[g] = length(low_data[which(low_data>cut_seq[g]|low_data==cut_seq[g])])/length(low_data)

      # high_high: high expr in high group
      # low_low: low expr in low group

      high_high[g] = length(high_data[which(high_data>cut_seq[g]|high_data==cut_seq[g])])/length(high_data)
      low_low[g] = length(low_data[which(low_data<cut_seq[g]|low_data==cut_seq[g])])/length(low_data)

    }
    threshold_point_final=cut_seq[which.min(low_high+high_low)]





  }

  ## create matrix binarization
  Group_Threshold_new=ifelse(as.numeric(t(data_intensity_cluster[,name]))>threshold_point_final, 1,0)


  return(list(Group_Threshold_new,threshold_point_final))
}




#' Assign Groups Based on Thresholds
#'
#' This function assigns group numbers to elements of an input vector based on a set of threshold values.
#' Each element of the input vector is categorized into a group based on which thresholds it falls between.
#' The function supports multiple thresholds, allowing for a range of group assignments.
#'
#' @param input_vector a numeric vector to which group numbers will be assigned.
#' @param thresholds a numeric vector of threshold values used for grouping.
#'        Each element in `input_vector` is assigned a group number based on these thresholds.
#'
#' @return A numeric vector of the same length as `input_vector`,
#'         with each element assigned a group number based on the specified thresholds.
#'
#' @examples
#' # Example usage:
#' input_vector <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' thresholds <- c(3, 6, 9)
#' groups <- threshold_groups(input_vector, thresholds)
#' print(groups)
#'
#' @export

threshold_groups <- function(input_vector, thresholds) {
  # Initialize an empty vector to hold the group assignments
  groups <- numeric(length(input_vector))

  # Loop through each threshold in the list
  for (i in seq_along(thresholds)) {
    # Create a logical vector indicating which elements in the input_vector are above the threshold
    group <- input_vector >= thresholds[i]

    # Assign group i to all elements in input_vector that are above threshold i and below threshold (i+1)
    if (i < length(thresholds)) {
      next_group <- input_vector >= thresholds[i+1]
      group <- group & !next_group
    }

    # Assign group i to all elements in input_vector that are above threshold i
    groups[group] <- i
  }

  # Return the vector of group assignments
  return(groups)
}




#' Replace Zero with Column Names in Binary Data
#'
#' This function takes a binary matrix (or dataframe) and replaces the zeros
#' with the corresponding column names. Elements with a value of 1 are replaced
#' with the name of the column they belong to. This is useful for converting
#' binary representations into a more descriptive format.
#'
#' @param binary_data a binary matrix or dataframe with numeric 0 and 1 values.
#'
#' @return A matrix of the same dimensions as `binary_data` where each element
#'         that is 1 in the input is replaced with the corresponding column name
#'         and each 0 remains as 0.
#'
#' @examples
#' # Example usage:
#' binary_data <- matrix(c(1, 0, 0, 1, 0, 1), nrow = 2, ncol = 3,
#'                      dimnames = list(NULL, c("A", "B", "C")))
#' result <- replace_zero_with_column_names(binary_data)
#' print(result)
#'
#' @export

replace_zero_with_column_names <- function(binary_data) {
  result_matrix <- matrix("", nrow = nrow(binary_data), ncol = ncol(binary_data))

  for (i in 1:nrow(binary_data)) {
    for (j in 1:ncol(binary_data)) {
      if (binary_data[i, j] == 1) {
        result_matrix[i, j] <- colnames(binary_data)[j]
      }else{
        result_matrix[i, j] <- 0
      }
    }
  }

  return(result_matrix)
}

#' Identify Unique Markers for Subsets of Cell Types
#'
#' This function analyzes a signature matrix to identify unique markers for subset cell types.
#' It compares each cell type to potential parent cell types, identifying markers that are unique
#' to the subset cell type compared to its immediate parent.
#'
#' @param signature a matrix with cell types as rows and markers as columns.
#'        The first column is expected to contain the cell type names.
#'
#' @return A list where the names are subset cell types and each element contains
#'         the unique marker for that cell type, if one exists. If a cell type does not
#'         have a unique marker or is not a subset, it will not appear in the list.
#'
#' @importFrom stats setdiff
#'
#' @export
find_unique_marker_for_subsets <- function(signature) {
  cell_types <- signature[,1]
  markers <- colnames(signature)[-1]  # Exclude the first column which is Cell_type

  # Initialize a list to store the unique marker for each subset cell type
  unique_markers <- list()

  # Loop over each cell type
  for (i in 1:length(cell_types)) {
    subset_type <- cell_types[i]
    subset_markers <- colnames(signature)[signature[i, ] == 1]

    # Find the immediate parent cell type
    parent_markers <- NULL
    for (j in 1:length(cell_types)) {
      if (i != j) {
        potential_parent_markers <- colnames(signature)[signature[j, ] == 1]
        if (all(potential_parent_markers %in% subset_markers) &&
            length(potential_parent_markers) < length(subset_markers)) {
          parent_markers <- potential_parent_markers
        }
      }
    }

    # Determine the unique marker
    if (!is.null(parent_markers)) {
      additional_markers <- setdiff(subset_markers, parent_markers)
      if (length(additional_markers) == 1) {
        unique_markers[[subset_type]] <- additional_markers
      }
    }
  }

  return(unique_markers)
}



#' Threshold-based Assignment of Cell Types from Multiplexed Imaging DaTa
#'
#' This function identifies cell type annotation using CELLxFEATURE matrix and Signature matrix in spatial omics data.
#' It uses segmented regression to find breakpoints in the median intensity values of cell type relevant score
#' across clusters.The function then assigns intensity groups based on these breakpoints and
#' calculates a threshold for each cell type relevant score. The KNN would be used in the feature subspace to gave final annotation.
#'
#' @param data_expression a matrix or dataframe representing cell in the column and marker in the row.
#' @param r resolution value for microclusters.
#' @param p number of dimension used for microclusters.
#' @param Signature a matrix or dataframe representing signature of the cell type.
#'
#' @return A list with four elements:
#' \itemize{
#'   \item{final_threshold}{The dataframe of threshold for each cell type relevant score.}
#'   \item{Group_threshold_data}{The binary matrix using the threshold.}
#'   \item{data_clean_final}{The dataframe of final annotation.}
#'   \item{clusters}{The microclusters index for each cells.}
#' }
#'
#'
#'
#' @importFrom stats setdiff colSums Seurat
#' @importFrom caret knn3
#' @examples
#' # Example usage:
#' # Assuming 'data_expression', and 'Signature' are already defined
#' annotation <- get_annotation(data_expression, r, Signature, p)
#' print(annotations[3])
#'
#' @export


TACIT <- function(data_expression, r, p, Signature) {
  
  data_anb=data_expression[,colnames(Signature)[-1]]
  
  ### Processing data
  #----
  ## Signature data (St)
  Signature[is.na(Signature) == T] <- 0 ## assign all NA as 0
  Signature <- as.data.frame(Signature)
  
  # Intensity data (A_ij)
  
  main_cell_types_with_markers <- find_main_cell_types_with_subsets_and_markers(Signature)
  ct_main <- names(main_cell_types_with_markers)
  marker_subset <- unique(unlist(main_cell_types_with_markers))
  
  
  
  orig_values <- as.matrix(data_anb)
  rownames(orig_values) <- 1:nrow(data_anb)
  orig_values <- t(orig_values)
  orig_values_metadata <- data.frame("CellID" = 1:nrow(data_anb))
  row.names(orig_values_metadata) <- orig_values_metadata$CellID
  scfp <- CreateSeuratObject(counts = orig_values, meta.data = orig_values_metadata)
  
  
  ### Create Seurat object
  #----
  print("-------------------REACHED HERE: NormalizeData---------------")
  # Set up parallel processing
  plan(sequential)
  scfp <- NormalizeData(scfp, normalization.method = "CLR", margin = 2) # do not normalize data
  print("-------------------REACHED HERE: FindVariableFeatures---------------")
  scfp <- FindVariableFeatures(scfp, selection.method = "disp", nfeatures = 500)
  print("-------------------REACHED HERE: ScaleData---------------")
  plan(multisession, workers = parallelly::availableCores() - 1)
  scfp <- ScaleData(scfp, features = rownames(scfp))
  
  
  
  
  
  #### Clustering
  #----
  print("-------------------REACHED HERE: RunPCA---------------")
  scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp))
  print("-------------------REACHED HERE: RunUMAP---------------")
  
  
  ####Clustering
  #----
  scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp))
  print("-------------------REACHED HERE: RunUMAP---------------")
  total_cells <- ncol(scfp)
  #### Clustering
  #----
  if (total_cells < 100000) {
    print("Less than 100000 cells, running UMAP on the entire dataset")
    scfp <- RunUMAP(scfp,
                    features = VariableFeatures(object = scfp), n.components = p, metric = "correlation", umap.method = "uwot",
                    n.neighbors = 15, learning.rate = 1, spread = 1, min.dist = 0.01, set.op.mix.ratio = 1, local.connectivity = 1
    )
  } else {
    print("More than 100000 cells, training UMAP using 30% of cells for training or max 100000 cells, whichever is smaller")
    # Use 30% of cells for training or max 100000 cells, whichever is smaller
    subset_size <- min(100000, ceiling(total_cells * 0.3))
    # For very small datasets, use all cells
    subset_size <- min(subset_size, total_cells)
    
    scfp <- FastRunUMAP(
      seurat_obj = scfp,
      subset_size = subset_size,
      n.components = p,
      n.neighbors = 15,
      metric = "correlation",
      min.dist = 0.01, # Same as original
      random_seed = 1, # For reproducibility
      verbose = TRUE
    )
  }
  
  print("-------------------REACHED HERE: FindNeighbors---------------")
  scfp <- FindNeighbors(scfp, dims = 1:p, reduction = "umap")
  
  
  
  
  # Run clusters
  print("-------------------REACHED HERE: FindClusters---------------")
  scfp <- FindClusters(scfp, resolution = r, random.seed = 1, algorithm = 4, n.iter = 2)
  print("---------------- now here ---------------")
  # Collected clusters ID for each cells
  clusters <- as.numeric(scfp@meta.data[["seurat_clusters"]])
  
  print(paste0("Number of clusters: ", length(unique(clusters)), " with ", round(mean(as.numeric(table(clusters))), 0), " cells per clusters"))
  
  ### Cell type score
  #----
  data_anb_matrix <- as.matrix((data_anb)) # scale data intensity and convert to matrix from data frame
  Signature_matrix <- as.matrix(Signature[, -1]) # convert signature to matrix from data frame
  ct <- data_anb_matrix %*% t(Signature_matrix)
  colnames(ct) <- Signature[, 1]
  
  
  #  data_anb=scale(data_anb)
  
  ### Threshold for cell type
  #----
  library(parallel)
  library(foreach)
  library(doParallel)
  library(caret)
  
  ct <- as.data.frame(ct)
  
  num_cores <- parallelly::availableCores() - 1 # Leave one core free for system processes
  registerDoParallel(cores = num_cores)
  
  # Parallelize first loop
  Group_threshold_data <- NULL
  final_threshold <- NULL
  
  cat("==== Phase 1: Processing Main Data Columns ====\n")
  cat("Processing", ncol(ct), "columns in parallel...\n")
  final_threshold <- foreach(
    k = 1:ncol(ct), .combine = c,
    .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  ) %dopar% {
    # Return the second element of the threshold_function result
    vector <- threshold_function(k, ct, clusters)
    return(vector[[2]])
  }
  
  # Process Group_threshold_data in parallel
  Group_threshold_data <- foreach(
    k = 1:ncol(ct), .combine = cbind,
    .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  ) %dopar% {
    # Return the first element of the threshold_function result
    vector <- threshold_function(k, ct, clusters)
    return(vector[[1]])
  }
  
  
  colnames(Group_threshold_data) <- colnames(ct)
  final_threshold <- data.frame(value = final_threshold, Name = colnames(ct)) ########### . Output 1
  
  # Conditionally handle marker subset processing
  if (length(marker_subset) > 1) {
    cat("==== Phase 2: Processing Marker Subset Data ====\n")
    cat("Detected", length(marker_subset), "marker subsets to process\n")
    
    ct2 <- data_anb[, marker_subset]
    ct2 <- as.data.frame(ct2)
    colnames(ct2) <- marker_subset
    
    cat("Processing marker subset columns in parallel...\n")
    
    # Parallelize second loop
    final_threshold2 <- foreach(
      k = 1:ncol(ct2), .combine = c,
      .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
    ) %dopar% {
      vector <- threshold_function(k, data_anb = ct2, clusters = clusters)
      return(vector[[2]])
    }
    
    Group_threshold_data2 <- foreach(
      k = 1:ncol(ct2), .combine = cbind,
      .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
    ) %dopar% {
      vector <- threshold_function(k, data_anb = ct2, clusters = clusters)
      return(vector[[1]])
    }
    
    final_threshold2 <- data.frame(value = final_threshold2, Name = marker_subset) ########### . Output 1
    colnames(Group_threshold_data2) <- marker_subset
    
    cat("==== Phase 3: Processing Cell Types and Markers ====\n")
    cat("Processing", length(ct_main), "main cell types...\n")
    
    # Vectorize operations where possible
    for (p in 1:length(ct_main)) {
      cat("  - Processing cell type", p, "of", length(ct_main), ":", ct_main[p], "\n")
      ct_score_current <- Group_threshold_data[, ct_main[p]]
      
      if (length(unlist(main_cell_types_with_markers[ct_main[p]])) == 1) {
        rSum <- (Group_threshold_data2[, unlist(main_cell_types_with_markers[ct_main[p]])])
      } else {
        rSum <- rowSums(Group_threshold_data2[, unlist(main_cell_types_with_markers[ct_main[p]])])
      }
      Group_threshold_data[, ct_main[p]] <- ifelse(rSum > 0, 0, ct_score_current)
    }
    
    unique_markers_for_subsets <- find_unique_marker_for_subsets(Signature)
    if (length(unique_markers_for_subsets) > 0) {
      cat("\nProcessing", length(unique_markers_for_subsets), "unique markers for subsets...\n")
      # Process unique markers more efficiently
      for (p in 1:length(unique_markers_for_subsets)) {
        ctype <- names(unique_markers_for_subsets[p])
        marker_id <- unique_markers_for_subsets[ctype][[1]]
        Group_threshold_data[, ctype] <- ifelse(Group_threshold_data2[, marker_id] == 0, 0, Group_threshold_data[, ctype])
      }
    }
  }
  
  cat("==== Phase 4: Post-Processing and Analysis ====\n")
  cat("Generating potential cell types...\n")
  potential_ct <- replace_zero_with_column_names_opt(Group_threshold_data)
  
  cat("Creating test data...\n")
  potential_ct_test <- replace_zero_with_column_names_opt(Group_threshold_data)
  combined_vector <- apply(potential_ct_test, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  
  cat("Calculating group mixtures...\n")
  group_mix <- rowSums(Group_threshold_data)
  Group_threshold_data <- data.frame(ID = 1:nrow(Group_threshold_data), clusters, Group_threshold_data)
  data_cluster <- data.frame(clusters, ct)
  colnames(data_cluster)[-1] <- Signature$cell_type
  
  cat("Calculating median values per cluster...\n")
  medians_per_cluster <- aggregate(. ~ clusters, data = data_cluster, FUN = median)
  medians_per_cluster_filter <- medians_per_cluster[, -1]
  
  cat("Applying thresholds and converting to binary values...\n")
  # Convert final_threshold to a named vector for easier lookup
  thresholds <- setNames(final_threshold$value, final_threshold$Name)
  # Apply thresholding to convert values to binary
  binary_df <- medians_per_cluster_filter
  binary_df[] <- sapply(names(medians_per_cluster_filter), function(column_name) {
    if (column_name %in% names(thresholds)) {
      as.integer(medians_per_cluster_filter[[column_name]] > thresholds[[column_name]])
    } else {
      medians_per_cluster_filter[[column_name]]
    }
  })
  
  cat("Generating potential cell type clusters...\n")
  potential_ct_cluster <- replace_zero_with_column_names_opt(binary_df)
  combined_vector <- apply(potential_ct_cluster, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  
  cat("Identifying clean and mixed clusters...\n")
  cluster_clean <- medians_per_cluster$clusters[rowSums(binary_df) == 1]
  cluster_mix <- medians_per_cluster$clusters[rowSums(binary_df) != 1]
  
  data_cluster_group <- data.frame(clusters = medians_per_cluster$clusters, combined_vector)
  data_merge_cluster <- merge(data_cluster_group, Group_threshold_data, "clusters")
  data_merge_cluster <- data_merge_cluster[order(data_merge_cluster$ID), ]
  
  
  cat("Separating clean and mixed groups...\n")
  data_clean_group <- data_merge_cluster[which(data_merge_cluster$clusters %in% cluster_clean), ]
  data_mix_group <- data_merge_cluster[which(data_merge_cluster$clusters %in% cluster_mix), ]
  
  cat("\nCluster Summary:\n")
  print(paste0("Number of clean MicroClusters: ", length(cluster_clean)))
  print(paste0("Number of mixed MicroClusters: ", length(cluster_mix)))
  cat("- Total MicroClusters:", length(cluster_clean) + length(cluster_mix), "\n\n")
  
  colnames(data_mix_group)[-c(1, 2, 3)] <- Signature$cell_type
  
  potential_ct <- replace_zero_with_column_names_opt(data_mix_group[, -c(1, 2, 3)])
  
  ### Identify clean and mix data
  
  cat("Processing cell type combinations...\n")
  ## vector of potential cell type
  combined_vector <- apply(potential_ct, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  ## adding the vector combination into data
  
  cat("Preparing final data matrices...\n")
  data_anb_matrix <- data_anb_matrix[, colnames(Signature)[-1]]
  data_check <- data.frame(combined_vector, data_anb_matrix[data_mix_group$ID, ], data_mix_group[, -c(1, 2, 3)], ID = data_mix_group$ID)
  
  
  colnames(data_check)[2:sum(ncol(data_anb_matrix), 1)] <- colnames(data_anb_matrix)
  colnames(data_check)[sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))] <- Signature[, 1]
  
  ## Count the number of instances where a cell's value exceeds the corresponding threshold.
  row_sums <- rowSums(data_check[, sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))])
  
  ## Selected clean cell types
  data_out <- data_check[row_sums %in% c(1, 0), ]
  data_clean <- data.frame(mem = c(data_out$combined_vector), ID = c(data_out$ID))
  
  ## Selected clean cell types
  data_mixed <- data_check[!(row_sums %in% c(1, 0)), ]
  data_mixed <- data_mixed[which(data_mixed$combined_vector %in% names(table(data_mixed$combined_vector))[as.numeric(table(data_mixed$combined_vector)) > 0]), ]
  data_mixed <- data_mixed[order(rowSums(data_mixed[sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))])), ]
  
  
  cat("==== Phase 5: Final Processing and KNN Analysis ====\n")
  ######## Solving for mixed group
  ## KNN
  #----
  library(parallel)
  library(foreach)
  library(doParallel)
  library(caret)
  
  data_out_test <- data_out
  data_clean_final <- data_clean
  data_clean_cluster <- data.frame(mem = data_clean_group$combined_vector, ID = data_clean_group$ID)
  data_clean_final <- rbind(data_clean_final, data_clean_cluster)
  
  
  print(paste0("Number of clean cells: ", nrow(data_clean_final), " with ", length(unique(data_clean_final$mem)), " cell types"))
  print(paste0("Number of mixed cells: ", nrow(data_mixed), " with ", length(unique(data_mixed$combined_vector)), " mixed cell types"))
  
  # Get unique mixed cell types
  mix_group <- unique(data_mixed$combined_vector)
  
  clean_cell_types <- unique(data_clean_final$mem)
  mixed_by_type <- split(data_mixed, data_mixed$combined_vector)
  
  # Get signature markers for each clean cell type
  signature_markers <- list()
  for (ct in clean_cell_types) {
    signature_markers[[ct]] <- names(which(colSums(Signature[Signature[, 1] == ct, ][-1]) > 0))
  }
  
  num_cores <- min(parallelly::availableCores() - 1, length(mix_group))
  registerDoParallel(cores = num_cores)
  
  results <- foreach(i = 1:length(mix_group), .combine = rbind, .packages = c("caret")) %dopar% {
    current_mix_type <- mix_group[i]
    
    # Get all cells with this specific mixed type
    data_choose_second <- mixed_by_type[[current_mix_type]]
    
    # Get potential component cell types
    mix_select <- strsplit(current_mix_type, "::")[[1]]
    mix_select <- intersect(mix_select, clean_cell_types)
    
    if (length(mix_select) == 0) {
      return(data.frame(ID = data_choose_second$ID, mem = "Others"))
    }
    
    # Selected markers that signature for component cell type
    anb_select <- names(which(colSums(Signature[Signature[, 1] %in% mix_select, ][-1]) > 0))
    
    if (length(anb_select) == 0) {
      return(data.frame(ID = data_choose_second$ID, mem = "Others"))
    }
    
    # Get the data intensity from the clean cells (all at once for these cell types)
    data_clean_work <- data_out_test[
      which(data_out_test$combined_vector %in% mix_select),
      c(anb_select, "combined_vector", "ID")
    ]
    
    if (nrow(data_clean_work) == 0) {
      return(data.frame(ID = data_choose_second$ID, mem = "Others"))
    } else {
      # Extract relevant data columns for mixed cells
      data_test <- data_choose_second[, c(anb_select, "ID")]
      
      # Add small random noise (vectorized operation)
      train_matrix <- as.matrix(data_clean_work[, anb_select]) +
        matrix(runif(nrow(data_clean_work) * length(anb_select), 0, 0.001),
               nrow = nrow(data_clean_work)
        )
      
      test_matrix <- as.matrix(data_test[, anb_select]) +
        matrix(runif(nrow(data_test) * length(anb_select), 0, 0.001),
               nrow = nrow(data_test)
        )
      
      knn_model <- caret::knn3(train_matrix,
                               as.factor(data_clean_work$combined_vector),
                               k = min(10, nrow(data_clean_work))
      )
      
      # Get predictions for all cells at once
      predictions <- predict(knn_model, test_matrix, type = "prob")
      
      max_indices <- max.col(predictions, ties.method = "first")
      predicted_types <- colnames(predictions)[max_indices]
      
      # Return batch results
      return(data.frame(ID = data_test$ID, mem = predicted_types))
    }
  }
  
  # Combine results with clean cells
  data_clean_final <- rbind(data_clean_final, results)
  
  # Final processing
  data_clean_final$mem <- ifelse(data_clean_final$mem == "", "Others", data_clean_final$mem)
  data_clean_final <- data_clean_final[order(as.numeric(data_clean_final$ID), decreasing = F), ]
  
  data_final=data.frame(TACIT=data_clean_final$mem,UMAP1=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,1]),
                        UMAP2=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,2]),data_expression)
  
  return(data_final)
}










#' Visualize TACIT Data with UMAP Coordinates
#'
#' This function applies UMAP (Uniform Manifold Approximation and Projection) to the provided data
#' and combines the resulting coordinates with existing spatial coordinates and cell annotations.
#' It's designed to prepare data for visualization in the context of TACIT analysis.
#'
#' @param data_anb a dataframe or matrix containing marker intensity data or other relevant data for UMAP.
#' @param coor_x a numeric vector representing the x-coordinates in a spatial context.
#' @param coor_y a numeric vector representing the y-coordinates in a spatial context.
#' @param cell_annotation a vector or factor providing cell annotations.
#'
#' @return A dataframe that includes original x and y coordinates, UMAP coordinates, cell annotations,
#'         and the original data used for UMAP. Column names are 'X', 'Y', 'UMAP1', 'UMAP2', 'TACIT',
#'         followed by the original columns from `data_anb`.
#'
#' @importFrom uwot umap
#' @examples
#' # Example usage:
#' # Assuming 'data_anb', 'coor_x', 'coor_y', and 'cell_annotation' are already defined
#' visualization_data <- TACIT_visualization(data_anb, coor_x, coor_y, cell_annotation)
#' print(visualization_data)
#'
#' @export
TACIT_visualization=function(data_anb,coor_x,coor_y,cell_annotation){
  data_umap=uwot::umap(data_anb)
  colnames(data_umap)=c("UMAP1","UMAP2")
  data_final=data.frame(X=coor_x,Y=coor_y,data_umap,TACIT=cell_annotation,data_anb)
}



#' Remove Far Connections in Triangulation
#'
#' This function removes edges from a triangulation object that are longer than a specified cutoff distance.
#' The cutoff is determined as a quantile of the distances between points in the triangulation.
#'
#' @param tri a triangulation object, typically from a Delaunay triangulation.
#' @param cutoff a numeric value between 0 and 1 indicating the quantile of edge distances to use as a cutoff.
#'
#' @return A modified triangulation object with edges longer than the cutoff distance removed.
#'
#' @examples
#' # Example usage:
#' # Assuming 'tri' is a triangulation object (e.g., from Delaunay triangulation)
#' cutoff <- 0.95
#' new_tri <- removeFarConnections(tri, cutoff)
#' print(new_tri)
#'
#' @export
removeFarConnections <- function(tri,cutoff) {


  # Extract the coordinates of the points from the triangulation
  point_coords <- tri$delsgs[, 1:4]

  # Compute the Euclidean distance between each pair of points
  distances <- sqrt((point_coords[, 1] - point_coords[, 3])^2+(point_coords[, 2] - point_coords[, 4])^2)
  threshold_distance <- as.numeric(quantile(as.numeric(distances),cutoff))

  # Get the indices of edges that are within the threshold distance
  valid_edges <- distances <= threshold_distance

  # Get the indices of points involved in the valid edges
  valid_points <- unique(c(tri$delsgs[valid_edges, 1], tri$delsgs[valid_edges, 2]))

  # Create a new triangulation object with the valid connections
  new_tri <- tri
  new_tri$delsgs <- new_tri$delsgs[valid_edges, ]
  new_tri$dirsgs <- new_tri$dirsgs[valid_edges, ]
  new_tri$ptlist <- new_tri$ptlist[valid_points, ]

  return(new_tri)
}








#' Create a Dot-Heatmap for Antibody Expressions
#'
#' This function generates a heatmap to visualize the expression levels and frequency of antibodies.
#' It uses median values for expression levels and calculates the percentage of cells expressing each
#' antibody in different groups.
#'
#' @param data_pos_neg a dataframe with binary (0/1) data indicating positive or negative expression.
#' @param data_anb a dataframe with quantitative data for antibody expressions.
#' @param group a vector indicating the group classification for each cell.
#' @param name_anb a vector of antibody names to order the heatmap.
#'
#' @return A ggplot object representing the dot-heatmap with points colored by mean value and sized by
#'         percentage of cells expressing each antibody.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_gradientn theme
#' @importFrom dplyr filter mutate
#' @importFrom cowplot theme_cowplot
#' @examples
#' # Example usage:
#' # Assuming 'data_pos_neg', 'data_anb', 'group', and 'name_anb' are already defined
#' heatmap <- createHeatmapForAntibodies(data_pos_neg, data_anb, group, name_anb)
#' print(heatmap)
#'
#' @export

heatmap_anb=function(data_pos_neg,data_anb,group,name_anb){



  dat100=data.frame(data_pos_neg,group)
  colnames(dat100)[-ncol(dat100)]=colnames(data_pos_neg)

  SLIDE_1_value_anb=data.frame(data_anb,group)
  mean_value=count=NULL


  for(i in 1:length(names(table(group)))){

    dat9=dat100[which(dat100$group==names(table(group))[i]),]
    per_cent_count=as.numeric(sapply(dat9[,c(1:ncol(data_pos_neg))], sum))/nrow(dat9)
    count=append(count,per_cent_count)

    data_final9=SLIDE_1_value_anb[which(group==names(table(group))[i]),-ncol(SLIDE_1_value_anb)]
    per_cent_mean=as.numeric(sapply(as.data.frame(data_final9), median))
    mean_value=append(mean_value,per_cent_mean)
  }

  antibody_name=rep(colnames(dat9[,c(1:ncol(data_pos_neg))]),length(names(table(group))))
  Cluster=rep(names(table(group)),each=ncol(dat9[,c(1:ncol(data_pos_neg))]))
  plot_dat=data.frame(mean_value,count,antibody_name,Cluster)




  markers <- plot_dat$antibody_name %>% unique()


  plot_dat2=plot_dat[,c(1,3,4)]







  Cluster_order=names(table(group))




  ############ NOTICE THE t() above)



  plot_dat2=data.frame(plot_dat2,count)

  Anb_order=name_anb
  plot_dat2=plot_dat2[which(plot_dat2$antibody_name%in%Anb_order),]


  my_colors <- c("black", "grey70", "firebrick")

  Anb_order=name_anb
  plot_dat2[is.na(plot_dat2)==T]=0
  p1=plot_dat2 %>% filter(antibody_name %in% markers) %>%
    mutate(`% Expressing` = count * 100,
           antibody_name = factor(antibody_name, levels = Anb_order),
           Cluster = factor(Cluster, levels = Cluster_order)) %>%
    ggplot(aes(x=antibody_name, y = Cluster, color = mean_value, size = `% Expressing`)) +
    geom_point() +
    cowplot::theme_cowplot() +
    #  theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 40, hjust=1)) +
    ylab('') +
    #  theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = my_colors, oob = scales::squish, name = 'Mean value', limits = c(-1,1))  +
    scale_y_discrete(position = "right")+xlab("")+ggtitle("")

  return(p1)
}





#' Faster UMAP implementation using subset-and-project approach
#'
#' @param seurat_obj The Seurat object
#' @param subset_size Number or proportion of cells to use for UMAP model training
#' @param n.components Number of UMAP dimensions
#' @param n.neighbors Number of neighbors for UMAP
#' @param random_seed Random seed for reproducibility
#' @param features Features to use (defaults to variable features)
#' @param metric Distance metric to use
#' @param min.dist Minimum distance parameter for UMAP
#' @param verbose Print progress messages
#' @return Updated Seurat object with UMAP embedding
FastRunUMAP <- function(seurat_obj,
                        subset_size = 50000,
                        n.components = 2,
                        n.neighbors = 15,
                        random_seed = 42,
                        features = NULL,
                        metric = "correlation",
                        min.dist = 0.1,
                        verbose = TRUE) {
  require(Seurat)
  require(uwot)
  require(Matrix)
  require(future)
  require(future.apply)
  
  # Set up parallel processing
  if (verbose) message("Preparing data for UMAP...")
  
  # Get features to use
  if (is.null(features)) {
    features <- VariableFeatures(seurat_obj)
  }
  
  # Get dimensionality reduction to use as input (PCA by default)
  if ("pca" %in% names(seurat_obj@reductions)) {
    use_dr <- "pca"
    input_data <- Embeddings(seurat_obj, reduction = "pca")
  } else {
    if (verbose) message("PCA not found, using scaled data...")
    input_data <- t(GetAssayData(seurat_obj, slot = "scale.data")[features, ])
  }
  
  # Determine subset size (if fraction, convert to number)
  cell_count <- ncol(seurat_obj)
  if (subset_size < 1) {
    subset_size <- ceiling(cell_count * subset_size)
  }
  subset_size <- min(subset_size, cell_count)
  
  if (verbose) {
    message(sprintf(
      "Using %d cells (out of %d) to train UMAP model...",
      subset_size, cell_count
    ))
  }
  
  # Randomly select subset of cells for UMAP model training
  set.seed(random_seed)
  subset_idx <- sample(cell_count, size = subset_size)
  
  # Train UMAP on subset
  if (verbose) message("Training UMAP model on subset...")
  
  subset_data <- input_data[subset_idx, ]
  
  # Train UMAP with uwot directly for more speed control
  umap_model <- uwot::umap(
    X = subset_data,
    n_components = n.components,
    n_neighbors = n.neighbors,
    metric = metric,
    min_dist = min.dist,
    ret_model = TRUE,
    n_epochs = 200, # Reduced for speed
    n_threads = future::availableCores() - 1,
    n_sgd_threads = "auto",
    batch = TRUE, # Use batched processing
    verbose = verbose
  )
  
  if (verbose) message("Projecting all cells onto UMAP model...")
  
  # Project all cells onto trained UMAP model
  if (cell_count <= 100000) {
    # For datasets that fit in memory, do it all at once
    all_umap <- uwot::umap_transform(
      X = input_data,
      model = umap_model,
      n_threads = future::availableCores() - 1,
      n_sgd_threads = "auto",
      verbose = verbose
    )
  } else {
    # For very large datasets, process in chunks
    chunk_size <- 100000
    n_chunks <- ceiling(cell_count / chunk_size)
    all_umap <- matrix(0, nrow = cell_count, ncol = n.components)
    
    for (i in 1:n_chunks) {
      if (verbose) message(sprintf("Processing chunk %d of %d...", i, n_chunks))
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, cell_count)
      chunk_data <- input_data[start_idx:end_idx, ]
      
      all_umap[start_idx:end_idx, ] <- uwot::umap_transform(
        X = chunk_data,
        model = umap_model,
        n_threads = future::availableCores() - 1,
        n_sgd_threads = "auto",
        verbose = TRUE
      )
    }
  }
  
  # Create new dimensional reduction object
  colnames(all_umap) <- paste0("UMAP_", 1:n.components)
  
  rownames(all_umap) <- colnames(seurat_obj)
  
  umap_dimred <- CreateDimReducObject(
    embeddings = all_umap,
    key = "UMAP_",
    assay = DefaultAssay(seurat_obj)
  )
  
  # Add UMAP to Seurat object
  seurat_obj[["umap"]] <- umap_dimred
  
  if (verbose) message("UMAP computation complete.")
  return(seurat_obj)
}








