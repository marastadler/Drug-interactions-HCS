

hiernet.lasso.weak <- function(x, y, q, l){
  
  if (!requireNamespace("hierNet", quietly = TRUE))
    stop("Package ", sQuote("hierNet"), " needed but not available")
  
  
  x.lin <- x[, 1:l]
  p_init <- ncol(x.lin)
  ## there might be zero columns (modifications) after removing observations
  ## (subsampling) which cumpute.interactions can't handle:
  ## Remove zero columns
  col_sub = apply(x.lin, 2, function(col) all(col == 0 ))
  x.lin = x.lin[, !col_sub]
  
  
  nlam = 30
  
  #library(tictoc)
  #tic("path")
  
  
  ## with passing zz we avoid that zz gets computed based on the scaled x
  ## which is not wanted for bunary features
  zz_binary <- compute.interactions.c(x.lin, diagonal = FALSE)
  fit <- hierNet::hierNet.path(x.lin, y, minlam = 1e-20, maxlam = 200, 
                               nlam = nlam, diagonal = FALSE, strong = FALSE,
                               stand.int = FALSE, stand.main = TRUE,
                               zz = zz_binary
                               )
  
  
  
  #-fitcv <- hierNet::hierNet.cv(fit, x.lin, y, nfolds = 5)
  #-lamhat = fitcv$lamhat
  # fitfinal <- hierNet(x = x.lin, 
  #                y = y,
  #                lam = lamhat, diagonal = FALSE, strong = TRUE, 
  #                stand.int = TRUE, stand.main = TRUE)
  #toc()
  p <- ncol(x.lin)
  coefmatrix_lam <- matrix(nrow = p * (p + 1)/2, ncol = nlam)
  
  linear_coef <- fit$bp - fit$bn ## linear coefficients
  
  for(lam in 1:length(fit$lamlist)){
    
    theta <- (fit$th[, , lam] + t(fit$th[, , lam]))/2 ## symmetrize theta
    interact_coef <- theta[lower.tri(theta)]
    
    coefmatrix_lam[, lam] <- c(linear_coef[, lam], interact_coef) ## combine linear effects and interactions
    
  }
  
  ## interaction feature names after potentially removing linear features
  A <- colnames(x.lin)
  n <- 0
  names_int <- c()
  for(i in 1:(ncol(x.lin) - 1)){
    for(j in (i + 1): ncol(x.lin)){
      n <- n + 1
      names_int[n] <- paste0(A[i], ":", A[j])
    }
  }
  
  rownames(coefmatrix_lam) <- c(colnames(x.lin), names_int)
  
  ## Empty matrix with all features (also the ones that where removed before)
  sequence <- matrix(nrow = p_init * (p_init + 1)/2, ncol = nlam)
  rownames(sequence) <- colnames(x)
  sequence[rownames(coefmatrix_lam), ] <- (coefmatrix_lam != 0)
  sequence[is.na(sequence)] <- 0
  
  ## select the lambda where number of selected features is <= q
  seq_q <- which(colSums(sequence) <= q)
  ret <- sequence[, seq_q[length(seq_q)]]
  
  ## return both
  return(list(selected = ret, path = sequence))
}


hiernet.stabsel.all <- function(X = X_filt, Y = Y_filt_sc, q = 50,
                                selection_probability = .6){
  
  fit.list.weak <- list()
  #l <- ncol(X)
  n <- 0
  
  X_interactions <- cbind(X, 
                          hierNet::compute.interactions.c(X, diagonal = F)) 
  
  for(i in 1:ncol(Y)){
    print(i)
    y1 = Y[, i]
    names(y1) = rownames(Y)
    n <- n + 1
    
    hiernet.lasso.weak_l <- function(x, y, q){
      return(hiernet.lasso.weak(x, y, q, l = ncol(X)))
    }
    stopifnot(all(rownames(X_interactions) == names(y1)))
    tryCatch(
      {
        fit.list.weak[[n]] = stabsel(x = X_interactions, y = y1,
                                     fitfun = hiernet.lasso.weak_l,
                                     cutoff = selection_probability,
                                     #PFER = 1,
                                     q = q, sampling.type = "SS" #folds = my_subsample(rep(1, nrow(X)), B = B)\
                                     , B = 50
        )
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  return(result = fit.list.weak)
  
}


refit.all <- function(X, Y, fit_list_weak){
  
  refit_list_weak <- list()
  coef_list_weak <- list()
  
  X_interactions <- cbind(X, 
                          hierNet::compute.interactions.c(X, diagonal = F))
  n <- 0 
  for(i in 1:ncol(Y)){
    n <- n + 1
    
    sel <- fit_list_weak[[n]]$selected
    if(length(sel)>0){
      y1 = Y[, i]
      names(y1) = rownames(Y)
      refit_list_weak[[n]] <- lm(y1 ~ X_interactions[, sel])
      coef <- refit_list_weak[[n]]$coefficients
      names(coef) <- substr(names(coef), 22, 
                            nchar(names(coef)))
      
      coef_list_weak[[n]] <- coef
    }
  }
  
  return(list(coef_list = coef_list_weak,
              refit_list = refit_list_weak))
}


## This function plots a matrix that shows which coefficients have been selected 
## by the model how often among all models

plot.count.coef <- function(coef_matrix_weak, X, cell_height = 14, cell_width = 14,
                            fontsize_row = 9, fontsize_col = 9){
  
  coef_matrix_weak <- t(coef_matrix_weak)
  ## check if matrix has the right dimension
  nrow(coef_matrix_weak) == ncol(X) * (ncol(X) + 1) /2
  
  df_overview <- matrix(ncol = 4, nrow = (ncol(X) * (ncol(X) + 1) /2 - ncol(X)))
  
  colnames(df_overview) <- c("Feat1", "Feat2", "BothLin", "Interaction")
  n = 0
  for(i in 1:(ncol(X)-1)
  ){
    for(j in (i + 1):ncol(X)
    ){
      #print(paste(i,":", j))
      n = n + 1
      sel_i <- (coef_matrix_weak[i, ] != 0)
      sel_j <- (coef_matrix_weak[j, ] != 0)
      bothlin_ij = sum(((sel_i + sel_j) == 2))
      int_ij <- paste0(rownames(coef_matrix_weak)[i], ":", rownames(coef_matrix_weak)[j])
      if(int_ij %in% rownames(coef_matrix_weak)){
        interaction_ij = sum(coef_matrix_weak[int_ij,] != 0)
      }
      else{interaction_ij = 0}
      
      df_overview[n, ] <- c(rownames(coef_matrix_weak)[i], rownames(coef_matrix_weak)[j],
                            bothlin_ij, interaction_ij)
      
    }
  }
  df_overview <- as.data.frame(df_overview)
  df_overview$BothLin <- as.numeric(df_overview$BothLin)
  df_overview$Interaction <- as.numeric(df_overview$Interaction)
  # plot(df_overview$BothLin, df_overview$Interaction)
  
  
  
  lin_mat = reshape2::acast(df_overview, Feat1 ~ Feat2, value.var = "BothLin")
  int_mat = reshape2::acast(df_overview, Feat1 ~ Feat2, value.var = "Interaction")
  lin_mat <- round(lin_mat)
  lin_mat[is.na(lin_mat)] <- 0
  int_mat[is.na(int_mat)] <- 0
  
  library(igraph)
  
  # Make undirected so that graph matrix will be symmetric
  g <- graph.data.frame(df_overview[,-4], directed=FALSE)
  # add value as a weight attribute
  lin_mat = get.adjacency(g, attr="BothLin", sparse=FALSE)
  
  g2 <- graph.data.frame(df_overview[,-3], directed=FALSE)
  # add value as a weight attribute
  int_mat = get.adjacency(g2, attr="Interaction", sparse=FALSE)
  diag(lin_mat) <- NA
  lin_int_mat <- lin_mat
  lin_int_mat[lower.tri(lin_int_mat)] <- int_mat[lower.tri(int_mat)]
  
  
  ret <- pheatmap(lin_int_mat, cluster_rows = F, cluster_cols = F,
           border_color = "white", cellheight = cell_height, cellwidth = cell_width,
           fontsize_row = fontsize_row, fontsize_col = fontsize_col,
           na_col = "white",
           color = colorRampPalette((brewer.pal(n = 7, name =
                                                  "Blues")))(100), 
           display_numbers = T, number_format =  "%.0f", number_color = "white",
           legend = F)
  return(plot = ret)
}


plot_interactions <- function(x = "DNA Meth. m5C", y = "H3K9 me3", 
                              # cutoff_labels = 0.15, 
                              coeff.matrix = coef_matrix, text.size = 2,
                              max.overlaps = 100
                              ){
  
  if(paste0(x, ":", y) %in% rownames(coeff.matrix)){
    xy <- paste0(x, ":", y)
  }
  else{xy <- paste0(y, ":", x)}
  df <- as.data.frame(t(coeff.matrix[c(x, y , xy), ]))
  df <- subset(df, abs(df[, xy]) > 0.05 )
  df$name <- rownames(df)
  
 
  dfs <- df
  
  hw_plot <- ggplot(df, aes(x = df[, x], y = df[, y],
                            colour = df[, xy]
  )) + 
    theme_minimal() +
    geom_vline(xintercept = 0, color = "grey") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_point(size = 3) +
    scale_colour_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      midpoint = 0,
      breaks = seq(100, from = -5, to = 4)
    ) +
    
    geom_text_repel(data = dfs,
                    aes(x = dfs[, x], y = dfs[, y],
                        label = name
                    ), size = text.size, color = "black", max.overlaps = max.overlaps)
  hw_plot <- hw_plot +  
    xlab(bquote(hat(beta)[.(colnames(df)[1])])) + 
    ylab(bquote(hat(beta)[.(colnames(df)[2])]))  + 
    labs(color=bquote(hat(theta)["interaction"]))  
   
  
  return(hw_plot)
}



scatterplt_all <- function(channel = "Ch1", C = coef_matrix_weak_10uM_21features,
                           textsize = 1.8, max.overlaps = 100, 
                           point.size = 1.8){
  
  channel_id <- substr(rownames(C), 
                       nchar(rownames(C)) - 2,
                       nchar(rownames(C)))
  ch_id <- which(channel_id == channel)
  coef_int_ch <- C[ch_id, 21:ncol(C)]
  
  #coef_int <- coef_matrix[, index_int]
  C_tch <- t(C[ch_id, ])
  
  coef_mat3ch <- matrix(nrow = 3, 
                        ncol = sum(
                          C_tch[
                            21:nrow(C_tch),] != 0))
  
  rownames(coef_mat3ch) <- c("comp. 1", "comp. 2", "Int") 
  # export_csv <- matrix(nrow = 54, ncol = 6)
  index_interaction_proteins <- which(rowSums(coef_int_ch)==0)
  n <- 0
  col_namesch <- c()
  for(prot in 1:ncol(C_tch)){
    
    for(x in 1:19){
      for(y in (x + 1):20){
        mod1 <- rownames(C_tch)[x]
        mod2 <- rownames(C_tch)[y]
        
        if(paste0(mod1, ":", mod2) %in% rownames(C_tch)){
          mod1mod2 <- paste0(mod1, ":", mod2)
        }
        else{mod1mod2 <- paste0(mod2, ":", mod1)}
        
        
        
        if(C_tch[mod1mod2, prot] != 0){
          n <- n + 1
          
          coef_mat3ch[, n] <- c(C_tch[mod1, prot], 
                                C_tch[mod2, prot], 
                                C_tch[mod1mod2, prot])
          
          col_namesch[n] <- colnames(C_tch)[prot]
          
          
        }
      }
      
    }
    
    
  }
  
  
  colnames(coef_mat3ch) <- col_namesch
  
  df <- as.data.frame(t(coef_mat3ch))
  col_names <- colnames(coef_mat3ch)
  df$label <- col_namesch
  library(scales)
  library(ggrepel)
  plot_int <- ggplot(df, aes(x = `comp. 1`, y = `comp. 2`,
                             colour = Int#,
                             #label = rownames(df)
  )) +
    
    theme_minimal() +
    theme(axis.title=element_text(size=15)
    ) +
    theme(
      legend.title=element_text(size=15)) +
    geom_vline(xintercept = 0, color = "grey") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_point(size = point.size) +
    scale_colour_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      midpoint = 0,
      breaks = seq(100, from = -5, to = 4)
    )  +
    
    geom_text_repel(data = df,
                    aes(x = `comp. 1`, y = `comp. 2`,
                        colour = "black",
                        label = label
                    ), size = textsize, color = "black", max.overlaps = max.overlaps)
  
  
  
  plot_int + 
    xlab(bquote(hat(beta)["comp. 1"])) + 
    ylab(bquote(hat(beta)["comp. 2"]))  + 
    labs(colour=bquote(hat(theta)["Interaction"])) 
}
