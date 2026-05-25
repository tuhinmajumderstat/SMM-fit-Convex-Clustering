#Load packages
library(clues)
library(combinat)
library(fossil)
library(geometry)
library(gstat)
library(LaplacesDemon)
library(plyr)
library(sets)
library(spatstat)
library(stringr)
library(mclust)




#### Set up SMM Object ####


#DEFINE SMM OBJECT AND FUNCTIONS TO MANIPULATE IT

SMM <- setRefClass("SMM",
                   fields = c("subset_belong", "subset_probs", "order", "sigma"),
                   methods = list(
                     
                     #Find probabilities from a history
                     find_probs = function(history){
                       history = history[(length(history) - order + 1):length(history)]
                       asstring = paste(history, sep="", collapse="")
                       subset = subset_belong[[asstring]]
                       return(subset_probs[subset,])
                     },
                     
                     #Generate next step given a history
                     next_step = function(history){
                       probs = find_probs(history)
                       index = sample.int(n=length(sigma), size=1, prob=probs)
                       return(sigma[index])
                     },
                     
                     #Generate n steps ahead from a history and return last observation
                     generate_n_steps = function(history, n){
                       history = strsplit(history, "")[[1]]
                       for(i in 1:n){
                         a = next_step(history)
                         history = c(history, a)
                       }
                       return(a)
                     },
                     
                     #Generate n steps ahead from a history and return all generated values
                     simult_n_steps = function(history, n){
                       history = strsplit(history, "")[[1]]
                       len = length(history)
                       output = rep("", n)
                       for(i in 1:n){
                         a = next_step(history)
                         history = c(history[2:length(history)], a)
                         output[i] = a
                       }
                       return(output)
                     },
                     
                     #Run bootstrap iterations times n steps ahead
                     bootstrap_sample = function(history, n, iterations, variety="n_ahead"){
                       retvec = rep("", iterations)
                       
                       if(variety == "n_ahead"){
                         for(i in 1:iterations){
                           retvec[i] = generate_n_steps(history, n)
                         }
                       }
                       
                       else if(variety == "simultaneous"){
                         for(i in 1:iterations){
                           retvec[i] = paste(simult_n_steps(history, n), sep="", collapse="")
                         }
                       }
                       
                       else{
                         stop("variety must be either n_ahead or simultaneous")
                       }
                       
                       return(retvec)
                     },
                     
                     #Calculate probability of a sequence
                     sequence_prob = function(sequence){
                       n = length(sequence)
                       histories = rep('', n - order)
                       probs = rep(0, n - order)
                       outcome_nums = rep(0, n-order)
                       for(i in 1:(n - order)){
                         histories[i] = paste(sequence[i:(i+order-1)], sep='', collapse='')
                         outcome_nums[i] = which(sigma == sequence[i+order])
                       }
                       for(i in 1:(n - order)){
                         probs[i] = subset_probs[subset_belong[[histories[i]]],outcome_nums[i]]
                       }
                       return(exp(mean(log(probs))))
                     }
                   ))

#Create a SMM of order with at most q clusters and alphabet size sigma
random_smm = function(order, q, sigma){
  hists = hcube(rep(sigma, order))
  assignments = sample(1:q, sigma^order, replace=TRUE)
  
  subset_belong = list()
  for(i in 1:(sigma^order)){
    subset_belong[[paste(as.character(hists[i,]), sep="", collapse="")]] = assignments[i]
  }
  rm(hists)
  
  subset_probs = matrix(rep(0, q*sigma), nrow=q)
  for(i in 1:q){
    rel_weight = sample(1:1000, sigma)
    subset_probs[i,] = rel_weight/sum(rel_weight)
  }
  
  
  return(SMM$new(subset_belong = subset_belong, subset_probs = subset_probs,
                 order=order, sigma = as.character(1:sigma)))
}

#Specify number, size, and probabilities associated with clusters and generate SMM
large_smm = function(order, q, sizes, sigma, subset_probs){
  hists = hcube(rep(sigma, order))
  
  #Generate random assignment of histories to classes
  numbers = rep(1, sizes[1])
  for(i in 2:q){
    numbers = c(numbers, rep(i, sizes[i]))
  }
  group_assignments = sample(numbers)
  
  subset_belong = list()
  
  for(i in 1:nrow(hists)){
    subset_belong[[paste(as.character(hists[i,]), sep="", collapse="")]] = group_assignments[i]
  }
  
  output_smm = SMM$new(subset_belong = subset_belong, subset_probs = subset_probs,
                       order=order, sigma = as.character(1:sigma))
}

#Reformat history into list with counts
reformat = function(history, order, sigma){
  
  #Find counts for each history/next step pair
  historylist = list()
  
  for(i in 1:(length(history)-order)){
    short_hist = history[i:(i+order-1)]
    short_hist = paste(short_hist, sep="", collapse="")
    next_char = history[i+order]
    index = which(sigma == next_char)
    
    if(short_hist %in% names(historylist)){
      historylist[[short_hist]][index] = (historylist[[short_hist]][index])+1
    }
    else{
      historylist[[short_hist]] = rep(0, length(sigma))
      historylist[[short_hist]][index] = 1
    }
  }
  return(list(counts=historylist))
}

#Obtain correct class labels
correct_class <- function(smm, histories){
  correct = rep(0, length(histories))
  for(zs in 1:length(histories)){
    correct[zs] = smm$subset_belong[[histories[zs]]]
  }
  return(correct)
}

#Create counts matrix
counts_matrix <- function(countslist, histories){
  sigmacard = length(countslist[[histories[1]]])
  numhists = length(countslist)
  counts = matrix(rep(0, sigmacard*numhists), nrow = numhists)
  for(zs in 1:numhists){
    counts[zs,] = countslist[[histories[zs]]]
    
  }
  return(counts)
}


#### GSDPMM ####

run_dpmm = function(histories, runs, order, sigma,
                    I, alpha, largemodel=FALSE){
  output_matrix = matrix(rep(0, (length(sigma)^order)*runs),
                         nrow=runs)
  V = length(sigma)
  
  for(run in 1:runs){
    history = strsplit(histories[run], "")[[1]]
    countslist = reformat(history, order, sigma)$counts
    
    histcounts = counts_matrix(countslist, names(countslist))
    Ns = rowSums(histcounts)
    Dnum = nrow(histcounts)
    
    #Create matrix to store assignments
    assignments_by_iter = matrix(rep(0, Dnum*I), nrow=I)
    
    #Initialize
    assigned = c(1, rep(0, Dnum-1))
    clusterdat = list(list(numhists = 1, counts = histcounts[1,]))
    for(i in 2:Dnum){
      history_counts = histcounts[i,]
      Nd = Ns[i]
      
      assignment_probs = rep(0, max(assigned)+1)
      term1s = rep(0, max(assigned)+1)
      term2s = rep(0, max(assigned)+1)
      term3s = rep(0, max(assigned)+1)
      
      correction = 0
      
      for(k in 1:max(assigned)){
        term1s[k] = clusterdat[[k]]$numhists/(i - 1 + alpha)
        term2s[k] = sum(lgamma(clusterdat[[k]]$counts + history_counts + 1)) -
          sum(lgamma(clusterdat[[k]]$counts + 1))
        term3s[k] = lgamma(sum(clusterdat[[k]]$counts) + V + Nd) -
          lgamma(sum(clusterdat[[k]]$counts) + V)
        while(exp(term2s[k] - term3s[k] + correction) == 0){
          correction = correction + 100
        }
      }
      
      #New cluster
      term1s[k+1] = alpha/(i - 1 + alpha)
      term2s[k+1] = sum(lgamma(history_counts + 1))
      term3s[k+1] = lgamma(V+Nd) - lgamma(V)
      while(exp(term2s[k+1] - term3s[k+1] + correction) == 0){
        correction = correction + 100
      }
      
      assignment_probs = term1s*exp(term2s-term3s+correction)
      
      #Assign to intial cluster
      new_cluster = which.max(assignment_probs)
      assigned[i] = new_cluster
      
      if(new_cluster <= k){
        clusterdat[[new_cluster]]$numhists = clusterdat[[new_cluster]]$numhists + 1
        clusterdat[[new_cluster]]$counts = clusterdat[[new_cluster]]$counts + histcounts[i,]
      }
      
      else{
        clusterdat[[new_cluster]] = list(numhists = 1, counts = histcounts[i,])
      }
    }
    
    #For I iterations, reassign all histories one by one
    for(i in 1:I){
      for(j in 1:Dnum){
        
        #Get data on history
        history_counts = histcounts[j,]
        Nd = Ns[j]
        
        #Remove history from its current cluster
        current_cluster = assigned[j]
        clusterdat[[current_cluster]]$numhists = clusterdat[[current_cluster]]$numhists - 1
        clusterdat[[current_cluster]]$counts = clusterdat[[current_cluster]]$counts - histcounts[j,]
        
        #If removing that element made the cluster empty, remove empty cluster and relabel
        if((clusterdat[[current_cluster]]$numhists == 0) & (current_cluster < max(assigned))){
          for(k in (current_cluster+1):max(assigned)){
            assigned[which(assigned == k)] = k-1
            clusterdat[[k-1]] = clusterdat[[k]]
          }
          clusterdat[k] = NULL
        }
        else if(clusterdat[[current_cluster]]$numhists == 0){
          clusterdat[k] = NULL
        }
        
        assigned[j] = 0
        
        #Find probabilities for each remaining cluster
        
        assignment_probs = rep(0, max(assigned)+1)
        term1s = rep(0, max(assigned)+1)
        term2s = rep(0, max(assigned)+1)
        term3s = rep(0, max(assigned)+1)
        
        correction = 0
        for(k in 1:max(assigned)){
          term1s[k] = clusterdat[[k]]$numhists/(Dnum - 1 + alpha)
          term2s[k] = sum(lgamma(clusterdat[[k]]$counts + history_counts + 1)) -
            sum(lgamma(clusterdat[[k]]$counts + 1))
          term3s[k] = lgamma(sum(clusterdat[[k]]$counts) + V + Nd) -
            lgamma(sum(clusterdat[[k]]$counts) + V)
          
          while(exp(term2s[k] - term3s[k] + correction) == 0){
            correction = correction + 100
          }
        }
        
        #Probability of new cluster
        term1s[k+1] = alpha/(Dnum - 1 + alpha)
        term2s[k+1] = sum(lgamma(history_counts + 1))
        term3s[k+1] = lgamma(V+Nd) - lgamma(V)
        while(exp(term2s[k+1]-term3s[k+1]+correction) == 0){
          correction = correction + 100
        }
        
        assignment_probs = term1s*exp(term2s-term3s+correction)
        
        #Find new cluster assignment
        new_cluster = sample.int(length(assignment_probs), size=1, prob = assignment_probs)
        assigned[j] = new_cluster
        
        if(new_cluster <= k){
          clusterdat[[new_cluster]]$numhists = clusterdat[[new_cluster]]$numhists + 1
          clusterdat[[new_cluster]]$counts = clusterdat[[new_cluster]]$counts + histcounts[j,]
        }
        
        else{
          clusterdat[[new_cluster]] = list(numhists = 1, counts = histcounts[j,])
        }
      }
      
      reformat_clusters = rep(0, Dnum)
      unique_clust = unique(assigned)
      for(l in 1:length(unique_clust)){
        reformat_clusters[which(assigned ==
                                  unique_clust[l])] = l
      }
      
      #Save clustering
      assignments_by_iter[i,] = reformat_clusters
    }
    
    #Find most common clustering
    frequencies = plyr::count(assignments_by_iter)
    clusters = unlist(unname(frequencies[which.max(frequencies$freq),1:Dnum]))
    
    if(length(clusters) < (V^order)){
      clusters = c(clusters, rep(0, V^order-
                                   length(clusters)))
    }
    output_matrix[run,] = clusters
    
  }
  return(list("output_matrix" = output_matrix))#, "histcounts" = histcounts))
}


#### Xiong et al. ####


#Bayes factor
bayes_factor = function(countsa, countsb, alpha, q){
  first_term = sum(lgamma(alpha*q)) - lgamma(alpha)
  second_term = sum(lgamma(alpha*q + countsa + countsb)) -
    lgamma(sum(alpha*q + countsa + countsb))
  third_term = lgamma(sum(countsa + alpha*q)) -
    sum(lgamma(countsa + alpha*q))
  fourth_term = lgamma(sum(countsb + alpha*q)) -
    sum(lgamma(countsb + alpha*q))
  return(first_term + second_term + third_term + fourth_term)
}


qhull.options <- function(options = NULL, output.options = NULL, supported_output.options = c("Fa", "Fn"), full = FALSE) {
  if (is.null(options)) options <- ""
  if (is.null(output.options)) output.options <- ""
  
  # Validate output options
  invalid <- setdiff(strsplit(output.options, " ")[[1]], supported_output.options)
  if (length(invalid) > 0) {
    stop("Unsupported output option(s): ", paste(invalid, collapse = ", "))
  }
  
  final <- paste(options, output.options)
  if (full) final <- paste(final, "Fd Fi Fv Fn")  # or other flags for full output
  return(final)
}

delaunayn_updated <- function (p, options = NULL, output.options = NULL, full = FALSE) 
{
  tmp_stdout <- tempfile("Rf")
  tmp_stderr <- tempfile("Rf")
  on.exit(unlink(c(tmp_stdout, tmp_stderr)))
  if (is.data.frame(p)) {
    p <- as.matrix(p)
  }
  storage.mode(p) <- "double"
  if (any(is.na(p))) {
    stop("The first argument should not contain any NAs")
  }
  default.options <- "Qt Qc Qx"
  if (ncol(p) < 4) {
    default.options <- "Qt Qc Qz"
  }
  if (is.null(options)) {
    options <- default.options
  }
  options <- tryCatch(qhull.options(options, output.options, 
                                    supported_output.options <- c("Fa", "Fn"), full = full), 
                      error = function(e) {
                        stop(e)
                      })
  if (!grepl("Qt", options) & !grepl("QJ", options)) {
    options <- paste(options, "Qt")
  }
  out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout, 
               tmp_stderr, PACKAGE = "geometry")
  if (nrow(out$tri) > 0) {
    missing.points <- length(setdiff(seq(1, nrow(p)), unique(as.vector(out$tri))))
    if (missing.points > 0) {
      warning(paste0(missing.points, " points missing from triangulation.\nIt is possible that setting the 'options' argument of delaunayn may help.\nFor example:\noptions = \"", 
                     default.options, " Qbb\"\noptions = \"", default.options, 
                     " QbB\"\nIf these options do not work, try shifting the centre of the points\nto the origin by subtracting the mean coordinates from every point."))
    }
  }
  out[which(sapply(out, is.null))] <- NULL
  if (is.null(out$areas) & is.null(out$neighbours)) {
    attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
    return(out$tri)
  }
  class(out) <- "delaunayn"
  out$p <- p
  return(out)
}


run_xiong = function(histories, runs, order, sigma){
  output_matrix = matrix(rep(0, (length(sigma)^order)*runs),
                         nrow=runs)
  for(run in 1:runs){
    #Get data from histories
    history = strsplit(histories[run], "")[[1]]
    a = reformat(history, order, sigma)
    count_data = a$counts
    
    #Reformat as matrix
    columns = length(count_data[[1]])
    num_hists = length(count_data)
    
    #Set uniform prior on transition probs
    alpha = columns
    q = rep(1/columns, columns)
    
    clust_counts = matrix(rep(0, num_hists*columns), nrow=num_hists)
    sample_probs = matrix(rep(0, num_hists*columns), nrow=num_hists)
    clustering = list()
    
    for(i in 1:num_hists){
      clust_counts[i,] = count_data[[i]]
      sample_probs[i,] = floor(1000000*count_data[[i]]/sum(count_data[[i]]))
      clustering[[i]] = i
    }
    
    #Perform Delaunay triangulation
    
    triangulation = delaunayn_updated(sample_probs, options = "QJ")
    
    #Reformat triangulation into adjacency matrix
    reformatted_tri = matrix(rep(0, num_hists^2),
                             nrow = num_hists)
    for(i in 1:nrow(triangulation)){
      rowvals = triangulation[i,]
      for(j in 1:(ncol(triangulation)-1)){
        for(k in (j+1):ncol(triangulation)){
          reformatted_tri[rowvals[j], rowvals[k]] = 1
          reformatted_tri[rowvals[k], rowvals[j]] = 1
        }
      }
    }
    
    #Calculate bfs between clusters/histories only if
    #they are adjacent in triangulation
    bf = matrix(rep(0, num_hists^2), nrow=num_hists)
    for(i in 1:(num_hists - 1)){
      for(j in (i+1): num_hists){
        if(reformatted_tri[i,j] == 1){
          bf[j,i] = bayes_factor(clust_counts[i,], clust_counts[j,], alpha, q)
        }
      }
    }
    
    #Find max (nonzero) value in bf matrix and combine clusters
    
    while(length(which(bf > 0)) > 0){
      maxbf = max(bf)
      to_combine = which(bf == maxbf)[1]
      column = to_combine%/%num_hists
      if(to_combine %% num_hists != 0){
        column = column + 1
      }
      row = to_combine%%num_hists
      if(row == 0){
        row = num_hists
      }
      
      #Combine these groups
      clustering[[column]] = c(clustering[[column]], clustering[[row]])
      clust_counts[column,] = clust_counts[row,] + clust_counts[column,]
      clustering[[row]] = 0
      clust_counts[row,] = rep(0, ncol(clust_counts))
      
      reformatted_tri[column,] = reformatted_tri[row,] + reformatted_tri[column,]
      reformatted_tri[column, column] = 0
      reformatted_tri[row,] = rep(0, num_hists)
      reformatted_tri[,row] = rep(0, num_hists)
      for(j in 1:num_hists){
        if(reformatted_tri[column,j] > 0){
          bf[j, column] = bayes_factor(clust_counts[column,], clust_counts[j,],
                                       alpha, q)
        }
      }
      bf[row, ] = rep(0, num_hists)
      bf[, row] = rep(0, num_hists)
      bf[column, column] = 0
    }
    
    #Reformat into list
    clusters = rep(0, num_hists)
    for(i in 1:num_hists){
      for(history in clustering[[i]]){
        if(history > 0){
          clusters[history] = i
        }
      }
    }
    if(length(clusters) < (length(sigma)^order)){
      clusters = c(clusters, rep(0, length(sigma)^order-
                                   length(clusters)))
    }
    output_matrix[run,] = clusters
  }
  return(output_matrix)
}



#### GGL (Garcia & Gonzalez-Lopez) ####


#Compute d(i,j)
distance_combined = function(i, j, L, counts){
  rows_i = unlist(L[[i]])
  rows_j = unlist(L[[j]])
  
  if(length(rows_i) > 1){
    i_counts = colSums(counts[rows_i,])
  }
  else{
    i_counts = counts[rows_i,]
  }
  
  if(length(rows_j) > 1){
    j_counts = colSums(counts[rows_j,])
  }
  else{
    j_counts = counts[rows_j,]
  }
  overall_counts = i_counts + j_counts
  
  n = sum(counts)
  
  dl_1 = sum(i_counts*log(i_counts/sum(i_counts)))
  dl_2 = sum(j_counts*log(j_counts/sum(j_counts)))
  dl_3 = sum(overall_counts*log(overall_counts/sum(overall_counts)))
  
  #Replaces NaNs caused by 0 counts with 0's
  dl_1[which(is.nan(dl_1))] = 0
  dl_2[which(is.nan(dl_2))] = 0
  dl_3[which(is.nan(dl_3))] = 0
  return((1/log(n))*(dl_1 + dl_2 - dl_3))
}

#Run BIC Method
ggl_method = function(L, counts, k){
  
  i = 0
  j = 1
  
  while(i < (k-1)){
    i = i+1
    j = i
    while(j < k){
      j = j+1
      d = distance_combined(i, j, L, counts)
      
      while(d < 1){
        L[[i]] = set_union(L[[i]], L[[j]])
        L[[j]] = NULL
        k = k-1
        
        if(j <= k){
          d = distance_combined(i, j, L, counts)
        }
        else{
          break
        }
        
      }
    }
    
  }
  return(L)
}

#Create list with partition
partition_list <- function(k){
  L = list(set(1))
  for(i in 2:k){
    L[[i]] = set(as.numeric(i))
  }
  return(L)
}

#Convert from partition to clustering
obtain_clustering <- function(L){
  clustering = rep(0, length(L))
  for(i in 1:length(lengths(L))){
    clustering[unlist(L[[i]])] = i
  }
  return(clustering)
}

#Run for multiple chains, output as matrix
run_ggl = function(histories, runs, order, sigma){
  
  output_matrix = matrix(rep(0, (length(sigma)^order)*runs),
                         nrow=runs)
  for(run in 1:runs){
    
    history = strsplit(histories[run], "")[[1]]
    countslist = reformat(history, order, sigma)$counts
    X = counts_matrix(countslist, names(countslist))
    
    
    #GGL
    L = partition_list(nrow(X))
    L = ggl_method(L, X, nrow(X))
    
    clustering = obtain_clustering(L)
    
    #adjust for chains where not all histories appear
    if(length(clustering) < (length(sigma)^order)){
      clustering = c(clustering, rep(0, length(sigma)^order-
                                       length(clustering)))
    }
    output_matrix[run,] = clustering
  }
  return(output_matrix)
}



#### K-means ####


#Calculate BIC for SMM
give_bic <- function(countslist, clustering, order){
  #Number of params estimated
  num_clusters = max(clustering)
  num_outputs = length(countslist[[1]])
  num_params = num_clusters*(num_outputs - 1)
  
  #Find modified maximum likelihood (Garcia and Gonzalez-Lopez)
  cluster_probs = matrix(rep(0,num_clusters*num_outputs), nrow = num_clusters)
  counts = matrix(rep(0,num_clusters*num_outputs), nrow = num_clusters)
  for(i in 1:num_clusters){
    histories = which(clustering == i)
    counts[i,] = Reduce('+', countslist[histories])
    cluster_probs[i,] = counts[i,]/sum(counts[i,])
  }
  
  
  ml = 0
  for(i in 1:num_clusters){
    if(0 %in% counts[i,]){
      rowcount = counts[i, which(counts[i,] != 0)]
      rowprobs = cluster_probs[i, which(counts[i,] != 0)]
      ml = ml + sum(rowcount*log(rowprobs))
    }
    else{
      ml = ml + sum(counts[i,]*log(cluster_probs[i,]))
    }
  }
  
  #Penalty term
  n = sum(counts) + order
  penalty_term = .5*num_params*log(n)
  
  return(ml - penalty_term)
}

run_kmeans = function(histories, runs, order, sigma){
  output_matrix = matrix(rep(0, (length(sigma)^order)*runs),
                         nrow=runs)
  for(run in 1:runs){
    history = strsplit(histories[run], "")[[1]]
    a = reformat(history, order, sigma)
    
    #Accepts as input data from reformat
    count_data = a$counts
    
    #Reformat as matrix
    columns = length(count_data[[1]])
    num_hists = length(count_data)
    coordinates = matrix(rep(0, num_hists*columns), nrow=num_hists)
    for(i in 1:num_hists){
      counts = count_data[[i]]
      raw_coords = counts/sum(counts)
      
      #0 correction
      zeroes = which(raw_coords == 0)
      raw_coords = raw_coords*
        (1 - (length(zeroes)/(sum(counts) + length(counts))))
      raw_coords[zeroes] = 1/(sum(counts) + length(counts))
      
      #CLR transformation
      geo_mean = prod(raw_coords)^(1/length(raw_coords))
      raw_coords = log(raw_coords/geo_mean)
      coordinates[i,] = raw_coords
    }
    
    
    
    #Run kmeans for all k, choose best using BIC
    bics = rep(0, nrow(unique(coordinates))-1)
    for(k in 1:(nrow(unique(coordinates))-1)){
      clustering=kmeans(coordinates, centers=k, iter.max=1000)$cluster
      bics[k] = give_bic(count_data, clustering, order)
    }
    k = which.max(bics)
    clustering = kmeans(coordinates, centers=k, iter.max=1000)$cluster
    
    #adjust for chains where not all histories appear
    if(length(clustering) < (length(sigma)^order)){
      clustering = c(clustering, rep(0, length(sigma)^order-
                                       length(clustering)))
    }
    output_matrix[run,] = clustering
  }
  return(output_matrix)
}




#### Comparison Metric ####

#Find correct clustering
correct_matrix = function(histories, models){
  output_matrix = matrix(rep(0, length(histories)*
                               (length(models$sigma)^models$order)),
                         nrow=length(histories))
  for(i in 1:length(histories)){
    history = strsplit(histories[i], "")[[1]]
    a = reformat(history, models$order, models$sigma)
    correct = correct_class(models, names(a$counts))
    if(length(correct) < (length(models$sigma)^models$order)){
      correct = c(correct, rep(0, length(models$sigma)^models$order -
                                 length(correct)))
    }
    output_matrix[i,] = correct
  }
  return(output_matrix)
}



#________________________________________________
#Compute adjusted rand indices for assigned clusterings vs correct model
compute_rand = function(assigned_matrix, correct_matrix){
  output = rep(0, nrow(assigned_matrix))
  for(i in 1:nrow(assigned_matrix)){
    #Remove histories that did not appear in simulation
    assigned = assigned_matrix[i,which(assigned_matrix[i,] > 0)]
    correct = correct_matrix[i,which(correct_matrix[i,] > 0)]
    if(length(correct) != length(assigned)){
      print(i)
    }
    output[i] = adjustedRand(assigned, correct, "HA")
  }
  return(output)
}

