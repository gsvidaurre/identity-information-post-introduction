# G. Smith-Vidaurre
# 20 January 2023

# Description: This function calculates Earth Mover's Distance (emdist::emd) for distributions of similarity values within or among categories (e.g. individuals or sites) at a given social scale. The function uses a histogram-based approach to create signatures compared between distributions (the analogy is to think of these signatures per distribution as piles of dirt). The function can iterate over different numbers of bins to calculate Earth Mover's Distance, or the minimum cost needed to turn one distribution into the other. The function also randomly samples an equivalent number of similarity values from each distribution, and this resampling can be performed over various iterations. This function was written for the datasets used in the accompanying manuscript, and will likely need to be updated before use with other datasets.

# Arguments:

# sim_df: A data frame representing pairwise comparisons within and among categories (e.g. individuals, sites) at the given social scale. This data frame should be the result of extract_simValues().

# social_scale: String. The given social scale represented by these pairwise comparisons. This argument is used to inform how the within and among pairwise comparisons are extracted. The possible values here are "Individual", "Site", and "Regional" to facilitate using the function for yellow-naped amazon contact calls as well. These social scales are used with the species to determine which pairwise comparisons will be returned 

# range: String. The range (native or introduced) of the given samples in the object sim_df. This argument is used solely as metadata in the output of this function

# dataset: String. The site scale dataset representing different methods of accounting for repeated sampling of unmarked individuals. This argument is used solely as metadata in the output of this function

# similarity_method: String. The similarity method used to generate the similarity matrix above. This argument is used solely as metadata in the output of this function

# city_year: String. The city and year of calls represented in the similarity matrix, relevant only for temporal analyses later. This argument is used solely as metadata in the output of this function

# distr_types: A character vector with two elements representing the two distributions for which Earth Mover's Distance will be calculated. The default value is c("within", "among")

# total_bins: A numeric vector with one or more elements specifying the total number of bins that will be used to bin each distribution between the bounds below. Set this to a single value if you are confident that the Earth Mover's Distance calculation doesn't change much across different numbers of total bins

# bounds: A numeric vector with two elements specifying the upper and lower bounds for the histogram bins. The default is c(0, 1), which represents the lower and upper limits of similarity values along a single dimension.

# Details: 
# In this histogram-based approach, each distribution is binned using the given number of total bins between the specified bounds. For each bin, the function calculates the center of each bin as the mean similarity value per bin, and the weight as the proportion of values assigned to each bin. The bin weights and means are used to calculate Earth Mover's Distance using the emdist package (https://cran.r-project.org/web/packages/emdist/index.html). To build this function I also used the colordistance package as a reference for how to calculate this metric through a histogram-based approach (https://cran.r-project.org/web/packages/colordistance/index.html). The resulting Earth Mover's Distance represents the amount of work needed to convert the distribution of similarity values within categories to the distribution of values among categories along a single dimension of similarity (bounded between 0 and 1). In the accompanying manuscript, we use these Earth Mover's Distance calculations as estimates of acoustic convergence at each social scale, such that greater Earth Mover's Distances represent more convergence (e.g. more work is need to turn one distribution into the other, indicating that calls converged more within than among a given category).

# The function uses emdist::emd with a maximum of 1000 iterations and Euclidean distance for calculating costs.

# Value: This function returns a data frame of similarity values representing all unique pairwise comparisons extracted from the original similarity matrix. The data frame also contains various metadata columns about each pairwise comparison.

resampled_earthMoversDistance <- function(sim_df, social_scale, range, dataset, similarity_method, city_year, distr_types = c("within", "among"), total_bins, bounds){
  
  # Find the number of values in the distribution of "within" comparisons and randomly sample that same number from the "among" distribution without replacement
  nw <- sim_df %>% 
    dplyr::filter(type == "within") %>% 
    nrow()
  
  rs_sim_df <- sim_df %>% 
    dplyr::filter(type == "within") %>% 
    bind_rows(
      sim_df %>% 
        dplyr::filter(type == "among") %>% 
        sample_n(., size = nw, replace = FALSE) %>% 
        ungroup()
    )
  
  # Iterate over bins
  res_df3 <- data.table::rbindlist(lapply(1:length(total_bins), function(b){
    
    # Create a vector of bins for the given number of total bins. From 0 to 1 or the absolute min and max possible in SPCC similarity space.
    # Using code in colordistance::getImageHist as a baseline (https://github.com/hiweller/colordistance/blob/master/R/02a_histogram_color_clustering.R)
    
    breaks <- lapply(total_bins[b] + 1, function(x){
      seq(bounds[1], bounds[2], length = x)
    })
    
    # Create a new data frame for each distribution, in which you get the average similarity value per bin (e.g. cluster) as well as the percentage of samples in each bin
    # Account for the fact that among all of the original bins, some of these may not have been filled. Also needs to return the bins in order
    bin_inds <- seq(1, total_bins[b], 1)
    
    res_df <- data.table::rbindlist(lapply(1:length(distr_types), function(j){
      
      dt <- distr_types[j]
      
      # Then cut each distribution into these bins
      binned_df <- data.frame(
        sim_vals = rs_sim_df %>% 
          dplyr::filter(type == dt) %>% 
          pull(similarity_value),
        binned = cut(
          rs_sim_df %>% 
            dplyr::filter(type == dt) %>% 
            pull(similarity_value),
          breaks = breaks[[1]],
          include.lowest = TRUE,
          labels = FALSE
        )
      )
      
      # For each bin, get the mean value and weight (e.g. the proportion of values assigned to that bin), and return these values in order across bins
      res_df2 <- data.table::rbindlist(lapply(1:length(bin_inds), function(i){
        
        bin <- bin_inds[i]
        
        bin_tmp_df <- binned_df %>% 
          dplyr::filter(binned == bin) 
        
        if(nrow(bin_tmp_df) > 0){
          
          bin_mean <- bin_tmp_df %>% 
            pull(sim_vals) %>% 
            mean(.)
          
          bin_weight <- nrow(bin_tmp_df)/nrow(binned_df)
          
        } else {
          
          bin_mean <- 0
          bin_weight <- 0
          
        }
        
        # Return the mean similarity value and weight per bin
        # Make sure the data frame is arranged by bin in descending order
        return(
          data.frame(
            distr_type = dt,
            bin = bin,
            bin_mean = bin_mean,
            bin_weight = bin_weight
          ) %>% 
            dplyr::arrange(-desc(bin))
        )
        
      }))
      
      return(res_df2)
      
    }))

    mat_A <- as.matrix(
      res_df %>% 
        dplyr::filter(distr_type == distr_types[1]) %>% 
        dplyr::select(bin_weight, bin_mean)
    )
    
    mat_B <- as.matrix(
      res_df %>% 
        dplyr::filter(distr_type == distr_types[2]) %>% 
        dplyr::select(bin_weight, bin_mean)
    )
    
    emd_res <- emdist::emd(
      A = mat_A, 
      B = mat_B,
      dist = "euclidean",
      max.iter = 5000
    )
    
    emd_df <- data.frame(
      social_scale = unique(rs_sim_df$social_scale),
      range = unique(rs_sim_df$range),
      total_bins = total_bins[b],
      n_within = rs_sim_df %>% 
        dplyr::filter(type == "within") %>% 
        nrow(),
      n_among = rs_sim_df %>% 
        dplyr::filter(type == "among") %>% 
        nrow(),
      EMD = emd_res,
      dataset = unique(rs_sim_df$dataset), 
      similarity_method = unique(rs_sim_df$similarity_method), 
      city_year = unique(rs_sim_df$city_year)
    )
    
    return(emd_df)
    
  }))
  
  return(res_df3)
  
}
