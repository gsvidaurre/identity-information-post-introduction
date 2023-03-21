# G. Smith-Vidaurre
# 20 January 2023

# Description: This function extracts similarity values from a symmetric matrix, and returns a data frame of unique pairwise comparisons within and among categories at a given social scale (e.g. individual, site, regional dialect) for a given set of calls. This function was written for the datasets used in the accompanying manuscript, and will likely need to be updated for other datasets. Note that this function should be run separately for each social scale, range, dataset (site scale calls), city-year (temporal analyses), and similarity method

# Arguments:

# sim_mat: A symmetric matrix filled with similarity values representing pairwise comparisons

# calls: A character vector of all the calls for which pairwise comparisons should be extracted 

# species: String. The name of the species, either "monk parakeet" or "yellow-naped amazon", which will determine how metadata is added to the data frame of pairwise comparisons

# social_scale: String. The given social scale represented by these pairwise comparisons. This argument is used to inform how the within and among pairwise comparisons are extracted. The possible values here are "Individual", "Site", and "Regional" to facilitate using the function for yellow-naped amazon contact calls as well. These social scales are used with the species to determine which pairwise comparisons will be returned 

# dataset: String. The site scale dataset representing different methods of accounting for repeated sampling of unmarked individuals. This argument is used solely as metadata in the output of this function

# similarity_method: String. The similarity method used to generate the similarity matrix above. This argument is used solely as metadata in the output of this function

# city_year: String. The city and year of calls represented in the similarity matrix, relevant only for temporal analyses later. This argument is used solely as metadata in the output of this function

# metadata_df: Data frame. A data frame object with metadata for the calls used to extract similarity values

# analysis_type: String. This can be either "EMD" for Earth Mover's Distance, or "bootstrap" for bootstrapping analyses. This argument determines how comparisons at the individual scale will be filtered for monk parakeets (more strictly for the bootstrapping analyses)

# Value: This function returns a data frame of similarity values representing all unique pairwise comparisons extracted from the original similarity matrix. The data frame also contains various metadata columns about each pairwise comparison.

extract_simValues <- function(sim_mat, calls, species, social_scale, dataset, similarity_method, city_year, metadata_df, analysis_type){
  
  tmp_sim_mat <- sim_mat[grep(paste(paste("^", calls, "$", sep = ""), collapse = "|"), dimnames(sim_mat)[[1]]), grep(paste(paste("^", calls, "$", sep = ""), collapse = "|"), dimnames(sim_mat)[[2]])]
  
  # Melt the pairwise symmetric matrix into a longer data frame such that each row is a pairwise comparison (upper and lower triangle comparisons are all here so comparisons are duplicated, diagonal is included)
  tmp_df <- reshape2::melt(tmp_sim_mat, varnames = c("target_call", "other_call"), na.rm = TRUE)
  # glimpse(tmp_df)
  
  # Checking, looks good
  # nrow(tmp_df) == dim(tmp_sim_mat)[1] * dim(tmp_sim_mat)[2]
  
  # Remove values of 1, as these are values along the original diagonal of the matrix (e.g. each call compared to itself)
  tmp_df <- tmp_df %>% 
    dplyr::filter(value < 1)
  
  if(grepl("parakeet", species)){
    
    # Add metadata
    tmp_df2 <- tmp_df %>% 
      # Add metadata for the target calls
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df$target_call) %>%
          as_tibble() %>% 
          dplyr::select(sound.files, Bird_ID, site_year, range, social_scale),
        by = c("target_call" = "sound.files")
      ) %>% 
      dplyr::rename(
        `target_bird_ID` = "Bird_ID",
        `target_site` = "site_year",
        `target_range` = "range",
        `target_ss` = "social_scale"
      ) %>% 
      # Add metadata for the other calls used for comparison
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df$other_call) %>%
          as_tibble() %>%  
          dplyr::select(sound.files, Bird_ID, site_year, range),
        by = c("other_call" = "sound.files")
      ) %>% 
      dplyr::rename(
        `other_bird_ID` = "Bird_ID",
        `other_site` = "site_year",
        `other_range` = "range"
      ) %>% 
      # Add columns indicating whether the comparison was done within or among birds or sites
      # Site comparisons were set to same or different by range, such that the 3 sites where invasive range individuals were sampled were considered a single site for the purposes of comparing a similar number of individuals in each range at similarly small geographic scales
      dplyr::mutate(
        indiv_comparison_type = ifelse(target_bird_ID == other_bird_ID, "same", "different"),
        site_comparison_type = ifelse(target_site == other_site, "same", "different"),
        range_comparison_type = ifelse(target_range == other_range, "same", "different")
      ) %>% 
      dplyr::rename(
        `similarity_value` = "value"
      ) %>% 
      dplyr::select(
        indiv_comparison_type,
        site_comparison_type,
        range_comparison_type,
        similarity_value,
        target_call,
        other_call,
        target_bird_ID,
        other_bird_ID,
        target_site,
        other_site,
        target_range,
        other_range,
        target_ss
      )
    
    # Check that the social scale is the expected value
    if(unique(tmp_df2$target_ss) != social_scale){
      stop('This dataset does not match the expected social scale')
    }
    
  } else if(grepl("amazon", species)){
    
    # Add metadata
    tmp_df2 <- tmp_df %>% 
      # Add metadata for the target calls
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df$target_call) %>%
          as_tibble() %>% 
          dplyr::select(sound.files, bird_ID, site, regional_dialect),
        by = c("target_call" = "sound.files")
      ) %>% 
      dplyr::rename(
        `target_bird_ID` = "bird_ID",
        `target_site` = "site",
        `target_dialect` = "regional_dialect"
      ) %>% 
      # Add metadata for the other calls used for comparison
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df$other_call) %>%
          as_tibble() %>%   
          dplyr::select(sound.files, bird_ID, site, regional_dialect),
        by = c("other_call" = "sound.files")
      ) %>% 
      dplyr::rename(
        `other_bird_ID` = "bird_ID",
        `other_site` = "site",
        `other_dialect` = "regional_dialect"
      ) %>% 
      # Add columns indicating whether the comparison was done within or among birds, sites, or dialects
      dplyr::mutate(
        indiv_comparison_type = ifelse(target_bird_ID == other_bird_ID, "same", "different"),
        site_comparison_type = ifelse(target_site == other_site, "same", "different"),
        dialect_comparison_type = ifelse(target_dialect == other_dialect, "same", "different")
      ) %>% 
      dplyr::rename(
        `similarity_value` = "value"
      ) %>% 
      dplyr::select(
        indiv_comparison_type,
        site_comparison_type,
        dialect_comparison_type,
        similarity_value,
        target_call,
        other_call,
        target_bird_ID,
        other_bird_ID,
        target_site,
        other_site,
        target_dialect,
        other_dialect
      )
    
  }
  
  # Then get the within and among comparisons
  # To do this, create a string of conditions for filtering that is conditional on the given species and social scale
  if(grepl("parakeet", species)){
    
    if(social_scale == "Individual" & analysis_type == "EMD"){
      
      withn_comp <- "indiv_comparison_type == 'same' & range_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & range_comparison_type == 'same'"
      
    } else if(social_scale == "Individual" & analysis_type == "bootstrap"){
      
      withn_comp <- "indiv_comparison_type == 'same' & site_comparison_type == 'same' & range_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'same' & range_comparison_type == 'same'"
      
    } else if(social_scale == "Site"){
      
      withn_comp <- "site_comparison_type == 'same' & range_comparison_type == 'same'"
      
      among_comp <- "site_comparison_type == 'different' & range_comparison_type == 'same'"
      
    } 
    
  } else if(grepl("amazon", species)){

    # Individual scale comparison: same or different individual, same site, same dialect
    if(social_scale == "Individual"){
      
      withn_comp <- "indiv_comparison_type == 'same' & site_comparison_type == 'same' & dialect_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'same' & dialect_comparison_type == 'same'"
     
    # Site scale comparison: different individuals, same or different site, same dialect 
    } else if(social_scale == "Site"){
      
      withn_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'same' & dialect_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'different' & dialect_comparison_type == 'same'"
      
    # Regional dialect scale comparison: different individuals, different sites, same or different dialects
    } else if(social_scale == "Regional Dialect"){
      
      withn_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'different' & dialect_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & site_comparison_type == 'different' & dialect_comparison_type == 'different'"
      
    }
    
  }
  
  sim_vals <- tmp_df2 %>%
    dplyr::mutate(
      type = ifelse(!!rlang::parse_expr(withn_comp), "within", "drop")
    ) %>% 
    dplyr::mutate(
      type = ifelse(!!rlang::parse_expr(among_comp), "among", type)
    ) %>% 
    # For monk parakeets, this line will drop the values that represent comparisons between ranges. For all datasets, this line will drop comparisons that do not represent the specified social scale for the given species
    dplyr::filter(type != "drop") %>% 
    # Retain only values that represent unique pairwise comparisons
    # I tested this and found that the number of duplicates and non-duplicates were the same, as expected for a symmetric matrix
    dplyr::filter(!duplicated(paste0(pmax(target_call, other_call), pmin(target_call, other_call))))

  if(grepl("parakeet", species)){
    
    simVals_df <- sim_vals %>% 
      dplyr::mutate(
        dataset = dataset,
        similarity_method = similarity_method,
        city_year = city_year
      ) %>% 
      dplyr::rename(
        `range` = "target_range",
        `social_scale` = "target_ss"
      )
    
  } else if(grepl("amazon", species)){
    
    simVals_df <- sim_vals %>% 
      dplyr::mutate(
        dataset = dataset,
        similarity_method = similarity_method,
        city_year = city_year
      )
    
  }

  return(simVals_df)
  
}