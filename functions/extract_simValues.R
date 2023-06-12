# G. Smith-Vidaurre
# 20 January 2023

# This code accompanies the following paper in PLOS Computational Biology:

# Smith-Vidaurre, G., Perez-Marrufo, V., Hobson, E.A., Salinas-Melgoza, A., and T.F. Wright. Individual identity information persists in learned calls of introduced parrot populations.

# Description: This function extracts similarity values from a symmetric matrix, and returns a data frame of unique pairwise comparisons within and among categories at a given social scale (e.g. individual, site, regional dialect) for a given set of calls. This function was written for the datasets used in the accompanying manuscript and has not been generalized for other datasets. Note that this function should be run separately for each social scale, range, dataset (site scale calls), city-year (temporal analyses), and similarity method

# Arguments:

# sim_mat: A symmetric matrix filled with similarity values representing pairwise comparisons (values along the diagonal should be 1)

# calls: A character vector of all the calls for which similarity values from pairwise comparisons should be extracted 

# species: String. The name of the species, either "monk parakeet" or "yellow-naped amazon", which will determine how metadata is added to the data frame of pairwise comparisons

# social_scale: String. The given social scale represented by these pairwise comparisons. This argument is used to inform how the within and among pairwise comparisons are extracted. The possible values here are "Individual", "Site", and "Regional" to facilitate using the function for yellow-naped amazon contact calls as well. These social scales are used with the species to determine which pairwise comparisons will be returned 
# dataset: String. The site scale dataset representing different methods of accounting for repeated sampling of unmarked individuals for monk parakeet site scale datasets. This argument is used solely as metadata in the output of this function

# similarity_method: String. The similarity method used to generate the similarity matrix above. This argument is used solely as metadata in the output of this function

# city_year: String. The city and year of calls represented in the similarity matrix, relevant only for temporal analyses later. This argument is used solely as metadata in the output of this function

# metadata_df: Data frame. A data frame object with metadata for the calls used to extract similarity values

# analysis_type: String. This can be either "EMD" for Earth Mover's Distance, or "bootstrap" for bootstrapping analyses. This argument determines how comparisons at the individual scale will be filtered for monk parakeets (more strictly for the bootstrapping analyses)

# Value: This function returns a data frame of similarity values representing all unique pairwise comparisons extracted from the original similarity matrix. The data frame also contains various metadata columns the pairwise comparisons.

extract_simValues <- function(sim_mat, calls, species, social_scale, dataset, similarity_method, city_year, metadata_df, analysis_type){
  
  tmp_sim_mat <- sim_mat[grep(paste(paste("^", calls, "$", sep = ""), collapse = "|"), dimnames(sim_mat)[[1]]), grep(paste(paste("^", calls, "$", sep = ""), collapse = "|"), dimnames(sim_mat)[[2]])]
  # dim(tmp_sim_mat) == length(calls)
  
  # Melt the pairwise symmetric matrix into a longer data frame such that each row is a pairwise comparison (the upper and lower triangle comparisons are all here so comparisons are duplicated, also the diagonal is included)
  tmp_df <- reshape2::melt(tmp_sim_mat, varnames = c("target_call", "other_call"), na.rm = TRUE) %>% 
    dplyr::mutate(
      target_call = as.character(target_call),
      other_call = as.character(other_call)
    )
  # glimpse(tmp_df)
  
  # Checking, looks good
  # nrow(tmp_df) == dim(tmp_sim_mat)[1] * dim(tmp_sim_mat)[2]
  
  # # Remove values of 1, as these are values along the original diagonal of the matrix (e.g. each call compared to itself)
  tmp_df2 <- tmp_df %>%
    dplyr::filter(value < 1)

  # glimpse(tmp_df2)

  if(grepl("parakeet", species)){
    
    # Add metadata
    tmp_df3 <- tmp_df2 %>% 
      # Add metadata for the target calls
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df2$target_call) %>%
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
          dplyr::filter(sound.files %in% tmp_df2$other_call) %>%
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
    if(unique(tmp_df3$target_ss) != social_scale){
      stop('This dataset does not match the expected social scale')
    }
    
  } else if(grepl("amazon", species)){
    
    # Add metadata
    tmp_df3 <- tmp_df2 %>% 
      # Add metadata for the target calls
      inner_join(
        metadata_df %>% 
          dplyr::filter(sound.files %in% tmp_df2$target_call) %>%
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
          dplyr::filter(sound.files %in% tmp_df2$other_call) %>%
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
  
  # glimpse(tmp_df3)
  
  # Then get the within and among comparisons
  # To do this, create a string of conditions for filtering that is tailored for the given species and social scale
  
  # For the monk parakeet individual scale Earth Mover's Distance and bootstrapping analyses, I filtered pairwise comparisons in different ways
  if(grepl("parakeet", species)){
    
    # For the Earth Mover's Distance analyses, I set the conditional to same individual and same range, or different individual and same range, and did not filter by site comparisons. This was because I did not place geographic restrictions on the individuals used for these analyses, but rather a sample size restriction. So all repeatedly sampled individuals could be selected in random sampling for this analysis, as long as they had 5 calls or more (which excluded 1 native range bird, NAT-RAW)
    if(social_scale == "Individual" & analysis_type == "EMD"){
      
      withn_comp <- "indiv_comparison_type == 'same' & range_comparison_type == 'same'"
      
      among_comp <- "indiv_comparison_type == 'different' & range_comparison_type == 'same'"
      
      # For the bootstrapping analyses, I did place a geographic restriction on the individuals used and assigned a single site ID to the 3 sites at which the introduced individuals were sampled. So here I did include the site comparison filter in the conditional
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
  
  sim_vals <- tmp_df3 %>%
    dplyr::mutate(
      type = ifelse(!!rlang::parse_expr(withn_comp), "within", "drop")
    ) %>%
    dplyr::mutate(
      type = ifelse(!!rlang::parse_expr(among_comp), "among", type)
    ) %>%
    # For monk parakeets, this line will drop the values that represent comparisons between ranges (if present). For all datasets, this line will drop comparisons that do not represent the specified social scale for the given species
    # For instance, for yellow-naped amazon calls at the individual scale, this line will drop comparisons of different individuals from different sites in the same dialect
    dplyr::filter(type != "drop")
  
  # For Earth Mover's Distance analyses, retain only values that represent unique pairwise comparisons
  # I tested the pmax and pmin filtering and found that the number of duplicates and non-duplicates were the same, as expected for a symmetric matrix. After this filtering, the number of pairwise comparisons obtained by filtering on the target ID or other ID column will not be equivalent, and this is the expected outcome for using these values to calculate Earth Mover's Distance in downstream analysis
  if(analysis_type == "EMD"){
    
    sim_vals2 <- sim_vals %>% 
      dplyr::filter(!duplicated(paste0(pmax(target_call, other_call), pmin(target_call, other_call))))
    
    # Checking: After this step, the number of similarity values should no longer be equivalent after filtering on the target or other ID column, since only unique comparisons remain. Looks good
    # sim_vals2 %>%
    #   dplyr::filter(target_bird_ID == "NAT-ZW8") %>%
    #   pull(indiv_comparison_type) %>%
    #   table()

    # sim_vals2 %>%
    #   dplyr::filter(other_bird_ID == "NAT-ZW8") %>%
    #   pull(indiv_comparison_type) %>%
    #   table()
    
    
  # For similarity value extraction for the comparative analyses with bootstrapping, do not drop duplicate comparisons. The code written for the comparative analyses filters similarity values based on unique IDs in a given target group column (e.g. target_bird_ID), and this is a different way of filtering out duplicate comparisons
  } else if(analysis_type == "bootstrap"){
    
    # Checking: The given data frame of similarity values should have all of the non-unique pairwise comparisons for the given target group. Looks good
    # sim_vals %>%
    #   dplyr::filter(target_bird_ID == "NAT-ZW8") %>%
    #   pull(indiv_comparison_type) %>%
    #   table()

    # sim_vals %>%
    #   dplyr::filter(other_bird_ID == "NAT-ZW8") %>%
    #   pull(indiv_comparison_type) %>%
    #   table()
    
    sim_vals2 <- sim_vals
    
  }
  
  if(grepl("parakeet", species)){
    
    simVals_df <- sim_vals2 %>% 
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
    
    simVals_df <- sim_vals2 %>% 
      dplyr::mutate(
        dataset = dataset,
        similarity_method = similarity_method,
        city_year = city_year
      )
    
  }
  
  return(simVals_df)
  
}
