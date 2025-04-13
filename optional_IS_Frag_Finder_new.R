library(tidyverse)

find_ppm <- function(mo, me){
  ppm <- (mo - me)/me * 10^6
  return(ppm)
}

read_mgf <- function(path){
  #Read in the mgf file
  l <- readLines(path)
  #Find the start of each feature
  ti <- which(str_detect(l, 'FEATURE_ID'))
  #Find the end of each feature
  tt <- which(str_detect(l, 'END IONS'))
  ck <- data.frame()
  #Pull out the necessary information
  for(i in 1:length(ti)){
    tmp <- tibble(Compound_ID = str_split(l[ti[i]], '=')[[1]][2],
                  ions = list(as.numeric(sapply(str_split(l[(ti[i]+8):(tt[i]-1)], ' '), function(x) x[1]))))
    ck <- rbind(ck, tmp)
  }
  # #Check to see if "FT_" is at the start of feature names
  # if(str_detect(ck$Compound_ID[1], 'FT_', negate = TRUE)){
  #   #If not, add it on
  # frags <- ck %>%
  #   modify_at('Compound_ID', ~paste0('FT_', .x))
  # return(frags)
  # }else{
  #   return(frags)
  # }
  return(ck)
}

setwd("~/PERMA")
tm <- read_mgf('PERMA_RPPOS.mgf')

duplicate_ids <- tm %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)


link_data <- link_josh
time_tol = 0.005
ppm_tol = 5
cor_cutoff = 0.9
mgf_file <- '/Users/benyang/Downloads/Josh/josh_HILIC_Neg.mgf'
tm <- read_mgf(mgf_file)


find_ISFrags <- function(link_data, mgf_file,
                         time_tol = 0.005, ppm_tol = 5, cor_cutoff = 0.9){
  
  link_data$Compound_ID <- paste0("FT_", link_data$RT..min., "_", link_data$m.z)
  
  intensity_data <- link_data %>%
    select(Compound_ID, starts_with("Group.Area"))
  
  #Check data integrity:
  if(sum(colnames(link_data) %in% c('Compound_ID', 'm.z', 'RT..min.')) != 3){
    stop('Check Link Data, it must contain the following columns named exactly: Compound_ID, RT [min], m/z')
  }
  
  #Format the fragment table to be correct
  fragtab <- link_data %>%
    dplyr::select(Compound_ID, m.z, RT..min.) %>%
    dplyr::rename('mz' =  m.z, 'rt' = RT..min.) 
  
  #Calculate the difference between each rt and all others:
  #Output the results into a list
  olist <- list()
  for(i in 1:nrow(fragtab)){
    ovec <- fragtab$rt - fragtab$rt[i]
    names(ovec) <- fragtab$Compound_ID
    olist[[i]] <- ovec
  }
  
  #Now trim each vector to only places where the rt dif <= 1 and pull the names
  olist_clean <- lapply(olist, function(x) names(x[abs(x) <= time_tol]))
  
  #Remove features that are a length of 1 - These have no matches
  olist_trimmed <- olist_clean[sapply(olist_clean, length) != 1]
  #And each set will show up twice so trim it down
  olist_final <- olist_trimmed[!duplicated(olist_trimmed)]
  
  #Now turn it into a dataframe
  #Making a unique group for each set
  dfo <- tibble()
  for(i in 1:length(olist_final)){
    df <- tibble(group = paste0('G', i),
                 matches = list(olist_final[[i]]))
    dfo <- rbind(dfo, df)
  }
  
  #Now make a table with the groups
  group_tab <- dfo %>%
    unnest('matches') %>%
    dplyr::rename('Compound_ID' = matches)
  
  #Make a table of intensities
  int_table <- intensity_data 
  
  #Make a table with a correlation matrix of each group
  big_group <- group_tab %>%
    #Join on teh intensity data
    left_join(int_table) %>%
    #Make long
    pivot_longer(cols = c(-group, -Compound_ID), names_to = 'sample') %>%
    #Group
    group_by(group) %>%
    nest() %>%
    #Make a correlation matrix for each set of candidates
    mutate(cor_tab = purrr::map(data, function(x){
      ct <- x %>%
        pivot_wider(names_from = Compound_ID) %>%
        column_to_rownames('sample') %>%
        dplyr::select_if(is.numeric) %>%
        cor()
      return(ct)
    })) 
  
  #Filter the candidates based on correlation info:
  IS_cand <- big_group %>%
    mutate(cands = purrr::map_chr(cor_tab, function(x){
      x %>%
        as.data.frame() %>%
        rownames_to_column('var1') %>%
        #Make the correlation table long for ease
        pivot_longer(starts_with('FT_'), names_to = 'var2') %>%
        #Remove the diag
        filter(var1 != var2) %>%
        #Grab only ones above the cutoff
        filter(value >= cor_cutoff) %>%
        #Pull out var1 since it will include both sides
        pull(var1) %>%
        #Remove duplicates
        unique() %>%
        paste0(collapse = '; ')
    })) %>%
    ungroup() %>%
    #Remove candidates where nothing passed
    filter(cands != '') %>%
    #Clear out the duplicates
    filter(!duplicated(cands))
  
  #Now using the grouping information prep for MS2 data
  grp <- IS_cand %>%
    #Clean up the fragment data
    mutate(split = str_split(cands, pattern = '; ')) %>%
    dplyr::select(group, split) %>%
    unnest('split') %>%
    #Rename for merging
    dplyr::rename('Compound_ID' = split) %>%
    #Merge in the fragment data
    left_join(fragtab) %>%
    group_by(group) %>%
    #Keep the highest mz as the parent
    mutate(type = ifelse(max(mz) == mz, 'parent', 'frag'))
  
  frags <- read_mgf(mgf_file)
  
  #Merge in the MS2 data
  grp_frag <- grp %>%
    left_join(frags) %>%
    #Grab the necessary information from the fragment tab
    reframe(parent_mz = mz[type == 'parent'],
            parent_ft = Compound_ID[type == 'parent'],
            frag_mz = list(mz[type == 'frag']),
            frag_fts = list(Compound_ID[type == 'frag']),
            ions_parent = ions[type == 'parent']) %>%
    #Calculate PPMs of each fragment candidate against the MS2 data of parent fragment
    mutate(ppm_vals = map2(frag_mz, ions_parent, function(x,y){
      purrr::map(x, function(a){
        find_ppm(a,y)
      })
    })) %>%
    mutate(match = purrr::map2(ppm_vals, frag_fts, function(x,y){
      lgls <- sapply(x, function(x) any(abs(x) <= ppm_tol))
      data.frame(frag_cand = y,
                 match = lgls)
    })) %>%
    dplyr::select(parent_ft, match) %>%
    unnest('match') %>%
    filter(match)
  
  final_output <- grp_frag %>%
    dplyr::select(-match) %>%
    group_by(parent_ft) %>%
    summarise(fragments = paste0(unique(frag_cand), collapse = '; '))
  
  return(final_output)
}

#link_data contains columns with m/z, RT [min.], and Compound_ID
#Intensity data in a numeric matrix of RAW (or adjusted for things like sample weight) of feature intensities
#Each column is a feature and each row is a sample with sample names as rownames ONLY
#mgf_file is the path of the mgf file you gave to SIRIUS
#time_tol is the tolerance for finding initial candidates - Defaults to 0.005
#cor_cutoff is the cutoff for correlations between two features to be considered a fragment - Default is 0.9
#ppm_tol is the ppm tolerance between candidate fragments (which pass the correlation test) and 
#the parent (heaviest) ion MS2 fragments - Default is 5

#Example use:

link_rp <- read_csv('PERMA_RPPOS.csv')

duplicate_ids <- link_rp %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)

intensity_RP <- read.csv('RPPOS_intensity.csv') 

duplicate_ids <- intensity_RP %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)

intensity_new <- intensity_RP %>%
  dplyr::select(Compound_ID, contains('.raw'), -contains('QC'), -contains('Blank'), -contains('Peak.Rating')) %>%
  column_to_rownames('Compound_ID') %>%
  rename_with(~(gsub('Area..', '', .x))) %>%
  t()


testout2 <- find_ISFrags(link_data = link_rp, intensity_data = intensity_new, 
                        mgf_file = 'PERMA_RPPOS.mgf',
                        time_tol = 0.005, cor_cutoff = 0.98)


write.csv(testout2, file = "in_source_RP.csv")

link_josh <- read.csv("/Users/benyang/Downloads/Josh/josh_RP_Pos.csv")
testout_josh <- find_ISFrags(link_data = link_josh, 
                         mgf_file = '/Users/benyang/Downloads/Josh/josh_RP_Pos.mgf',
                         time_tol = 0.005, cor_cutoff = 0.98)
testout_josh_RP <- testout_josh %>%
  mutate(
    parent_ft = str_c(parent_ft, "_RP"),
    fragments = str_replace_all(fragments, "(FT_[^;]+)", "\\1_RP")
  )

link_josh <- read.csv("/Users/benyang/Downloads/Josh/josh_HILIC_Neg.csv")
testout_josh <- find_ISFrags(link_data = link_josh, 
                             mgf_file = '/Users/benyang/Downloads/Josh/josh_HILIC_Neg.mgf',
                             time_tol = 0.005, cor_cutoff = 0.98)
testout_josh_HN <- testout_josh %>%
  mutate(
    parent_ft = str_c(parent_ft, "_HN"),
    fragments = str_replace_all(fragments, "(FT_[^;]+)", "\\1_HN")
  )

write.csv(testout_josh_RP, file = "/Users/benyang/Downloads/Josh/in_source_RP.csv")
write.csv(testout_josh_HN, file = "/Users/benyang/Downloads/Josh/in_source_HN.csv")

##hneg

link_hn <- read_csv('PERMA_HNEG.csv')

duplicate_ids <- link_hn %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)

intensity_HN <- read.csv('HNEG_intensity.csv') 

duplicate_ids <- intensity_HN %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)

intensity_new <- intensity_HN %>%
  dplyr::select(Compound_ID, contains('.raw'), -contains('QC'), -contains('Blank'), -contains('Peak.Rating')) %>%
  column_to_rownames('Compound_ID') %>%
  rename_with(~(gsub('Area..', '', .x))) %>%
  t()


testout2_hn <- find_ISFrags(link_data = link_hn, intensity_data = intensity_new, 
                         mgf_file = 'PERMA_HNEG.mgf',
                         time_tol = 0.005, cor_cutoff = 0.98)


write.csv(testout2_hn, file = "in_source_HN.csv")





make_comparison <- function(int_df, frag_df, row, small = NULL){
  int_df = intensity_new
  frag_df = testout2
  row = 2
  small = c(1,10)
  
  parent <- frag_df$parent_ft[row] 
  frags <- gsub(' ', '', unlist(str_split(frag_df$fragments[[row]], ';')))
  if(!is.null(small)){
    frags <- frags[small[1]:small[2]]
  }
  
  plot.new()
  dplyr::select(data.frame(int_df), all_of(c(parent, frags))) %>%
    psych::pairs.panels(smooth = F, rug = F, ellipses = F)
  #return(frags)
}

make_comparison(int_df = intensity_new, frag_df = testout2, row = 2, small = c(1,10)) 



# make_comparison2 <- function(int_df, frag_df, row){
#   parent <- frag_df$parent_ft[row] 
#   frags <- gsub(' ', '', unlist(str_split(frag_df$fragments[[row]], ';')))
#   return(frags)
#   dplyr::select(data.frame(int_df), all_of(c(parent, frags))) %>%
#   cor()
# }

#test <- make_comparison2(int_df = intensity_new, frag_df = testout2, row = 2)

make_comparison(int_df = intensity_new, frag_df = testout2, row = 2, small = c(1,10)) 





