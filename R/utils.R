# Helper functions

# function to clean up col names
cleanup_colname <- function(raw, minlen=16){
  raw <- ifelse(is.na(raw), " ", raw)
  proc <- stringr::str_trim(raw, side = "left")
  proc <- gsub("%", "pct", proc)
  proc <- gsub("Schwarz-Weiss", "sw", proc)
  proc <- gsub("[^[:alnum:] ]", "", proc)
  proc <- stringr::str_replace_all(proc, c('ü' = 'ue', 'ï' = 'ie', 'ë' = 'ee', 'ä' = 'ae','ö' = 'oe'))
  proc <- suppressWarnings(
    do.call(rbind, lapply(proc, function(x) abbreviate(x, minlength = minlen)))[,1])
  proc <- gsub(" ", "", proc)
  proc <- tolower(proc)
  return(list(raw, proc))
}

# function to remove zero
RemoveZeroVar <- function(xdf){
  dfvar <- apply(xdf, 2, function(x) var(x, na.rm=TRUE)) 
  dfvar_idx <- ifelse(is.na(dfvar) | dfvar == 0, FALSE, TRUE)
  
  return(dfvar_idx)
}

# function to get NA columns
GetNAcols <- function(xdf){
  dfvar <- apply(xdf, 2, function(x) sum(is.na(x))) 
  dfvar_idx <- ifelse(dfvar > 0, TRUE, FALSE)
  
  return(dfvar_idx)
}

# function to replace NA with its mean
FunReplaceNAmean <- function(x){ 
  tidyr::replace_na(x, mean(x, na.rm = TRUE)) 
  }


# Extract features: column aggregation
# reference: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
# 1) taking into account NA, tiergruppe, and number of symptoms found

# Please make sure you consider missing values adequately (e.g. 3x1 and 7x0 (farm 1) should be 
# rated better than 3x1, 2x0 and 5xmissing value (farm 2, where 5 symptoms were not rated during 
# the farm visit)). The maximum number of symptoms to rate is identical for breeding farms (n=14) and
# for fattening farms (n=12), but depends on the animal category (see upper screenshot for breeding and 
# lower screenshot for fattening farms). The grey fields are not rated (since they are not relevant 
# in the respective animal category). In the white fields you can have missing data, either if 
# an entire animal category was not present during the farm visit or if the reason for the visit 
# was very specific (e.g. only concentrated on diarrhea, so that the other symptoms were not considered). 
# That’s why I would suggest to work with means (e.g. sum of symptom scores in “Saugferkel” 
# and divide it by 12).

count_symptoms2 <- function(data, piggrp = '_eber', ptrn = symt_ptrn){
  # column naming
  rated <- paste0('smrated', piggrp) # max symptoms - NA (symptoms not rated)
  found <- paste0('smfound', piggrp) # nr. of symptoms found
  score <- paste0('smscore', piggrp) # nr. of  symptoms rated - nr. of symptoms found 
  percentage <- paste0('smpct', piggrp) # score / max symptoms * 100
  
  # total symptoms as negative sign
  symt_cnt <- data %>%
    select(ends_with(piggrp), betriebsid, ttlpig_info) %>%
    # Record how many symptoms available for this particular group
    mutate( smmax = dim(.)[2] - 2 ) %>% # take ncol - [non-symp-related-col]
    rowwise() %>%
    # symptoms found as a relative to total pig in the farm
    mutate( sm_found = ifelse( ttlpig_info == 0, sum(across(matches(ptrn)), na.rm = TRUE ), 
              sum(across(matches(ptrn)) / ttlpig_info, na.rm = TRUE )) ) %>%
    # 1) how many symptoms are actually checked? = max - sym_not_checked
    mutate( sm_rated = smmax - sum(is.na(across(matches(ptrn)))) ) %>%
    ungroup %>%
    # 2) of all checked symptoms, how many symptoms detected?
    mutate( score =  sm_rated - sm_found ) %>%
    # 3) score as percentage
    mutate( sm_pct = ifelse(score == 0, 0, score/smmax * 100) ) %>%
    # 4) dynamic col names
    rename( !! rated := sm_rated, 
            !! found := sm_found, 
            !! score := score, 
            !! percentage := sm_pct ) %>%
    select( betriebsid, !! rated, !! found, !! score, !! percentage)
  
  return(symt_cnt)
}

# Function : check empty data
ShowNARows <- function(db=db) {
  # get table dimension
  dim.db <- getDim(db)
  
  # calculate percent empty rows
  db.narows <- db %>% 
    summarise(across(everything(), ~sum(is.na(.x), na.rm = TRUE))) %>%
    collect() %>%
    gather("varname", "na_rows") %>% 
    arrange(desc(na_rows)) %>%
    mutate(pct = na_rows/dim.db[1])
  
  return(db.narows)
}

# Function to extract y coefficient from multiblock model
# Extract coefficient from bootstrapped multiblock object

GetMBCoef <- function(mbbootstrap, yname = "lagedesbetriebes"){
  coef_y <- as.data.frame(cbind(
    mbbootstrap$XYcoef[[yname]]$obs,
    mbbootstrap$XYcoe[[yname]]$stats
    ))
  names(coef_y) <- paste(yname, names(coef_y), sep = "_")
  return(coef_y)
}

# Function to extract block importance

GetBlockImp <- function(mbbootstrap){
  bloimp <- data.frame(
    varname = names(mbbootstrap$bipc$obs),
    bloimp = mbbootstrap$bipc$obs,
    cilow = mbbootstrap$bipc$stats[,1],
    cihigh = mbbootstrap$bipc$stats[,2]
  )
  return(bloimp)
}

# Function to extract variable importance

GetVarImp <- function(mbbootstrap){
  varimp <- data.frame(
    varname = names(mbbootstrap$vipc$obs),
    varimp = mbbootstrap$vipc$obs,
    cilow = mbbootstrap$vipc$stats[,1],
    cihigh = mbbootstrap$vipc$stats[,2]
  )
  return(varimp)
}




# Function : calculate variance
CalcVariance <- function(db=db, dropcols=dropcols){
  db.var <- db %>% 
    select(-one_of(dropcols)) %>% 
    summarise_if(is.numeric,  ~ var(., na.rm = TRUE)) %>%
    collect %>%
    gather() %>%
    arrange(desc(value))
  
  return(db.var)
  dbDisconnect()
}


# Function : calculate pca
CalcPCA <- function(db=db, dropcols=dropcols){
  # remove na rows
  temp.db <- db %>% 
    mutate(rowidx = dplyr::row_number()) %>%
    select(-one_of(dropcols)) %>% 
    collect %>% 
    drop_na # remove na rows
  
  pca.fit <- temp.db %>%
    select(-rowidx) %>%
    select_if(is.numeric) %>% # select only numeric cols
    select_if(colSums(.) != 0) %>% # remove cols that contains only 0
    scale() %>% 
    as_tibble() %>% 
    prcomp
  
  # get row index of cleaned data
  rowidx <- temp.db %>% select(rowidx) %>% collect %>% pull
  
  return( list(pca.fit, rowidx) )
  DBI::dbDisconnect()
}

# Function: draw biplot
PCbiplot <- function(PC, x="PC1", y="PC2", xxmin=-10, xxmax=15, alphapoint = .007) {
  # PC being a prcomp object
  data <- data.frame(obsnames=seq(1:nrow(PC$x)), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_point(size = .1, alpha = alphapoint) +
    scale_x_continuous(limits = c(xxmin, xxmax))
  plot <- plot + geom_hline(aes(yintercept = 0), size=.2) + geom_vline(aes(xintercept = 0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = 1.1 * mult * (get(x)),
                      v2 = 1.1 * mult * (get(y))
  )
  # arrow
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="lightblue") 
  # text
  plot <- plot + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3, vjust=-1.5, color="blue", check_overlap = TRUE) +
    theme_minimal() + coord_equal()
  plot
}

# Function : PC Screeplot
PCIscreeplot <- function(PCobj, pcmax=7){
  var_explained_df <- data.frame(PC= seq.int(1:pcmax), # paste0("PC",1:pcmax),
                                 var_explained=(PCobj$sdev)^2/sum((PCobj$sdev)^2))
  
  splot <- ggplot(var_explained_df, aes(x=as.factor(PC),y=var_explained, group=1)) +
    geom_point(size=2) +
    scale_y_continuous(limits = c(0,1)) +
    geom_line() +
    labs(title="Scree plot: PCA on scaled data") +
    xlab("PC") +
    ylab("Variation explained") +
    #theme(panel.background = element_rect(fill = NA),
    #      panel.grid.minor = element_line(colour = "grey"),
    #      axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) #+
    theme_minimal()
  
  splot
}

# function to plot sample points based on PC
PCscatterplot <- function(ispca=pcafit, isdb=mergedtbl.db, marker=c("qtl_gesamtergebnis"), x="PC1", y="PC2"){
  # get data
  df <- data.frame(obsnames=seq(1:nrow(ispca[[1]]$x)), ispca[[1]]$x)
  
  # get the rank
  mkr <- isdb %>% select(one_of(marker)) %>% collect %>% filter(row(.) %in% ispca[[2]]) %>% pull
  
  # add rank to data
  df$mkr <- mkr
  names(df$mkr) <- marker[1]
  
  # plot
  plot <- ggplot(df, aes_string(x=x, y=y, col="mkr")) + 
    geom_point(size = .1) + 
    theme_minimal() 
  plot
}


