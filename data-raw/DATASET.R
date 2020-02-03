# Assembles the melt.cells data structure from XML
addOI <- function(cl.dat, melt.cells){
  df.data <- lapply(cl.dat, function(i){
    oi <- i$'derived-from'$'cv-term'$'.attrs'['accession']
    if(is.null(oi)) NA else as.character(oi)
  })
  melt.df <- data.frame("CVCL"=names(cl.dat),
                        "OI"=unlist(df.data))
  melt.cells <- merge(melt.cells, melt.df, by="CVCL", all.x=TRUE)
  melt.cells$OI <- as.character(melt.cells$OI)
  return(melt.cells)
}

add99Problems <- function(cl.dat, melt.cells){
  bitch <- 0
  pcl.data <- sapply(cl.dat, function(i){
    pcl <- grepl("Problematic cell line", sapply(i$'comment-list', function(j) j$'.attrs'))
    msi <- grepl("Microsatellite instability", sapply(i$'comment-list', function(j) j$'.attrs'))
    return(data.frame(MSI=any(msi), PCL=any(pcl)))
  })
  pcl.data <- as.data.frame(t(pcl.data))
  pcl.data$CVCL <- names(cl.dat)

  pcl.idx <- which(as.logical(pcl.data$PCL))
  pcl.data$comment <- NA
  pcl.data$comment[pcl.idx] <- sapply(cl.dat[pcl.idx], function(i){
    idx <- which(sapply(i$'comment-list', function(comment) comment$'.attrs' == 'Problematic cell line'))
    i$'comment-list'[[idx]]$'text'
  })

  melt.cells <- merge(melt.cells, pcl.data, by="CVCL", all.x=TRUE)
  melt.cells$MSI <- as.logical(melt.cells$MSI)
  melt.cells$PCL <- as.logical(melt.cells$PCL)
  melt.cells$comment <- as.character(melt.cells$comment)
  return(melt.cells)
}

removeFrogs <- function(cl.dat, melt.cells){
  species <- sapply(cl.dat, function(i) i$'species-list'$'cv-term'[['text']])
  species.df <- data.frame("CVCL"=names(cl.dat),
                           "species"=unlist(species))
  species.cells <- merge(melt.cells, species.df, by="CVCL", all.x=TRUE)
  melt.cells <- melt.cells[-which(!species.cells$species=='Homo sapiens'),]
  return(melt.cells)
}

createMeltCells <- function(cpath=NULL){
  require(reshape)
  require(XML)

  #cpath <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cellosaurus/cellosaurus.xml'
  data <- xmlParse(cpath)
  cellosaurus <- xmlToList(data)

  # Split cellosaurus into unique cellosaurus ids
  cl.dat <- cellosaurus$`cell-line-list`
  names(cl.dat) <- sapply(cl.dat, function(i) i$`accession-list`$accession$text)

  # Melt into data-frame of CEllosaurus IDs with all synonymous IDs
  melt.cells <- lapply(cl.dat, function(i) {
    sapply(i$`name-list`, function(j) j$text)
  })
  melt.cells <- reshape::melt(melt.cells)
  colnames(melt.cells) <- c("UID", "CVCL")
  melt.cells$ID <- gsub(" \\[.*", "", melt.cells$UID)

  melt.cells <- addOI(cl.dat, melt.cells)
  melt.cells <- add99Problems(cl.dat, melt.cells)
  melt.cells <- removeFrogs(cl.dat, melt.cells)

  # saveRDS(melt.cells, file="cellosaurus.RDS")
  # save(cl.dat, file="cellosaurus_raw.rda")
  return(melt.cells)
}

melt.cells <- createMeltCells('/mnt/work1/users/pughlab/projects/cancer_cell_lines/cellosaurus/cellosaurus.xml')
usethis::use_data(melt.cells, overwrite = T)
