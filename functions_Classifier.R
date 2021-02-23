psmFDR <- function(protPep, pepScore="hyperscore", decoy_id="DECOY", concat_decoy=0, prot_id="Label")
{
	protPep_decoy <- protPep[order(as.numeric(paste(protPep[,pepScore])), decreasing=TRUE),]
	decoyVec <- vector(mode = "numeric", nrow(protPep_decoy))
	decoyVec[grep(decoy_id, protPep_decoy[,prot_id])] <- 1
	if(concat_decoy==1)
	{
		protPep_decoy$psmFDR <- (2*cumsum(decoyVec))/(1:nrow(protPep_decoy))
	}else
	{
		protPep_decoy$psmFDR <- cumsum(decoyVec)/((1:nrow(protPep_decoy))-(cumsum(decoyVec)))
	}
	# todo lo mayor a 1, lo ponemos a 1, ya que no debería poder darse
	protPep_decoy[protPep_decoy$psmFDR > 1, "psmFDR"] <- 1
	return(protPep_decoy)

}

parseMascotDATFileSimple_G<-function(fileName) {

  cat("Reading file ", fileName, " ... \n");

  conn=file(fileName,open="r")
  linn=readLines(conn)

  indQuery <- grep("^q[0-9]", linn)

  queries <- linn[indQuery]

  cat("Processing queries ... \n");

  queries.lst <- strsplit(queries, "=")

  queries.lst_1 <- lapply(queries.lst, FUN = function(x) x[1])
  queries.lst_2 <- lapply(queries.lst, FUN = function(x) x[2])

  queries.lst_1.df <- strsplit(unlist(queries.lst_1), "_")

  queries.lst_1.df2 <- data.frame("query" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[1])), "rank" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[2])), "ob1" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[3])), "ob2" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[4])))

  queries.lst_2_str <- unlist(queries.lst_2)
  queries.lst_2_strlst <- strsplit(queries.lst_2_str, ";")

  queries.lst_2_df <- data.frame("assignPep" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[1])),  "assignProt" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[2])))

  queries.df <- data.frame(queries.lst_1.df2, queries.lst_2_df)
  queries.assign.df <- queries.df[!is.na(queries.df$assignProt),]


  queries.assign.df.p1 <- as.data.frame(strsplit(paste(queries.assign.df$assignPep), ","))
  queries.assign.df.p1 <- t(queries.assign.df.p1)
  queries.assign.df.p1 <- data.frame(queries.assign.df[, c(1,2)], queries.assign.df.p1)
  queries.assign.df.p2 <- strsplit(paste(queries.assign.df$assignProt), ",")

  queries.assign.prot.df <- unlist(strsplit(unlist(queries.assign.df.p2), ";"))
  pepProt <- unlist(lapply(queries.assign.df.p2, length))
  queries.assign.prot.df2 <- strsplit(queries.assign.prot.df, ":")
  queries.assign.prot.df3 <- data.frame(unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[1])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[2])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[3])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[4])), unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[5])))
  colnames(queries.assign.prot.df3) <- NULL
  rownames(queries.assign.prot.df3) <- NULL

  indexPep <- list()

  for (i in 1:length(pepProt)) {
	  indexPep[i] <- list(rep(i, pepProt[i]))
  }

  indexPep.vec <- unlist(indexPep)

  queriesPepProt <- data.frame(queries.assign.df.p1[indexPep.vec,], queries.assign.prot.df3)

  colnames(queriesPepProt) <- c("Query", "Rank", "MissedCleavages", "PeptideMr", "Delta", "NumberOfIonsMatched", "PeptideSeq", "PeaksFromIons1", "VarMod", "ion_score", "IonSeries", "PeaksFromIons2", "PeaksFromIons3", "ProteinAccession", "FrameNumber", "Start", "End", "Multiplicity")
  rownames(queriesPepProt) <- NULL
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])

  # Merged to add the Enzime cuts
  enzimeCuts_df <- setNames(queries.df[,grep(paste(c("query", "rank", "ob1", "assignPep"), collapse="|"), colnames(queries.df))], c("Query", "Rank", "ob1", "enzCut"))
  enzimeCuts_df <- enzimeCuts_df[c(which(enzimeCuts_df$ob1 == "terms")),]
  enzimeCuts_df <- enzimeCuts_df[,grep(paste(c("Query", "Rank", "enzCut"), collapse="|"), colnames(enzimeCuts_df))]

  queriesPepProt <- merge(queriesPepProt, enzimeCuts_df, by=c("Query", "Rank"), all=T)

  return(queriesPepProt)

}

parseMascotStandardOutput <- function(fileName) {

  cat("Reading file ", fileName, " ... \n");
  conn=file(fileName,open="r")
  linn=readLines(conn)
  indQuery <- grep("^q[0-9]", linn)
  queries <- linn[indQuery]
  queryIndexes <- grep("Content-Type: application/x-Mascot; name=\"query", linn)
  allIndexIndexes <- grep("index=", linn, fixed = TRUE)
  allTitleIndexes <- grep("title=", linn, fixed = TRUE)
  allIndexIndexes <- allIndexIndexes[allIndexIndexes > queryIndexes[1]]
  allTitleIndexes <- allTitleIndexes[allTitleIndexes > queryIndexes[1]]
  allIndexIndexes <- linn[allIndexIndexes]
  allTitleIndexes <- linn[allTitleIndexes]
  indexMGF <- as.numeric(paste(sapply(allIndexIndexes, FUN = function(x) unlist(strsplit(x, "=")[[1]][2]))))
  titlesDat <- paste(sapply(allTitleIndexes, FUN = function(x) unlist(strsplit(x, "=")[[1]][2])))
  querynames <- linn[queryIndexes]

  querynames <- strsplit(querynames, " ")
  querynames <- unlist(sapply(1:length(querynames), FUN = function(x) querynames[x][[1]][3]))
  querynames <- strsplit(querynames, "=")
  querynames <- unlist(lapply(querynames, FUN = function(x) x[2]))
  querynames <- substr(querynames, 2, nchar(querynames) - 1)
  querynames <- gsub("query", "q", querynames)

  titlesDat <- gsub("%2e", ".", titlesDat)
  titlesDat <- gsub("%3a", ":", titlesDat)
  titlesDat <- gsub("%20", " ", titlesDat)
  titlesDat <- gsub("%22", "\"", titlesDat)

   # this index is 0-based, and R is not. so we translate it.
  indexMGF <- as.numeric(indexMGF) + 1

  queryTitleDF <- data.frame("queryname" =  querynames, "ScanString" = titlesDat, "mgfIndex" =  indexMGF)

  ## ver como analizar y añadir al objeto resultado

  cat("Processing queries ... \n");

  if (length(grep("terms", queries)) != 0) {
  	queries <- queries[-(grep("terms", queries))]
  }
  if (length(grep("primary_nl", queries)) != 0) {
  	queries <- queries[-(grep("primary_nl", queries))]
  }
  if (length(grep("et_mods", queries)) != 0) {
  	queries <- queries[-(grep("et_mods", queries))]
  }

  queries.lst <- strsplit(queries, "=")

  queries.lst_1 <- lapply(queries.lst, FUN = function(x) x[1])
  queries.lst_2 <- lapply(queries.lst, FUN = function(x) x[2])

  queries.lst_1.df <- strsplit(unlist(queries.lst_1), "_")

  queries.lst_1.df2 <- data.frame("query" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[1])), "rank" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[2])), "ob1" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[3])), "ob2" = unlist(lapply(queries.lst_1.df, FUN = function(x) x[4])))

  queries.lst_2_str <- unlist(queries.lst_2)
  queries.lst_2_strlst <- strsplit(queries.lst_2_str, ";")

  queries.lst_2_df <- data.frame("assignPep" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[1])),  "assignProt" = unlist(lapply(queries.lst_2_strlst, FUN = function(x) x[2])))

  queries.df <- data.frame(queries.lst_1.df2, queries.lst_2_df)
  queries.assign.df <- queries.df[!is.na(queries.df$assignProt),]
  queries.assign.df.p1 <- as.data.frame(strsplit(paste(queries.assign.df$assignPep), ","))
  queries.assign.df.p1 <- t(queries.assign.df.p1)
  queries.assign.df.p1 <- data.frame(queries.assign.df[, c(1,2)], queries.assign.df.p1)
  queries.assign.df.p2 <- strsplit(paste(queries.assign.df$assignProt), ",")
  queries.assign.prot.df <- unlist(strsplit(unlist(queries.assign.df.p2), ";"))
  pepProt <- unlist(lapply(queries.assign.df.p2, length))
  queries.assign.prot.df2 <- strsplit(queries.assign.prot.df, ":")
  queries.assign.prot.df3 <- data.frame(unlist(lapply(queries.assign.prot.df2, FUN = function(x) x[1])))
  colnames(queries.assign.prot.df3) <- NULL
  rownames(queries.assign.prot.df3) <- NULL

  indexPep <- list()

  for (i in 1:length(pepProt)) {
	  indexPep[i] <- list(rep(i, pepProt[i]))
  }

  indexPep.vec <- unlist(indexPep)

  queriesPepProt <- data.frame(queries.assign.df.p1[indexPep.vec,], queries.assign.prot.df3)

  colnames(queriesPepProt) <- c("Query", "Rank", "MissedCleavages", "PeptideMass", "PeptideDelta", "NumberOfIonsMatched", "PeptideSeq", "PeaksFromIons1", "Modifications", "score", "IonSeries", "PeaksFromIons2", "PeaksFromIons3", "ProteinAccession")
  queriesPepProt <- queriesPepProt[, c("Query", "Rank", "MissedCleavages", "PeptideMass", "PeptideDelta", "PeptideSeq", "Modifications", "score", "ProteinAccession")]
  rownames(queriesPepProt) <- NULL
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])
  queriesPepProt[, "ProteinAccession"] <- sub("\\\"", "", queriesPepProt[, "ProteinAccession"])
  queriesPepProt <- merge(queriesPepProt, queryTitleDF, by.x = "Query", by.y = "queryname")
  queriesPepProt$ScanNumber <- NA
  queriesPepProt$charge <- NA
  return(queriesPepProt)

}
### scan is the index of the spectra in the mgf
parseCometStandardOutput <- function(fileName) {
# xcorr the higher the better
  cat("Reading file ", fileName, " ... \n");
  cometSearchResultFile <- read.csv2(fileName, header = TRUE, sep = "\t", fill = TRUE, skip = 1)
  if ((colnames(cometSearchResultFile)[1] != "plain_peptide") & (colnames(cometSearchResultFile)[1] != "scan")) {
 	 cometSearchResultFile <- read.csv2(fileName, header = TRUE, sep = "\t", fill = TRUE)
  }
  # if(length(colnames(cometSearchResultFile)) == 18) {
  # 	cometSearchResultFile$duplicate_protein_count <- cometSearchResultFile$protein
  # 	cometSearchResultFile$protein <- cometSearchResultFile$next_aa
  # }
  # they only report one hit per query, so we set the rownumber as
  queries <- paste("q_", rownames(cometSearchResultFile), sep = '')
  # allowed missed cleavages are 0 or 1
  cometSearchResult <- data.frame(Query = queries, Rank = NA, MissedCleavages = "0/1", PeptideMass = as.numeric(paste(cometSearchResultFile$exp_neutral_mass)), PeptideDelta = as.numeric(paste(cometSearchResultFile$delta_cn)), PeptideSeq = paste(cometSearchResultFile$plain_peptide), Modifications = paste(cometSearchResultFile$modifications), score = as.numeric(paste(cometSearchResultFile$xcorr)), ProteinAccession = paste(cometSearchResultFile$protein), ScanString = paste(fileName, ":scan", cometSearchResultFile$scan, sep = ""), mgfIndex = NA, ScanNumber = paste(cometSearchResultFile$scan), charge = as.numeric(paste(cometSearchResultFile$charge)))
  return(cometSearchResult)

}

parseXTandemStandardOutput <- function(fileName) {

	cat("Reading file ", fileName, " ... \n");
	tandemSearchResultFile <- read.csv2(fileName, header = TRUE, sep = "\t", fill = TRUE)
	spectrum <- as.numeric(paste(unlist(lapply(strsplit(paste(tandemSearchResultFile$spectrum), "@"), FUN = function(x) x[1]))))
	# x.tandem.expect the lower the better.
	# convertimos para que sea the higher the better igual que comet y mascot
	score = -log(as.numeric(paste(tandemSearchResultFile$X.Tandem.expect)))
	# in order to avoid, Inf values. Inf arises when expect is 0.
	score[score == Inf] <- 1000
	queries <- paste("q_", spectrum, sep = '')
	tandemSearchResult <- data.frame(Query = queries, Rank = tandemSearchResultFile$rank, MissedCleavages = "0/1", PeptideMass = as.numeric(paste(tandemSearchResultFile$expMass)), PeptideDelta = NA, PeptideSeq = paste(tandemSearchResultFile$peptideSequence), Modifications = paste(tandemSearchResultFile$peptideMods), score = score, ProteinAccession = paste(tandemSearchResultFile$protein), ScanString = paste(tandemSearchResultFile$psm), mgfIndex = spectrum, ScanNumber = NA, charge = as.numeric(paste(tandemSearchResultFile$charge)))
	return(tandemSearchResult)

}

parseXTandemStandardOutput_v_NEW <- function(fileName) {

	cat("Reading file ", fileName, " ... \n");
	tandemSearchResultFile <- read.csv2(fileName, header = TRUE, sep = "\t", fill = TRUE)
	# x.tandem.expect the lower the better.
	# convertimos para que sea the higher the better igual que comet y mascot
	score = -log(as.numeric(paste(tandemSearchResultFile$expect)))
	# in order to avoid, Inf values. Inf arises when expect is 0.
	score[score == Inf] <- 1000
	queries <- paste("q_", tandemSearchResultFile$index_mgf, sep = '')

    scan_tmp<-lapply(strsplit(paste(tandemSearchResultFile$scan)," "),"[",5)
    scan_tmp2<-str_extract(scan_tmp, "[0-9]+")
	tandemSearchResult <- data.frame("Query" = queries, "Rank" = NA, "MissedCleavages" = tandemSearchResultFile$missed_cleavages, "PeptideMass" = as.numeric(paste(tandemSearchResultFile$PrecursorMZ)), "PeptideDelta" = NA, "PeptideSeq" = paste(tandemSearchResultFile$PeptideSeq), "Modifications" = paste(tandemSearchResultFile$modification), "score" = score, "ProteinAccession" = paste(tandemSearchResultFile$Protein), "ScanString" = paste(tandemSearchResultFile$scan),
    "mgfIndex" = NA, "ScanNumber" = paste(scan_tmp2), "charge" = as.numeric(paste(tandemSearchResultFile$charge)))
	return(tandemSearchResult)

}
