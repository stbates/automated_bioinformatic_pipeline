# automated_bioinformatic_pipeline_for_conservation_purposes.R
# v1.2
# 01 Mar. 2025

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(httr)
library(jsonlite)
library(measurements)
library(pracma)
library(RCurl)
library(rvest)
library(stringdist)
library(tcltk)
library(XML)

# Select the state to survey by removing the hashtag

#theState <- "Alabama"
#theState <- "Alaska"
##theState <- "Arizona"
#theState <- "Arkansas"
#theState <- "California"
theState <- "Colorado"
#theState <- "Connecticut"
#theState <- "Delaware"
#theState <- "Florida"
#theState <- "Georgia"
#theState <- "Hawaii"
#theState <- "Idaho"
#theState <- "Illinois"
#theState <- "Indiana"
#theState <- "Iowa"
#theState <- "Kansas"
#theState <- "Kentucky"
#theState <- "Louisiana"
#theState <- "Maine"
#theState <- "Maryland"
#theState <- "Massachusetts"
#theState <- "Michigan"
#theState <- "Minnesota"
#theState <- "Mississippi"
#theState <- "Missouri"
#theState <- "Montana"
#theState <- "Nebraska"
#theState <- "Nevada"
#theState <- "New Hampshire"
#theState <- "New Jersey"
#theState <- "New Mexico"
#theState <- "New York"
#theState <- "North Carolina"
#theState <- "North Dakota"
#theState <- "Ohio"
#theState <- "Oklahoma"
#theState <- "Oregon"
#theState <- "Pennsylvania"
#theState <- "Rhode Island"
#theState <- "South Carolina"
#theState <- "South Dakota"
#theState <- "Tennessee"
#theState <- "Texas"
#theState <- "Utah"
#theState <- "Vermont"
#theState <- "Virginia"
#theState <- "Washington"
#theState <- "West Virginia"
#theState <- "Wisconsin"
#theState <- "Wyoming"

# Read Fungal Red List csv file into the R environment

pth_to_csv <- file.choose()
redListFungi <- fread(pth_to_csv, header = TRUE)

resultsTable <- data.table(occid=character(), collid=character(), dbpk=character(), catalogNumber=character(), basisOfRecord=character(), sciname=character(), scientificNameAuthorship=character(), recordedBy=character(), eventDate=character(), year=character(), month=character(), day=character(), verbatimEventDate=character(), country=character(), stateProvince=character(), county=character(), municipality=character(), locality=character(), decimalLatitude=character(), decimalLongitude=character(), verbatimCoordinates=character(), geodeticDatum=character(), coordinateUncertaintyInMeters=character(), locationRemarks=character(), georeferenceProtocol=character(), georeferenceSources=character(), georeferenceVerificationStatus=character(), georeferenceRemarks=character(), recordIDcharacter=character(), habitat=character(), substrate=character(), associatedTaxa=character(), fieldnumber=character(), fieldNotes=character(), occurrenceRemarks=character())
taxonHoldTable <- data.frame(occid=character(), collid=character(), originalName=character(), originalNameAuthorship=character(), element1=character(), element2=character(), element3=character(), element4=character(), elementsLength=character(), validGenusName=character(), validSpeciesName=character(), validInfraspecificRank=character(), validInfraspecificName=character(), genusValid=character(), speciesValid=character(), infraspecificValid=character(), currentName=character(), currentAuthors=character(), currentMbNumber=character(), type=character())
finalDataFrame <- data.frame(occid=character(), collid=character(), dbpk=character(), catalogNumber=character(), basisOfRecord=character(), originalName=character(), originalNameAuthorship=character(), validGenusName=character(), validSpeciesName=character(), validInfraspecificRank=character(), validInfraspecificName=character(), currentName=character(), currentAuthors=character(), currentMbNumber=character(), recordedBy=character(), eventDate=character(), year=character(), month=character(), day=character(), verbatimEventDate=character(), country=character(), stateProvince=character(), county=character(), municipality=character(), locality=character(), decimalLatitude=character(), decimalLongitude=character(), verbatimCoordinates=character(), geodeticDatum=character(), coordinateUncertaintyInMeters=character(), locationRemarks=character(), georeferenceProtocol=character(), georeferenceSources=character(), georeferenceVerificationStatus=character(), georeferenceRemarks=character(), recordIDcharacter=character(), habitat=character(), substrate=character(), associatedTaxa=character(), fieldnumber=character(), fieldNotes=character(), occurrenceRemarks=character())
finalCondensedDataFrame <- data.frame(occid=character(), currentName=character(), currentAuthors=character(), currentMbNumber=character(), recordedBy=character(), eventDate=character(), year=character(), month=character(), day=character(), verbatimEventDate=character(), country=character(), stateProvince=character(), county=character(), municipality=character(), locality=character(), decimalLatitude=character(), decimalLongitude=character(), verbatimCoordinates=character(), geodeticDatum=character(), coordinateUncertaintyInMeters=character())
taxaCountTable <- data.frame(currentName = character(), tally = numeric(), state = character())
redListCountTable <- data.frame(currentName = character(), tally = numeric())
finalRareTaxaTable <- data.frame(currentName=character(), currentAuthors=character(), currentMbNumber=character(), tally = numeric(), distanceKm=character(), status=character(), occid_1=character(), decimalLatitude_1=character(), decimalLongitude_1=character(), county_1=character(), occid_2=character(), decimalLatitude_2=character(), decimalLongitude_2=character(), county_2=character())

mp_api_url <- paste("https://www.mycoportal.org/portal/api/v2/occurrence/search?stateProvince=", theState, "&limit=300&offset=", sep = "")
fd_api_url <- "https://www.mycoportal.org/fdex/services/api/query.php?"
if_api_url <- "https://www.indexfungorum.org/ixfwebservice/fungus.asmx/NameSearch?AnywhereInText=FALSE"
fg_api_url <- "https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText="
nonFungi <- "This generic name in not considered to apply to an organism within the fungal clade"
formInfraRank <- c("f.", "f ", "f. sp.", "form", "forma")
subspInfraRank <- c("ssp.", "ssp ", "subsp.", "sub species")
varInfraRank <- c("v.", "var.", "var ", "variety")

mp_api_resultCountURL <- paste("https://www.mycoportal.org/portal/api/v2/occurrence/search?stateProvince=", theState, sep = "")
mp_api_resultCount <- try(fromJSON(mp_api_resultCountURL), silent = TRUE)
mp_api_resultCount <- mp_api_resultCount$count

mp_api_get <- function(v){
  # get data from MyCoPortal API
  resultLength <- 0
  mp_api_url_pg <- paste(mp_api_url, v, sep = "")
  mpPage <- try(fromJSON(mp_api_url_pg), silent = TRUE)
  mpPageLength <- length(mpPage)
  if(mpPageLength == 5){
    result <- mpPage$results
    resultLength <- dim(result)[1]
    if(is.null(resultLength) == TRUE){
      result <- NULL
    }
  }else{
    result <- NULL
  }
  return(result)
}

fd_api_get <- function(w, x){
  resultLength <- 0
  fd_api_url_pg <- paste(fd_api_url, "qField=", w, "&qText=", x, sep = "")
  fdPage <- try(fromJSON(fd_api_url_pg), silent = TRUE)
  fdPageLength <- length(fdPage)
  if(fdPageLength == 16){
    result <- subset(fdPage, fdPage$recordSource == "Index Fungorum")
    result <- subset(result, result$taxonomicStatus == "Assumed legitimate" | result$taxonomicStatus == "Legitimate" | result$taxonomicStatus == "Conserved" | result$taxonomicStatus == "Orthographic variant")
    resultLength <- dim(result)[1]
    if(resultLength == 0){
      result <- NULL
    }
  }else{
    result <- NULL
  }
  return(result)
}

if_api_get <- function(y, z){
  if_api_url_pg <- paste(if_api_url, "&MaxNumber=", y, "&SearchText=", z, sep = "")
  ifPage <- RETRY("GET", if_api_url_pg)
  ifPage <- xmlInternalTreeParse(ifPage, encoding = "utf-8")
  result <- xmlToList(ifPage)
  return(result)
}

fg_api_get <- function(z){
  resultLength <- 0
  fg_api_url_pg <- paste(fg_api_url, z, sep = "")
  fgPage <- try(fromJSON(fg_api_url_pg), silent = TRUE)
  fgPageLength <- length(fgPage)
  if(fgPageLength == 11){
    result <- fgPage
    resultLength <- length(result)
    if(resultLength == 0){
      result <- NULL
    }
  }else{
    result <- NULL
  }
  return(result)
}

correct_record_get <- function(records){
  recordsCheckHold <- subset(records, records$taxonomicAgreement == "Synchronous")
  recordsCheckHoldLength <- dim(recordsCheckHold)[1]
  if(recordsCheckHoldLength == 0){
    recordsCheckHold <- subset(records, records$mbNumber == records$currentMbNumber)
    recordsCheckHoldLength <- dim(recordsCheckHold)[1]
    if(recordsCheckHoldLength == 0){
      recordsCheckHold <- subset(records, records$currentStatus == "Stable")
      recordsCheckHold <- recordsCheckHold[1,]
      recordsCheckHoldLength <- dim(recordsCheckHold)[1]
    }
  }
  if(recordsCheckHoldLength > 1){
    recordsCheckHold <- subset(records, records$mbNumber == records$currentMbNumber)
    recordsCheckHoldLength <- dim(recordsCheckHold)[1]
    if(recordsCheckHoldLength == 0){
      recordsCheckHold <- subset(records, records$currentStatus == "Stable")
      recordsCheckHold <- recordsCheckHold[1,]
    }
    if(recordsCheckHoldLength > 1){
      recordsCheckHold <- subset(records, records$currentStatus == "Stable")
      recordsCheckHold <- recordsCheckHold[1,]
    }
    if(recordsCheckHoldLength > 1){
      recordsCheckHold <- records[1,]
    }
  }
  result <- recordsCheckHold
  return(result)
}

slime_mold_check <- function(aGenus, aGenusCheck){
  if(aGenus != "Lycoperdon"){
    genusCheck <- aGenusCheck$IndexFungorum
    if(aGenus == genusCheck$NAME_x0020_OF_x0020_FUNGUS){
      slimeMoldCheck <- genusCheck$EDITORIAL_x0020_COMMENT
      slimeMoldCheckLength <- length(slimeMoldCheck)
      if(slimeMoldCheckLength != 0){
        if(slimeMoldCheck == nonFungi){
          result <- "Slime Mold"
        }else{
          result <- "Valid Genus"
        }
      }else{
        result <- "Valid Genus"
      }
    }else{
      result <- "Not Valid Genus"
    }
  }else{
    result <- "Valid Genus"
  }
  return(result)
}

macrofungi_check <- function(aCheck){
  result <- "Macrofungus"
  lichenCheck <- grep("Lichenized", aCheck$guild)
  lichenCheckLength <- length(lichenCheck)
  if(lichenCheckLength == 1){
    result <- "Lichen"
  }
  microfungiCheck <- grep("Microfungus", aCheck$growthForm)
  microfungiCheckLength <- length(microfungiCheck)
  if(microfungiCheckLength == 1){
    result <- "Microfungus"
  }
  rustCheck <- grep("Rust", aCheck$growthForm)
  rustCheckLength <- length(rustCheck)
  if(rustCheckLength == 1){
    result <- "Microfungus"
  }
  smutCheck <- grep("Smut", aCheck$growthForm)
  smutCheckLength <- length(smutCheck)
  if(smutCheckLength == 1){
    result <- "Microfungus"
  }
  return(result)
}

closest_match <- function(string, stringVector){
  result <- stringVector[amatch(string, stringVector, maxDist=Inf)]
  return(result)
}

county_check <- function(aCounty){
  countyFound <- "No"
  countyCheck <- strsplit(aCounty, " ")
  countyCheck <- countyCheck[[1]]
  countyCheckLength <- length(countyCheck)
  if(countyCheckLength >= 2){
    theCheck <- tolower(countyCheck[countyCheckLength])
    if(theCheck == "county"){
      countyFound <- "Yes"
    }
    if(theCheck == "co."){
      countyFound <- "Yes"
    }
    if(theCheck == "co"){
      countyFound <- "Yes"
    }
  }
  if(countyFound == "Yes"){
    pasteLength <- countyCheckLength - 1
    result <- trimws(paste(countyCheck[1:pasteLength], collapse = " "))
  }else{
    result <- trimws(aCounty)
  }
  return(result)
}

numberResults <- 0
n <- 0

#gather the records from MyCoPortal API

while(numberResults < mp_api_resultCount){
  theData <- mp_api_get(n)
  theDataLength <- dim(theData)[1]
  if(is.null(theData) != TRUE){
    for(j in 1:theDataLength){
      theRow <- theData[j,]
      occidSearch <- paste("^", theRow$occid, "$", sep = "")
      occidCheck <- grep(occidSearch, resultsTable$occid)
      occidCheckLength <- length(occidCheck)
      if(occidCheckLength == 0){
        addRow <- list(theRow$occid, theRow$collid, theRow$dbpk, theRow$catalogNumber, theRow$basisOfRecord, theRow$sciname, theRow$scientificNameAuthorship, theRow$recordedBy, theRow$eventDate, theRow$year, theRow$month, theRow$day, theRow$verbatimEventDate, theRow$country, theRow$stateProvince, theRow$county, theRow$municipality, theRow$locality, theRow$decimalLatitude, theRow$decimalLongitude, theRow$verbatimCoordinates, theRow$geodeticDatum, theRow$coordinateUncertaintyInMeters, theRow$locationRemarks, theRow$georeferenceProtocol, theRow$georeferenceSources, theRow$georeferenceVerificationStatus, theRow$georeferenceRemarks, theRow$recordIDcharacter, theRow$habitat, theRow$substrate, theRow$associatedTaxa, theRow$fieldnumber, theRow$fieldNotes, theRow$occurrenceRemarks)
        resultsTable <- rbind(resultsTable, addRow)
        numberResults <- numberResults + 1
        print(numberResults)
        n <- n + 1
      }
    }
  }
}

#gather the records from MyCoPortal API with over 400K records (e.g., California)

if(numberResults != 0){
  while(numberResults < mp_api_resultCount){
    if(n == 299999){
      n <- 0
    }else{
      for(n in 1:299999){
        theData <- mp_api_get(n)
        theDataLength <- dim(theData)[1]
        if(is.null(theData) != TRUE){
          for(j in 1:theDataLength){
            theRow <- theData[j,]
            occidSearch <- paste("^", theRow$occid, "$", sep = "")
            occidCheck <- grep(occidSearch, resultsTable$occid)
            occidCheckLength <- length(occidCheck)
            if(occidCheckLength == 0){
              addRow <- list(theRow$occid, theRow$collid, theRow$dbpk, theRow$catalogNumber, theRow$basisOfRecord, theRow$sciname, theRow$scientificNameAuthorship, theRow$recordedBy, theRow$eventDate, theRow$year, theRow$month, theRow$day, theRow$verbatimEventDate, theRow$country, theRow$stateProvince, theRow$county, theRow$municipality, theRow$locality, theRow$decimalLatitude, theRow$decimalLongitude, theRow$verbatimCoordinates, theRow$geodeticDatum, theRow$coordinateUncertaintyInMeters, theRow$locationRemarks, theRow$georeferenceProtocol, theRow$georeferenceSources, theRow$georeferenceVerificationStatus, theRow$georeferenceRemarks, theRow$recordIDcharacter, theRow$habitat, theRow$substrate, theRow$associatedTaxa, theRow$fieldnumber, theRow$fieldNotes, theRow$occurrenceRemarks)
              resultsTable <- rbind(resultsTable, addRow)
              numberResults <- numberResults + 1
              print(numberResults)
            }
          }
        }
      }
    }
  }
}


# create taxon elements for checking each element for validity
#numberResults <- dim(resultsTable)[1]

for(g in 1:numberResults){
  print(g)
  resultsRow <- resultsTable[g]
  theOccid <- resultsRow$occid
  theCollid <- resultsRow$collid
  theOriginalName <- resultsRow$sciname
  theOriginalNameAuthorship <- resultsRow$scientificNameAuthorship
  theElements <- strsplit(theOriginalName, " ")
  theElements <- theElements[[1]]
  theElementsLength <- length(theElements)
  if(theElementsLength < 3){
    theElement1 <- theElements[1]
    theElement1 <- tolower(theElement1)
    theElement1 <- paste(toupper(substr(theElement1, 1, 1)), substr(theElement1, 2, nchar(theElement1)), sep="")
    theElement2 <- theElements[2]
    theElement2 <- tolower(theElement2)
    theElement3 <- NA
    theElement4 <- NA
    infraRank <- NA
  }else{
    theElement1 <- theElements[1]
    theElement1 <- tolower(theElement1)
    theElement1 <-paste(toupper(substr(theElement1, 1, 1)), substr(theElement1, 2, nchar(theElement1)), sep="")
    theElement2 <- theElements[2]
    theElement2 <- tolower(theElement2)
    theElement3 <- theElements[3]
    theElement3 <- tolower(theElement3)
    theElement4 <- theElements[4]
    theElement4 <- tolower(theElement4)
    infraCheck <- theElement3
    infraCheck <- gsub("\\(", "", infraCheck)
    infraCheck <- gsub("\\)", "", infraCheck)
    infraCheck <- gsub("\\[", "", infraCheck)
    infraCheck <- gsub("\\]", "", infraCheck)
    infraFormCheck <- grep(infraCheck, formInfraRank)
    infraFormCheckLength <- length(infraFormCheck)
    if(infraFormCheckLength != 0){
      infraRank <- "f."
    }
    infraSubSpCheck <- grep(infraCheck, subspInfraRank)
    infraSubSpCheckLength <- length(infraSubSpCheck)
    if(infraSubSpCheckLength != 0){
      infraRank <- "subsp."
    }
    infraVarCheck <- grep(infraCheck, varInfraRank)
    infraVarCheckLength <- length(infraVarCheck)
    if(infraVarCheckLength != 0){
      infraRank <- "var."
    }
  }
  addRow <- list(occid=theOccid, collid=theCollid, originalName=theOriginalName, originalNameAuthorship=theOriginalNameAuthorship, element1=theElement1, element2=theElement2, element3=theElement3, element4=theElement4, elementsLength=theElementsLength, validGenusName=NA, validSpeciesName=NA, validInfraspecificRank=infraRank, validInfraspecificName=NA, genusValid=NA, speciesValid=NA, infraspecificValid=NA, currentName=NA, currentAuthors=NA, currentMbNumber=NA, type=NA)
  taxonHoldTable <- rbind(taxonHoldTable, addRow)
}


nonSpeciesTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength == 1)
taxonHoldTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength != 1)
genericNames <- unique(taxonHoldTable$element1)
genericNamesLength <- length(genericNames)

for(h in 1:genericNamesLength){
  print(h)
  theGenus <- genericNames[h]
  genusCheck <- if_api_get("1", theGenus)
  genusCheckLength <- length(genusCheck)
  if(genusCheckLength == 1){
    theType <- slime_mold_check(theGenus, genusCheck)
    exactGenusSearch <- paste("^", theGenus, "$", sep = "")
    foundRows <- grep(exactGenusSearch, taxonHoldTable$element1)
    foundRowsLength <- length(foundRows)
    if(theType == "Slime Mold"){
      for(i in 1:foundRowsLength){
        foundRowsNumber <- foundRows[i]
        taxonHoldTable[foundRowsNumber, 10] <- theGenus
        taxonHoldTable[foundRowsNumber, 14] <- "yes"
        taxonHoldTable[foundRowsNumber, 20] <- theType
      }
    }else{
      if(theType == "Valid Genus"){
        for(i in 1:foundRowsLength){
          foundRowsNumber <- foundRows[i]
          taxonHoldTable[foundRowsNumber, 10] <- theGenus
          taxonHoldTable[foundRowsNumber, 14] <- "yes"
        }
      }
    }
  }
}


badTaxaTable <- subset(taxonHoldTable, is.na(taxonHoldTable$validGenusName) == TRUE)
badTaxaNames <- unique(badTaxaTable$element1)
badTaxaNamesLength <- length(badTaxaNames)

for(j in 1:badTaxaNamesLength){
  print(j)
  badGenus <- badTaxaNames[j]
  badGenusBegins <- substring(badGenus, 1, 3)
  badGenusBegins <- paste(badGenusBegins, "%", sep = "")
  if(badGenus != "_"){
    genusSuggestions <- fd_api_get("taxon", badGenusBegins)
    genusSuggestions <- subset(genusSuggestions, genusSuggestions$rank == "Genus")
    genusSuggestions <- unique(genusSuggestions$taxon)
    genusSuggestionsLength <- length(genusSuggestions)
  }else{
    genusSuggestionsLength <- 0
  }
  if(genusSuggestionsLength > 0){
    genusMatch <- closest_match(badGenus, genusSuggestions)
    theStringSim <- stringsim(badGenus, genusMatch, method = "jaccard")
    if(theStringSim > 0.90){
      genusCheck <- if_api_get("1", genusMatch)
      genusCheckLength <- length(genusCheck)
      if(genusCheckLength == 1){
        theType <- slime_mold_check(genusMatch, genusCheck)
        exactGenusSearch <- paste("^", badGenus, "$", sep = "")
        foundRows <- grep(exactGenusSearch, taxonHoldTable$element1)
        foundRowsLength <- length(foundRows)
        if(theType == "Slime Mold"){
          for(i in 1:foundRowsLength){
            foundRowsNumber <- foundRows[i]
            taxonHoldTable[foundRowsNumber, 10] <- genusMatch
            taxonHoldTable[foundRowsNumber, 14] <- "yes"
            taxonHoldTable[foundRowsNumber, 20] <- "Slime Mold"
          }
        }else{
          if(theType == "Valid Genus"){
            for(i in 1:foundRowsLength){
              foundRowsNumber <- foundRows[i]
              taxonHoldTable[foundRowsNumber, 10] <- genusMatch
              taxonHoldTable[foundRowsNumber, 14] <- "yes"
            }
          }
        }
      }
    }
  }
}

slimeMoldTable <- subset(taxonHoldTable, taxonHoldTable$type == "Slime Mold")
taxonHoldTable <- subset(taxonHoldTable, is.na(taxonHoldTable$type) == TRUE)
valdidGenericNames <- unique(taxonHoldTable$validGenusName)
valdidGenericNamesLength <- length(valdidGenericNames)

for(k in 1:valdidGenericNamesLength){
  print(k)
  validGenus <- valdidGenericNames[k]
  fgCheck <- fg_api_get(validGenus)
  if(is.null(fgCheck) != TRUE){
    theType <- macrofungi_check(fgCheck)
    if(theType != "Macrofungus"){
      exactValidGenusSearch <- paste("^", validGenus, "$", sep = "")
      foundRows <- grep(exactValidGenusSearch, taxonHoldTable$validGenusName)
      foundRowsLength <- length(foundRows)
      for(l in 1:foundRowsLength){
        foundRowsNumber <- foundRows[l]
        taxonHoldTable[foundRowsNumber, 20] <- theType
      }
    }
  }
}

lichenTable <- subset(taxonHoldTable, taxonHoldTable$type == "Lichen")
microfungusTable <- subset(taxonHoldTable, taxonHoldTable$type == "Microfungus")
taxonHoldTable <- subset(taxonHoldTable, is.na(taxonHoldTable$type) == TRUE)
taxonHoldTableLength <- dim(taxonHoldTable)[1]

name_count <- unique(taxonHoldTable$originalName)
name_count_length <- length(name_count)

speciesTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength >= 2)
speciesNames <- unique(speciesTable$originalName)
speciesNamesLength <- length(speciesNames)

for(m in 1:speciesNamesLength){
  print(m)
  theSpecies <- speciesNames[m]
  theSpecies <- gsub("\\(", "\\\\(", theSpecies)
  theSpecies <- gsub("\\)", "\\\\)", theSpecies)
  exactSpeciesSearch <- paste("^", theSpecies, "$", sep = "")
  foundRows <- grep(exactSpeciesSearch, taxonHoldTable$originalName)
  foundRowsLength <- length(foundRows)
  theGenusRow <- taxonHoldTable[foundRows[1], ]
  theGenus <- theGenusRow$validGenusName
  theSpecificEpithet <- theGenusRow$element2
  checkSpecies <- paste(theGenus, theSpecificEpithet, sep = "%20")
  if(is.na(theGenus) == FALSE){
    speciesCheck <- fd_api_get("taxon", checkSpecies)
    speciesCheckLength <- dim(speciesCheck)[1]
    if(is.null(speciesCheck) == TRUE){
      speciesCheckLength <- 0
    }
    if(speciesCheckLength > 1){
      speciesCheck <- correct_record_get(speciesCheck)
      speciesCheckLength <- dim(speciesCheck)[1]
    }
    if(is.null(speciesCheck) == TRUE){
      speciesCheckLength <- 0
    }
    if(speciesCheckLength == 1){
      for(l in 1:foundRowsLength){
        foundRowsNumber <- foundRows[l]
        taxonHoldTable[foundRowsNumber, 11] <- theSpecificEpithet
        taxonHoldTable[foundRowsNumber, 15] <- "yes"
        if(speciesCheck$mbNumber == speciesCheck$currentMbNumber){
          taxonHoldTable[foundRowsNumber, 17] <- speciesCheck$taxon
          taxonHoldTable[foundRowsNumber, 18] <- speciesCheck$authors
          taxonHoldTable[foundRowsNumber, 19] <- speciesCheck$mbNumber
        }else{
          currentSpeciesCheck <- fd_api_get("mbNumber", speciesCheck$currentMbNumber)
          currentSpeciesCheckLength <- dim(currentSpeciesCheck)[1]
          if(is.null(currentSpeciesCheck) == TRUE){
            currentSpeciesCheckLength <- 0
          }
          if(currentSpeciesCheckLength == 1){
            taxonHoldTable[foundRowsNumber, 17] <- currentSpeciesCheck$taxon
            taxonHoldTable[foundRowsNumber, 18] <- currentSpeciesCheck$authors
            taxonHoldTable[foundRowsNumber, 19] <- currentSpeciesCheck$mbNumber
          }
        }
      }
    }else{
      badSpeciesBegins <- substring(theSpecificEpithet, 1, 3)
      badSpeciesBegins <- paste(badSpeciesBegins, "%", sep = "")
      badSpeciesCheck <- paste(theGenus, badSpeciesBegins, sep = "%20")
      badSpecies <- paste(theGenus, theSpecificEpithet, sep = " ")
      speciesSuggestionRecords <- fd_api_get("taxon", badSpeciesCheck)
      speciesSuggestionRecords <- subset(speciesSuggestionRecords, speciesSuggestionRecords$rank == "Species")
      speciesSuggestions <- unique(speciesSuggestionRecords$taxon)
      speciesSuggestionsLength <- length(speciesSuggestions)
      if(speciesSuggestionsLength > 0){
        speciesMatch <- closest_match(badSpecies, speciesSuggestions)
        theStringSim <- stringsim(badSpecies, speciesMatch, method = "jaccard")
        if(theStringSim > 0.83){
          speciesMatchRecord <- subset(speciesSuggestionRecords, speciesSuggestionRecords$taxon == speciesMatch)
          speciesMatchRecordLength <- dim(speciesMatchRecord)[1]
          if(speciesMatchRecordLength > 1){
            speciesMatchRecord <- correct_record_get(speciesMatchRecord)
            speciesMatchRecordLength <- dim(speciesMatchRecord)[1]
          }
          if(is.null(speciesMatchRecord) == TRUE){
            speciesMatchRecordLength <- 0
          }
          if(speciesMatchRecordLength == 1){
            for(l in 1:foundRowsLength){
              foundRowsNumber <- foundRows[l]
              taxonHoldTable[foundRowsNumber, 11] <- theSpecificEpithet
              taxonHoldTable[foundRowsNumber, 15] <- "yes"
              if(speciesMatchRecord$mbNumber == speciesMatchRecord$currentMbNumber){
                taxonHoldTable[foundRowsNumber, 17] <- speciesMatchRecord$taxon
                taxonHoldTable[foundRowsNumber, 18] <- speciesMatchRecord$authors
                taxonHoldTable[foundRowsNumber, 19] <- speciesMatchRecord$mbNumber
              }else{
                currentSpeciesMatchCheck <- fd_api_get("mbNumber", speciesMatchRecord$currentMbNumber)
                currentSpeciesMatchCheckLength <- dim(currentSpeciesMatchCheck)[1]
                if(is.null(currentSpeciesMatchCheck) == TRUE){
                  currentSpeciesMatchCheckLength <- 0
                }
                if(currentSpeciesMatchCheckLength == 1){
                  taxonHoldTable[foundRowsNumber, 17] <- currentSpeciesMatchCheck$taxon
                  taxonHoldTable[foundRowsNumber, 18] <- currentSpeciesMatchCheck$authors
                  taxonHoldTable[foundRowsNumber, 19] <- currentSpeciesMatchCheck$mbNumber
                }
              }
            }
          }
        }
      }
    }
  }
}

badSpeciesTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength >= 2)
badSpeciesTable <- subset(badSpeciesTable, is.na(badSpeciesTable$speciesValid) == TRUE)
badSpeciesNames <- unique(badSpeciesTable$originalName)
badSpeciesNamesLength <- length(badSpeciesNames)

for(n in 1:badSpeciesNamesLength){
  print(n)
  badSpeciesName <- badSpeciesNames[n]
  badSpeciesNameCharacterCheck <- grepl("[^ -~]", badSpeciesName)
  if(badSpeciesNameCharacterCheck == TRUE){
    searchFindLength <- 0
  }else{
    bingSearchName <- strsplit(badSpeciesName, " ")
    bingSearchName <- bingSearchName[[1]]
    bingSearchName <- paste(bingSearchName, collapse = "+")
    bingURL <- paste("https://www.bing.com/search?q=", bingSearchName, sep = "")
    possibleNamePage <- try(read_html(bingURL), silent = TRUE)
    possibleNamePageHold <- html_nodes(possibleNamePage, "div")
    searchFind <- grep('href="https://en.wikipedia.org/wiki/', possibleNamePageHold)
    searchFindLength <- length(searchFind)
  }
  if(searchFindLength >= 1){
    rando <- round(runif(1, min = 0, max = 15), 0)
    Sys.sleep(rando)
    searchFind <- searchFind[1]
    wikiHold <- as.character(possibleNamePageHold[searchFind])
    wikiHold <- strsplit(wikiHold, "/wiki/")
    wikiHold <- wikiHold[[1]][2]
    wikiHold <- strsplit(wikiHold, '"')
    wikiHold <- wikiHold[[1]][1]
    wikiHold <- strsplit(wikiHold, "_")
    wikiHold <- wikiHold[[1]]
    possibleName <- paste(wikiHold, collapse = " ")
    theStringSim <- stringsim(badSpeciesName, possibleName, method = "jaccard")
    if(theStringSim > 0.83){
      foundNameCheck <- strsplit(possibleName, " ")
      foundNameCheck <- foundNameCheck[[1]]
      foundNameCheck <- paste(foundNameCheck, collapse = "%20")
      foundNameCheckRecords <- fd_api_get("taxon", foundNameCheck)
      foundNameCheckRecords <- subset(foundNameCheckRecords, foundNameCheckRecords$rank == "Species")
      foundNameCheckRecordsLength <- dim(foundNameCheckRecords)[1]
      if(is.null(foundNameCheckRecords) == TRUE){
        foundNameCheckRecordsLength <- 0
      }
      if(foundNameCheckRecordsLength >= 1){
        if(foundNameCheckRecordsLength > 1){
          foundNameCheckRecords <- correct_record_get(foundNameCheckRecords)
          foundNameCheckRecordsLength <- dim(foundNameCheckRecords)[1]
          if(is.null(foundNameCheckRecords) == TRUE){
            foundNameCheckRecordsLength <- 0
          }
        }
        if(foundNameCheckRecordsLength == 1){
          exactSpeciesSearch <- paste("^", badSpeciesName, "$", sep = "")
          foundRows <- grep(exactSpeciesSearch, taxonHoldTable$originalName)
          foundRowsLength <- length(foundRows)
          insertName <- strsplit(possibleName, " ")
          insertName <- insertName[[1]]
          insertGenus <- insertName[1]
          insertSpecies <- insertName[2]
          genusCheck <- if_api_get("1", insertGenus)
          genusCheckLength <- length(genusCheck)
          if(genusCheckLength == 1){
            theType <- slime_mold_check(insertGenus, genusCheck)
            if(theType != "Slime Mold"){
              fgCheck <- fg_api_get(insertGenus)
              theType <- macrofungi_check(fgCheck)
            }
          }
          for(l in 1:foundRowsLength){
            foundRowsNumber <- foundRows[l]
            theFoundRow <- taxonHoldTable[foundRowsNumber,]
            taxonHoldTable[foundRowsNumber, 10] <- insertGenus
            taxonHoldTable[foundRowsNumber, 11] <- insertSpecies
            taxonHoldTable[foundRowsNumber, 14] <- "yes"
            taxonHoldTable[foundRowsNumber, 15] <- "yes"
            if(foundNameCheckRecords$mbNumber == foundNameCheckRecords$currentMbNumber){
              taxonHoldTable[foundRowsNumber, 17] <- foundNameCheckRecords$taxon
              taxonHoldTable[foundRowsNumber, 18] <- foundNameCheckRecords$authors
              taxonHoldTable[foundRowsNumber, 19] <- foundNameCheckRecords$mbNumber
              if(theType != "Macrofungus"){
                taxonHoldTable[foundRowsNumber, 20] <- theType
                addRow <- list(occid=theFoundRow$occid, collid=theFoundRow$collid, originalName=theFoundRow$originalName, originalNameAuthorship=theFoundRow$originalNameAuthorship, element1=theFoundRow$element1, element2=theFoundRow$element2, element3=theFoundRow$element3, element4=theFoundRow$element4, elementsLength=theFoundRow$elementsLength, validGenusName=theFoundRow$validGenusName, validSpeciesName=theFoundRow$validSpeciesName, validInfraspecificRank=theFoundRow$validInfraspecificRank, validInfraspecificName=theFoundRow$validInfraspecificName, genusValid=theFoundRow$genusValid, speciesValid=theFoundRow$speciesValid, infraspecificValid=theFoundRow$infraspecificValid, currentName=theFoundRow$currentName, currentAuthors=theFoundRow$currentAuthors, currentMbNumber=theFoundRow$currentMbNumber, type=theFoundRow$type)
                if(theType == "Lichen"){
                  lichenTable <- rbind(lichenTable, addRow)
                }
                if(theType == "Microfungus"){
                  microfungusTable <- rbind(microfungusTable, addRow)
                }
                if(theType == "Slime Mold"){
                  slimeMoldTable <- rbind(slimeMoldTable, addRow)
                }
              }
            }else{
              currentSpeciesCheck <- fd_api_get("mbNumber", foundNameCheckRecords$currentMbNumber)
              currentSpeciesCheckLength <- dim(currentSpeciesCheck)[1]
              if(is.null(currentSpeciesCheck) == TRUE){
                currentSpeciesCheckLength <- 0
              }
              if(currentSpeciesCheckLength == 1){
                taxonHoldTable[foundRowsNumber, 17] <- currentSpeciesCheck$taxon
                taxonHoldTable[foundRowsNumber, 18] <- currentSpeciesCheck$authors
                taxonHoldTable[foundRowsNumber, 19] <- currentSpeciesCheck$mbNumber
                if(theType != "Macrofungus"){
                  taxonHoldTable[foundRowsNumber, 20] <- theType
                  addRow <- list(occid=theFoundRow$occid, collid=theFoundRow$collid, originalName=theFoundRow$originalName, originalNameAuthorship=theFoundRow$originalNameAuthorship, element1=theFoundRow$element1, element2=theFoundRow$element2, element3=theFoundRow$element3, element4=theFoundRow$element4, elementsLength=theFoundRow$elementsLength, validGenusName=theFoundRow$validGenusName, validSpeciesName=theFoundRow$validSpeciesName, validInfraspecificRank=theFoundRow$validInfraspecificRank, validInfraspecificName=theFoundRow$validInfraspecificName, genusValid=theFoundRow$genusValid, speciesValid=theFoundRow$speciesValid, infraspecificValid=theFoundRow$infraspecificValid, currentName=theFoundRow$currentName, currentAuthors=theFoundRow$currentAuthors, currentMbNumber=theFoundRow$currentMbNumber, type=theFoundRow$type)
                  if(theType == "Lichen"){
                    lichenTable <- rbind(lichenTable, addRow)
                  }
                  if(theType == "Microfungus"){
                    microfungusTable <- rbind(microfungusTable, addRow)
                  }
                  if(theType == "Slime Mold"){
                    slimeMoldTable <- rbind(slimeMoldTable, addRow)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

taxonHoldTable <- subset(taxonHoldTable, is.na(taxonHoldTable$type) == TRUE)
nonInfraspecificTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength == 3)
taxonHoldTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength != 3)
infraspecificTable <- subset(taxonHoldTable, taxonHoldTable$elementsLength > 3)
infraspecificNames <- unique(infraspecificTable$originalName)
infraspecificNamesLength <- length(infraspecificNames)

for(o in 1:infraspecificNamesLength){
  print(o)
  theInfraspecificName <- infraspecificNames[o]
  theInfraspecificName <- gsub("\\(", "\\\\(", theInfraspecificName)
  theInfraspecificName <- gsub("\\)", "\\\\)", theInfraspecificName)
  exactInfraspecificSearch <- paste("^", theInfraspecificName, "$", sep = "")
  foundRows <- grep(exactInfraspecificSearch, taxonHoldTable$originalName)
  foundRowsLength <- length(foundRows)
  theInfraspecificName <- gsub("\\[", "", theInfraspecificName)
  theInfraspecificName <- gsub("\\]", "", theInfraspecificName)
  theGenusRow <- taxonHoldTable[foundRows[1], ]
  theGenus <- theGenusRow$validGenusName
  theSpecificEpithet <- theGenusRow$validSpeciesName
  if(is.na(theSpecificEpithet) == TRUE){
    theSpecificEpithet <- theGenusRow$element2
  }
  if(is.na(theSpecificEpithet) == TRUE){
    theSpecificEpithet <- "NOTHING"
  }
  theInfraspecificRank <- theGenusRow$element3
  theInfraspecificName <- theGenusRow$element4
  if(is.na(theInfraspecificName) == TRUE){
    theInfraspecificName <- "NOTHING"
  }
  if(theSpecificEpithet == theInfraspecificName){
    checkInfraspecific <- paste(theGenus, theSpecificEpithet, sep = "%20")
    rankLevel <- "species"
  }else{
    checkInfraspecific <- paste(theGenus, theSpecificEpithet, theInfraspecificRank, theInfraspecificName, sep = "%20")
    rankLevel <- "intraspecific"
  }
  infraspecificCheck <- fd_api_get("taxon", checkInfraspecific)
  infraspecificCheckLength <- dim(infraspecificCheck)[1]
  if(is.null(infraspecificCheck) == TRUE){
    infraspecificCheckLength <- 0
  }
  if(infraspecificCheckLength != 0){
    if(infraspecificCheckLength > 1){
      infraspecificCheck <- correct_record_get(infraspecificCheck)
      infraspecificCheckLength <- dim(infraspecificCheck)[1]
    }
    if(is.null(infraspecificCheck) == TRUE){
      infraspecificCheckLength <- 0
    }
    if(infraspecificCheckLength == 1){
      for(l in 1:foundRowsLength){
        foundRowsNumber <- foundRows[l]
        if(rankLevel == "intraspecific"){
          taxonHoldTable[foundRowsNumber, 13] <- theInfraspecificName
          taxonHoldTable[foundRowsNumber, 16] <- "yes"
          if(infraspecificCheck$mbNumber == infraspecificCheck$currentMbNumber){
            taxonHoldTable[foundRowsNumber, 17] <- infraspecificCheck$taxon
            taxonHoldTable[foundRowsNumber, 18] <- infraspecificCheck$authors
            taxonHoldTable[foundRowsNumber, 19] <- infraspecificCheck$mbNumber
          }else{
            currentInfraspecificCheck <- fd_api_get("mbNumber", infraspecificCheck$currentMbNumber)
            currentInfraspecificCheckLength <- dim(currentInfraspecificCheck)[1]
            if(is.null(currentInfraspecificCheck) == TRUE){
              currentInfraspecificCheckLength <- 0
            }
            if(currentInfraspecificCheckLength == 1){
              taxonHoldTable[foundRowsNumber, 17] <- currentInfraspecificCheck$taxon
              taxonHoldTable[foundRowsNumber, 18] <- currentInfraspecificCheck$authors
              taxonHoldTable[foundRowsNumber, 19] <- currentInfraspecificCheck$mbNumber
            }
          }
        }else{
          taxonHoldTable[foundRowsNumber, 7] <- NA
          taxonHoldTable[foundRowsNumber, 8] <- NA
          taxonHoldTable[foundRowsNumber, 9] <- 2
          taxonHoldTable[foundRowsNumber, 12] <- NA
          taxonHoldTable[foundRowsNumber, 13] <- NA
        }
      }
    }
  }else{
    badInfraspecific <- paste(theGenus, theSpecificEpithet, theInfraspecificRank, theInfraspecificName, sep = " ")
    badInfraspecificNameBegins <- substring(theInfraspecificName, 1, 3)
    badInfraspecificNameBegins <- paste(badInfraspecificNameBegins, "%", sep = "")
    badInfraspecificCheck <- paste(theGenus, theSpecificEpithet, theInfraspecificRank, badInfraspecificNameBegins, sep = "%20")
    infraspecificSuggestionRecords <- fd_api_get("taxon", badInfraspecificCheck)
    infraspecificSuggestions <- unique(infraspecificSuggestionRecords$taxon)
    infraspecificSuggestionsLength <- length(infraspecificSuggestions)
    if(is.null(infraspecificSuggestions) == TRUE){
      infraspecificSuggestionsLength <- 0
    }
    if(infraspecificSuggestionsLength > 0){
      infraspecificMatch <- closest_match(badInfraspecific, infraspecificSuggestions)
      theStringSim <- stringsim(badInfraspecific, infraspecificMatch, method = "jaccard")
      if(theStringSim > 0.83){
        infraspecificMatchRecord <- subset(infraspecificSuggestionRecords, infraspecificSuggestionRecords$taxon == infraspecificMatch)
        infraspecificMatchRecordLength <- dim(infraspecificMatchRecord)[1]
        if(infraspecificMatchRecordLength > 1){
          infraspecificMatchRecord <- correct_record_get(infraspecificMatchRecord)
          infraspecificMatchRecordLength <- dim(infraspecificMatchRecord)[1]
        }
        if(is.null(infraspecificMatchRecord) == TRUE){
          infraspecificMatchRecordLength <- 0
        }
        if(infraspecificMatchRecordLength == 1){
          for(l in 1:foundRowsLength){
            foundRowsNumber <- foundRows[l]
            taxonHoldTable[foundRowsNumber, 13] <- theInfraspecificName
            taxonHoldTable[foundRowsNumber, 16] <- "yes"
            if(infraspecificMatchRecord$mbNumber == infraspecificMatchRecord$currentMbNumber){
              taxonHoldTable[foundRowsNumber, 17] <- speciesMatchRecord$taxon
              taxonHoldTable[foundRowsNumber, 18] <- speciesMatchRecord$authors
              taxonHoldTable[foundRowsNumber, 19] <- speciesMatchRecord$mbNumber
            }else{
              currentInfraspecificMatchCheck <- fd_api_get("mbNumber", infraspecificMatchRecord$currentMbNumber)
              currentInfraspecificMatchCheckLength <- dim(currentInfraspecificMatchCheck)[1]
              if(is.null(currentInfraspecificMatchCheck) == TRUE){
                currentInfraspecificMatchCheckLength <- 0
              }
              if(currentInfraspecificMatchCheckLength == 1){
                taxonHoldTable[foundRowsNumber, 17] <- currentInfraspecificMatchCheck$taxon
                taxonHoldTable[foundRowsNumber, 18] <- currentInfraspecificMatchCheck$authors
                taxonHoldTable[foundRowsNumber, 19] <- currentInfraspecificMatchCheck$mbNumber
              }
            }
          }
        }
      }
    }
  }
}

taxonHoldTable <- subset(taxonHoldTable, is.na(taxonHoldTable$type) == TRUE)
badTaxaTableFinal <- subset(taxonHoldTable, is.na(taxonHoldTable$currentMbNumber) == TRUE)
finalBadSpeicesTable <- subset(badTaxaTableFinal, badTaxaTableFinal$elementsLength == 2)
finalBadInfraspecificTable <- subset(badTaxaTableFinal, badTaxaTableFinal$elementsLength == 4)
finalTaxonHoldTable <- subset(taxonHoldTable, is.na(taxonHoldTable$currentMbNumber) == FALSE)
finalTaxonHoldTableLength <- dim(finalTaxonHoldTable)[1]

#Create the final data.frame from taxonHoldTable and resultsTable

for(p in 1:finalTaxonHoldTableLength){
  print(p)
  theTaxonHoldRow <- taxonHoldTable[p, ]
  theTaxonHoldOccid <- theTaxonHoldRow$occid
  theResultsRow <- subset(resultsTable, resultsTable$occid == theTaxonHoldOccid)
  addRow1 <- list(occid=theTaxonHoldRow$occid, collid=theTaxonHoldRow$collid, dbpk=theResultsRow$dbpk, catalogNumber=theResultsRow$catalogNumber, basisOfRecord=theResultsRow$basisOfRecord, originalName=theTaxonHoldRow$originalName, originalNameAuthorship=theTaxonHoldRow$originalNameAuthorship, validGenusName=theTaxonHoldRow$validGenusName, validSpeciesName=theTaxonHoldRow$validSpeciesName, validInfraspecificRank=theTaxonHoldRow$validInfraspecificRank, validInfraspecificName=theTaxonHoldRow$validInfraspecificName, currentName=theTaxonHoldRow$currentName, currentAuthors=theTaxonHoldRow$currentAuthors, currentMbNumber=theTaxonHoldRow$currentMbNumber, recordedBy=theResultsRow$recordedBy, eventDate=theResultsRow$eventDate, year=theResultsRow$year, month=theResultsRow$month, day=theResultsRow$day, verbatimEventDate=theResultsRow$verbatimEventDate, country=theResultsRow$country, stateProvince=theResultsRow$stateProvince, county=theResultsRow$county, municipality=theResultsRow$municipality, locality=theResultsRow$locality, decimalLatitude=theResultsRow$decimalLatitude, decimalLongitude=theResultsRow$decimalLongitude, verbatimCoordinates=theResultsRow$verbatimCoordinates, geodeticDatum=theResultsRow$geodeticDatum, coordinateUncertaintyInMeters=theResultsRow$coordinateUncertaintyInMeters, locationRemarks=theResultsRow$locationRemarks, georeferenceProtocol=theResultsRow$georeferenceProtocol, georeferenceSources=theResultsRow$georeferenceSources, georeferenceVerificationStatus=theResultsRow$georeferenceVerificationStatus, georeferenceRemarks=theResultsRow$georeferenceRemarks, recordIDcharacter=theResultsRow$recordIDcharacter, habitat=theResultsRow$habitat, substrate=theResultsRow$substrate, associatedTaxa=theResultsRow$associatedTaxa, fieldnumber=theResultsRow$fieldnumber, fieldNotes=theResultsRow$fieldNotes, occurrenceRemarks=theResultsRow$occurrenceRemarks)
  addRow2 <- list(occid=theTaxonHoldRow$occid, currentName=theTaxonHoldRow$currentName, currentAuthors=theTaxonHoldRow$currentAuthors, currentMbNumber=theTaxonHoldRow$currentMbNumber, recordedBy=theResultsRow$recordedBy, eventDate=theResultsRow$eventDate, year=theResultsRow$year, month=theResultsRow$month, day=theResultsRow$day, verbatimEventDate=theResultsRow$verbatimEventDate, country=theResultsRow$country, stateProvince=theResultsRow$stateProvince, county=theResultsRow$county, municipality=theResultsRow$municipality, locality=theResultsRow$locality, decimalLatitude=theResultsRow$decimalLatitude, decimalLongitude=theResultsRow$decimalLongitude, verbatimCoordinates=theResultsRow$verbatimCoordinates, geodeticDatum=theResultsRow$geodeticDatum, coordinateUncertaintyInMeters=theResultsRow$coordinateUncertaintyInMeters)
  finalDataFrame <- rbind(finalDataFrame, addRow1)
  finalCondensedDataFrame <- rbind(finalCondensedDataFrame, addRow2)
}

#Create the tally table for all macrofungi taxa

finalTaxa <- unique(finalDataFrame$currentName)
finalTaxaLength <- length(finalTaxa)

for(q in 1:finalTaxaLength){
  print(q)
  countTaxon <- finalTaxa[q]
  exactCountTaxon <- paste("^", countTaxon, "$", sep = "")
  taxonCount <- grep(exactCountTaxon, taxonHoldTable$currentName)
  taxonCountLength <- length(taxonCount)
  addRow <- list(currentName = countTaxon, tally = taxonCountLength, state = theState)
  taxaCountTable <- rbind(taxaCountTable, addRow)
}

# Plot top-100 most abundant fungal taxa from tally

topValues <- taxaCountTable %>% slice_max(tally, n = 100)
plotTitle <- paste(theState, "Macrofungi: Top-100 Most Abundant Taxa", sep = " ")
ggplot(topValues, aes(x=reorder(currentName, -tally), y=tally)) + geom_col(fill = "darkgreen") + xlab("") + ylab("Record Numbers") + ggtitle(plotTitle) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, hjust=1))

topValues <- taxaCountTable %>% slice_max(tally, n = 50)
plotTitle <- paste(theState, "Macrofungi: Top-50 Most Abundant Taxa", sep = " ")
ggplot(topValues, aes(x=reorder(currentName, -tally), y=tally)) + geom_col(fill = "darkgreen") + xlab("") + ylab("Record Numbers") + ggtitle(plotTitle) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, hjust=1))

# Create the NatureServe 'Element Occurrence' table for rare taxa (currently 1 or 2 reports in state)

elementTable <- subset(taxaCountTable, taxaCountTable$tally <= 2)
elementTableLength <- dim(elementTable)[1]

for(r in 1:elementTableLength){
  print(r)
  theDistance <- "Undetermined"
  elementTaxonRow <- elementTable[r, ]
  elementTaxon <- elementTaxonRow$currentName
  elementTally <- elementTaxonRow$tally
  exactElementTaxonSearch <- paste("^", elementTaxon, "$", sep = "")
  exactElementTaxonRows <- grep(exactElementTaxonSearch, finalDataFrame$currentName)
  exactElementTaxonRowsLength <- length(exactElementTaxonRows)
  if(exactElementTaxonRowsLength == 1){
    elementFinalDataRow1 <- finalDataFrame[exactElementTaxonRows[1], ]
    occid1 <- elementFinalDataRow1$occid
    latitude1 <- elementFinalDataRow1$decimalLatitude
    longitude1 <- elementFinalDataRow1$decimalLongitude
    county1 <- elementFinalDataRow1$county
    county1 <- county_check(county1)
    occid2 <- NA
    latitude2 <- NA
    longitude2 <- NA
    county2 <- NA
    theStatus <- "One distinct element"
  }else{
    doHaversine <- "No"
    elementFinalDataRow1 <- finalDataFrame[exactElementTaxonRows[1], ]
    elementFinalDataRow2 <- finalDataFrame[exactElementTaxonRows[2], ]
    occid1 <- elementFinalDataRow1$occid
    latitude1 <- elementFinalDataRow1$decimalLatitude
    longitude1 <- elementFinalDataRow1$decimalLongitude
    county1 <- elementFinalDataRow1$county
    county1 <- county_check(county1)
    occid2 <- elementFinalDataRow2$occid
    latitude2 <- elementFinalDataRow2$decimalLatitude
    longitude2 <- elementFinalDataRow2$decimalLongitude
    county2 <- elementFinalDataRow2$county
    county2 <- county_check(county2)
    if(is.na(latitude1) == FALSE && is.na(longitude1) == FALSE){
      if(is.na(latitude2) == FALSE && is.na(longitude2) == FALSE){
        doHaversine <- "Yes"
      }
    }
    if(doHaversine == "Yes"){
      lat1 <- latitude1
      lat1 <- conv_unit(lat1, "dec_deg", "deg_min_sec")
      lat1 <- strsplit(lat1, " ")
      lat1 <- lat1[[1]]
      lat1Deg <- as.numeric(lat1[1])
      if(lat1Deg > 180){
        doHaversine <- "No"
      }
      lat1Min <- as.numeric(lat1[2])
      if(lat1Min == 0){
        doHaversine <- "No"
      }
      lat1Sec <- round(as.numeric(lat1[3]), 0)
      if(lat1Sec == 60){
        lat1Sec <- 59
      }
      long1 <- as.numeric(longitude1)
      if(long1 < 0){
        long1 <- long1 * -1
      }
      long1 <- conv_unit(long1, "dec_deg", "deg_min_sec")
      long1 <- strsplit(long1, " ")
      long1 <- long1[[1]]
      long1Deg <- as.numeric(long1[1])
      if(long1Deg > 180){
        doHaversine <- "No"
      }
      long1Min <- as.numeric(long1[2])
      if(long1Min == 0){
        doHaversine <- "No"
      }
      long1Sec <- round(as.numeric(long1[3]), 0)
      if(long1Sec == 60){
        long1Sec <- 59
      }
      lat_long1 <- paste(lat1Deg, " ", lat1Min, " ", lat1Sec, "N, ", long1Deg, " ", long1Min, " ", long1Sec, "W", sep = "")
      lat2 <- latitude2
      lat2 <- conv_unit(lat2, "dec_deg", "deg_min_sec")
      lat2 <- strsplit(lat2, " ")
      lat2 <- lat2[[1]]
      lat2Deg <- as.numeric(lat2[1])
      if(lat2Deg > 180){
        doHaversine <- "No"
      }
      lat2Min <- as.numeric(lat2[2])
      if(lat2Min == 0){
        doHaversine <- "No"
      }
      lat2Sec <- round(as.numeric(lat2[3]), 0)
      if(lat2Sec == 60){
        lat2Sec <- 59
      }
      long2 <- as.numeric(longitude2)
      if(long2 < 0){
        long2 <- long2 * -1
      }
      long2 <- conv_unit(long2, "dec_deg", "deg_min_sec")
      long2 <- strsplit(long2, " ")
      long2 <- long2[[1]]
      long2Deg <- as.numeric(long2[1])
      if(long2Deg > 180){
        doHaversine <- "No"
      }
      long2Min <- as.numeric(long2[2])
      if(long2Min == 0){
        doHaversine <- "No"
      }
      long2Sec <- round(as.numeric(long2[3]), 0)
      if(long2Sec == 60){
        long2Sec <- 59
      }
      lat_long2 <- paste(lat2Deg, " ", lat2Min, " ", lat2Sec, "N, ", long2Deg, " ", long2Min, " ", long2Sec, "W", sep = "")
      if(doHaversine == "Yes"){
        theDistance <- haversine(lat_long1, lat_long2)
        #theDistance <- round(conv_unit(theDistance, "km", "mi"), 2)
        theDistance <- round(theDistance, 2)
        if(theDistance > 1.0){
          theStatus <- "Two distinct elements"
        }else{
          theStatus <- "Two non-distinct elements"
          theDistance <- "Undetermined"
        }
      }else{
        theDistance <- "Undetermined"
      }
    }else{
      if(is.na(county1) == TRUE || is.na(county2) == TRUE){
        theStatus <- "Two non-distinct elements"
      }else{
        if(county1 == county2){
          theStatus <- "Two non-distinct elements"
        }else{
          theStatus <- "Two distinct elements"
        }
      }
    }
  }
  addRow <- list(currentName=elementTaxon, currentAuthors=elementFinalDataRow1$currentAuthors, currentMbNumber=elementFinalDataRow1$currentMbNumber, tally=elementTally, distanceKm=theDistance, status=theStatus, occid_1=occid1, decimalLatitude_1=latitude1, decimalLongitude_1=longitude1, county_1=county1, occid_2=occid2, decimalLatitude_2=latitude2, decimalLongitude_2=longitude2, county_2=county2)
  finalRareTaxaTable <- rbind(finalRareTaxaTable, addRow)
}

# Create the Red List table

redListFungiNames <- unique(redListFungi$scientificName)
redListFungiNamesLength <- length(redListFungiNames)

for(s in 1:redListFungiNamesLength){
  print(s)
  redListName <- redListFungiNames[s]
  exactRedListNameSearch <- paste("^", redListName, "$", sep = "")
  foundRowRedList <- grep(exactRedListNameSearch, taxaCountTable$currentName)
  foundRowRedListLength <- length(foundRowRedList)
  if(foundRowRedListLength !=0){
    foundRowFinalData <- grep(exactRedListNameSearch, finalDataFrame$currentName)
    theRowFinalData <- finalDataFrame[foundRowFinalData[1], ]
    theRowRedList <- taxaCountTable[foundRowRedList, ]
    addRowRedList <- list(currentName = redListName, tally = theRowRedList$tally)
    redListCountTable <- rbind(redListCountTable, addRowRedList)
    foundRowFinalRareTaxaTable <- grep(exactRedListNameSearch, finalRareTaxaTable$currentName)
    foundRowFinalRareTaxaTableLength <- length(foundRowFinalRareTaxaTable)
    if(foundRowFinalRareTaxaTableLength > 0){
      theRowFinalRareTaxaTable <- finalRareTaxaTable[foundRowFinalRareTaxaTable, ]
      rowStatus <- theRowFinalRareTaxaTable$status
      newStatus <- paste(rowStatus, "Fungal Red List", sep = "; ")
      finalRareTaxaTable[foundRowFinalRareTaxaTable, 6] <- newStatus
    }else{
      addRowRareTaxa <- list(currentName=redListName, currentAuthors=theRowFinalData$currentAuthors, currentMbNumber=theRowFinalData$currentMbNumber, tally=theRowRedList$tally, distanceKm="Undetermined", status="Fungal Red List", occid_1=NA, decimalLatitude_1=NA, decimalLongitude_1=NA, county_1=NA, occid_2=NA, decimalLatitude_2=NA, decimalLongitude_2=NA, county_2=NA)
      finalRareTaxaTable <- rbind(finalRareTaxaTable, addRowRareTaxa)
    }
  }
}

# Output all the results tables

setwd('~')
directoryPath <- getwd()
thePath <- paste(directoryPath, '/Desktop/', theState, sep = "")
setwd(thePath)

finalFileName <- paste("Macrofungi_of_", theState, ".csv", sep = "")
finalCondensedFileName <- paste("Macrofungi_of_", theState, "_Condensed", ".csv", sep = "")
finalResultsFileName <- paste(theState, "_raw_data.csv")
finalLichenFileName <- paste(theState, "_lichens.csv", sep = "")
finalMicrofungiFileName <- paste(theState, "_microfungi.csv", sep = "")
finalSlimeMoldFileName <- paste(theState, "_slime_molds.csv", sep = "")
finalNonSpeciesFileName <- paste(theState, "_non_species.csv")
finalBadInfraspecificFileName <- paste(theState, "_bad_Infraspecific.csv", sep = "")
finalBadSpeciesFileName <- paste(theState, "_bad_species.csv", sep = "")
finalBadTaxaFileName <- paste(theState, "_bad_taxa.csv", sep = "")
finalTaxaCountFileName <- paste(theState, "_taxa_tally.csv", sep = "")
finalTaxonHoldFileName <- paste(theState, "_taxa_hold.csv", sep = "")
finalRedListCountFileName <- paste(theState, "_red_list_taxa.csv", sep = "")
finalRareTaxaTableFileName <- paste(theState, "_rare_and_red_list_taxa.csv", sep = "")

finalDataTable <- setDT(finalDataFrame)
finalCondensedDataTable <- setDT(finalCondensedDataFrame)
finalLichenTable <- setDT(lichenTable)
finalMicrofungiTable <- setDT(microfungusTable)
finalSlimeMoldTable <- setDT(slimeMoldTable)
finalNonSpecies <- setDT(nonSpeciesTable)
finalBadInfraspecific <- setDT(finalBadInfraspecificTable)
finalBadSpeicesTable <- setDT(finalBadSpeicesTable)
finalBadTaxaTable <- setDT(badTaxaTableFinal)
finalTaxonHoldTable <- setDT(finalTaxonHoldTable)
finalTaxaCountTable <- setDT(taxaCountTable)
finalRedListCountTable <- setDT(redListCountTable)
finalRareTaxaTable <- setDT(finalRareTaxaTable)

fwrite(finalDataTable, file = finalFileName, row.names = TRUE)
fwrite(finalCondensedDataTable, file = finalCondensedFileName, row.names = TRUE)
fwrite(resultsTable, file = finalResultsFileName, row.names = TRUE)
fwrite(finalLichenTable, file = finalLichenFileName, row.names = TRUE)
fwrite(finalMicrofungiTable, file = finalMicrofungiFileName, row.names = TRUE)
fwrite(finalSlimeMoldTable, file = finalSlimeMoldFileName, row.names = TRUE)
fwrite(finalNonSpecies, file = finalNonSpeciesFileName, row.names = TRUE)
fwrite(finalBadInfraspecific, file = finalBadInfraspecificFileName, row.names = TRUE)
fwrite(finalBadSpeicesTable, file = finalBadSpeciesFileName, row.names = TRUE)
fwrite(finalBadTaxaTable, file = finalBadTaxaFileName, row.names = TRUE)
fwrite(finalTaxonHoldTable, file = finalTaxonHoldFileName, row.names = TRUE)
fwrite(finalTaxaCountTable, file = finalTaxaCountFileName, row.names = TRUE)
fwrite(finalRedListCountTable, file = finalRedListCountFileName, row.names = TRUE)
fwrite(finalRareTaxaTable, file = finalRareTaxaTableFileName, row.names = TRUE)

