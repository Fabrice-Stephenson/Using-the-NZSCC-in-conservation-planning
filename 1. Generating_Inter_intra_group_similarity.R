##====== USING THE NZSCC FOR SPATIAL PLANNING   ==============================##

# SCRIPT 1. WITHIN (INTRA) AND BETWEEN (INTER) GROUP SIMILARITY 

##------ Authors: Fabrice Stephenson & Jordi Tablada based on code by John Leathwick
##------ Start date : 28/02/2024
##------ End date : 01/01/2025

##============================================================================##
# This code was produced as part of a DOC funded project to assess the use
# of the the New Zealand Seafloor Community Classification (NZSCC) and 
# associated biodiversity layers in Systematic Conservation Planning described
# in:
# 1.  Stephenson et al. (2021) Species composition and turnover models provide 
#     robust approximations of biodiversity in marine conservation planning. 
#     Ocean and Coastal Management.
# 2.  Tablada, J.,Geange, S. & Stephenson, F. (in review) Evaluation of current 
#     spatial management areas using the New Zealand Seafloor Community 
#     Classification. Aquatic Conservation: Marine And Freshwater Ecosystems
# 3.  Tablada, J.,Geange, S. Hiddink, JG. & Stephenson, F. (draft). Enhancing 
#     representativity of protected seafloor communities in Aotearoa New Zealand

##============================================================================##

# DESCRIPTION: Output spatial layers of within (intra) and between (inter) group 
# similarity for the NZSCC.

# 1.  Load files and packages
# 2.  Loop through the NZSCC and calculate within (intra) and between (inter)
#     Group similarity / dissimilarity for use in Zonation (Spatial prioritisation)

####==========    1. LOAD  PACKAGES & DATA  ================================####
rm(list = ls())
require(raster); require(cluster)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))
mask <- raster::raster("Template_1km.tif") ### Getting FID from raster to sum sightings
plot(r)
mask[mask > 1] <- 1
plot(mask)

MPIproj <- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
imp.vars <- c("Bathy","BedDist", "BotOxy", "BotNi", "BotPhos","BotSal",
              "BotSil", "BotTemp","BPI_broad", "BPI_fine", "ChlAGrad",
              "DET", "PB555nm","SeasTDiff", "Slope", "SSTGrad", "sed.class", 
              "TC", "POCFlux","Ebed")

load("Turnover_CMB.source")
load('EEZClaraClassification_fullEEZ.source')

n.grp <- 75
EEZSeventyFiveMeans <- matrix(0, nrow = n.grp, ncol = length(imp.vars))
dimnames(EEZSeventyFiveMeans)[[1]] <- paste('Grp_',c(1:n.grp),sep='')
dimnames(EEZSeventyFiveMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) EEZSeventyFiveMeans[,i] <- tapply(Turnover_CMB[,imp.vars[[i]]],
                                                                   EEZClaraClassification$SeventyFiveGroups,mean)

####==========    2. LOOP: CALCULATE INTRA AND INTER GROUP SIM / DISSIM  ===####
# Cell values indicate the distance from individual cells to the group centroid (mediod)
#
# input data are as follows:
#  ...transformed spatial data are contained in Turnover_CMB[,im.vars]
#  ...classification vector is contained in EEZClaraClassification[11] or EEZClaraClassification$SeventyFiveGroups
#  ...group centroids are contained in EEZSeventyFiveMeans

# load the FINAL group names used for the NZSCC
grp.conv <- read.csv2("Group_conversion_model_finalNZSCC.csv", sep = ",")

for (i in 1:nrow(EEZSeventyFiveMeans)) {
  # calculate the Manhattan distance between each cell in the model
  Distances <- apply(abs(matrix(c(t(Turnover_CMB[,imp.vars])) - c(t(EEZSeventyFiveMeans[i,])),
                                nrow(Turnover_CMB), byrow = TRUE)),1,sum)
  
  norm.distance <- Distances/max(Distances) # normalise distance 
  Similarities <- (1 - norm.distance)*100   # convert to similarities
  
  Similarities[EEZClaraClassification$SeventyFiveGroups != i] <- 0  # now set non group cells to zero
  # print(summary(Similarities[Similarities > 0]))
  
  # combine objects with inshore at 250m resolution and offshore at 1km
  CMB_SeventyFiveGroups.DF <- cbind(Turnover_CMB[,c(1:3)], Similarities)
  
  # save as a TIFF - 1 km resolution - $res = 1
  CMB_SeventyFiveGroups.DF.1km <- CMB_SeventyFiveGroups.DF[CMB_SeventyFiveGroups.DF$res ==1,]
  
  CMB_SeventyFiveGroups_1km <- rasterFromXYZ(data.frame(x = CMB_SeventyFiveGroups.DF.1km[,1], 
                                                        y = CMB_SeventyFiveGroups.DF.1km[,2], 
                                                        z = CMB_SeventyFiveGroups.DF.1km[,4]),
                                             crs = MPIproj)
  # save as a TIFF - 275 m resolution - $res = 0
  CMB_SeventyFiveGroups.DF.250m <- CMB_SeventyFiveGroups.DF[CMB_SeventyFiveGroups.DF$res ==0,]
  # export as raster
  CMB_SeventyFiveGroups_250m <- rasterFromXYZ(data.frame(x = CMB_SeventyFiveGroups.DF.250m[,1], 
                                                         y = CMB_SeventyFiveGroups.DF.250m[,2], 
                                                         z = CMB_SeventyFiveGroups.DF.250m[,4]),
                                              crs = MPIproj)

  # average out the 250m resolution and upscale to 1km
  CMB_SeventyFiveGroups_250mAgg <- round(aggregate(CMB_SeventyFiveGroups_250m, fact = c(4,4), fun=median, 
                                                   expand=F, na.rm=TRUE))
  
  NZSCC_Intra <- merge(CMB_SeventyFiveGroups_250mAgg, CMB_SeventyFiveGroups_1km, overlap = T)
  
  # remove artefact in kermadecs and rename to match the groups of NZSCC 
  NZSCC_Intra[is.na( NZSCC_Intra)] <- 0;  NZSCC_Intra <-  NZSCC_Intra * mask
  names(NZSCC_Intra) <- paste0("Intra_",grp.conv[i,2])
  # plot(NZSCC_Intra)
  
  # save raster with correct NZSCC name
  setwd(paste0(dir, "/Intra"))
  writeRaster(NZSCC_Intra, filename =  paste0("Intra_", grp.conv[i,2],".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
  
  # generate group mask for manipulation of intra sim. 
  Grp.mask <- NZSCC_Intra
  Grp.mask[Grp.mask == 0] <- -1
  Grp.mask[Grp.mask > 0] <- 0
  Grp.mask[Grp.mask == -1] <- 1
  # plot(Grp.mask)

  # INTER GROUP DISSIM
  Similarities[EEZClaraClassification$SeventyFiveGroups != i] <- 100 * norm.distance[EEZClaraClassification$SeventyFiveGroups != i]
  
  CMB_SeventyFiveGroups.DF <- cbind(Turnover_CMB[,c(1:3)], Similarities)
  
  # save as a TIFF - 1 km resolution - $res = 1
  CMB_SeventyFiveGroups.DF.1km <- CMB_SeventyFiveGroups.DF[CMB_SeventyFiveGroups.DF$res ==1,]
  
  CMB_SeventyFiveGroups_1km <- rasterFromXYZ(data.frame(x = CMB_SeventyFiveGroups.DF.1km[,1], 
                                                        y = CMB_SeventyFiveGroups.DF.1km[,2], 
                                                        z = CMB_SeventyFiveGroups.DF.1km[,4]),
                                             crs = MPIproj)
  # plot(CMB_SeventyFiveGroups_1km)
  # save as a TIFF - 250 m resolution - $res = 0
  CMB_SeventyFiveGroups.DF.250m <- CMB_SeventyFiveGroups.DF[CMB_SeventyFiveGroups.DF$res ==0,]
  # export as raster
  CMB_SeventyFiveGroups_250m <- rasterFromXYZ(data.frame(x = CMB_SeventyFiveGroups.DF.250m[,1], 
                                                         y = CMB_SeventyFiveGroups.DF.250m[,2], 
                                                         z = CMB_SeventyFiveGroups.DF.250m[,4]),
                                              crs = MPIproj)
  
  CMB_SeventyFiveGroups_250mAgg <- round(aggregate(CMB_SeventyFiveGroups_250m, fact = c(4,4), fun=median, 
                                                   expand=F, na.rm=TRUE))
  
  NZSCC_Inter <- merge(CMB_SeventyFiveGroups_250mAgg, CMB_SeventyFiveGroups_1km, overlap = T)
  NZSCC_Inter[is.na(NZSCC_Inter)] <- 100; NZSCC_Inter <- NZSCC_Inter * mask
  names(NZSCC_Inter) <- paste0("Inter_",grp.conv[i,2])
  # plot(NZSCC_Inter)
  
  # save raster with correct NZSCC name
  setwd(paste0(dir, "/Inter"))
  writeRaster(NZSCC_Inter, filename = paste0(home, "Inter_",grp.conv[i,2],".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
  
  # inter without values for group - inter disimilarties (used for zonation)
  NZSCC_Inter0 <- NZSCC_Inter
  NZSCC_Inter0 <- NZSCC_Inter0 * Grp.mask
  NZSCC_Inter0[is.na(NZSCC_Inter0)] <- 0; NZSCC_Inter0 <- NZSCC_Inter0 * mask
  names(NZSCC_Inter0) <- paste0("Inter0_",grp.conv[i,2])
  # plot(NZSCC_Inter0)
  
  setwd(paste0(dir, "/Inter0_"))
  writeRaster(NZSCC_Inter0, filename = paste0(home, "Inter0_",grp.conv[i,2],".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
  
  # inter similarties (rather than dissimilarities - used for reporting)
  SimInter0 <- NZSCC_Inter0 
  SimInter0 <- 100 - NZSCC_Inter0 
  names(SimInter0) <- paste0("SimInter0_",grp.conv[i,2])
  # plot(SimInter0)
  
  setwd(paste0(dir, "/SimInter"))
  writeRaster(SimInter0, filename = paste0(home, "SimInter0_",grp.conv[i,2],".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
  
  print(paste("Completed NZSCC Group: ", i))
}

rm(Distances, Similarities,GF_SeventyFive_Groups,CMB_SeventyFiveGroups_250mAgg,
   CMB_SeventyFiveGroups_250m, CMB_SeventyFiveGroups.DF.250m, CMB_SeventyFiveGroups_1km,
   CMB_SeventyFiveGroups.DF.1km, CMB_SeventyFiveGroups.DF)