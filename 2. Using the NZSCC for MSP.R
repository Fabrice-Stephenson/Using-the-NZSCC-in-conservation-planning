##====== USING THE NZSCC FOR SPATIAL PLANNING   ==============================##

# SCRIPT 2. USING THE NZSCC FOR REPORTING AND CONSERVATION PLANNING 

##------ Authors: Fabrice Stephenson & Jordi Tablada
##------ Start date : 28/02/2024
##------ End date : 01/01/2025

##============================================================================##
# This code was produced as part of a DOC funded project to assess the use
# of the the New Zealand Seafloor Community Classification (NZSCC) and 
# associated biodiversity layers in Systematic Conservation Planning described
# in:
# 1.  Tablada, J.,Geange, S. & Stephenson, F. (in review) Evaluation of current 
#     spatial management areas using the New Zealand Seafloor Community 
#     Classification. Aquatic Conservation: Marine And Freshwater Ecosystems
# 2.  Tablada, J.,Geange, S. Hiddink, JG. & Stephenson, F. (draft). Enhancing 
#     representativity of protected seafloor communities in Aotearoa New Zealand

##============================================================================##

# DESCRIPTION: Calculate the representativity of Spatial Management Areas as 
# assessed by NZSCC group extent, within and between group proportion and 
# species richness

# 1.  Load files and packages
# 2.  NZSCC group extents within Spatial Management Areas
# 3.  Within and Between group proportion across Spatial Management Areas
# 4.  Species richness within Spatial Management Areas
# 5.  Within and Between group proportion within Spatial Management Areas

####==========    1. LOAD  PACKAGES & DATA  ================================####
rm(list = ls())
require(raster)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(dir, "/Test data"))
mask <- raster::raster("Template_1km.tif") ### Getting FID from raster to sum sightings
mask[mask > 1] <- 1
plot(mask)

MPIproj <- CRS("+proj=aea +lat_1=-30 +lat_2=-50 +lat_0=-40 +lon_0=175 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# BATHY
Bathy <- raster("bathy.tif")*-1
Bathy[Bathy <0] <- 0
plot(Bathy)

# FISHABLE AND UNFISHABLE DEPTHS
Fish.D <- Bathy
Fish.D[Fish.D < 1600] <- 1; Fish.D[Fish.D >= 1600] <- 0
plot(Fish.D)

Un.Fish.D <- Bathy
Un.Fish.D[Un.Fish.D < 1600] <- 0; Un.Fish.D[Un.Fish.D >= 1600] <- 1
plot(Un.Fish.D)

# MASK
mask <- Bathy
mask[mask > 0] <- 1
mask.pt <- as.data.frame(mask, xy = T)
plot(mask)
  
# DEMERSAL FISH RICHNESS
DF.rich <- raster("DF_rich_proj.tif")
plot(DF.rich)

# MANAGEMENT POLYGONS
BPA <- raster("BPA.tif")
BPA[BPA > 0] <- 1; BPA[is.na(BPA)] <- 0; BPA <- BPA * mask
plot(BPA)

CSA <- raster("CSA.tif")
CSA[CSA > 0] <- 1; CSA[is.na(CSA)] <- 0; CSA <- CSA * mask
plot(CSA)

MR <- raster("MR.tif")
MR[MR > 0] <- 1; MR[is.na(MR)] <- 0; MR <- MR * mask
plot(MR)

SMCPP <- raster("SMCPP.tif")
SMCPP[SMCPP > 0] <- 1; SMCPP[is.na(SMCPP)] <- 0; SMCPP <- SMCPP * mask
plot(SMCPP)

All.PA <- BPA + CSA + MR + SMCPP
All.PA[All.PA >= 1] <- 1
plot(All.PA)
sum(na.omit(values(All.PA)))
sum(na.omit(values(mask)))
round((sum(na.omit(values(All.PA)))/sum(na.omit(values(mask)))*100),2)
# plot(All.PA, xlim = c(450000,600000), ylim = c(500000,1500000))

# combine into list of looping 
PA.list <- list(BPA = BPA, CSA= CSA, MR = MR, SMCPP = SMCPP, All.PA = All.PA)
PA.list <- stack(PA.list)
# plot(PA.list[["BPA"]])

# load NZSCC 
NZSCC <- raster("CMB.GF75_EEZ_TS-1km_Final.tif")
plot(NZSCC)

# Test extracting group 20
# grp20 <- NZSCC
# grp20[!grp20 == 31] <- 0; grp20[grp20 == 31] <- 1; grp20 <- grp20 * mask
# plot(grp20)

####==========    2. NZSCC GROUP EXTENT WITHIN SMAs   ======================####
# SPLIT PA INTO FISHABLE AND UNFISHABLE DEPTHS
PA.list.fish <- PA.list * Fish.D
plot(PA.list.fish[[5]])

PA.list.unfish <- PA.list * Un.Fish.D
plot(PA.list.unfish[[5]])

# Calcualte the extent (and percent) of fishable and unfishable PAs 
DF <- setNames(data.frame(matrix(ncol = 4,
                                 nrow = (5))),
               c("PA","Extent_km2", "Fishable","Unfishable"))
# i = 1
for(i in 1:nlayers(PA.list)){
  DF[i,1] <- names(PA.list[[i]])  
  DF[i,2] <- sum(na.omit(values(PA.list[[i]])))
  DF[i,3] <- round((sum(na.omit(values(PA.list.fish[[i]]))) / sum(na.omit(values(PA.list[[i]])))) * 100, 1)
  DF[i,4] <- round((sum(na.omit(values(PA.list.unfish[[i]]))) / sum(na.omit(values(PA.list[[i]])))) *100, 1)
}

# EXTENT AND PERCENT OF NZSCC GROUPS WITHIN PAs, OVERALL, AT FISHABLE AND UNFISHABLE DEPTHS
# EXTENT dataframte
DF.Ext <- setNames(data.frame(matrix(ncol = 16,
                                 nrow = (75))),
               c("GF.grp",
                 "BPA.Tot", "CSA.Tot","MR.Tot","SMCPP.Tot","ALL.Tot",
                 "BPA.Fish", "CSA.Fish","MR.Fish","SMCPP.Fish","ALL.Fish",
                 "BPA.UnFish", "CSA.UnFish","MR.UnFish","SMCPP.UnFish","ALL.UnFish"))
DF.Ext$GF.grp <- seq(1:75)

# PROPORTION dataframe
DF.Prop <- setNames(data.frame(matrix(ncol = 16,
                                     nrow = (75))),
                    c("GF.grp",
                      "BPA.Tot", "CSA.Tot","MR.Tot","SMCPP.Tot","ALL.Tot",
                      "BPA.Fish", "CSA.Fish","MR.Fish","SMCPP.Fish","ALL.Fish",
                      "BPA.UnFish", "CSA.UnFish","MR.UnFish","SMCPP.UnFish","ALL.UnFish"))
DF.Prop$GF.grp <- seq(1:75)

# LOOP THROUGH THE NZSCC GROUPS
# i = 20
# j = 1
for (i in 1:nrow(DF.Ext)){
  # Isolate NZSCC group
  Grp.R <- NZSCC
  Grp.R[!Grp.R == i] <- 0; Grp.R[Grp.R == i] <- 1
  # plot(Grp.R)
  
  for (j in 1:nlayers(PA.list)){ # all PA
    P <- Grp.R * PA.list[[j]]
    DF.Ext[i,j+1] <- sum(na.omit(values(P)))
    DF.Prop[i,j+1] <- round((sum(na.omit(values(P)))/sum(na.omit(values(Grp.R))))*100, 1)
  }
  for (j in 1:nlayers(PA.list)){ # fishable depths in PA
    P <- Grp.R * PA.list.fish[[j]]
    DF.Ext[i,j+6] <- sum(na.omit(values(P)))
    DF.Prop[i,j+6] <- round((sum(na.omit(values(P)))/sum(na.omit(values(Grp.R))))*100, 1)
  }
  for (j in 1:nlayers(PA.list)){# unfishable depths in PA
    P <- Grp.R * PA.list.unfish[[j]]
    DF.Ext[i,j+11] <- sum(na.omit(values(P)))
    DF.Prop[i,j+11] <- round((sum(na.omit(values(P)))/sum(na.omit(values(Grp.R))))*100, 1)
  }
  
  print(paste("Finished calculating for NZSCC Group: ", i))
}

setwd(paste0(dir, "/Outputs"))
write.csv2(DF.Ext, file = "NZSCC_GRP_PA_Ext.csv")
write.csv(DF.Prop, sep = ";", file = "NZSCC_GRP_PA_Prop.csv")

####==========    3. INTER AND INTRA GROUP PROP ACROSS SMAs   ==============####
# Mean and max intra for each NZSCC group
Intra.mat <- (data.frame(matrix(ncol = 3,nrow = 75)))
colnames(Intra.mat) <- c("Min", "Mean", "Max"); rownames(Intra.mat) <- paste("Grp_", seq(1:75))

setwd(dir)
# i = 10
for (i in 1:75){
  Intra <- raster(paste0("Intra/Intra_",i,".tif"))
  Intra[Intra == 0] <- NA
  # plot(Intra)
  q <- quantile(na.omit(values(Intra)), probs = c(0.05, 0.5, 0.95))
  Intra.mat[i,1] <- q[1]
  Intra.mat[i,2] <- q[2]
  Intra.mat[i,3] <- q[3]
  
  print(paste("Finished calculating for NZSCC Group: ", i))
}

setwd(paste0(dir, "/Outputs"))
Intra.mat <- round(Intra.mat, 2)
write.csv(Intra.mat, file = "Intra.mat.csv")

# Mean inter for each GF group 
Sim.mat <- (data.frame(matrix(ncol = 75,nrow = (75))))
colnames(Sim.mat) <- paste("Grp_", seq(1:75)); rownames(Sim.mat) <- paste("Grp_", seq(1:75))

Sim.max <- (data.frame(matrix(ncol = 75,nrow = (75))))
colnames(Sim.max) <- paste("Grp_", seq(1:75)); rownames(Sim.max) <- paste("Grp_", seq(1:75))

setwd(dir)
# i = 20; j = 34
for (i in 1:75){
  Inter0 <- raster(paste0("Sim_Inter/SimInter0_",i,".tif"))
  plot(Inter0)
 
 for (j in 1:75){
   grp2 <- NZSCC; grp2[!grp2 == j] <- NA; grp2[grp2 == j] <- 1;
   grp <- raster(paste0("Sim_Inter/SimInter0_",j,".tif")); grp[!grp == 100] <- NA
   grp[grp == 100] <- 1; grp <- grp * grp2
   # plot(grp)
   # plot(grp2)
   
   Inter <- Inter0 * grp 
   Inter[Inter == 100] <- NA # remove artifacts
   # plot(Inter)
   q <- quantile(na.omit(values(Inter)), probs = c(0.5, 0.975))
   Sim.mat[j,i] <- q[1]
   Sim.max[j,i] <- q[2]
 }
 print(paste("Finished calculating for NZSCC Group: ", i))
}

setwd(paste0(dir, "/Outputs"))
write.csv(Sim.mat, file = "Inter.mean.mat.csv")
write.csv(Sim.max, file = "Inter.max.mat.csv")

####==========    4. SPECIES RICHNESS WITHIN SMAs   ========================####
# Demersal fish richness
df.rich <- raster("DF_rich_proj.tif")
plot(df.rich)

# Benthic Invert richness
bi.rich <- raster("BI_rich_proj.tif")
plot(bi.rich)

# Reef fish richness
rf.rich <- raster("rf_rich_proj.tif")
plot(rf.rich)

# Macroalgae richness
ma.rich <- raster("ma_rich_proj.tif")
plot(ma.rich)

rich.mat <- (data.frame(matrix(ncol = 7,nrow = (75))))
colnames(rich.mat) <- c("min.GRP","mean.GRP","max.GRP","min.PA","mean.PA","max.PA","All.PA.prop.rep"); 
rownames(rich.mat) <- paste("Grp_", seq(1:75))

# i = 20
rich <- ma.rich # CHANGE THIS FOR EACH RICHNESS LAYER
for (i in 1:nrow(rich.mat)){
  # Isolate NZSCC group
  Grp.R <- NZSCC
  Grp.R[!Grp.R == i] <- NA; Grp.R[Grp.R == i] <- 1
  # plot(Grp.R)
  
  # richness in group
  r <- rich * Grp.R
  # plot(r)
  q <- quantile(na.omit(values(r)), probs = c(0.025, 0.5, 0.975))
  
  rich.mat[i,1] <- q[1]
  rich.mat[i,2] <- q[2]
  rich.mat[i,3] <- q[3]
  
  # richness in PA within groups
  P <- Grp.R * PA.list[[5]]
  P[P == 0] <- NA
  # plot(P)
  # richness in PA
  r.PA <- r * P
  # plot(r.PA)
  
  q2 <- quantile(na.omit(values(r.PA)), probs = c(0.025, 0.5, 0.975))
  
  rich.mat[i,4] <- q2[1]
  rich.mat[i,5] <- q2[2]
  rich.mat[i,6] <- q2[3]
  
  # prop of group in PA
  # proportion of sum of richness in r.PA compared to sum of richness in group 
  rich.mat[i,7] <- round((sum(na.omit(values(r.PA)))/sum(na.omit(values(r))))*100, 1)/
    round((sum(na.omit(values(P)))/sum(na.omit(values(Grp.R))))*100, 1)
  # plot(P)
  print(paste("Finished calculating for NZSCC Group: ", i))
}

setwd(paste0(dir, "/Outputs"))
rich.mat <- round(rich.mat, 2)
# write.csv(rich.mat, file = "DF.rich.mat.csv")
# write.csv(rich.mat, file = "BI.rich.mat.csv")
# write.csv(rich.mat, file = "RF.rich.mat.csv")
write.csv(rich.mat, file = "MA.rich.mat.csv")

####==========    5. INTER AND INTRA GROUP PROP WITHIN SMAs    =============####
Rep.Intra.Inter <- (data.frame(matrix(ncol = 4,nrow = (75))))
colnames(Rep.Intra.Inter) <- c("All.PA.Intra","All.PA.prop.Intra","All.PA.Inter", "All.PA.prop.Inter"); 
rownames(Rep.Intra.Inter) <- paste("Grp_", seq(1:75))

# i = 20
for (i in 1:75){
  # group for reporting
  grp2 <- NZSCC; grp2[!grp2 == i] <- NA; grp2[grp2 == i] <- 1;
  grp <- raster(paste0("Sim_Inter/SimInter0_",i,".tif")); grp[!grp == 100] <- NA; grp[grp == 100] <- 1; grp <- grp * grp2
  # plot(grp)
  
  # prop of intra sim within PAs
  Intra <- raster(paste0("Intra/Intra_",i,".tif"))
  Intra[Intra == 0] <- NA; Intra <- Intra * grp
  # plot(Intra)
  
  # prop of inter sim within PAs
  Inter0 <- raster(paste0("Sim_Inter/SimInter0_",i,".tif"))
  grp.inv <- grp;grp.inv[is.na(grp.inv)] <-2; grp.inv[grp.inv == 1] <- NA; grp.inv[grp.inv == 2] <- 1; grp.inv <- grp.inv * mask
  Inter0 <- Inter0 * grp.inv
  # plot(Inter0)
  
  # PA within group for intra
  P <- grp * PA.list[[5]]
  P[P == 0] <- NA
  Intra.P <- Intra * P
  # plot(Intra.P)
  Rep.Intra.Inter[i,1] <- round((sum(na.omit(values(Intra.P)))/sum(na.omit(values(Intra))))*100, 1)
  Rep.Intra.Inter[i,2] <- round((sum(na.omit(values(Intra.P)))/sum(na.omit(values(Intra))))*100, 1)/
    round((sum(na.omit(values(P)))/sum(na.omit(values(grp))))*100, 1)
  
  # PA across area for inter
  PA <- PA.list[[5]]
  PA[PA == 0] <- NA
  Inter.P <- Inter0 * PA
  # plot(Inter.P)
  Rep.Intra.Inter[i,3] <- round((sum(na.omit(values(Inter.P)))/sum(na.omit(values(Inter0))))*100, 1)
  Rep.Intra.Inter[i,4] <- round((sum(na.omit(values(Inter.P)))/sum(na.omit(values(Inter0))))*100, 1)/
    round((sum(na.omit(values(PA)))/sum(na.omit(values(mask))))*100, 1)
  
  print(paste("Finished calculating for NZSCC Group: ", i))
}

setwd(paste0(dir, "/Outputs"))
write.csv(Rep.Intra.Inter, file = "Rep.Intra.Inter.csv")