#### To Start ####

library(ggplot2) # plots etc.
library(nlme) # lme()
library(MuMIn) # r.squaredGLMM()
library(ggfortify) # plotting PCAs
library(plyr) # revalue()
library(lawstat) # levene.test()
library(PMCMRplus) # kwAllPairsDunnTest()
library(forcats) # fct_relevel()
library(ggradar) # radar plots
library(dplyr) # general
library(geomorph) # morphometric analyses
library(Evomorph) # Getting deviation from consensus mean morphology


setwd("/Users/helen/Desktop/BaldwinGitHub/")

Theme <- theme(axis.line=element_line(colour="black", linewidth=0.5), panel.background=element_blank())
TypeCols <- scale_colour_manual(values=c("Anadromous"="red","Young Freshwater"="Green", "Old Freshwater"="Blue" ))
TypeFill <- scale_fill_manual(values=c("Anadromous"="red","Young Freshwater"="Green", "Old Freshwater"="Blue" ))
SalCols <- scale_colour_manual(values=c("0"="blue", "10"="green","20"="yellow", "30"="red", "40"="purple"))
SalFill <- scale_fill_manual(values=c("0"="blue", "10"="green","20"="yellow", "30"="red", "40"="purple"))
PopCols <- scale_colour_manual(values=c("A_S"="#990033", "A_A"="#FF3333","YF1_S"="#FFFF00", "YF_A"="#66FF00", "YF2_S"="#339900", "YF3_S"="#006600", "OF1_A"="#3399FF", "OF2_A"="#0000CC", "OF_S"="#330066"))
PopFill <- scale_fill_manual(values=c("A_S"="#990033", "A_A"="#FF3333", "YF1_S"="#FFFF00", "YF_A"="#66FF00", "YF2_S"="#339900", "YF3_S"="#006600", "OF1_A"="#3399FF", "OF2_A"="#0000CC", "OF_S"="#330066"))

#### Salinity Tolerance Breadth ####

#### Table 2 ####
# Hatching against salinity
Hatching <- read.table("./Clutch_HatchingSurvival.csv", header=TRUE, sep=",")
Hatching <- subset(Hatching, Salinity != 40)
HatchingModel <- lme(PercentSurvivalToHatch ~ poly(Salinity, 2)*PopulationType, random=~1|ClutchID, data=Hatching)
anova(HatchingModel)
summary(HatchingModel)

#### Figure 1 ####

## Plot 1.1 - survival against salinity, split by population type # colour se by pop

OverallSurvival <- read.table("./Tank_11MonthSurvival.csv", header=TRUE, sep=",")
OverallSurvival$PopulationType <- fct_relevel(OverallSurvival$PopulationType, c("Anadromous", "Young Freshwater", "Old Freshwater"))

AvgSurvivalPlot <- ggplot(data=OverallSurvival, mapping=aes(Salinity, StandardisedSurvival)) + geom_line(mapping=aes(Salinity, StandardisedSurvival, col=PopulationType, group=Population), stat="smooth", alpha=0.6, linetype=2) + geom_smooth(mapping=aes(Salinity, StandardisedSurvival, col=PopulationType, fill=PopulationType), se=T, alpha=0.1, linewidth=2) + labs(x="Salinity (ppt)", y="Survival (%)", col= "Population Type", fill="Population Type") + Theme + TypeCols + TypeFill + coord_cartesian(y=c(0,80))
AvgSurvivalPlot

## Plot 1.2 - plasticity against population

PopPlasticity <- read.table("./Population_PlasticityIndeces.csv", header=TRUE, sep=",")
PopPlasticity$Population <- factor(PopPlasticity$Population, levels=c("KB", "RS", "FK", "LB", "CR", "HW", "BB", "CL", "DM"))
PopPlasticity$Population <- revalue(PopPlasticity$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
PopPlasticity$Population <- fct_relevel(PopPlasticity$Population, c("A_A", "A_S", "YF1_S","YF_A", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))
PopPlasticity$Type <- fct_relevel(PopPlasticity$Type, c("Anadromous", "Young Freshwater", "Old Freshwater"))

SalinityToleranceBreadthBar <- ggplot(data=PopPlasticity, mapping=aes(Population, SalinityToleranceBreadth, fill=Type)) + geom_col(mapping=aes(Population, SalinityToleranceBreadth)) + labs(x= "Population", y="Cross-Salinity Tolerance (%)") + coord_cartesian(ylim=c(0,100)) + TypeFill+ Theme
SalinityToleranceBreadthBar


#### Physiological Plasticity ####

PopPlasticity <- read.table("./Population_PlasticityIndeces.csv", header=TRUE, sep=",")

PopPlasticity$Population <- factor(PopPlasticity$Population, levels=c("KB", "RS", "FK", "LB", "CR", "HW", "BB", "CL", "DM"))
PopPlasticity$Population <- revalue(PopPlasticity$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))

PCARawPI <- prcomp(PopPlasticity[, c(4:10)], center=T, scale=F)
PCARawPI # loadings for each term in each component
summary(PCARawPI) # proportion of variance for each component


#### Figure 2 ####

### Figure 2A ###

PlasticityPhysiologyPIRadar <- PopPlasticity %>% select(Population, FilamentWidthPI, InterLamellarWidthPI, LamellarWidthPI, PercentLamellarPI, IonocyteCoveragePI, IonocyteAreaPI, IonocyteRoundnessPI)

A_Radar <- subset(PlasticityPhysiologyPIRadar, Population=="A_A"|Population=="A_S") %>% ggradar(group.point.size=0, group.line.width=1, fill=T, fill.alpha=0.25) + PopCols + PopFill + labs(fill="Population", col="Population")
A_Radar
YF_Radar <- subset(PlasticityPhysiologyPIRadar, Population=="YF_A"|Population=="YF1_S"|Population=="YF2_S"|Population=="YF3_S") %>% ggradar(group.point.size=0, group.line.width=1, fill=T, fill.alpha=0.25) + PopCols + PopFill + labs(fill="Population", col="Population")
YF_Radar
OF_Radar <- subset(PlasticityPhysiologyPIRadar, Population=="OF1_A"|Population=="OF2_A"|Population=="OF_S") %>% ggradar(group.point.size=0, group.line.width=1, fill=T, fill.alpha=0.25) + PopCols + PopFill + labs(fill="Population", col="Population")
OF_Radar

### Figure 2B ###
# PCA Plot
autoplot(PCARawPI, scale=0, label=FALSE, loadings=TRUE, loadings.label=TRUE, loadings.colour="black", loadings.label.size=3, loadings.label.vjust = -1, loadings.label.colour="black", loadings.label.hjust = -0.01) + Theme + geom_point(aes(colour=PopPlasticity$Type, size=PopPlasticity$SalinityToleranceBreadth)) + TypeCols + scale_size(range=c(2,8))  + geom_text(label=PopPlasticity$Population, size=2,hjust=-0.7, vjust=1) + labs(x="PC1 (64.1%)", y="PC2 (21.8%)",col="Population", size="Population Salinity Tolerance (%)")


#### Trade-offs with growth and survival ####

# density-corrected growth in freshwater with pop salinity tolerance breadth
FreshwaterLength <- read.table("./Individ_FreshwaterLength.csv", header=TRUE, sep=",")

DensityCorrectSize <- lm(Length ~ AvgFishInTank, data=FreshwaterLength)
FreshwaterLength$DensityCorrectedLength <- resid(DensityCorrectSize)

aggregate(FreshwaterLength$DensityCorrectedLength, list(FreshwaterLength$Population), FUN=mean)

PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=TRUE, sep=",")
cor.test(PopFreshwaterMetrics$FreshwaterDensityCorrectedLength, PopFreshwaterMetrics$SalinityToleranceBreadth, method="spearman")

# freshwater survival with pop salinity tolerance breadth
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=TRUE, sep=",")
cor.test(PopFreshwaterMetrics$FreshwaterSurvival, PopFreshwaterMetrics$SalinityToleranceBreadth, method="spearman")

  # Missing fraction check
    cor.test(PopFreshwaterMetrics$FreshwaterDensityCorrectedLength, PopFreshwaterMetrics$FreshwaterSurvival, method="spearman")

# individual within-individual ionocyte area with density-corrected freshwater length
IndivFreshwaterMetrics <- read.table("./Individ_FreshwaterMetrics.csv", header=T, sep=",")
cor.test(IndivFreshwaterMetrics$DensityCorrectedLength, IndivFreshwaterMetrics$CoeffVarIonocyteArea)
# population average coeff of variation of within-individual ionocyte area in freshwater, population plasticity index of ionocyte area
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")
cor.test(PopFreshwaterMetrics$WithinIndivCoeffVarIonocyteArea, PopFreshwaterMetrics$IonocyteAreaPI, method="spearman", alternative="greater")

#### Processing Morphological Landmark Data ####

#read in all coordinates, turn 1st column (photoID) into row labels, exclude first 5 columns (classifiers) 
# Remember to invert X coordinates for R side fish first!
# Missing R-side fish: CR2c_7, DM0b_13, DM2c_1, KB3b_4, LB1a_5, RS0b_7 are not included
AllLandmarkCoords <- (as.matrix(read.table("./Individ_LRLandmarks.csv", header=TRUE, row.names=1, sep=",")[,-(1:6)]))

# read in first 5 columns as the classifiers, with row labels as 1st column, convert all to factors:
Classifiers <- (read.table("./Individ_LRLandmarks.csv", header=TRUE, row.names=1, sep=",")[,1:6])
Classifiers[sapply(Classifiers, is.character)] <- lapply(Classifiers[sapply(Classifiers, is.character)], as.factor)
summary(Classifiers)

# Convert 2D array into 3D array: (landmarks x number of landmarks x dimensions)
AllLandmarkCoords3D <- arrayspecs(AllLandmarkCoords, 24, 2)

# Check for missing coordinates
any(is.na(AllLandmarkCoords3D))
which(is.na(AllLandmarkCoords3D))

# Estimate missing landmarks using thin-plate spline interpolation
AllLandmarkCoords3D <- estimate.missing(AllLandmarkCoords3D, method="TPS")
any(is.na(AllLandmarkCoords3D))

# geomorph dataframe
FryCoords <- geomorph.data.frame(shape=AllLandmarkCoords3D, ind=Classifiers$FishID, tankID=Classifiers$TankID, side=Classifiers$Side, population=Classifiers$Population, salinity=Classifiers$Salinity, length=Classifiers$Length)

## Procrustes analysis with matching symmetry
BaldwinProcCoords <- bilat.symmetry(A=shape, ind=ind, side=side, object.sym = FALSE, RRPP = TRUE,iter = 499, data = FryCoords, print.progress = TRUE)
summary(BaldwinProcCoords)

# Symmetric component of shape variation
BaldwinProcCoordsInd <- BaldwinProcCoords$symm.shape

# Individual L-R asymmetry (fluctuating asymmetry, unsigned asymmetry index)
BaldwinLRAsymm <- BaldwinProcCoords$unsigned.AI

# Create a geomorph dataframe with symmetric shape Procrustes coordinates
MorphSamples <- read.table("./Individ_LandmarkSamples.csv",header=TRUE, sep=",")
MorphSamples[sapply(MorphSamples, is.character)] <- lapply(MorphSamples[sapply(MorphSamples, is.character)], as.factor)

BaldwinProcDataframeIndiv <- geomorph.data.frame(shape=BaldwinProcCoords$symm.shape, ind=MorphSamples$FishID, tank=MorphSamples$TankID, population=MorphSamples$Population, salinity=MorphSamples$Salinity, length=MorphSamples$Length)

# Run a factorial Procrustes MANOVA e.g. y ~ a * b, controlling for allometry:

MorphMANOVA <- procD.lm(shape~length + population*salinity, data=BaldwinProcDataframeIndiv)
summary(MorphMANOVA)

PCA <- gm.prcomp(BaldwinProcDataframeIndiv$shape)
summary(PCA)
plot(PCA)

PCAData <- data.frame(Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity, PC1=PCA$x[,1], PC2=PCA$x[,2])
PCAData$Salinity <- as.factor(PCAData$Salinity)
ggplot(PCAData, aes(PC1, PC2, col=Population, shape=Salinity)) + geom_point()


#### Table 3 ####
# Correlations between tolerance breadth and variation in physiological and morphological metrics in freshwater

### Inter-Individual Physiological Variation ###
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")

cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarFilamentWidth, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarLamellarWidth, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarInterLamellarWidth, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarIonocyteAre, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarIonocyteRoundness, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarIonocyteCoverage, method="spearman")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$BetweenIndivCoefVarIonocyteDistribution, method="spearman")

### Inter-Individual Morphological Variation ###
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$MorphologicalDisparity, method="spearman")

### Within-Individual Physiological Variation ###

PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$WithinIndivCoeffVarLamellarWidth, method="spearman", alternative="greater")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$WithinIndivCoeffVarInterLamellarWidth, method="spearman", alternative="greater")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$WithinIndivCoeffVarIonocyteArea, method="spearman", alternative="greater")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$WithinIndivCoeffVarIonocyteRoundness, method="spearman", alternative="greater")

### Left-Right Morphological Asymmetry ###
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")
cor.test(PopFreshwaterMetrics$SalinityToleranceBreadth, PopFreshwaterMetrics$AvgOfAsymmetry, method="spearman")


#### Figure 3 ####

### Figure 3a ###
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=TRUE, sep=",")
PopFreshwaterMetrics$Population <- revalue(PopFreshwaterMetrics$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
PopFreshwaterMetrics$Population <- fct_relevel(PopFreshwaterMetrics$Population, c("A_A", "A_S", "YF_A", "YF1_S", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))

LengthPlot <- ggplot(PopFreshwaterMetrics, aes(SalinityToleranceBreadth, FreshwaterDensityCorrectedLength))  + geom_smooth(method=lm, linetype=2, col="black") + geom_errorbar(aes(ymin=FreshwaterDensityCorrectedLength - StErrCorrectedLength, ymax=FreshwaterDensityCorrectedLength + StErrCorrectedLength)) + geom_point(aes(col=Population), size=5) + Theme + PopCols + labs(x="Population Salinity Tolerance Breadth (%)", y="Density-Corrected Length in Freshwater (mm)") + coord_cartesian(xlim=c(15,100))
LengthPlot

### Figure 3b ###
PopFreshwaterMetrics <- read.table("./Population_FreshwaterMetrics.csv", header=T, sep=",")

PopFreshwaterMetrics$Population <- revalue(PopFreshwaterMetrics$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
PopFreshwaterMetrics$Population <- fct_relevel(PopFreshwaterMetrics$Population, c("A_A", "A_S", "YF_A", "YF1_S", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))

AreaVarPlot <- ggplot(PopFreshwaterMetrics, aes(SalinityToleranceBreadth, WithinIndivCoeffVarIonocyteArea)) + geom_smooth(method=lm, linetype=2, col="black") + geom_errorbar(aes(ymin=WithinIndivCoeffVarIonocyteArea - StErrWithinIndivCoeffVarIonocyteArea, ymax=WithinIndivCoeffVarIonocyteArea + StErrWithinIndivCoeffVarIonocyteArea)) + geom_point(aes(col=Population), size=5) + Theme + PopCols + labs(x="Population Salinity Tolerance Breadth (%)", y="Coefficient of Variance of Ionocyte Area in Freshwater")
AreaVarPlot

### Figure 3c ###

ProcrustesVariances <- read.table("./PopSal_ProcrustesVariances.csv", header=T, sep=",")
ProcrustesVariances$Population <- factor(ProcrustesVariances$Population, levels=c("KB", "RS", "FK", "LB", "CR", "HW", "BB", "CL", "DM"))
ProcrustesVariances$Population <- revalue(ProcrustesVariances$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
ProcrustesVariances$Population <- fct_relevel(ProcrustesVariances$Population, c("A_A", "A_S", "YF1_S", "YF_A", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))

MorphVarPlot <- ggplot(ProcrustesVariances, aes(Salinity, MorphologicalDisparity, col=Population)) + geom_line(size=1) + labs(x="Salinity (ppt)", y="Inter-individual Morphological Variation (Procrustes Variance)", col="Population") + PopCols + Theme + geom_vline(aes(xintercept=30), linetype=2, alpha=0.4)
MorphVarPlot

#### Trade-offs with salinity tolerance breadth and change in morphological variation with salinity ####
PopPlasticity <- read.table("./Population_PlasticityIndeces.csv", header=T, sep=",")

# overall pop salinity tol breadth and plasticity index of inter-individ variation
cor.test(PopPlasticity$SalinityToleranceBreadth, PopPlasticity$MorphInterIndivVarPI, method="spearman")
# overall pop salinity tol breadth and plasticity index of within-individ variation
cor.test(PopPlasticity$SalinityToleranceBreadth, PopPlasticity$MorphWithinIndivVarPI, method="spearman")


#### Table 4 ####
# Across-population comparison of inter-individual morphological variation, and correlation of L-R asymmetry within individual deviation from average

### Pairwise differences in morphological disparity
MorphDisparity <- morphol.disparity(shape~length + population*salinity, groups=~population*salinity, data=BaldwinProcDataframeIndiv)
summary(MorphDisparity)

MorphDisparity$PV.dist.Pval
MorphDisparity$PV.dist

### Associations between individual L-R asymmetry and distance from group mean shape 

## Calculating within-individual variation (fluctuating asymmetry)

BaldwinLRAsymm <- BaldwinProcCoords$unsigned.AI
MorphSamples <- read.table("./Individ_LandmarkSamples.csv",header=TRUE, sep=",")
MorphSamples[sapply(MorphSamples, is.character)] <- lapply(MorphSamples[sapply(MorphSamples, is.character)], as.factor)

BaldwinProcDataframeAsymm <- data.frame(fluctasymm=BaldwinProcCoords$unsigned.AI, ind=MorphSamples$FishID, tank=MorphSamples$TankID, population=MorphSamples$Population, salinity=MorphSamples$Salinity, length=MorphSamples$Length)
BaldwinProcDataframeAsymm$salinity <- as.factor(BaldwinProcDataframeAsymm$salinity)

## Calculating Procrustes Distances between each individual and their population/salinity group mean

# This uses the earlier-calculated 'BaldwinProcDataframeIndiv' which has all symmetric shape component landmarks for all individuals
# And mshape(BaldwinProcDataframeIndiv$shape) to get the mean (consensus configuration) for each salinity/population

# Subset whole dataset by salinity/population groups
group <- factor(paste(BaldwinProcDataframeIndiv$population, BaldwinProcDataframeIndiv$salinity))
levels(group)
GroupMorphologies <- coords.subset(A = BaldwinProcDataframeIndiv$shape, group = group)
names(GroupMorphologies) # see the list levels
# group shape means
GroupMeanMorphologies <- lapply(GroupMorphologies, mshape)
# Calculate the Procrustes distance between each individual's coordinates and the mean shape for all groups:
PDLB0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`LB 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="LB" & Salinity=="0"))
PDLB1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`LB 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="LB" & Salinity=="1"))
PDLB2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`LB 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="LB" & Salinity=="2"))
PDLB3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`LB 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="LB" & Salinity=="3"))

PDBB0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`BB 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="BB" & Salinity=="0"))
PDBB1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`BB 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="BB" & Salinity=="1"))
PDBB2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`BB 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="BB" & Salinity=="2"))
PDBB3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`BB 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="BB" & Salinity=="3"))

PDCL0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CL 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CL" & Salinity=="0"))
PDCL1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CL 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CL" & Salinity=="1"))
PDCL2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CL 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CL" & Salinity=="2"))
PDCL3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CL 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CL" & Salinity=="3"))

PDRS0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`RS 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="RS" & Salinity=="0"))
PDRS1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`RS 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="RS" & Salinity=="1"))
PDRS2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`RS 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="RS" & Salinity=="2"))
PDRS3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`RS 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="RS" & Salinity=="3"))
PDRS4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`RS 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="RS" & Salinity=="4"))

PDKB0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`KB 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="KB" & Salinity=="0"))
PDKB1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`KB 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="KB" & Salinity=="1"))
PDKB2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`KB 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="KB" & Salinity=="2"))
PDKB3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`KB 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="KB" & Salinity=="3"))
PDKB4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`KB 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="KB" & Salinity=="4"))

PDCR0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CR 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CR" & Salinity=="0"))
PDCR1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CR 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CR" & Salinity=="1"))
PDCR2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CR 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CR" & Salinity=="2"))
PDCR3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CR 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CR" & Salinity=="3"))
PDCR4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`CR 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="CR" & Salinity=="4"))

PDDM0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`DM 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="DM" & Salinity=="0"))
PDDM1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`DM 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="DM" & Salinity=="1"))
PDDM2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`DM 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="DM" & Salinity=="2"))
PDDM3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`DM 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="DM" & Salinity=="3"))
PDDM4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`DM 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="DM" & Salinity=="4"))

PDFK0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`FK 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="FK" & Salinity=="0"))
PDFK1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`FK 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="FK" & Salinity=="1"))
PDFK2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`FK 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="FK" & Salinity=="2"))
PDFK3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`FK 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="FK" & Salinity=="3"))
PDFK4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`FK 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="FK" & Salinity=="4"))

PDHW0 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`HW 0`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="HW" & Salinity=="0"))
PDHW1 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`HW 1`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="HW" & Salinity=="1"))
PDHW2 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`HW 2`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="HW" & Salinity=="2"))
PDHW3 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`HW 3`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="HW" & Salinity=="3"))
PDHW4 <- (data.frame(ProcDist = ShapeDist(shapes=BaldwinProcDataframeIndiv$shape, reference=GroupMeanMorphologies$`HW 4`), Individual=BaldwinProcDataframeIndiv$ind, Population=BaldwinProcDataframeIndiv$population, Salinity=BaldwinProcDataframeIndiv$salinity) %>% subset(Population=="HW" & Salinity=="4"))

# Combine into a single dataset
ProcrustesDistancesAll <- bind_rows(PDLB0, PDLB1, PDLB2, PDLB3, PDBB0, PDBB1, PDBB2, PDBB3, PDCL0, PDCL1, PDCL2, PDCL3, PDRS0, PDRS1, PDRS2, PDRS3, PDRS4, PDKB0, PDKB1, PDKB2, PDKB3, PDKB4, PDFK0, PDFK1, PDFK2, PDFK3, PDFK4,PDCR0, PDCR1, PDCR2, PDCR3, PDCR4, PDHW0, PDHW1, PDHW2, PDHW3, PDHW4, PDDM0, PDDM1, PDDM2, PDDM3, PDDM4)

ggplot(ProcrustesDistancesAll, aes(Salinity,ProcDist, col=Population)) + geom_point() + geom_smooth() + Theme

# Calculate mean distances from group average for each population/salinity combination (just to check)
aggregate(ProcrustesDistancesAll$ProcDist, list(ProcrustesDistancesAll$Population, ProcrustesDistancesAll$Salinity), FUN=mean)

# Bring up the individual L-R asymmetry
BaldwinProcDataframeAsymm <- data.frame(fluctasymm=BaldwinProcCoords$unsigned.AI, Individual=MorphSamples$FishID, tank=MorphSamples$TankID, population=MorphSamples$Population, salinity=MorphSamples$Salinity, length=MorphSamples$Length)
BaldwinProcDataframeAsymm$salinity <- as.factor(BaldwinProcDataframeAsymm$salinity)
# Combine individual L-R symmetry with individual procrustes distance from mean morphology
BaldwinIndivVariationMetrics <- full_join(ProcrustesDistancesAll, BaldwinProcDataframeAsymm, by="Individual")
  # write.csv(BaldwinIndivVariationMetrics, file="Individ_MorphologyMetrics.csv")

# And *finally* run the correlation test!
cor.test(BaldwinIndivVariationMetrics$ProcDist, BaldwinIndivVariationMetrics$fluctasymm)

# For all populations (use cor_test for pipelining, I am too lazy)
RSSub <- subset(BaldwinIndivVariationMetrics, population=="RS") # A_A
cor.test(RSSub$ProcDist,RSSub$fluctasymm)
KBSub <- subset(BaldwinIndivVariationMetrics, population=="KB") # A_S
cor.test(KBSub$ProcDist,KBSub$fluctasymm)

FKSub <- subset(BaldwinIndivVariationMetrics, population=="FK") # YF1_S
cor.test(FKSub$ProcDist,FKSub$fluctasymm)
LBSub <- subset(BaldwinIndivVariationMetrics, population=="LB") # YF_A
cor.test(LBSub$ProcDist,LBSub$fluctasymm)
CRSub <- subset(BaldwinIndivVariationMetrics, population=="CR") # YF2_S
cor.test(CRSub$ProcDist,CRSub$fluctasymm)
HWSub <- subset(BaldwinIndivVariationMetrics, population=="HW") # YF3_S
cor.test(HWSub$ProcDist,HWSub$fluctasymm)

BBSub <- subset(BaldwinIndivVariationMetrics, population=="BB") # OF1_A
cor.test(BBSub$ProcDist,BBSub$fluctasymm)
CLSub <- subset(BaldwinIndivVariationMetrics, population=="CL") # OF2_A
cor.test(CLSub$ProcDist,CLSub$fluctasymm)
DMSub <- subset(BaldwinIndivVariationMetrics, population=="DM") # OF_S
cor.test(DMSub$ProcDist,DMSub$fluctasymm)


#### Extended Data Table 4 ####
# PCA loadings of variation in plasticity indeces of physiological metrics between populations (from Fig 2)
PopPlasticity <- read.table("./Population_PlasticityIndeces.csv", header=TRUE, sep=",")

PopPlasticity$Population <- factor(PopPlasticity$Population, levels=c("KB", "RS", "FK", "LB", "CR", "HW", "BB", "CL", "DM"))
PopPlasticity$Population <- revalue(PopPlasticity$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))

PCARawPI <- prcomp(PopPlasticity[, c(4:10)], center=T, scale=F)
PCARawPI # loadings for each term in each component
summary(PCARawPI) # proportion of variance for each component



#### Appendix ####

### Figure A1 ###

ClutchCoV <- read.table("./Clutch_DailySurvivalCoefVar.csv", header=TRUE, sep=",")
ClutchCoV$Population <- revalue(ClutchCoV$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
ClutchCoV$Population <- fct_relevel(ClutchCoV$Population, c("A_A", "A_S", "YF_A", "YF1_S", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))

ggplot(ClutchCoV, aes(Day,CoV, col=Population, fill=Population))+ geom_smooth(alpha=0.2) + labs(x="Day Since Fertilisation", y="Between-treatment Coefficient of Variation in Survival") + Theme + PopCols + PopFill + coord_cartesian(ylim=c(0,0.5)) + geom_vline(aes(xintercept=21.7), linetype=2)

### Figure A2 ###

TankCoV <- read.table("./Tank_MonthlySurvivalCoefVar.csv", header=TRUE, sep=",")
TankCoV$Population <- revalue(TankCoV$Population, c("KB"="A_S", "RS"="A_A", "FK"="YF1_S", "LB"="YF_A", "CR"="YF2_S", "HW"="YF3_S", "BB"="OF1_A", "CL"="OF2_A", "DM"="OF_S"))
TankCoV$Population <- fct_relevel(TankCoV$Population, c("A_A", "A_S", "YF_A", "YF1_S", "YF2_S", "YF3_S", "OF1_A", "OF2_A", "OF_S" ))

ggplot(TankCoV, aes(MeasureMonth,CoV, col=Population, fill=Population))+ geom_smooth(alpha=0.2) + labs(x="Month", y="Between-treatment Coefficient of Variation in Survival") + Theme + coord_cartesian(xlim=c(0,11), ylim=c(0,0.6)) + PopCols + PopFill

