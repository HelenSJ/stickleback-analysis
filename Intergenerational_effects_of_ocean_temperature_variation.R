#### Intergenerational effects of ocean temperature variation: early life benefits are short-lived in threespine stickleback ####
#### Spence-Jones, Pein, Shama 

#### To Start ####

setwd()

library(nlme) # lme()
library(PMCMRplus) # kwAllPairsDunnTest()
library(lawstat) # Levene's test
library(ggsurvfit) # Kaplan-Meyer survival functions
library(survival) # Cox survival analysis
library(MASS) # Linear Discriminant analysis (Morphometrics data)
library(geomorph) # Morphometric analyses
library(plyr) # revalue() function
library(dplyr) # Always useful for wrangling datasets
library(forcats) # fct_relevel() function


library(ggplot2) # Plots
Theme <- theme(axis.line=element_line(colour="black", linewidth=0.5), panel.background=element_blank())
# revalue(X$Y,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))
# fct_relevel(X$Y, c("Constant", "Natural Variation", "Increased Variation"))
IncVarCols <- scale_colour_manual(values=c("Constant"="#104E8B", "Natural Variation"="#008B00","Increased Variation"="#FFB90F"), drop=TRUE)
IncVarFill <- scale_fill_manual(values=c("Constant"="#104E8B", "Natural Variation"="#008B00","Increased Variation"="#FFB90F"), drop=TRUE)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#### Fecundity ####

# Adult Length Comparison #

Adults <- read.table("./AdultLength.csv", header=TRUE, sep=",")
summary(aov(Length.mm.~Sex*TreatmentGroup, Adults))
AdultMale <- subset(Adults, Sex=="M")
mean(AdultMale$Length.mm.)
sd(AdultMale$Length.mm.)
kruskal.test(AdultMale$Length.mm.~AdultMale$TreatmentGroup)
AdultFemale <- subset(Adults, Sex=="F")
mean(AdultFemale$Length.mm.)
sd(AdultFemale$Length.mm.)
kruskal.test(AdultFemale$Length.mm.~AdultFemale$TreatmentGroup)

ClutchTiming <- read.table("./ClutchTiming.csv", header=TRUE, sep=",")
ggplot(ClutchTiming, aes(Day, col=TreatmentGroup)) + stat_ecdf(geom="step")
TimingModel <- lme(Day~TreatmentGroup + Length.mm., random=~1|DamID, data=ClutchTiming)
summary(TimingModel)
kruskal.test(ClutchTiming$Day~ClutchTiming$TreatmentGroup)


# Clutch Size #

Clutches <- read.table("./ClutchSize.csv", header=TRUE, sep=",")
Clutches$TreatmentGroup <- as.factor(Clutches$TreatmentGroup)
hist(Clutches$FertilisationPercent, breaks=10)
hist(Clutches$ClutchSize)
kruskal.test(Clutches$FertilisationPercent~Clutches$TreatmentGroup)

ggplot(Clutches, aes(Length.mm., ClutchSize, col=FertilisationPercent)) + geom_point(size=5) +Theme

# subset to only clutches with >10% fertilisation
Clutches <- subset(Clutches, FertilisationPercent > 10)
shapiro.test(Clutches$ClutchSize) # normally distributed
summary(Clutches$ClutchSize)
sd(Clutches$ClutchSize)
ClutchSizeModel <- lme(ClutchSize~TreatmentGroup*Length.mm., random=~1|DamID, data=Clutches)
anova(ClutchSizeModel)
summary(ClutchSizeModel)

ggplot(Clutches, aes(Length.mm., ClutchSize, col=TreatmentGroup)) + geom_point() + geom_smooth(se=FALSE, method=lm) + Theme

kruskal.test(Clutches$ClutchSize~Clutches$TreatmentGroup)


# Time to hatching #

HatchTime <- read.table("./HatchTiming.csv", header=TRUE, sep=",")

kruskal.test(HatchTime$Length~HatchTime$OffspringTreatment)
kruskal.test(HatchTime$Length~HatchTime$ParentTreatment)

HatchTime$OffspringTreatment <- as.factor(HatchTime$OffspringTreatment)
HatchTime$OffspringTreatment <- relevel(HatchTime$OffspringTreatment, "Var")

HatchLengthModel <- lme(Length~OffspringTreatment*ParentTreatment, random=~1|ClutchID, data=HatchTime, na.action=na.omit)
summary(HatchLengthModel)
anova(HatchLengthModel)

# Egg Size #

EggSize <- read.table("./EggSize.csv", header=TRUE, sep=",")
EggSize$ClutchID <- as.factor(EggSize$ClutchID)
EggSize$TreatmentGroup <- as.factor(EggSize$TreatmentGroup)
EggSize <- subset(EggSize, ClutchID !="PV005" | ClutchID !="PC055") # Issues with measuring for PV005; PC055 only had 8 eggs measured out of 30
EggSize <- subset(EggSize, ClutchSize >=60 & ClutchSize <=100) # Excluding unusually small and unusually large clutches which are not represented across all treatments
EggSize <- subset(EggSize, FertilisationPercent >10)

EggSizeModel <- lme(Length.um.~TreatmentGroup*ClutchSize, data=EggSize, random=~1|DamID/ClutchID)
anova(EggSizeModel)
summary(EggSizeModel)

mean(EggSize$Length.um.)
sd(EggSize$Length.um.)

EggVolumeModel <- lme(Volume.mm3.~TreatmentGroup*ClutchSize, data=EggSize, random=~1|DamID/ClutchID)
anova(EggVolumeModel)
summary(EggVolumeModel)

mean(EggSize$Volume.mm3.)
sd(EggSize$Volume.mm3.)

EggSizeModel <- lme(Length.um.~TreatmentGroup, data=EggSize, random=~1|DamID/ClutchID)
anova(EggSizeModel)
summary(EggSizeModel)

shapiro.test(EggSize$Length.um.) # non-normally distributed
kruskal.test(EggSize$Length.um. ~ EggSize$TreatmentGroup)
kwAllPairsDunnTest(Length.um.~TreatmentGroup, data=EggSize)
ggplot(EggSize, aes(Length.um., col=TreatmentGroup)) + geom_freqpoly() + Theme

ggplot(EggSize, aes(ClutchSize, Length.um., col=TreatmentGroup)) + geom_point() + geom_smooth(method=lm) + Theme

Clutches <- read.table("./Clutches.csv", header=TRUE, sep=",")
Clutches <- subset(Clutches, ClutchID !="PV005" | ClutchID != "PC055") # Issues with measuring for PV005
Clutches$TreatmentGroup <- as.factor(Clutches$TreatmentGroup)
ClutchesEggSize <- subset(Clutches, ClutchSize >=60 & ClutchSize<=100)

kruskal.test(ClutchesEggSize$AvgEggLength~ClutchesEggSize$TreatmentGroup)

# Variation Within/Between Clutches #

Clutches <- read.table("./Clutches.csv", header=TRUE, sep=",")
Clutches <- subset(Clutches, ClutchID !="PV005"| ClutchID != "PC055") # Issues with measuring for PV005
Clutches$TreatmentGroup <- as.factor(Clutches$TreatmentGroup)
Clutches <- subset(Clutches, FertilisationPercent >10)

# Variation in size within clutch
kruskal.test(Clutches$CoefVar~Clutches$TreatmentGroup)
kruskal.test(Clutches$EggStDev~Clutches$TreatmentGroup)

# Variation in size between clutches
levene.test(Clutches$AvgEggLength, Clutches$TreatmentGroup, location="mean")
levene.test(Clutches$AvgEggLength, Clutches$Treatment, location="median")


# Clutch Survival #

ClutchSurvival <- read.table("./ClutchSurvival.csv", header=TRUE, sep=",")
ClutchSurvival$OffspringTreatment <- as.factor(ClutchSurvival$OffspringTreatment)
ClutchSurvival$ParentTreatment <- as.factor(ClutchSurvival$ParentTreatment)
ClutchSurvival$ComboTreatment <- as.factor(ClutchSurvival$ComboTreatment)
ClutchSurvival <- subset(ClutchSurvival, FertilisationPercent >10)
shapiro.test(ClutchSurvival$HatchPercent) # non-normal
hist(ClutchSurvival$HatchPercent)
kruskal.test(ClutchSurvival$HatchPercent~ClutchSurvival$ComboTreatment)

ClutchSurvivalModel <- lme(HatchPercent~OffspringTreatment*ParentTreatment, data=ClutchSurvival, random=~1|DamID/ClutchID, na.action=na.omit)
anova(ClutchSurvivalModel)
summary(ClutchSurvivalModel)


#### Offspring Survival and Growth ####

# Post-hatch Survival #

KMSurvival <- read.table("./KMFrySurvival.csv", header=TRUE, sep=",")
FrySurv <- Surv(KMSurvival$Time, KMSurvival$Status)
ggsurvfit(survfit2(Surv(Time, Status) ~ 1, data=KMSurvival))
survdiff(Surv(Time, Status) ~ ComboTreatment, KMSurvival)
coxph(Surv(Time, Status)~OffspringTreatment*ParentTreatment, KMSurvival)

FrySurvival <- read.table("./FrySurvival.csv", header=TRUE, sep=",")
FrySurvival90 <- subset(FrySurvival, Day=="90")
mean(FrySurvival90$PercentSurvival)
sd(FrySurvival90$PercentSurvival)

# Growth #

Growth <- read.table("./FryLengths.csv", head=TRUE, sep=",")
Growth <- subset(Growth, Fish>5)
Growth$BatchID <- as.factor(Growth$BatchID)

ggplot(Growth, aes(Day, Length, col=ComboTreatment)) + geom_jitter(height=0, width=5) + geom_smooth(se=FALSE) + Theme

# 30d Growth #

Growth30 <- subset(Growth, Day==30)
levene.test(Growth30$Length,Growth30$ComboTreatment)
print(games_howell_test(Growth30, Length~ComboTreatment), n=36)
ggplot(Growth30, aes(Fish, Length)) + geom_point() + geom_smooth() + Theme
Growth30Model <- lme(Length~Fish + ParentTreatment*OffspringTreatment, random=~1|DamID/ClutchID, data=Growth30)
summary(Growth30Model)

# 60d Growth #

Growth60 <- subset(Growth, Day==60)
ggplot(Growth60, aes(Fish, Length)) + geom_point() + geom_smooth() + Theme
Growth60Model <- lme(Length~Fish + ParentTreatment*OffspringTreatment, random=~1|DamID/ClutchID, data=Growth60)
summary(Growth60Model)

# 90 Growth #

Growth90 <- subset(Growth, Day==90)
summary(Growth90$BatchID)
ggplot(Growth90, aes(Fish, Length)) + geom_point() + geom_smooth() + Theme
Growth90Model <- lme(Length~Fish + ParentTreatment*OffspringTreatment, random=~1|DamID/ClutchID, data=Growth90)
summary(Growth90Model)
anova(Growth90Model)
ranef(Growth90Model)


# Length Variation #

AvgGrowth <- read.table("./FryAvgGrowth.csv", header=TRUE, sep=",")
AvgGrowth <- subset(AvgGrowth, Day<100 & N_start > 5)
AvgGrowth <- subset(AvgGrowth, BatchID!="C028V")

ggplot(AvgGrowth, aes(Day, CoefVar, colour=ParentTreatment)) + geom_jitter(height=0, width=1) + geom_smooth() + Theme

LengthVariationModel <- lme(CoefVar~ParentTreatment*OffspringTreatment + Day + N_start, data=AvgGrowth, random=~1|ClutchID, na.action=na.omit)
anova(LengthVariationModel)
summary(LengthVariationModel)

AvgGrowth30 <- subset(AvgGrowth, Day==30)
kruskal.test(AvgGrowth30$CoefVar~AvgGrowth30$ComboTreatment)

AvgGrowth60 <- subset(AvgGrowth, Day==60)
kruskal.test(AvgGrowth60$CoefVar~AvgGrowth60$ComboTreatment)

AvgGrowth90 <- subset(AvgGrowth, Day==90)
kruskal.test(AvgGrowth90$CoefVar~AvgGrowth90$ComboTreatment)


#### Offspring Morphology ####

#read in all coordinates, turn 1st column (photoID) into row labels, exclude first 5 columns (classifiers) 
# Remember to invert X coordinates for R side fish first!
AllLandmarkCoords <- (as.matrix(read.table("./AllMorphologyLandmarksLR.csv", header=TRUE, row.names=1, sep=",")[,-(1:6)]))
# read in first 6 columns as the classifiers, with row labels as 1st column, convert all to factors:
Classifiers <- (read.table("./AllMorphologyLandmarksLR.csv", header=TRUE, row.names=1, sep=",")[,1:6])
Classifiers[sapply(Classifiers, is.character)] <- lapply(Classifiers[sapply(Classifiers, is.character)], as.factor)
summary(Classifiers)

# Convert 2D array into 3D array: (landmarks x number of landmarks x dimensions)
AllLandmarkCoords3D <- arrayspecs(AllLandmarkCoords, 22, 2)

# Check for missing coordinates
any(is.na(AllLandmarkCoords3D))

# Estimate missing landmarks on three specimens (eyes (landmarks 7, 8, 9) on 211019_I049C_05L.jpg and 211019_I049C_07R.jpg, spine (landmark 22) on 211110_V063V_01R.jpg) using thin-plate spline interpolation
AllLandmarkCoords3D <- estimate.missing(AllLandmarkCoords3D, method="TPS")
any(is.na(AllLandmarkCoords3D))

# geomorph dataframe
IncVarCoords <- geomorph.data.frame(shape=AllLandmarkCoords3D, ind=Classifiers$FishID, side=Classifiers$Side, tank=Classifiers$BatchID, parent=Classifiers$ParentTreatment, offspring=Classifiers$OffspringTreatment, length=Classifiers$Length)


#Procrustes transformation
IncVarProc <- gpagen(IncVarCoords$shape)
plotAllSpecimens(IncVarProc$coords)
IncVarProcDataframe <- geomorph.data.frame(shape=IncVarProc$coords, ind=IncVarCoords$ind, side=IncVarCoords$side, tank=IncVarCoords$length, parent=IncVarCoords$parent, offspring=IncVarCoords$offspring, length=IncVarCoords$length)

IncVarShapeModel <- procD.lm(shape~length + parent*offspring, data=IncVarProcDataframe, iter=499)
summary(IncVarShapeModel)


## Procrustes analysis with matching symmetry
IncVarProcCoords <- bilat.symmetry(A=shape, ind=ind, side=side, object.sym = FALSE, RRPP = TRUE,iter = 499, data = IncVarCoords, print.progress = TRUE)
summary(IncVarProcCoords)

# Symmetric component of shape variation
IncVarProcCoordsInd <- IncVarProcCoords$symm.shape
plotAllSpecimens(IncVarProcCoordsInd)

#Asymmetric component of shape variation
IncVarProcCoordsAsym <- IncVarProcCoords$asymm.shape
plotAllSpecimens(IncVarProcCoordsAsym)

# Create a geomorph dataframe with symmetric shape Procrustes coordinates
MorphSamples <- read.table("./AllMorphology_SampleDetails.csv",header=TRUE, sep=",")
MorphSamples[sapply(MorphSamples, is.character)] <- lapply(MorphSamples[sapply(MorphSamples, is.character)], as.factor)

IncVarProcDataframeIndiv <- geomorph.data.frame(shape=IncVarProcCoords$symm.shape, ind=MorphSamples$FishID, tank=MorphSamples$BatchID, parent=MorphSamples$ParentTreatment, offspring=MorphSamples$OffspringTreatment, length=MorphSamples$Length)

# Run a factorial Procrustes MANOVA e.g. y ~ a * b, controlling for allometry:

MorphMANOVA <- procD.lm(shape~length + offspring*parent, data=IncVarProcDataframeIndiv)
summary(MorphMANOVA)

PCA <- gm.prcomp(IncVarProcDataframeIndiv$shape)
summary(PCA)
plot(PCA)

PCAData <- data.frame(ParentTreatment=IncVarProcDataframeIndiv$parent, OffspringTreatment=IncVarProcDataframeIndiv$offspring, PC1=PCA$x[,1], PC2=PCA$x[,2])

ggplot(PCAData, aes(PC1, PC2, col=ParentTreatment, shape=OffspringTreatment)) + geom_point()


### LDA and jackknifed assignment accuracy
# Put everything into a nice table:
IncVarProcDataframeIndiv <- geomorph.data.frame(shape=IncVarProcCoords$symm.shape, ind=MorphSamples$FishID, tank=MorphSamples$BatchID, parent=MorphSamples$ParentTreatment, offspring=MorphSamples$OffspringTreatment, length=MorphSamples$Length)

IncVarLDATable <- data.frame(two.d.array(IncVarProcCoords$symm.shape))
IncVarLDATable <- tibble::rownames_to_column(IncVarLDATable, "Individual")
IncVarLDATable$OffspringTreatment <- substr(IncVarLDATable$Individual, 4,4) 
IncVarLDATable$ParentTreatment <- substr(IncVarLDATable$Individual, 8,8) 
IncVarLDATable$Combo <- paste(IncVarLDATable$ParentTreatment, IncVarLDATable$OffspringTreatment, sep="")
IncVarLDATable <- IncVarLDATable %>% select(-ParentTreatment)
IncVarLDATable <- IncVarLDATable %>% select(-OffspringTreatment)
IncVarLDATable <- IncVarLDATable %>% select(-Individual)

# Jackknifed assignment accuracies
IncVarLDA <- lda(Combo ~., CV=TRUE, data=IncVarLDATable)
IncVarLDA
diag(prop.table(table(IncVarLDATable$Combo, IncVarLDA$class),1))
sum(diag(prop.table(table(IncVarLDATable$Combo, IncVarLDA$class))))


# Looking at between-individual variation
MorphDisparity <- morphol.disparity(shape~length + parent*offspring, groups=~parent*offspring, data=IncVarProcDataframeIndiv)
summary(MorphDisparity)


ProcrustesVariances <- data.frame(TreatmentGroup=c("Con.Con","Con.Inc","Con.Var","Inc.Con","Inc.Inc", "Inc.Var", "Var.Con", "Var.Inc", "Var.Var"), ProcVar=c(0.0008859612, 0.0008149803, 0.0008472204, 0.0006822980, 0.0006916555, 0.0006372386, 0.0007628143, 0.0008848111, 0.0008727115))
ggplot(ProcrustesVariances, aes(TreatmentGroup, ProcVar)) + geom_col() + labs(x="Treatment Group (Parent:Offspring)", y="Procrustes Variance")

## Now on to fluctuating asymmetry between groups 
FAData <- data.frame(FA=IncVarProcCoords$unsigned.AI, ind=MorphSamples$FishID, tank=MorphSamples$BatchID, parent=MorphSamples$ParentTreatment, offspring=MorphSamples$OffspringTreatment, length=MorphSamples$Length)
FAModel <- lm(FA~length + parent*offspring, FAData)
summary(FAModel)
anova(FAModel)

ggplot(FAData, aes(offspring, FA, fill=parent)) + geom_boxplot()


#### Plots ####

### Figure 1 - Density-corrected Growth

Growth <- read.table("./FryLengths.csv", head=TRUE, sep=",")
Growth <- subset(Growth, Fish>5)
Growth$BatchID <- as.factor(Growth$BatchID)
Growth$ParentTreatment <- revalue(Growth$ParentTreatment,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))
Growth$OffspringTreatment <- revalue(Growth$OffspringTreatment,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))
Growth$ParentTreatment <- fct_relevel(Growth$ParentTreatment, c("Constant", "Natural Variation", "Increased Variation"))
Growth$OffspringTreatment <- fct_relevel(Growth$OffspringTreatment, c("Constant", "Natural Variation", "Increased Variation"))

Growth30 <- subset(Growth, Day==30)
Growth30DensityModel <- lme(Length~Fish, random=~1|DamID/ClutchID, data=Growth30)
Growth30$LengthResidual <- resid(Growth30DensityModel)

Growth30Plot <- ggplot(Growth30, aes(OffspringTreatment, LengthResidual)) + geom_boxplot(aes(fill=ParentTreatment)) + Theme + IncVarFill + theme(legend.position=c(0.5,0.9)) + labs(x= "Offspring Temperature Treatment", y="30-Day Length Residual (mm)", fill="Parent Temperature Treatment")


Growth60 <- subset(Growth, Day==60)
Growth60DensityModel <- lme(Length~Fish, data=Growth60,  random=~1|DamID/ClutchID)
Growth60$LengthResidual <- resid(Growth60DensityModel)

Growth60Plot <- ggplot(Growth60, aes(OffspringTreatment, LengthResidual)) + geom_boxplot(aes(fill=ParentTreatment)) + Theme + IncVarFill + guides(col='none', fill='none') + labs(x= "Offspring Temperature Treatment", y="60-Day Length Residual (mm)", fill="Parent Temperature Treatment")
Growth60Plot

Growth90 <- subset(Growth, Day==90)
Growth90DensityModel <- lme(Length~Fish, data=Growth90,  random=~1|DamID/ClutchID)
Growth90$LengthResidual <- resid(Growth90DensityModel)

Growth90Plot <- ggplot(Growth90, aes(OffspringTreatment, LengthResidual)) + geom_boxplot(aes(fill=ParentTreatment)) + Theme + IncVarFill + guides(col='none', fill='none') + labs(x= "Offspring Temperature Treatment", y="90-Day Length Residual (mm)", fill="Parent Temperature Treatment")

multiplot(cols=1, Growth30Plot, Growth60Plot, Growth90Plot)

### Figure 2 - Morphological Variation

## PCA Plot
PCAData <- data.frame(ParentTreatment=IncVarProcDataframeIndiv$parent, OffspringTreatment=IncVarProcDataframeIndiv$offspring, PC1=PCA$x[,1], PC2=PCA$x[,2])

ggplot(PCAData, aes(PC1, PC2, col=ParentTreatment, shape=OffspringTreatment)) + geom_point()

PCAData$ParentTreatment <- revalue(PCAData$ParentTreatment,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))
PCAData$OffspringTreatment <- revalue(PCAData$OffspringTreatment,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))
PCAData$ParentTreatment <- fct_relevel(PCAData$ParentTreatment, c("Constant", "Natural Variation", "Increased Variation"))
PCAData$OffspringTreatment <- fct_relevel(PCAData$OffspringTreatment, c("Constant", "Natural Variation", "Increased Variation"))

IncVarPCAPlot <- ggplot(aes(PC1, PC2), data=PCAData) + stat_ellipse(geom="polygon", aes(fill=OffspringTreatment), alpha=0.3, level=0.95) + stat_ellipse(geom="polygon", aes(linetype=ParentTreatment), col="black", alpha=0, level=0.95) +  geom_point(aes(PC1, PC2, col=OffspringTreatment, shape=ParentTreatment), size=1) + Theme  + labs(x="PC1 (23.8%)", y="PC2 (17.6%)", shape="Parent Treatment", linetype="Parent Treatment", col="Offspring Treatment", fill="Offspring Treatment") + IncVarFill + IncVarCols 
IncVarPCAPlot


## PC deformation plots

IncVarProcCoords$symm.shape

IVProWire <- mshape(IncVarProcCoords$symm.shape)
plotRefToTarget(PCA$shapes$shapes.comp1$min, IVProWire)
plotRefToTarget(PCA$shapes$shapes.comp1$max, IVProWire)

PC1Lolli <- plotRefToTarget(PCA$shapes$shapes.comp1$min, PCA$shapes$shapes.comp1$max, method = "vector", mag = 2, gridPars=gridPar(pt.bg = "green", pt.size = 1))
PC2Lolli <- plotRefToTarget(PCA$shapes$shapes.comp2$min, PCA$shapes$shapes.comp2$max, method = "vector", mag = 2, gridPars=gridPar(pt.bg = "green", pt.size = 1))



### Figure 3 - Morphological variation

ProcrustesVariances <- data.frame(ParentTreatment=c("Constant","Constant","Constant","Increased Variation","Increased Variation", "Increased Variation", "Natural Variation", "Natural Variation", "Natural Variation"), OffspringTreatment=c("Constant","Increased Variation","Natural Variation","Constant","Increased Variation", "Natural Variation", "Constant", "Increased Variation", "Natural Variation"), ProcVar=c(0.0008859612, 0.0008149803, 0.0008472204, 0.0006822980, 0.0006916555, 0.0006372386, 0.0007628143, 0.0008848111, 0.0008727115))
ProcrustesVariances$ParentTreatment <- fct_relevel(ProcrustesVariances$ParentTreatment, c("Constant", "Natural Variation", "Increased Variation"))
ProcrustesVariances$OffspringTreatment <- fct_relevel(ProcrustesVariances$OffspringTreatment, c("Constant", "Natural Variation", "Increased Variation"))

ggplot(data=ProcrustesVariances, mapping=aes(OffspringTreatment, ProcVar, fill=ParentTreatment)) + geom_col(position="dodge") + Theme + labs(x="Offspring Treatment", y="Procrustes Variance", fill="Parental Treatment") + IncVarFill


# Figure S1 - experimental plan

Experiment <- read.table("./TemperatureTimeline.csv", header=TRUE, sep=",")
Experiment$Treatment <- revalue(Experiment$Treatment,c("Con"="Constant", "Inc"="Increased Variation", "Var"="Natural Variation"))

ggplot(Experiment, aes(Day, AvgOfTemperature, col=Treatment)) + geom_line() + geom_vline(xintercept=24, linetype=2) + geom_vline(xintercept=81, linetype=2) + Theme + labs(y="Temperature (\u00B0C)", col="Temperature Treatment") + theme(legend.position=c(0.85,0.15)) + scale_y_continuous(breaks = seq(14, 22, len = 17))

# Figure S2 - density-dependence of growth

Growth <- read.table("./FryLengths.csv", head=TRUE, sep=",")
Growth$Day <- as.factor(Growth$Day)
Growth6 <- subset(Growth, Growth$Fish<6)
GrowthHigh <- subset(Growth, Growth$Fish>5)

ggplot(Growth, aes(Fish, Length, col=Day)) + geom_point() + geom_smooth(data=Growth6, aes(Fish, Length, col=Day), method=lm) + geom_smooth(data=GrowthHigh,aes(Fish, Length, col=Day), method=lm)+ geom_vline(xintercept=5.5, linetype=2) + Theme + labs(x="Number of Fish in Tank", y="Length (mm)", col="Fish Age (days)")


### Supplementary figure - Shape change between treatment groups

IncVarProcDataframeIndiv <- geomorph.data.frame(shape=IncVarProcCoords$symm.shape, ind=MorphSamples$FishID, tank=MorphSamples$BatchID, parent=MorphSamples$ParentTreatment, offspring=MorphSamples$OffspringTreatment, length=MorphSamples$Length)

# Subset whole dataset by parent/offspring groups
group <- factor(paste(IncVarProcDataframeIndiv$parent, IncVarProcDataframeIndiv$offspring))
levels(group)
GroupMorphologies <- coords.subset(A = IncVarProcDataframeIndiv$shape, group = group)

# group shape means
GroupMeanMorphologies <- lapply(GroupMorphologies, mshape)

CCWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Con Con`, method="TPS", mag=6)
CVWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Con Var`, method="TPS", mag=6)
CIWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Con Inc`, method="TPS", mag=6)

VCWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Var Con`, method="TPS", mag=6)
VVWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Var Var`, method="TPS", mag=6)
VIWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Var Inc`, method="TPS", mag=6)

ICWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Inc Con`, method="TPS", mag=6)
IVWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Inc Var`, method="TPS", mag=6)
IIWireframe <- plotRefToTarget(GroupMeanMorphologies$`Con Con`, GroupMeanMorphologies$`Inc Inc`, method="TPS", mag=6)






