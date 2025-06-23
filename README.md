# Codici-Tesi-Magistrale-Rolfi
here below the codes for each experiment performed for my thesis 
library(readxl)
library(chron)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

###### figure parameters #######
fntSizeAxisLabel = 17
fntSizeTitle     = 20
fntSizeLegend    = 15
lnWidthAxes      = 1.5
lnWidthPlot      = 1.5
mrkSize          = 3

##### input ######
ThisFolder       = dirname(rstudioapi::getSourceEditorContext()$path)
InputFolder      = ThisFolder
DateEXP          = '240314'
ExperimentName  = paste0(DateEXP,'_GA30C_WT_YPD')
GAFileName       = paste0(ExperimentName,'.xlsx')
GASheetRange     = 'A58:DB248'
GASheetName      = 'Sheet2'
# plate descriptors
ListWellNames   = list(paste0(LETTERS[1:8],rep( 2,each = 1)))
ListWellStrains = c('5063')
ListWellDescrs  = c('WT')
BlankWells      = paste0(LETTERS[1:8],rep(1,each = 1))


# experiment descriptors
Temperature      = 30
Timepoints       = 1:105
TimepointsForFit = 6:Timepoints[length(Timepoints)]
#### analysis parameters ####
# how to read excel
UseMean = T
MeasurementsPerTimePoint = 5

#### output ####
OutputMainFolder  = file.path(ThisFolder)
OutputImgsFolder  = file.path(OutputMainFolder,'imgs')
AnalysisPrefix    = ExperimentName
ImgsPrefix        = AnalysisPrefix

dir.create(OutputImgsFolder, showWarnings = FALSE)
# plot switches
PlotFit = T


#### create database ####
# read excel
DBtemp = as.data.frame(read_xlsx(file.path(InputFolder,GAFileName),sheet = GASheetName,range = GASheetRange,col_names = F))
# take out actual time 
ActualTime = as.numeric(as.matrix(DBtemp[which(DBtemp[,1] == 'Time [s]')[1],2:dim(DBtemp)[2]]))
# convert seconds to hours
ActualTime = ActualTime/60/60 
colnames(ActualTime) = NULL
# take out actual Temperature
ActualTemperature = as.matrix(DBtemp[which(DBtemp[,1] == "Temp. [°C]")[1],2:dim(DBtemp)[2]])
colnames(ActualTemperature) = NULL
## start building DB ##
# read wells
AllWells = DBtemp[seq(2,dim(DBtemp)[1],MeasurementsPerTimePoint+6+1),1]
if (UseMean) {  Wells = rep(AllWells,each=length(Timepoints))  }
# read ODraw and SD
if (UseMean) {
  ODs = DBtemp[DBtemp[,1] == 'Mean'  & !is.na(DBtemp[,1]),2:dim(DBtemp)[2]]
  SDs = DBtemp[DBtemp[,1] == 'StDev' & !is.na(DBtemp[,1]),2:dim(DBtemp)[2]]
}
# create DB
DB = data.frame(Well        = Wells, 
                Timepoint   = Timepoints, 
                ActualTime  = as.vector(t(ActualTime)), 
                Temperature = as.vector(t(ActualTemperature)),
                ODraw       = as.numeric(as.vector(t(ODs))), 
                SD          = as.numeric(as.vector(t(SDs))))
## new variables
#  descriptive 
DB$Well = as.character(DB$Well)
DB$ColNum      = unlist(lapply(DB$Well,function (x) {as.numeric(substr(x,2,nchar(x)))}))
DB$Row         = unlist(lapply(DB$Well,function (x) {substr(x,1,1)}))
DB$RowNum      = unlist(lapply(DB$Row, function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}))
DB$White       = DB$Well %in% BlankWells
DB$Target      = !DB$Well %in% BlankWells
#  measures
DB$CV = DB$SD/DB$ODraw
# assign iStrain / Strain / Descr
DB$Strain = 'blank'
DB$Descr  = 'blank'
for (iStrain in 1:length(ListWellStrains)) {
  Index = DB$Well %in% ListWellNames[[iStrain]]
  DB$iStrain[Index] = iStrain
  DB$Strain[Index]  = ListWellStrains[iStrain]
  DB$Descr[Index]   = ListWellDescrs[iStrain]
}
# read Wells
TargetWells = unique(DB$Well[DB$Target])

#### clean wells/timepoint ####
DB$keep = 1

#### measure net OD ####
bkg = aggregate(DB$ODraw[DB$White],list(DB$Timepoint[DB$White]),mean)
colnames(bkg)[colnames(bkg) == 'Group.1'] = 'Timepoint'
colnames(bkg)[colnames(bkg) == 'x'] = 'bkg'
DB = merge(DB,bkg,all=T)
# 
DB$ODnet = DB$ODraw - DB$bkg

### measure net initial OD
# ODnet0 = initial ODnet for target. one per Well
ODnet0 = DB[DB$Timepoint == 0,c('Well','ODnet')]
colnames(ODnet0)[2] = 'ODnet0'
DB = merge(DB,ODnet0,all=T)



##### fit and plot fitting ####
# define a specific OD for fitting
DB$OD = DB$ODnet
# define a specific measure of time for fitting
DB$Time = DB$ActualTime

RESULTS = data.frame()
if (PlotFit) {
  pdfScale = 0.8
  pdf(file.path(OutputImgsFolder,paste0(ImgsPrefix,'_fit.pdf')), height = 9*pdfScale, width = 16*pdfScale)
}
yLIM = c(log(min(DB$OD[DB$Target])),log(max(DB$OD[DB$Target])))
for (Well in TargetWells) {
  Indexes       = DB$Well == Well
  #  descriptors
  Strain  = unique(DB$Strain[Indexes])
  iStrain = unique(DB$iStrain[Indexes])
  Descr   = unique(DB$Descr[Indexes])
  Row     = unique(DB$Row[Indexes])
  RowNum  = unique(DB$RowNum[Indexes])
  ColNum  = unique(DB$ColNum[Indexes])
  #  at time 0
  ODraw0 = DB$ODraw[DB$Well == Well & DB$Timepoint == Timepoints[1]]
  ODnet0 = DB$ODnet[DB$Well == Well & DB$Timepoint == Timepoints[1]]

  # fit
  IndexesForFit = DB$Well == Well & DB$keep ==T & DB$Timepoint %in% TimepointsForFit
  if (sum(IndexesForFit)) {
    LinearFit = lm(log(OD) ~ Time, data=DB[IndexesForFit,])
    CI        = confint(LinearFit, level = 0.95)
    #
    Slope        = LinearFit$coefficients[2]
    SlopeMin     = CI[2,1]
    SlopeMax     = CI[2,2]
    Intercept    = LinearFit$coefficients[1]
    InterceptMin = CI[1,1]
    InterceptMax = CI[1,2]
    R2        = summary(LinearFit)$r.squared
    Y         = log(DB$OD[IndexesForFit])
    YFit      = LinearFit$fitted.values
    RMSD      = sqrt(mean((Y-YFit)^2))
    NFit      = length(unique(DB$Time[IndexesForFit]))
  } else {
    Slope        = NaN
    SlopeMin     = NaN
    SlopeMax     = NaN
    Intercept    = NaN
    InterceptMin = NaN
    InterceptMax = NaN
    R2           = NaN
    RMSD         = NaN
    NFit         = 0
  }
  # 
  temp = data.frame(ExpDate = DateEXP, Strain, iStrain, Descr, 
                    Well, Row, RowNum, ColNum,
                    ODraw0, ODnet0,
                    Slope, SlopeMin, SlopeMax,
                    Intercept, InterceptMin, InterceptMax,
                    R2, RMSD, NFit)
  rownames(temp) = ''
  
  if (PlotFit) {
    TITLE = paste(Strain,Descr,Well,'\nSlope:',format(Slope,digits =3),'\nR2:',format(R2,digits =3),' RMSD:',format(RMSD,digits =3))
    xLAB = 'Time (hours)'
    
    LogFIT =  ggplot() +
      geom_point(data=DB[Indexes,],      aes(x=Time,y=log(OD),color=Descr),size=1) +
      geom_point(data=DB[IndexesForFit,],aes(x=Time,y=log(OD),color=Descr),size=3) +
      geom_abline(data=temp, aes(slope = Slope, intercept = Intercept)) +
      scale_color_manual(values = ColorzStrain,name='') +
      ggtitle(TITLE) + ylab(paste0('log(ODnet)')) + xlab(xLAB) +
      coord_cartesian(ylim = yLIM) +
      theme_bw() +
      # set axes 
      theme(
        legend.position = 'none',
        axis.title   = element_text(size = fntSizeTitle),
        axis.text    = element_text(size = fntSizeAxisLabel),
        axis.line    = element_line(size = lnWidthAxes),
        plot.title   = element_text(hjust = 0.5, size = fntSizeTitle)
      )
    
    RawFIT =  ggplot() +
      geom_point(data=DB[Indexes,],      aes(x=Time,y=OD,color=Descr),size=1) +
      geom_point(data=DB[IndexesForFit,],aes(x=Time,y=OD,color=Descr),size=3) +
      scale_color_manual(values = ColorzStrain,name='') +
      ggtitle(TITLE) + ylab(WhichODforFit) + xlab(xLAB) +
      coord_cartesian(ylim = c(0,max(DB$OD,na.rm = T))) +
      theme_bw() +
      # set axes 
      theme(
        legend.position = 'none',
        axis.title   = element_text(size = fntSizeTitle),
        axis.text    = element_text(size = fntSizeAxisLabel),
        axis.line    = element_line(size = lnWidthAxes),
        plot.title   = element_text(hjust = 0.5, size = fntSizeTitle)
      )      
    print(grid.arrange(LogFIT,RawFIT,nrow = 1))
  }
  RESULTS = rbind(RESULTS,temp)
}
if(PlotFit) {
  dev.off()
}

RESULTS$Descr = factor(RESULTS$Descr, levels = names(ColorzStrain), ordered = T)
RESULTS$N     = 1:dim(RESULTS)[1]

#### clean RESULTS ####
RESULTS$keep = 1


#### plot final results ####
if (T) {
  PLOT = RESULTS[RESULTS$keep == 1,]
  
  yLIM = c(min(RESULTS$Slope),max(RESULTS$Slope))
  yLAB = 'Fitted slope'
  
  pdfScale = 0.8
  pdf(file.path(OutputImgsFolder,paste0(ImgsPrefix,'_RESULTS_boxplot.pdf')),
      width=16*pdfScale,height =9*pdfScale)
  print(
    ggplot(PLOT,aes(x = Descr, y=Slope)) +
      geom_boxplot(aes(fill=Descr),outlier.shape = NA) +
      geom_jitter(width = 0.3,height = 0,size = .5) +
      scale_fill_manual(values = ColorzStrain, name = "") +
      ylab(yLAB) + xlab('') +
      coord_cartesian(ylim=yLIM) +
      theme_bw()+
      theme(
        legend.position='none',
        legend.text  = element_text(size = fntSizeAxisLabel),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.0,"cm"),
        plot.margin = unit(c(4, 1, 4, 1), "lines"),
        axis.title   = element_text(size = fntSizeTitle),
        axis.text    = element_text(size = fntSizeAxisLabel),
        axis.text.x  = element_text(angle = 45,hjust=0.95,vjust=.9), 
        axis.line    = element_line(size = lnWidthAxes),
        plot.title   = element_text(hjust = 0.5, size = fntSizeTitle),
        strip.text.x = element_text(size = fntSizeAxisLabel),
        strip.text.y = element_text(size = fntSizeAxisLabel)
      )
  )
  dev.off()
  
  
  
  pdfScale = 0.8
  pdf(file.path(OutputImgsFolder,paste0(ImgsPrefix,'_RESULTS_scatter.pdf')),
      width=16*pdfScale,height =9*pdfScale)
  print(
    ggplot(data = PLOT, aes(x = N, y = Slope, color = Descr)) + 
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = SlopeMin, ymax = SlopeMax),width = 0) +
      scale_color_manual(values = ColorzStrain, name = "") +
      ylab(yLAB) + xlab('') +
      coord_cartesian(ylim=yLIM) +
      theme_bw() +
      # set axes 
      theme(
        legend.text  = element_text(size = fntSizeAxisLabel),
        axis.title   = element_text(size = fntSizeTitle),
        axis.text    = element_text(size = fntSizeAxisLabel),
        axis.text.x  = element_blank(),
        axis.line    = element_line(size = lnWidthAxes),
        plot.title   = element_text(hjust = 0.5, size = fntSizeTitle)
      )
  )
  dev.off()
}




