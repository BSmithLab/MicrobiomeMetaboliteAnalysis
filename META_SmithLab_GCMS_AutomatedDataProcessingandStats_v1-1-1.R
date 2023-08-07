
################################################# META ANALYSIS #########################################################################################################
#This is a companion script designed to facilitate comparison of GC-MS results across multiple independent experiments

#Automated Alignment and Analysis of Smith Lab GC-MS data
#Written By: Jarrett Eshima
#Fall 2022 - Summer 2023
#SBHSE, Arizona State University

#Version History
#V1-0-0: Base script, statistic functionality not included
#V1-1-0: Allows functional group analysis to be performed at the experimental level
#V1-1-1: Quick fix to address extra column bug in MatrixtoHeatmap2() when customsamps = T
#vDev: Statistics and bar plot functionality added

#Required File:
#UsefulFunctions_v2.R

#########################################################################################################################################################################
###############################################   READ ME: CODE USAGE   #################################################################################################
#########################################################################################################################################################################

#This code performs the following:
#1) Not Optional: Re-aligns VOCs across multiple experiments. No additional filtering or normalization is applied. 
#2) Optional: Plots data as heatmap - this requires the file UsefulFunctions_v2.R (see below)
#3) Optional: Performs principal component analysis
#4) Optional: Performs analysis of functional groups using established suffix conventions 
#5) Optional: Runs parametrics statistics 
#6) Optional: Plots the data as bar plots, automatically adds standard error bars and significance


#From the UsefulFunctions_v2.R script, you MUST load the following functions into your workspace:
#1) zscore()
#2) MatrixtoHeatmap2() - V5 and above requires the upgraded function: includes row/column clustering functionality
#3) FunctionalGroups()

#Important Notes:
#1) This code will only function properly if Prep4Meta is set to TRUE for all individual experiments.
#See Previous Script: SmithLab_GCMS_AutomatedDataProcessingandStats_v6-6-0 and later

######################################################################################################################################################################
################################################              ########################################################################################################
################################################  USER INPUT  ########################################################################################################
################################################              ########################################################################################################
######################################################################################################################################################################

######################################  FILE INFORMATION  ##############################################

#Get .csv files in a specific directory (this should match the MetaDirectory object from the previous script, the code will automatically navigate to the correct folder)
MetaDirectory = "C:/Users/jeshima/Documents/Smith Lab/Spring 2023/BMC Microbiology Paper/Data"

#Provide the order for displaying individual experiments
#Note: This order may change depending on which experiment has the greatest number of VOCs observed. This cannot be adjusted due to the script logic.
ExperimentOrder = c("EFaecalis","CSporogenes","ELenta") #These must be the same as the "organism" object in the previous script (CASE SENSITIVE)


#Provide the order for experimental conditions/groups; each list of conditions should correspond to the order provided in "ExperimentOrder"
#Allows up to 4; UNIQUE FILE HANDLE MUST BE EXACTLY 3 CHARACTERS LONG - CASE SENSITIVE
#Only accepted strings (in any order): Pos, Neg, Ctr, Ext
#Ctr is used if the number of conditions is greater than or equal to 3 in the previous script; Ext is used if the number of conditions is equal to 4 in the previous script
conditions = c("Pos","Ext","Neg","Ctr") 


#Experiment Name (without spaces)
AnalysisName = "Meta"

#Write Files?
writef = FALSE

###################################  FILTERING AND ALIGNMENT #################################################

#Retention time tolerance for aligning samples across experiments
RTtol = 0.011 # ~1 second variability between samples


#Remove features/VOCs seen in a reference level?
FilterUsingControls = F
#If TRUE, specify the reference level (typically "Ctr")
reflevel = "Ctr"
#If TRUE, do you want to remove the control samples from downstream analysis?
removecontrols = T #If this is set to false, control samples are retained (they will have missing / 0 abundance observations for all remaining VOCs)

#Do you want to re-order the data matrix so that the peaks are ordered by chromatographic elution / retention time across ALL experiments?
#If false, then the peaks are provided in order of retention time, but separated by experimental group
PeaksInOrder = F

##########################  REMOVE KNOWN CONTAMINANTS AND MANUALLY FILTER  ####################################

#Do you want known contaminants (siloxanes) to be removed?
removecontam = T

##Remove VOCs in a specific retention time window 
#(v6.4.1 allows multiple filter windows, with high and low RTs paired together in the order they are provided)
ClearNoise = F #If false, the VOCs are not removed.
Remove_LowRT = c(2,14.14) #Not needed, if ClearNoise = False
Remove_HighRT = c(3,15.6) #Not needed, if ClearNoise = False


###################################################################################################################
########################################  PLOTTING  ###############################################################
###################################################################################################################

################################### Heatmap Parameters ############################################################

#Plot the Heatmap?
PlotHeatmap = T

#Optional arguments (defaults are provided below):
heatcolorlow = "#17bac7"
heatcolormid = "white"
heatcolorhigh = "#ff6600"
limits = c(-3,3)
clusterrows = T
clustercols = F
      clustermethod = "ward.D2"
      distancemethod = "manhattan"
      #A custom row/column order can be provided. This will overwrite clustering results if used. 
      #With 3 groups
      #customsamps = c("Replicate1_Pos","Replicate2_Pos","Replicate3_Pos","Replicate1_Neg","Replicate2_Neg","Replicate3_Neg","Replicate1_Ctr","Replicate2_Ctr","Replicate3_Ctr") #Column order for plotting, starting with the left side
      #With 4 groups
      usecustomsamps = T
      customsamps = c("CSporogenesReplicate1_Neg","CSporogenesReplicate2_Neg","CSporogenesReplicate3_Neg","EFaecalisReplicate1_Neg","EFaecalisReplicate2_Neg","EFaecalisReplicate3_Neg","ELentaReplicate1_Neg","ELentaReplicate2_Neg","ELentaReplicate3_Neg")
      #customsamps = c("CSporogenesReplicate1_Pos","CSporogenesReplicate2_Pos","CSporogenesReplicate3_Pos","CSporogenesReplicate1_Neg","CSporogenesReplicate2_Neg","CSporogenesReplicate3_Neg","CSporogenesReplicate1_Ctr","CSporogenesReplicate2_Ctr","CSporogenesReplicate3_Ctr","EFaecalisReplicate1_Pos","EFaecalisReplicate2_Pos","EFaecalisReplicate3_Pos","EFaecalisReplicate1_Ext","EFaecalisReplicate2_Ext","EFaecalisReplicate3_Ext","EFaecalisReplicate1_Neg","EFaecalisReplicate2_Neg","EFaecalisReplicate3_Neg","EFaecalisReplicate1_Ctr","EFaecalisReplicate2_Ctr","EFaecalisReplicate3_Ctr","ELentaReplicate1_Pos","ELentaReplicate2_Pos","ELentaReplicate3_Pos","ELentaReplicate1_Ext","ELentaReplicate2_Ext","ELentaReplicate3_Ext","ELentaReplicate1_Neg","ELentaReplicate2_Neg","ELentaReplicate3_Neg","ELentaReplicate1_Ctr","ELentaReplicate2_Ctr","ELentaReplicate3_Ctr") #Column order for plotting with 4 conditions, starting with the left side
      usecustomfeats = F
      customfeats = c("VOC1","VOC2","VOC3","etc.") #This pipes well into external clustering methods
      #Adjust z-score values greater than 3 to 3 and less than -3 to -3 (for plotting purposes)
      maxz3 = T
      
      ############################### Principal Component Analysis #######################################################
      RunPCA = T
      
      #Which principal components do you want to plot? Typically, the first and second components capture the greatest amount of variability in your data.
      comp1 = 1
      comp2 = 2
      
      #Plot 3-axis (3D) PCA?
      pca3d = F
      #If TRUE:
      comp3 = 3
      
      #Colors for plotting groups?
      pcacolor1 = "#61bd7e"
        pcacolor2 = "#fc8435"
          pcacolor3 = "#07bdd4"
            pcacolor4 = "mediumaquamarine" #This is only needed if 4 experiments are provided
              
            #Use a specific experimental condition? (for example, only consider the negative or positive experimental conditions in the PCA space)
            specificExpCond = T
            PCAfocus = "Neg"  #THIS MUST BE ONE OF THE FOUR VALID STRINGS PROVIDED IN "CONDITIONS" OBJECT (ABOVE)
            
            
############################ Analysis of Functional Groups ##########################################################
#install.packages("grDevices")
library(grDevices)

#Should the VOC functional groups be assessed and plotted as a pie chart?
AnalyzeFunctionalGroups = T

#Do you want to split the pie chart according to EXPERIMENT?
SplitPie = T

#If you're plotting the nested pie chart, do you want to use default ggplot2 colors?
usedefaultcolors = F

#Should VOCs with an ambiguous functional group (using substring searching) be assigned as unknown? If yes, set the variable below to TRUE
flagambiguous = F

#Do you want to overwrite the functional group names assigned? 
customFGlist = F

#Provide custom list
customlist = T

#IMPORTANT: Specify the functional group names below as a list. You MUST provide a list that is equal to the number of unique funtional groups. You can check this easily after running the code once with: length(names(table(FunctionalGroups(rownames(data)))))

List1 = c("Heteroaromatic","Heteroaromatic","Unknown","Heteroaromatic","Aldehyde","Heteroaromatic","Aromatic","Aldehyde","Other","Heteroaromatic","Aldehyde","Aldehyde","Nitrile","Hydrocarbon","Other","Acid","Heteroaromatic","Heteroaromatic","Aromatic","Hydrocarbon","Other")
List2 = c("Amine","Unknown","Heteroaromatic","Aromatic","Heteroaromatic","Heteroaromatic","Heteroaromatic","Aldehyde","Aldehyde","Ketone","Other","Aldehyde","Hydrocarbon","Hydrocarbon","Other","Other")
List3 = c("Unknown","Unknown","Unknown","Unknown","Unknown","Unknown","Amine","Unknown","Unknown","Alcohol","Heteroaromatic","Aldehyde","Ester","Heteroaromatic","Other","Ester","Ester","Ester","Other","Ester","Ester","Ester","Ester","Alcohol","Nitrile","Ester","Hydrocarbon","Ester","Ester","Alcohol","Ester","Unknown","Unknown","Hydrocarbon","Ester","Ester","Ester","Other","Unknown","Ester","Ester","Hydrocarbon","Ester")

FGlabels = list(List1,List2,List3)

############# Colors for the piechart
color1 = "tomato"
color2 = "skyblue"
  
userainbow = T #If this is true, color1 and color2 are ignored
usecustomrainbow = T #If this is true, the color palette specified below is used
piecolorfnc = colorRampPalette(c("#ff4040","#ff754e","#fd9a49","#ffed73","#6dda6d","#00d1eb","#12b4f5","#7791ff","#ce6fb5","#f75582")) #You can select as many or as few colors as you would like. The color order can be in any order you like.
#DEFAULT: piecolorfnc = colorRampPalette(c("red","yellow","green","blue","violet"))

#Later versions...

# ################################  Target VOCs to Plot  ###################################
# 
# #Do you want a subset of VOCs plotted? If True, only the names provided in 'VOClist' are considered during plotting
# SpecificVOCs = F
# 
# #List of VOCs to plot. IMPORTANT: Names are CASE SENSITIVE. Copy the name exactly, including symbols
# #MINIMUM 2 VOCs
# VOClist = c("Isopentyl hexanoate","Hexanoic acid, 2-methylbutyl ester")
# 
# 
# ##############################################  STATISTICS  ########################################################
# 
# #Hypothesis test tails (Options - case sensitive - "less","greater","two.sided")
# #Less and greater options are one-tailed t-tests
# tails = "two.sided"
# #Do you want to assume equal variance? (unequal variance is often a safer assumption)
# eqvar = F
# #Is the data paired? (TRUE if same sample/individual/experiment at two different time points)
# paired = F
# 
# #Multiple hypothesis test correction method (see ?p.adjust.methods):
# MHCmethod = "BH"
# 
# 
# #User specified statistics. If uselevels = F, then a student's t-test is applied for all pairwise comparisons. The FDR/Benjamini-Hochberg method is used for multiple hypothesis test corrections.
# uselevels = F
# 
# #If uselevels = TRUE:
# #Provide one or more pairs for statistical testing. If you have two or more pairs, specify them in order. 
# #Important: These should be the exact strings provided in "conditions" (line 76)
# #Example: c("Pos","Neg","Pos","Ctr") will compare Pos vs Neg AND Pos vs Control.
# pairs = c("Pos","Neg")
# #pairs = c("Pos","Neg","Pos","Ctr")
# 
# ################################################################################################
# 

# ## Parameters for Bar Plots ##
# 
# PlotSigOnly = T #Should only the significant genes be plotted (after p-value adjustment for multiple hypothesis testing)
# 
# #Note if PlotSigOnly = F, you will need to specify up to 30 rows to plot at a time.
# startrow = 1
# stoprow = 5
# 
# AllowUnadjusted = F #Should the unadjusted p-values be considered/used? This might be useful to run if VOCs are weakly significant.
# sigthresh = 0.1 #Significance threshold
# addpvalues = T #Should p-values be shown for VOCs that do NOT pass the significance threshold?
# #Note: For 3 conditions, this works as expected. For 4 conditions, this provides the necessary plotting space to manually add significance in a graphic editor (Illustrator/Inkscape)
# givemepspace = T
# 
# #Barplot Colors:
# poscolor = "coral"
# negcolor = "deepskyblue2"
# ctrcolor = "grey60"
# excolor = "mediumaquamarine" #This is only needed if 4 conditions are provided
# lowerlim = -1 #manually adjust lower limit

#######################################################################################################################################
#######################################  END OF USER INPUT ############################################################################
#######################################                    ############################################################################
#######################################   **HIT SOURCE**   ############################################################################
#######################################################################################################################################



#########################################################  Load Data  #################################################################

wd = paste(MetaDirectory,"/Meta",sep="")
setwd(wd)
files = list.files(wd,".csv")
refs = list.files(wd,".RData")
nExperiments = length(ExperimentOrder)
Norms = expsamp = rep(NA,nExperiments)


for(i in 1:nExperiments){
  
  #Load reference file to the workspace
  load(refs[i])
  fillindex = which(substr(refs,1,nchar(gsub("_ref.RData","",refs))) == ExperimentOrder[i])
  assign(paste(ExperimentOrder[fillindex],"_ref",sep=""),reffile)
  
  assign(paste("n",fillindex,sep=""),sum(reffile$nsamples))
  
  Norms[fillindex] = reffile$NormMethod
  expsamp[fillindex] = sum(reffile$nsamples)
  nconditions = length(reffile$conditions)
  assign(paste(ExperimentOrder[fillindex],"conditions",sep=""),reffile$conditions)
  
  #Defines the column order that is always output by the previous script.
  if(nconditions == 2){
    setorder = c("Pos","Neg")
  }else if(nconditions == 3){
    setorder = c("Pos","Neg","Ctr")
  }else if(nconditions == 4){
    setorder = c("Pos","Neg","Ctr","Ext")
  }
  
  nnames = paste(ExperimentOrder[fillindex],"_n",setorder,sep="")
  
  #Assign the number of samples for all groups
  
  for(z in 1:length(setorder)){
    ref = get(paste(ExperimentOrder[fillindex],"_ref",sep=""))
    assign(nnames[z],ref$nsamples[z])
  }
  
  if(i == 1){
    totsamp = sum(ref$nsamples)
  }else{
    totsamp = totsamp + sum(ref$nsamples)
  }
  
  
  #Read in expression
  
  tind = which(gsub(".csv","",files) == ExperimentOrder[fillindex])
  tmpVOC = read.csv(files[tind])
  rownames(tmpVOC) = tmpVOC$X
  tmpVOC = tmpVOC[,-1]
  
  #Re-order according to user input
  index = rep(NA,ncol(tmpVOC))
  for(j in 1:ncol(tmpVOC)){
    
    tmpstr = colnames(tmpVOC)[j]
    
    for(k in 1:length(conditions)){
      for(l in 1:nchar(tmpstr)){
        
        if(substr(tmpstr,l,l+2) == conditions[k]){
          index[j] = k
        }
        
      }
    }
    
  }
  
  index = index[!is.na(index)]
  findex = order(index)
  
  fVOC = tmpVOC[,findex]
  fVOC$RT = tmpVOC$RT
  
  #Write file to workspace
  assign(ExperimentOrder[fillindex],fVOC)
  
}

if(length(table(Norms))>1){
  cat("MAJOR CONCERN: THE SAME NORMALIZATION METHOD WAS NOT USED FOR ALL EXPERIMENTS.\n RESULTS ARE LIKELY INVALID. \n")
}


################################################# Experimental Alignment ##############################################################

#Define reference experiment for matrix filling
maxVOCs = rep(NA,length(ExperimentOrder))
for(i in 1:length(maxVOCs)){
  tmp = get(ExperimentOrder[i])
  maxVOCs[i] = nrow(tmp)
}

nVOCs = max(maxVOCs)
exref = which(maxVOCs == nVOCs)
myseq = seq(1,nExperiments)

ExpLeft = ExperimentOrder[-exref]
trueexp = c(ExperimentOrder[exref],ExpLeft)

#Build aligned matrix
tmp = get(ExperimentOrder[exref])
colnames(tmp) = paste(ExperimentOrder[exref],colnames(tmp),sep="")
colnames(tmp)[which(colnames(tmp) == paste(ExperimentOrder[exref],"RT",sep=""))] = "RT"
rest = myseq[!myseq %in% exref]
fillmat = data.frame(matrix(NA,nrow(tmp),sum(expsamp[rest])))
MetaData = cbind(tmp,fillmat)
MetaData = MetaData[,-which(colnames(MetaData) == "RT")]
MetaData$RT = tmp$RT

blank = rownames(MetaData)
MetaData = cbind(blank,MetaData)
colnames(MetaData)[1] = "VOC"


#Fill in missing column names


fillcols = as.numeric(which(is.na(colSums(MetaData[,-1]))))+1
mycols = rep(NA,length(fillcols))
for(i in 1:length(ExpLeft)){
  
  for(j in 1:length(conditions)){
    
    if(exists(paste(ExpLeft[i],"_n",conditions[j],sep=""))){
      
      nsamp = get(paste(ExpLeft[i],"_n",conditions[j],sep=""))
      
      if(i == 1 && j == 1){
        
        mycols = paste(ExpLeft[i],"Replicate",seq(1,nsamp),"_",conditions[j],sep="")
        
      }else{
        
        mycols = c(mycols,paste(ExpLeft[i],"Replicate",seq(1,nsamp),"_",conditions[j],sep=""))
        
      }
      
    }
    
  }
  
}

colnames(MetaData)[fillcols] = mycols


#Fill the matrix, flag VOCs that don't fall within the RT tolerance
MetaData$RT = as.numeric(MetaData$RT)
nfill = rep(NA,length(ExpLeft))

for(i in 1:length(ExpLeft)){
  
  dat = get(ExpLeft[i])
  colnames(dat) = paste(ExpLeft[i],colnames(dat),sep="")
  colnames(dat)[which(colnames(dat) == paste(ExpLeft[i],"RT",sep=""))] = "RT"
  nfill[i] = expsamp[which(ExpLeft[i] == ExperimentOrder)]
  
  if(i == 1){
    overflowind = seq(1,nfill[i])
    overflowind = fillcols[overflowind]
  }else{
    overflowind = seq(nfill[i-1]+1,nfill[i-1]+nfill[i])
    overflowind = fillcols[overflowind]
  }
  
  MetaData$RT = as.numeric(MetaData$RT)
  MetaRTs = MetaData$RT
  ExpRTs = dat$RT
  
  LowerB = MetaRTs - RTtol
  UpperB = MetaRTs + RTtol
  keep = remove = rep(NA,length(ExpRTs))
  
  for(j in 1:length(MetaRTs)){
    
    for(k in 1:length(ExpRTs)){
      
      if(ExpRTs[k] > LowerB[j] && ExpRTs[k] < UpperB[j]){
        keep[k] = k
      }
      
      
    }
  }
  
  
  tmpseq = seq(1,length(ExpRTs))
  keep = keep[!is.na(keep)]
  remove = tmpseq[! tmpseq %in% keep]
  
  filldat = dat[keep,]
  adddat = dat[remove,]
  
  adddat2 = adddat[,-which(colnames(adddat) == "RT")]
  
  for(j in fillcols){
    
    for(k in 1:ncol(filldat)){
      
      if(colnames(MetaData)[j] == colnames(filldat)[k] && colnames(filldat)[k] != "RT"){  #Column names must match before alignment proceeds
        
        cat(paste("Aligning Sample:",colnames(filldat)[k],"\n"))
        
        #Alignment
        for(m in 1:nrow(MetaData)){
          
          for(n in 1:nrow(filldat)){
            
            if(filldat$RT[n] > (as.numeric(MetaData$RT[m])-RTtol) && filldat$RT[n] < (as.numeric(MetaData$RT[m])+RTtol)){
              MetaData[m,j] = filldat[n,k]
            }
            
          }
          
        }
      }
      
    }
    
  }
  
  for(z in 1:nrow(adddat)){
    
    if(z == 1){
      blankvector = rep(NA,ncol(MetaData))
      blankvector[overflowind] = as.numeric(adddat2[z,])
      blankvector[length(blankvector)] = adddat$RT[z]
      blankvector[1] = rownames(adddat)[z]
      flagged = blankvector
    }else{
      blankvector = rep(NA,ncol(MetaData))
      blankvector[overflowind] = as.numeric(adddat2[z,])
      blankvector[length(blankvector)] = adddat$RT[z]
      blankvector[1] = rownames(adddat)[z]
      flagged = rbind(flagged,blankvector)
    }
    
  }
  flagged = data.frame(flagged)
  colnames(flagged) = colnames(MetaData)
  MetaData = rbind(MetaData,flagged)
  
}

rownames(MetaData) = make.unique(MetaData$VOC)
MetaData = MetaData[,-1]
MetaData$RT = as.numeric(MetaData$RT)

if(PeaksInOrder == T){
  MetaData = MetaData[order(MetaData$RT,decreasing = F),]
}


if(FilterUsingControls == T){
  
  FilterData = MetaData
  
  nstr = nchar(reflevel)
  keep = rep(NA,ncol(FilterData))
  
  for(i in 1:ncol(FilterData)){
    
    currentcol = colnames(FilterData)[i]
    
    for(j in 1:nchar(currentcol)){
      
      if(substr(currentcol,j,j+nstr-1) == reflevel || currentcol == "RT"){
        
        keep[i] = i
        
      } 
      
    }
    
  }
  keep = keep[!is.na(keep)]
  FilterData = FilterData[,keep]
  keep = keep[-length(keep)]
  
  tmp = matrix(as.numeric(unlist(FilterData)),nrow(FilterData))
  colnames(tmp) = colnames(FilterData)
  rownames(tmp) = rownames(FilterData)
  FilterData = tmp
  
  for(i in 1:nrow(FilterData)){
    for(j in 1:ncol(FilterData)){
      if(is.na(FilterData[i,j])){
        FilterData[i,j] = 0
      }
    }
  }
  
  
  remove2 = as.numeric(which(rowSums(FilterData[,-ncol(FilterData)],na.rm = T) == 0))
  FilterData = FilterData[-remove2,]
  
  FilterRTs = as.numeric(FilterData[,ncol(FilterData)])
  
  MetaData = MetaData[! round(MetaData$RT,2) %in% round(FilterRTs,2),]
  
  if(removecontrols == T){
    MetaData = MetaData[,-keep]
  }
  
  
}


#Set missing observations to 0 for plotting purposes
for(i in 1:nrow(MetaData)){
  for(j in 1:ncol(MetaData)){
    if(is.na(MetaData[i,j])){
      MetaData[i,j] = 0
    }
  }
}

#Converts to numeric matrix
tmp = matrix(as.numeric(unlist(MetaData)),nrow(MetaData))
colnames(tmp) = colnames(MetaData)
rownames(tmp) = rownames(MetaData)
FinalData = tmp




################################################################################################

#This is the aligned and normalized data containing all experimental conditions/data
View(FinalData)

################################################################################################

############################################  POST-PROCESSING  ####################################################

### Strongly Recommended: Remove Known Contaminants ("silane")
#Note: "Sil" compounds have been retained in case there is noise driven mis-naming (73,144 ions). These can be confirmed by running purified standards using the same method. 
if(removecontam == T){
  contamindex = rep(NA,nrow(FinalData))
  count = 1
  for(i in 1:nrow(FinalData)){
    
    VOC = rownames(FinalData)[i]
    
    ncharVOC = nchar(VOC)
    
    for(j in 1:ncharVOC-5){
      
      if(substr(VOC,j,j+5) == "Silane" || substr(VOC,j,j+5) == "silane"){
        contamindex[count] = i
        count = count+1
      }
      
      if(substr(VOC,j,j+7) == "Siloxane" || substr(VOC,j,j+7) == "siloxane"){
        contamindex[count] = i
        count = count+1
      }
      
    }
    
  }
  contamindex = contamindex[!is.na(contamindex)]
  rownames(FinalData)[contamindex]
  if(length(contamindex)>0){
    FinalData = FinalData[-contamindex,]
  }
}


### Optional: Manually remove VOCs in a specific retention time window
if(ClearNoise == T){
  tmpdata = data.frame(FinalData)
  
  if(length(Remove_LowRT) == length(Remove_HighRT)){
    for(a in 1:length(Remove_LowRT)){
      
      manualrem = rep(NA,nrow(FinalData))
      count = 1
      
      for(i in 1:nrow(FinalData)){
        
        if(tmpdata$RT[i] <= Remove_HighRT[a] && tmpdata$RT[i] >= Remove_LowRT[a]){
          manualrem[count] = i
          count = count+1
        }
        
      }
      manualrem = manualrem[!is.na(manualrem)]
      FinalData = FinalData[-manualrem,]
    }
  }else{
    cat("ERROR: You have provided an unequal number of retention time pairs for filtering. \n")
    break
  }
}

##################################  Write Aligned File  ###############################################

#Write out the aligned file?
if(writef == T){
  
  setwd(MetaDirectory)
  filename = paste("AlignedExperiments_MetaData",Sys.Date(),".csv",sep="")
  write.csv(FinalData,filename)
  
}


####################################################  HEATMAP  #############################################################
mydataforstats = FinalData

#Heatmap (remember to z-score normalize your data before plotting or else the scale will not come out in the correct range (-3,3))

if(PlotHeatmap == T){
  NormMethod = Norms[1]
  
  if(NormMethod != 9 && NormMethod != 17 && NormMethod != 18 && NormMethod != 19){
    mydataforstats2 = mydataforstats[,-which(colnames(mydataforstats) == "RT")]
    finaldata_z = zscore(mydataforstats2,featuresarerows = T,removeNA = T)
    colnames(finaldata_z) = colnames(mydataforstats2)
    rownames(finaldata_z) = rownames(mydataforstats2)
    
    if(maxz3 == T){
      for(i in 1:nrow(finaldata_z)){
        for(j in 1:ncol(finaldata_z)){
          if(finaldata_z[i,j] > 3){
            finaldata_z[i,j] = 3
          }else if(finaldata_z[i,j] < -3){
            finaldata_z[i,j] = -3
          }
        }
      }
    }
    
    if(usecustomsamps == T){
      if(usecustomfeats == T){
        MatrixtoHeatmap2(finaldata_z,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customsamps = customsamps,customfeats = customfeats)
      }else{
        MatrixtoHeatmap2(finaldata_z,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customsamps = customsamps)
      }
    }else{
      if(usecustomfeats == T){
        MatrixtoHeatmap2(finaldata_z,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customfeats = customfeats)
      }else{
        MatrixtoHeatmap2(finaldata_z,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T)
      }
    }
    
  }else{
    TransformedData2 = TransformedData[,-which(colnames(TransformedData) == "RT")]
    
    
    if(maxz3 == T){
      for(i in 1:nrow(TransformedData2)){
        for(j in 1:ncol(TransformedData2)){
          if(TransformedData2[i,j] > 3){
            TransformedData2[i,j] = 3
          }else if(TransformedData2[i,j] < -3){
            TransformedData2[i,j] = -3
          }
        }
      }
    }
    
    if(usecustomsamps == T){
      if(usecustomfeats == T){
        MatrixtoHeatmap2(TransformedData2,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customsamps = customsamps,customfeats = customfeats)
      }else{
        MatrixtoHeatmap2(TransformedData2,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customsamps = customsamps)
      }
    }else{
      if(usecustomfeats == T){
        MatrixtoHeatmap2(TransformedData2,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T,customfeats = customfeats)
      }else{
        MatrixtoHeatmap2(TransformedData2,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T)
      }
    }
  }
  
}


############################  Principal Component Analysis  #############################################################

if(RunPCA == T){
  
  tmpseq = seq(1,nExperiments)
  mypal = paste("pcacolor",tmpseq,sep="")
  
  for(i in 1:nExperiments){
    if(i == 1){
      tmp = get(mypal[i])
      pcacolors = rep(tmp,n1)
    }else{
      tmp = get(mypal[i])
      ntmp = get(paste("n",i,sep=""))
      pcacolors = c(pcacolors,rep(tmp,ntmp))
    }
  }
  
  pcadata = mydataforstats[,-which(colnames(mydataforstats) == "RT")]
  
  if(specificExpCond == T){
    
    nstr = nchar(PCAfocus)
    keep = rep(NA,ncol(pcadata))
    
    for(i in 1:ncol(pcadata)){
      
      currentcol = colnames(pcadata)[i]
      
      for(j in 1:nchar(currentcol)){
        
        if(substr(currentcol,j,j+nstr-1) == PCAfocus){
          
          keep[i] = i
          
        } 
        
      }
      
    }
    keep = keep[!is.na(keep)]
    pcadata = pcadata[,keep]
    
    for(i in 1:nExperiments){
      if(i == 1){
        tmp = paste(trueexp[i],"_n",PCAfocus,sep="")
        nrep = get(tmp)
        pcacolors = rep(pcacolor1,nrep)
      }else if(i == 2){
        tmp = paste(trueexp[i],"_n",PCAfocus,sep="")
        nrep = get(tmp)
        pcacolors = c(pcacolors,rep(pcacolor2,nrep))
      }else if(i == 3){
        tmp = paste(trueexp[i],"_n",PCAfocus,sep="")
        nrep = get(tmp)
        pcacolors = c(pcacolors,rep(pcacolor3,nrep))
      }else if(i == 4){
        tmp = paste(trueexp[i],"_n",PCAfocus,sep="")
        nrep = get(tmp)
        pcacolors = c(pcacolors,rep(pcacolor4,nrep))
      }
    }
    
    
    
  }
  
  #Run PCA
  pca = prcomp(t(pcadata))
  PC1 = c(pca$x[,comp1])
  PC2 = c(pca$x[,comp2])
  
  #Get axis bounds
  min1 = min(PC1)
  max1 = max(PC1)
  min2 = min(PC2)
  max2 = max(PC2)
  
  if(nExperiments == 2){
    legcolors = c(pcacolor1,pcacolor2)
  }else if(nExperiments == 3){
    legcolors = c(pcacolor1,pcacolor2,pcacolor3)
  }else if(nExperiments == 4){
    legcolors = c(pcacolor1,pcacolor2,pcacolor3,pcacolor4)
  }
  
  if(pca3d == T){
    require(rgl)
    PC3 = c(pca$x[,comp3])
    min3 = min(PC3)
    max3 = max(PC3)
    rgl::plot3d(PC1,PC2,PC3,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",AnalysisName,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5,size=12)
    legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = trueexp,col = legcolors,pch=20,pt.cex = 1.5)
  }else{
    plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",AnalysisName,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
    legend(max1 + (0.05*max1),max2 + (0.55*max2),legend = trueexp,col = legcolors,pch=20,pt.cex = 1.5)
  }
  
  
}

#######################################  Functional Group Analysis  #############################################################

if(AnalyzeFunctionalGroups == T){
  
  if(userainbow == T){
    if(usecustomrainbow == T){
      piecolorfnc = piecolorfnc
    }else{
      piecolorfnc = colorRampPalette(c("red","yellow","green","blue","violet"))
    }
  }else{
    piecolorfnc = colorRampPalette(c(color1,color2))
  }
  
  
  if(SplitPie == T){
    
    for(i in 1:nExperiments){
      if(customlist == F){
        handle = ExperimentOrder[i]
        nc = nchar(handle)
        ind = which(substr(colnames(mydataforstats),1,nc) == handle)
        tmp = mydataforstats[,ind]
        tmp = tmp[-which(rowSums(tmp) == 0),]
        assign(paste(handle,"_Filt",sep=""),tmp)
        tmp2 = FunctionalGroups(rownames(tmp),flagambiguous)
        assign(paste(handle,"_FG",sep=""),tmp2)
      }else{
        handle = ExperimentOrder[i]
        nc = nchar(handle)
        tmp2 = FGlabels[i]
        assign(paste(handle,"_FG",sep=""),tmp2)
      }
      
      
      FGs = names(table(tmp2))
      assign(paste(handle,"_Groups",sep=""),FGs)
      
      if(i == 1){
        runningGroups = FGs
      }else{
        runningGroups = c(runningGroups,FGs)
      }
    }
    nFG = length(names(table(runningGroups)))
    piecolors = piecolorfnc(nFG)
    
    CleanGroups = names(table(runningGroups))
    for(j in 1:length(CleanGroups)){
      if(j == 1){
        FunctionalGroupList = rep(CleanGroups[j],nExperiments)
        customcolors = rep(piecolors[j],nExperiments)
      }else{
        FunctionalGroupList = c(FunctionalGroupList,rep(CleanGroups[j],nExperiments))
        customcolors = c(customcolors,rep(piecolors[j],nExperiments))
      }
    }
    CleanExp = rep(trueexp,length(CleanGroups))
    
    MetaFG = data.frame(matrix(NA,nrow=length(FunctionalGroupList),ncol = 4))
    colnames(MetaFG) = c("Experiment","FunctionalGroup","Percent","Color")
    MetaFG$Experiment = CleanExp
    MetaFG$FunctionalGroup = FunctionalGroupList
    MetaFG$Color = customcolors
    
    i=1
    for(i in 1:nExperiments){
      handle = ExperimentOrder[i]
      tmp = get(paste(handle,"_FG",sep=""))
      tmp = table(tmp)
      tmp = tmp/sum(as.numeric(tmp))
      tmp = tmp*100
      
      for(m in 1:nrow(MetaFG)){
        for(n in 1:length(names(tmp))){
          if(handle == MetaFG$Experiment[m] && names(tmp)[n] == MetaFG$FunctionalGroup[m]){
            MetaFG$Percent[m] = as.numeric(tmp[n])
          }
        }
      }
    }
    
    for(i in 1:nrow(MetaFG)){
      if(is.na(MetaFG$Percent[i])){
        MetaFG$Color[i] = NA
      }
    }
    
    customcolors = MetaFG$Color[!is.na(MetaFG$Percent)]
    
    cat("Note: The inner cicle is from Experiment:",trueexp[1],"\n")
    cat("Note: The outer cicle is from Experiment:",trueexp[length(trueexp)],"\n")
    
    if(usedefaultcolors == F){
      p = ggplot(MetaFG, aes(x = Experiment, y = Percent, fill = FunctionalGroup)) +
        geom_bar(stat = "identity",position = "stack") +
        scale_x_discrete(limits = trueexp)
      p = p+theme_bw()
      p = p +ggtitle("Bacterial Functional Group Analysis") + xlab("Strain") + ylab("Percent")
      p = p+theme(axis.text = element_text(size=12), axis.title = element_text(size=14),plot.title = element_text(size=22))
      p = p+theme(axis.title.y=element_text(angle=90, vjust=2,size = 18))
      p = p+theme(plot.title = element_text(hjust = 0.5))
      p = p+theme(text = element_text(size=14))
      p = p+theme(axis.text.x = element_text(size = 13))
      p = p+theme(axis.text.y = element_text(size = 16))
      p = p+theme(legend.text = element_text(size = 16))
      p = p+theme(legend.title = element_text(size = 16))
      p = p+scale_fill_manual(values = piecolors)
      p = p+theme(legend.title = element_text(size=16),legend.text = element_text(size=14))
      
      print(p) #Print Stacked Bar Chart
      #p = p +theme_void() #NO LABELS AT ALL
      p = p + coord_polar(theta = "y")
      print(p) #Print Nested Pie Chart
    }else{
      p = ggplot(MetaFG, aes(x = Experiment, y = Percent, fill = FunctionalGroup)) +
        geom_bar(stat = "identity",position = "stack") +
        scale_x_discrete(limits = trueexp) +
        p = p+theme_bw()
        print(p)
        p = p+coord_polar(theta = "y")
        print(p)
    }
    
    
    
    
  }else{
    FG = FunctionalGroups(rownames(mydataforstats),flagambiguous)
    
    nFG = length(names(table(FG)))
    
    piecolors = piecolorfnc(nFG)
    
    if(customFGlist == T){
      pie(table(FG),labels = FGlabels,main = paste("Functional Groups for: ",AnalysisName,sep=""),col = piecolors)
    }else{
      pie(table(FG),main = paste("Functional Groups for: ",AnalysisName,sep=""),col = piecolors)
    }
  }
}


###### Clean Heatmaps (v1.1.1)

CleanHM = finaldata_z[,colnames(finaldata_z) %in% customsamps]
MatrixtoHeatmap2(CleanHM,samplesarecols = T,featuresasrows = T,limits = limits,colors = c(heatcolorlow,heatcolormid,heatcolorhigh),clusterrows = clusterrows,clustercols = clustercols,clustermethod = clustermethod,distancemethod = distancemethod,reverserows = T)