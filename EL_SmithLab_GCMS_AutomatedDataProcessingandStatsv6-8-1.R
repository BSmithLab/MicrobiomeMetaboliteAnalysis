
#Automated Alignment and Analysis of Smith Lab GC-MS data
#Written By: Jarrett Eshima
#Fall 2022 - Summer 2023
#SBHSE, Arizona State University

#Raw and derived data files are publicly available at: https://www.ebi.ac.uk/metabolights/editor/study/MTBLS8293

#Version History
#V2 - Includes p-value results in barplots and allows barplots with fewer than 3 controls, minor bugs if conditions = 2
#V3 - PCA functionality added
#V4 - Allows variable number of positive and negative samples
#V4.1.X *dev version*
#V5 - Allows for a fourth experimental condition. Adds optional functionality to apply scale factors to the third and fourth experimental groups. Added CtrOnlyData and ExOnlyData for use with transformation and stats. Greater control given to scale factor normalization.
#V5.1 - Minor bugs addressed when two groups/conditions are provided. Bug where retention time column is plotted in the heatmap has been fixed. 
#V6 - Greater control given to custom row/column order in heatmap. 3D PCA functionality added. Y-axis label on barplots is adjusted depending on normalization method. 
#V6.1 - Allows user to specify the group pair(s) for t.test(s), no other pairs are tested. Gives user greater control over statistical assumptions and p-value adjustment methods.Two condition analysis has been fully tested and minor bugs addressed. 
#V6.2 - Allows user to provide a subset of VOCs for barplots
#V6.3 - Allows an unlimited number of samples in each experimental condition
#v6.3.1 - Bug fixed related to the combination of broad peaks
#v6.4 - Allows functional group analysis
#v6.4.1 - Allows vector input into 'ClearNoise' start and stop retention times (allows filtering of multiple regions)
#v6.5 - Allows baseline correction
#v6.6 - Allows preparation for companion script (SmithLab_GCMS_AutomatedDataProcessingandStats_Meta_vX-X-X.R). Companion script compares results across independent EXPERIMENTS
#v6.7 - Allows removal / filtering of VOCs found in a control/reference sample
#v6.8 - Fixed minor bugs related to peak filtering using the minimum abundance threshold and peak filtering using a retention time window. Added heatmap functionality to adjust z-score values outside [-3,3] to -3 and 3. Fixed minor bugs related to filtering using a retention time window. Upgraded Filtering of known contaminants. Fixed minor bugs related to the printing of p-values when PlotSigOnly = F.
#v.6.8.1 - Fixed p.adjust bug

#########################################################################################################################################################################
###############################################   READ ME: DATA EXPORT   ################################################################################################
#########################################################################################################################################################################

#Raw and derived data files are publicly available at: https://www.ebi.ac.uk/metabolights/editor/study/MTBLS8293

#Instructions to access data in "Axon Data Analysis" Software:
#1) Open "Axon Data Analysis"
#2) Click "Integrate" --> "Select MS Integrator" --> ChemStation Integrator --> click OK to all text boxes
#3) Click "Integrate" --> "Integrate and label peak areas" --> Wait for integration to finish
#4) Click "Export" --> "Multiple sheets to XLS" --> Leave only the "Library Hit Table" box checked --> Click OK
#5) Wait for excel to finish exporting all peaks and hits - THIS CAN TAKE SEVERAL MINUTES, wait until the excel software closes and reopens automatically
#6) Open the "Library" tab at the bottom of the excel sheet, if step 4 was followed this will be open already
#7) NOTE: You do not need to edit the header or any aspects of the file, this code expects the raw GC-MS output.
#8) Save that sheet as a .CSV file - the file type is IMPORTANT
#9) Click OK to all text boxes and close excel
#10) Repeat for all desired samples

#Instructions to set up a directory to allow proper code usage
#1) Create a CLEAN directory/project folder (no existing/current files)
#2) Copy ALL .CSV files (library hit tables) from your desired experiment to this folder
#3) Create a new folder within the current project folder called "RTs" - THIS IS CASE SENSITIVE
#4) Create a new .CSV file following the format detailed in the "Input" section below. This file will be used to combine broad peaks. 
#5) Read this file in using the "Input" section below. If this function is not desired, set combinebroadpeaks = FALSE

#########################################################################################################################################################################
###############################################   READ ME: CODE USAGE   #################################################################################################
#########################################################################################################################################################################

#This code performs the following:
#1) Optional: Combines broad peaks using manual input file
#2) Not Optional: Filters noisy peaks
#3) Optional, Strongly Recommended : Normalizes by internal standard and other metrics
#4) Not Optional: Aligns data by experimental condition and across all samples
#5) Optional, Strongly Recommended: Performs final cleanup (removes known contaminants)
#6) Optional, Strongly Recommended: Normalizes and transforms data
#7) Not Optional: Runs parametrics statistics (advanced users can manually change the statistical test; lines 885-901)
#8) Not Optional: Plots the data as bar plots, automatically adds standard error bars and significance
#9) Optional: Plots data as heatmap - this requires the file UsefulFunctions_v2.R (see below)
#10) Optional: Performs principal component analysis
#11) Optional: Performs analysis of functional groups using established suffix conventions (v6-5-0 and later)
#12) Optional: Prepares workspace for Meta Analysis. Meta analysis is defined as the re-alignment of peaks across multiple experiments. (v6-6-0 and later)

#Meta Analysis:
#This is performed in a separate script. The companion script name is included below:
#(SmithLab_GCMS_AutomatedDataProcessingandStats_Meta_vX-X-X.R)

#From the UsefulFunctions_v2.R script, you MUST load the following functions into your workspace:
#1) zscore()
#2) MatrixtoHeatmap2() - V5 and above requires the upgraded function: includes row/column clustering functionality
#3) PREP2.0()
#4) PQN()
#5 - version 6-4-0 and later) FunctionalGroups()

#Important Notes:
#1) This code will function properly if you have less than 3 control samples, including if no control samples are available.
#1a) This is not the case for the positive, negative, and fourth experimental conditions. You MUST have three OR MORE replicates of each (code will technically function "correctly" with 2).
#1b) This code does not currently support more than 4 experimental conditions
#2) The "Experimental Conditions" input below are file handles that are repeated in each file name. There should be three.
#3a) The length of the bar plot code will make it difficult to manually edit the plots for asthetics. It is recommended to download the plot as a .SVG image, then edit the plot legend/titles in inkscape or Adobe Illustrator.
#3b) VERY IMPORTANT: If four experimental groups are specified, p-values and significance bars are NOT added to the barplots, however space is provided to add them manually in a graphic editor (Adobe Illustrator or Inkscape)
#4) .CSV files should contain the following somewhere in the name: 1) one of the three file handles, specified in line 61; 2) RX, where X is the replicate number (capital R)
#5) Unused arguments below will not cause issues if loaded into the environment/workspace. Once parameters are set, hit "Source".

######################################################################################################################################################################
################################################  USER INPUT  ########################################################################################################
######################################################################################################################################################################

########################  FILE INFORMATION ####################################

#Get .csv files in a specific directory
wd2 = "C:/Users/jeshima/Documents/Smith Lab/Spring 2023/BMC Microbiology Paper/Data/EL"
setwd(wd2)
files = list.files(wd2,".csv") #make sure ALL files are shown correctly -->

#Experimental Conditions; Allows up to 4 (see note above for more detail - CASE SENSITIVE - For best usage, provide experimental conditions in the following order: Positive, Negative, Control, Optional)
conditions = c("Pos","Neg","Ctr","AFM") #UNIQUE FILE HANDLE MUST BE EXACTLY 3 CHARACTERS LONG
conditions2 = c("Pos","Neg","Ctr","ExD") #You can rearrange these but the string cannot change
#Note: conditions and conditions2 should have the same number of elements / experimental condition names

#Default conditions2 with 4 experimental groups
#conditions2 = c("Pos","Neg","Ctr","ExD") #You can rearrange these but the string cannot change


## Number of Replicates (V6.2.0 max=3; V6.2.1 max=5; unlimited on 6.3.0+)
npos=3 #Min 2 
nneg=3 #Min 2 
#If 3 conditions are provided:
ncontrol=3 #No minimum
#If 4 conditions are provided:
nextra=3 #If specified, Min 2

#Organism/Experiment Name (without spaces)
organism = "ELenta"

#Write Files?
writef = FALSE

########################  FILTERING AND ALIGNMENT ####################################

#Minimum integrated peak signal
minthreshold = 1000000

#Filter before or after coming peaks:
FilterBefore = T

#Retention time tolerance for aligning samples within an experimental condition (i.e. all positive samples)
RTtol = 0.01 # 0.6 second variability between samples

#Retention time tolerance for aligning experimental conditions (positive, negative, and control)
RTtol2 = 0.01

#Apply baseline correction?
Baselinecorrection = F

#If TRUE, specify the threshold (peaks above this threshold are moved and the remaining values are averaged to define the baseline)
baselinethreshold = minthreshold*0.75 #it is recommended to select a value between 0.50 and 0.90 

#Remove features/VOCs seen in a reference level?
FilterUsingControls = F
#If TRUE, specify the reference level (typically "Ctr")
reflevel = "Ctr"

########################  COMBINE BROAD PEAKS  ####################################
combinebroadpeaks = T

##If combinebroadpeaks is TRUE:

#Directory where the Retention Time file is located
if(combinebroadpeaks == T){
  wd = "C:/Users/jeshima/Documents/Smith Lab/Spring 2023/BMC Microbiology Paper/Data/EL/RTs"
  setwd(wd)
  RTind = read.csv("EL_RTs.csv")
}


#RTind should be a file with 2 or 3 columns - "Start" and "Stop" must be there, case sensitive 
# VOC    Start    Stop
# Unknown   1.4    1.55
# Toluene   2.2    2.44
# Etc.

########################  REMOVE KNOWN CONTAMINANTS AND MANUALLY FILTER  ####################################
#Do you want known contaminants (siloxanes) to be removed?
removecontam = T

##Remove VOCs in a specific retention time window
#(v6.4.1 allows multiple filter windows, with high and low RTs paired together in the order they are provided)
ClearNoise = T #If false, the VOCs are not removed.
Remove_LowRT = c(0,3.31,4.48,5.01,5.1,6.02,6.27,6.45,6.727,6.77,7.33,7.9,8.27,8.434,8.46,8.58,8.67,9.305,9.4,9.58,10.1,10.565,10.67,10.83,10.95,11.03,11.144,11.22,11.3,11.43,11.65,12.01) #Not needed, if ClearNoise = False
Remove_HighRT = c(3.115,3.53,4.95,5.07,5.31,6.14,6.41,6.66,6.75,6.78,7.61,7.915,8.406,8.436,8.49,8.595,8.68,9.308,9.41,9.75,10.39,10.635,10.68,10.902,10.99,11.083,11.148,11.253,11.33,11.58,11.94,15.6) #Not needed, if ClearNoise = False

########################  Use Normalization Factors  ####################################
#Factors can be any values recorded for each sample/experiment to normalize (leave the name as OD600, as this is arbitrary)

usenormfactors = T

#If usenormfactors is TRUE (groups 1 and 2):
OD600_pos = c(0.040,0.056,0.023)
OD600_neg = c(0.033,0.147,0.031)

#If usecontrolfactors is TRUE (group 3): 
usecontrolfactors = F
OD600_ctr = c(0,0,0)

#If useexfactors is TRUE (group 4):
useexfactors = T
OD600_ex = c(0.021,0.028,0.021)

#Should you divide through by the factors or use them to scale by group? If true, each sample is divided by its factor.
dividefactors = F #If true, this divides the data by the normalization factor directly. Not recommended unless normalization factors exist for all groups being considered.

#If dividebyfactors is FALSE you must specify the scalebygroup argument
#Do you want to use the normalization factors to scale within each group? 
#For example, if you had two groups A and B with norm factors (1,1,0.2) and (2,2,2), and scalebygroup is true, 
#then factors become (1,1,5) and (1,1,1). If scalebygroup is false, then the factors become (2,2,10) and (1,1,1) - i.e. scaled globally
scalebygroup = T #True reflects the original coding in versions 4 and previous

########################  Apply Internal Standard Correction  ####################################
#Internal Standard Abundances

useIS = F

#If useIS is TRUE:
IS_Abundance_pos = c(743104,743104,743104)
IS_Abundance_neg = c(743104,743104,743104)
IS_Abundance_ctr = c(743104,743104,743104)

#If the number of conditions = 4:
IS_Abundance_ex = c(743104,743104,743104)

#Should you divide through by the abundance of the internal standard peak? If true, each peak abundance is divided by the abundance of the IS.
divideIS = F

#If dividebyfactors is FALSE you must specify the ISscalebygroup argument. See (norm factor) description above for more information.
ISscalebygroup = F #True reflects the original coding in versions 4 and previous


########################  Would you like to filter out VOCs with a low number of observations?  ####################################
filterlowobs = F
Xmissing = 3 #This is total missing observations, not missing observation per experimental condition.
###################################################################################

########################  SELECT DATA FOR NORMALIZATION & STATISTICS  ####################################

#Possible Dataframe Names:
#1) FinalData - raw data after filtering and normalization
#2) PosNegData - VOCs found in all positive (condition 1) and negative (condition 2) samples; no transformation/PQN
#3) PosOnlyData - VOCs found in ALL positive (condition 1) samples and ZERO negative (condition 2) samples
#4) NegOnlyData - VOCs found in ALL negative (condition 2) samples and ZERO positive (condition 1) samples
#5) AlwaysObservedData - VOCs observed in ALL SAMPLES.
#6) ThreshObservedData - VOCs found in all samples, up to a threshold X. This is total missing observations, not missing observation per experimental condition.
#7) TransformedData - data after performing the desired transformation above

#### V5-0-0 and later (make sure to read lines 188-190 as this is an important argument in later versions) ####
#8) CtrOnlyData - VOCs found in ALL control (condition 3) samples and ZERO positive,negative,and group 4 samples
#9) ExOnlyData (number of conditions must = 4) - VOCs found in ALL condition 4 samples and ZERO positive,negative,and control samples
#3) PosOnlyData - Redefined if maintainposnegonly = F: VOCs found in ALL positive (condition 1) samples and ZERO negative, control, and group 4 samples
#4) NegOnlyData - Redefined if maintainposnegonly = F: VOCs found in ALL negative (condition 2) samples and ZERO positive, control, and group 4 samples

#Dataframe for normalization
normdata = "FinalData"
#Dataframe for statistics
statdata = "TransformedData"

#IMPORTANT NOTE: statdata should always equal "TransformedData" if RunNorm = T and you wish to use the normalized data for statistics. 
#If RunNorm = F, then 1-6,8,9 SHOULD be selected for statdata
#IMPORTANT NOTE: the name should be entered as a string (quotation marks)

#If you would only like to compare the positive and negative conditions when identifying unique VOCs, set the argument below to TRUE
#If maintainposnegonly = T, then this runs the original code from V3 and before (this ignores VOCs that may be present in group 4)
maintainposnegonly = T #Setting this to false allows identification of control and group 4 only VOCs. THIS MUST BE TRUE IF YOU HAVE 2 CONDITIONS/GROUPS.

#Should the V4 and earlier versions of PosNegData be built (comparing only the positive and negative conditions)? This can be ignored if you are using 3 conditions.
#If this is set to false, the PosNegData matrix will contain VOCs seen in all Positive, Negative, and Extra (group 4) samples
originalposnegdata = T #Only relevant if 4 conditions are specified

#You can also manually change this on line 1680

###################################################################################

########################  NORMALIZATION METHOD  ####################################
#Should normalization be performed?
RunNorm = T

#Which normalization method? Select only one (current functionality does not allow more than one)
#Note: transformations/normalization performed sequentially, starting from the left.

#Options:
#1) No Transformation (same as RunNorm = F)
#2) PQN
#3) PQN + log10 transform
#4) PQN + log10 transform + centering
#5) PQN + log10 transform + centering + scaling
#6) log10 transformation only
#7) centering only
#8) scaling only
#9) zscore normalize
#10) log10 + centering
#11) log10 + centering + scaling
#12) centering + scaling only
#13) Impute missing values 
#14) Impute missing values + log10 transformation (careful with 0s)
#15) Impute missing values + log10 transformation + centering
#16) Impute missing values + log10 transformation + centering + scaling
#17) PQN + zscore normalize
#18) PQN + log10 transform + zscore normalize
#19) PQN + log10 transform + centering + zscore normalize
#20) PQN + log10 transform + centering + scaling + zscore normalize
#21) log10 transform + PQN

NormMethod = 6

#If methods 13-16 are selected, please specify the desired imputation value:

#imputevalueto = 0

################################################################################################

#################################  Target VOCs to Plot  #################################################

#Do you want a subset of VOCs plotted? If True, only the names provided in 'VOClist' are considered during plotting
SpecificVOCs = F

#List of VOCs to plot. IMPORTANT: Names are CASE SENSITIVE. Copy the name exactly, including symbols
#MINIMUM 2 VOCs
VOClist = c("Isopentyl hexanoate","Hexanoic acid, 2-methylbutyl ester")


#################################  STATISTICS  #################################################

#Hypothesis test tails (Options - case sensitive - "less","greater","two.sided")
#Less and greater options are one-tailed t-tests
tails = "two.sided"
#Do you want to assume equal variance? (unequal variance is often a safer assumption)
eqvar = F
#Is the data paired? (TRUE if same sample/individual/experiment at two different time points)
paired = F

#Multiple hypothesis test correction method (see ?p.adjust.methods):
MHCmethod = "BH"


#User specified statistics. If uselevels = F, then a student's t-test is applied for all pairwise comparisons. The FDR/Benjamini-Hochberg method is used for multiple hypothesis test corrections.
uselevels = F

#If uselevels = FALSE, do you want to write out the adjusted pvalue table? v6.8.1 and later
writeadjustedp = T

#If uselevels = TRUE:
#Provide one or more pairs for statistical testing. If you have two or more pairs, specify them in order. 
#Important: These should be the exact strings provided in "conditions" (line 76)
#Example: c("Pos","Neg","Pos","Ctr") will compare Pos vs Neg AND Pos vs Control.
pairs = c("Pos","Neg")
#pairs = c("Pos","Neg","Pos","Ctr")

################################################################################################

#################################  PLOTTING  #################################################

## Parameters for Bar Plots ##

PlotSigOnly = T #Should only the significant genes be plotted (after p-value adjustment for multiple hypothesis testing)

#Note if PlotSigOnly = F, you will need to specify up to 30 rows to plot at a time.
startrow = 1
stoprow = 5

AllowUnadjusted = F #Should the unadjusted p-values be considered/used? This might be useful to run if VOCs are weakly significant.
sigthresh = 0.1 #Significance threshold
addpvalues = T #Should p-values be shown for VOCs that do NOT pass the significance threshold?
#Note: For 3 conditions, this works as expected. For 4 conditions, this provides the necessary plotting space to manually add significance in a graphic editor (Illustrator/Inkscape)
givemepspace = T

#Barplot Colors:
poscolor = "coral"
negcolor = "deepskyblue2"
ctrcolor = "grey60"
excolor = "mediumaquamarine" #This is only needed if 4 conditions are provided
lowerlim = -1 #manually adjust lower limit


## Heatmap Parameters ##

#Plot the Heatmap?
PlotHeatmap = T

#Optional arguments (defaults are provided below):
heatcolorlow = "deepskyblue2"
heatcolormid = "white"
heatcolorhigh = "darkorange2"
limits = c(-3,3)
clusterrows = T
clustercols = T
clustermethod = "ward.D2"
distancemethod = "manhattan"
#A custom row/column order can be provided. This will overwrite clustering results if used. 
#With 3 groups
#customsamps = c("Replicate1_Pos","Replicate2_Pos","Replicate3_Pos","Replicate1_Neg","Replicate2_Neg","Replicate3_Neg","Replicate1_Ctr","Replicate2_Ctr","Replicate3_Ctr") #Column order for plotting, starting with the left side
#With 4 groups
usecustomsamps = T
customsamps = c("Replicate1_Pos","Replicate2_Pos","Replicate3_Pos","Replicate1_Ext","Replicate2_Ext","Replicate3_Ext","Replicate1_Neg","Replicate2_Neg","Replicate3_Neg","Replicate1_Ctr","Replicate2_Ctr","Replicate3_Ctr") #Column order for plotting with 4 conditions, starting with the left side
#customsamps = c("Pos1","Pos2","Pos3","Inhib1","Inhib2","Inhib3","Neg1","Neg2","Neg3","Ctr1","Ctr2","Ctr3")
usecustomfeats = F
customfeats = c("VOC1","VOC2","VOC3","etc.") #This pipes well into external clustering methods, although the updated code 

#Adjust z-score values greater than 3 to 3 and less than -3 to -3 (for plotting purposes)
maxz3 = T

#Defaults:
# customsamps = c("Replicate1_Pos","Replicate2_Pos","Replicate3_Pos","Replicate1_Neg","Replicate2_Neg","Replicate3_Neg","Replicate1_Ctr","Replicate2_Ctr","Replicate3_Ctr") #Column order for plotting, starting with the left side
# customsamps = c("Replicate1_Pos","Replicate2_Pos","Replicate3_Pos","Replicate1_Ext","Replicate2_Ext","Replicate3_Ext","Replicate1_Neg","Replicate2_Neg","Replicate3_Neg","Replicate1_Ctr","Replicate2_Ctr","Replicate3_Ctr") #Column order for plotting with 4 conditions, starting with the left side
# heatcolorlow = "blue"
# heatcolormid = "white"
# heatcolorhigh = "red"
# limits = c(-3,3)
# clusterrows = F
# clustercols = F
# clustermethod = "ward.D"
# distancemethod = "euclidean"


## Principal Component Analysis ##
RunPCA = T

#Which principal components do you want to plot? Typically, the first and second components capture the greatest amount of variability in your data.
comp1 = 1
comp2 = 2

#Plot 3-axis (3D) PCA?
pca3d = F
#If TRUE:
comp3 = 3

#Colors for plotting groups?
pcacolorpos = "coral"
pcacolorneg = "deepskyblue"
pcacolorctr = "grey60"
pcacolorex = "mediumaquamarine" #This is only needed if 4 conditions are provided
  


## Analysis of Functional Groups ##

library(grDevices)

#Should the VOC functional groups be assessed and plotted as a pie chart?
AnalyzeFunctionalGroups = T

#If TRUE,
#Should VOCs with an ambiguous functional group (using substring searching) be assigned as unknown? If yes, set the variable below to TRUE
flagambiguous = F

#If TRUE,
#Do you want to overwrite the functional group names assigned? 
customFGlist = F

#IMPORTANT: Specify the functional group names below as a list. You MUST provide a list that is equal to the number of unique functional groups. You can check this easily after running the code once with: length(names(table(FunctionalGroups(rownames(data)))))
FGlabels = c("aromatic hydrocarbon","nitrogenated ring","etc.","etc.")

#Colors for the piechart
color1 = "deepskyblue1"
color2 = "coral"
piecolorfnc = colorRampPalette(c(color1,color2))

########################################  PREP FOR META ANALYSIS  ######################################################

#Do you want to prepare the workspace and directory for meta analysis?
#This will output the matrix used for statistics to a new working directory. 
#Note: This analysis is only relevant if you would like to compare results across different experiments (multiple sets of positives,negatives,and controls)

Prep4Meta = T

MetaDirectory = "C:/Users/jeshima/Documents/Smith Lab/Spring 2023/BMC Microbiology Paper/Data"
#Important: This directory should not change as you execute multiple instances of this script (one script execution for each experiment)

#######################################################################################################################################
#######################################  END OF USER INPUT ############################################################################
#######################################                    ############################################################################
#######################################   **HIT SOURCE**   ############################################################################
#######################################################################################################################################

#Raw and derived data files are publicly available at: https://www.ebi.ac.uk/metabolights/editor/study/MTBLS8293

#########################################################  Load Data  #################################################################

setwd(wd2)

#Read the files in
npeaks = rep(NA,length(files))
for(i in 1:length(files)){
  
  tmpname = paste(files[i])
  tmp = read.csv(files[i])[-1:-4,]
  colnames(tmp) = tmp[1,]
  tmp = tmp[-1,]
  assign(tmpname,tmp)
  
  #Get the peak number
  tmpvec = unlist(tmp$`Compound number (#)`)
  remind = rep(NA,length(tmpvec))
  for(j in 1:length(tmpvec)){
    
    if(tmpvec[j] == ""){
      remind[j] = j
    }
    
  }
  remind = remind[!is.na(remind)]
  npeaks[i] = max(as.numeric(tmpvec[-remind]))
  
}

#Condense the hit table and build data frames
cleannames = gsub("\\.csv","",files)
for(i in 1:length(files)){
  
  data = get(files[i])
  
  cleanHT = data[which(data$`Hit Number` == 1),]
  
  tmp = cbind(cleanHT$`Hit Name`,cleanHT$`RT (min)`,cleanHT$`Area (Ab*s)`)
  tmp = data.frame(tmp)
  colnames(tmp) = c("VOC","RT","Abundance")
  
  assign(cleannames[i],tmp)
}

files2 = cleannames

#Reproduce the chromatograms - you may have to adjust the y axis limits
for(i in 1:length(files)){
  
  tmp = get(cleannames[i])
  plot.ts(seq(1:nrow(tmp)),tmp$Abundance,ylim = c(0,1E8),type = "l",main = paste(cleannames[i],"Chromatogram",sep=" "))
  
  #Define the baseline
  if(Baselinecorrection == T){
    tmp$Abundance = as.numeric(tmp$Abundance)
    tmp2 = tmp[-which(tmp$Abundance > baselinethreshold),]
    baseline = mean(tmp2$Abundance,na.rm=T)
    abline(h=baseline,col="red",lty=2)
    cat(paste("Baseline for ",cleannames[i]," = ",baseline,sep=""))
    
    tmp$Abundance = tmp$Abundance - baseline
    assign(cleannames[i],tmp)
  }
  
}


#####################################################  Pre-Processing  ################################################################

pnames = paste(files2,"_Processed",sep="")
for(i in 1:length(files2)){
  
  myfile = get(files2[i])
  
  #Filter peaks less than 1 million abundance
  if(FilterBefore == T){
    count = 1
    index = rep(NA,nrow(myfile))
    for(k in 1:nrow(myfile)){
      
      if(as.numeric(myfile$Abundance[k]) > minthreshold){
        index[count] = k
        count = count+1
      }
      
    }
    
    index = index[!is.na(index)]
    
    cleandata = myfile[index,]
    
    if(combinebroadpeaks == T){
      #Empty data frame 
      cleandata2 = data.frame(matrix(NA,nrow=nrow(RTind),ncol = ncol(cleandata)))
      colnames(cleandata2) = colnames(cleandata)
      cleandata2$VOC = RTind$VOC
      for(z in 1:nrow(RTind)){
        cleandata2$RT[z] = mean(c(RTind$Start[z],RTind$Stop[z]))
      }
      
      
      #Sum up peaks in a certain retention time window
      for(j in 1:nrow(RTind)){
        
        index1 = which(cleandata$RT >= RTind$Start[j])
        index2 = which(cleandata$RT <= RTind$Stop[j])
        finalindex = index1[index1 %in% index2]
        
        SumAb = sum(as.numeric(cleandata$Abundance[finalindex]),na.rm = TRUE)
        
        cleandata2$Abundance[j] = SumAb
        
        #cleandata$Start.Time..min.[i] = RTind$Start[i]
        #cleandata$End.Time..min.[i] = RTind$Stop[i]
        
        if(j == 1){
          runningindex = finalindex
        }else{
          runningindex = c(runningindex,finalindex)
        }
        
      }
      
      runningindex = as.numeric(names(table(runningindex)))
      cleandata = cleandata[-runningindex,]
      
    }
    
    if(combinebroadpeaks == T){
      ProcessedData = rbind(cleandata2,cleandata)
    }else{
      ProcessedData = cleandata
    }
    
    rownames(ProcessedData) = ProcessedData$Hit.Name
    
    assign(pnames[i],ProcessedData)
    
    if((i %% 1) == 0) cat("Broad Peaks Combined and Peak Filtering Completed for file:",files2[i],"\n")
  }else{
    myfile = get(files2[i])
    
    if(combinebroadpeaks == T){
      #Empty data frame 
      cleandata = data.frame(matrix(NA,nrow=nrow(RTind),ncol = ncol(myfile)))
      colnames(cleandata) = colnames(myfile)
      cleandata$VOC = RTind$VOC
      for(z in 1:nrow(RTind)){
        cleandata$RT[z] = mean(c(RTind$Start[z],RTind$Stop[z]))
      }
      
      
      #Sum up peaks in a certain retention time window
      for(j in 1:nrow(RTind)){
        
        index1 = which(myfile$RT >= RTind$Start[j])
        index2 = which(myfile$RT <= RTind$Stop[j])
        finalindex = index1[index1 %in% index2]
        
        SumAb = sum(as.numeric(myfile$Abundance[finalindex]),na.rm = TRUE)
        
        cleandata$Abundance[j] = SumAb
        
        #cleandata$Start.Time..min.[i] = RTind$Start[i]
        #cleandata$End.Time..min.[i] = RTind$Stop[i]
        
        if(j == 1){
          runningindex = finalindex
        }else{
          runningindex = c(runningindex,finalindex)
        }
        
      }
      
      runningindex = as.numeric(names(table(runningindex)))
      filterfile = myfile[-runningindex,]
      
    }else{
      filterfile = myfile
    }
    
    #Filter peaks less than 1 million abundance
    count = 1
    index = rep(NA,nrow(filterfile))
    for(k in 1:nrow(filterfile)){
      
      if(as.numeric(filterfile$Abundance[k]) > minthreshold){
        index[count] = k
        count = count+1
      }
      
    }
    
    index = index[!is.na(index)]
    
    cleandata2 = filterfile[index,]
    
    if(combinebroadpeaks == T){
      ProcessedData = rbind(cleandata,cleandata2)
    }else{
      ProcessedData = cleandata2
    }
    
    rownames(ProcessedData) = ProcessedData$Hit.Name
    
    assign(pnames[i],ProcessedData)
    
    if((i %% 1) == 0) cat("Broad Peaks Combined and Peak Filtering Completed for file:",files2[i],"\n")
  }
}


################################################# Replicate Alignment ##############################################################

#Indices

pnames = paste(files2,"_Processed",sep="")
count = count2 = count3 = count4 = 1
PosInd = NegInd = CtrInd = ExInd = rep(NA,length(pnames))
for(i in 1:length(pnames)){
  
  filen = pnames[i]
  
  for(j in 1:nchar(filen)){
    
    if(substr(filen,j,j+2) == conditions[1]){
      
      PosInd[count] = i
      count = count+1
      
    }
    
    if(substr(filen,j,j+2) == conditions[2]){
      
      NegInd[count2] = i
      count2 = count2+1
      
    }
    
    if(length(conditions) == 3){
      if(substr(filen,j,j+2) == conditions[3]){
        
        CtrInd[count3] = i
        count3 = count3+1
        
      }
    }
    
    
    if(length(conditions) == 4){
      
      if(substr(filen,j,j+2) == conditions[3]){
        
        CtrInd[count3] = i
        count3 = count3+1
        
      }
      
      if(substr(filen,j,j+2) == conditions[4]){
        
        ExInd[count4] = i
        count4 = count4+1
        
      }
      
    }
    
    
  }
  
}
PosInd = PosInd[!is.na(PosInd)]
NegInd = NegInd[!is.na(NegInd)]
if(length(conditions) == 3){
  CtrInd = CtrInd[!is.na(CtrInd)]
}
if(length(conditions) == 4){
  CtrInd = CtrInd[!is.na(CtrInd)]
  ExInd = ExInd[!is.na(ExInd)]
}


#Positive Dataframe
Posfiles = pnames[PosInd]
nVOCs_p = rep(NA,length(Posfiles))
for(i in 1:length(Posfiles)){
  
  nVOCs_p[i] = nrow(get(Posfiles[i]))
  
}
maxind = which(nVOCs_p == max(nVOCs_p))
if(length(maxind)>1){
  maxind = maxind[1]
}
tmp = get(Posfiles[maxind])
PosRTs = as.numeric(tmp$RT)
PosAb = as.numeric(tmp$Abundance)
PosData = data.frame(matrix(NA,nrow=max(nVOCs_p),ncol=npos))
PosData[,maxind] = PosAb
PosData = cbind(PosData,PosRTs)
rownames(PosData) = make.unique(get(Posfiles[maxind])$VOC)
colnames(PosData) = c(paste("Replicate",seq(1,npos),"_Pos",sep=""),"RT")
fillPosfiles = Posfiles[-maxind]



#Negative Dataframe
Negfiles = pnames[NegInd]
nVOCs_n = rep(NA,length(Negfiles))
for(i in 1:length(Negfiles)){
  
  nVOCs_n[i] = nrow(get(Negfiles[i]))
  
}
maxind = which(nVOCs_n == max(nVOCs_n))
if(length(maxind)>1){
  maxind = maxind[1]
}
tmp = get(Negfiles[maxind])
NegRTs = as.numeric(tmp$RT)
NegAb = as.numeric(tmp$Abundance)
NegData = data.frame(matrix(NA,nrow=max(nVOCs_n),ncol=nneg))
NegData[,maxind] = NegAb
NegData = cbind(NegData,NegRTs)
rownames(NegData) = make.unique(get(Negfiles[maxind])$VOC)
colnames(NegData) = c(paste("Replicate",seq(1,nneg),"_Neg",sep=""),"RT")
fillNegfiles = Negfiles[-maxind]



#Control Dataframe
if(length(conditions) == 3){
  Ctrfiles = pnames[CtrInd]
  nVOCs_c = rep(NA,length(Ctrfiles))
  for(i in 1:length(Ctrfiles)){
    
    nVOCs_c[i] = nrow(get(Ctrfiles[i]))
    
  }
  maxind = which(nVOCs_c == max(nVOCs_c))
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  tmp = get(Ctrfiles[maxind])
  CtrRTs = as.numeric(tmp$RT)
  CtrAb = as.numeric(tmp$Abundance)
  CtrData = data.frame(matrix(NA,nrow=max(nVOCs_c),ncol=ncontrol))
  CtrData[,maxind] = CtrAb
  CtrData = cbind(CtrData,CtrRTs)
  rownames(CtrData) = make.unique(get(Ctrfiles[maxind])$VOC)
  colnames(CtrData) = c(paste("Replicate",seq(1,ncontrol),"_Ctr",sep=""),"RT")
  fillCtrfiles = Ctrfiles[-maxind]
}


#Optional Experimental Condition Dataframe
if(length(conditions) == 4){
  
  Ctrfiles = pnames[CtrInd]
  nVOCs_c = rep(NA,length(Ctrfiles))
  for(i in 1:length(Ctrfiles)){
    
    nVOCs_c[i] = nrow(get(Ctrfiles[i]))
    
  }
  maxind = which(nVOCs_c == max(nVOCs_c))
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  tmp = get(Ctrfiles[maxind])
  CtrRTs = as.numeric(tmp$RT)
  CtrAb = as.numeric(tmp$Abundance)
  CtrData = data.frame(matrix(NA,nrow=max(nVOCs_c),ncol=ncontrol))
  CtrData[,maxind] = CtrAb
  CtrData = cbind(CtrData,CtrRTs)
  rownames(CtrData) = make.unique(get(Ctrfiles[maxind])$VOC)
  colnames(CtrData) = c(paste("Replicate",seq(1,ncontrol),"_Ctr",sep=""),"RT")
  fillCtrfiles = Ctrfiles[-maxind]
  
  
  Exfiles = pnames[ExInd]
  nVOCs_e = rep(NA,length(Exfiles))
  for(i in 1:length(Exfiles)){
    
    nVOCs_e[i] = nrow(get(Exfiles[i]))
    
  }
  maxind = which(nVOCs_e == max(nVOCs_e))
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  tmp = get(Exfiles[maxind])
  ExRTs = as.numeric(tmp$RT)
  ExAb = as.numeric(tmp$Abundance)
  ExData = data.frame(matrix(NA,nrow=max(nVOCs_e),ncol=nextra))
  ExData[,maxind] = ExAb
  ExData = cbind(ExData,ExRTs)
  rownames(ExData) = make.unique(get(Exfiles[maxind])$VOC)
  colnames(ExData) = c(paste("Replicate",seq(1,nextra),"_Ext",sep=""),"RT")
  fillExfiles = Exfiles[-maxind]
  
}


#Fill Positive Data Frame
for(i in 1:length(fillPosfiles)){
  
  data = get(fillPosfiles[i])
  
  fillcolumn = as.numeric(gsub(".*?([0-9]+).*", "\\1", fillPosfiles[i]))
  
  for(j in 1:nrow(data)){
    for(k in 1:nrow(PosData)){
     
      if(as.numeric(data$RT[j]) > PosData$RT[k] - RTtol && as.numeric(data$RT[j]) < PosData$RT[k] + RTtol ){
        
        PosData[k,fillcolumn] = data$Abundance[j]
        
      } 
      
    }
    
  }
  
  if((i %% 1) == 0) cat("File Alignment Completed:",fillPosfiles[i],"\n")
  
}


#Fill Negative Data Frame
for(i in 1:length(fillNegfiles)){
  
  data = get(fillNegfiles[i])
  fillcolumn = as.numeric(gsub(".*?([0-9]+).*", "\\1", fillNegfiles[i]))
  
  for(j in 1:nrow(data)){
    for(k in 1:nrow(NegData)){
      
      if(as.numeric(data$RT[j]) > NegData$RT[k] - RTtol && as.numeric(data$RT[j]) < NegData$RT[k] + RTtol ){
        
        NegData[k,fillcolumn] = data$Abundance[j]
        
      } 
      
    }
    
  }
  
  if((i %% 1) == 0) cat("File Alignment Completed:",fillNegfiles[i],"\n")
  
}


#Fill Control Data Frame
if(length(conditions) == 3){
  for(i in 1:length(fillCtrfiles)){
    
    if(length(fillCtrfiles) != 0){
      data = get(fillCtrfiles[i])
      fillcolumn = as.numeric(gsub(".*?([0-9]+).*", "\\1", fillCtrfiles[i]))
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(CtrData)){
          
          if(data$RT[j] > CtrData$RT[k] - RTtol && data$RT[j] < CtrData$RT[k] + RTtol ){
            
            CtrData[k,fillcolumn] = data$Abundance[j]
            
          } 
          
        }
        
      }
      
      if((i %% 1) == 0) cat("File Alignment Completed:",fillCtrfiles[i],"\n")
    }
  }
}



#Fill the optional experimental condition Data Frame
if(length(conditions) == 4){
  
  for(i in 1:length(fillCtrfiles)){
    
    if(length(fillCtrfiles) != 0){
      data = get(fillCtrfiles[i])
      fillcolumn = as.numeric(gsub(".*?([0-9]+).*", "\\1", fillCtrfiles[i]))
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(CtrData)){
          
          if(data$RT[j] > CtrData$RT[k] - RTtol && data$RT[j] < CtrData$RT[k] + RTtol ){
            
            CtrData[k,fillcolumn] = data$Abundance[j]
            
          } 
          
        }
        
      }
      
      if((i %% 1) == 0) cat("File Alignment Completed:",fillCtrfiles[i],"\n")
    }
  }
  
  
  for(i in 1:length(fillExfiles)){
    
    if(length(fillExfiles) != 0){
      data = get(fillExfiles[i])
      fillcolumn = as.numeric(gsub(".*?([0-9]+).*", "\\1", fillExfiles[i]))
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(ExData)){
          
          if(data$RT[j] > ExData$RT[k] - RTtol && data$RT[j] < ExData$RT[k] + RTtol ){
            
            ExData[k,fillcolumn] = data$Abundance[j]
            
          } 
          
        }
        
      }
      
      if((i %% 1) == 0) cat("File Alignment Completed:",fillExfiles[i],"\n")
    }
  }
}


################################################  Abundance Normalization  ###########################################################
#Levodopa study note: This could also be tried after z-score normalization, but then the scale factors should also be dependent on the z-scores... I think

############################  OD600  ###############################

if(usenormfactors == T){
  
  if(dividefactors == T){
    #Positive OD600 normalization
    for(i in 1:nrow(PosData)){
      ind = which(substr(colnames(PosData),1,3)=="Rep")
      PosData[i,ind] = as.numeric(PosData[i,ind]) / OD600_pos[ind] #This works correctly
      
    }
    
    #Negative OD600 normalization
    for(i in 1:nrow(NegData)){
      ind = which(substr(colnames(NegData),1,3)=="Rep")
      NegData[i,ind] = as.numeric(NegData[i,ind]) / OD600_neg[ind]
      
    }
    
    #Control OD600 normalization
    if(usecontrolfactors == T){
      for(i in 1:nrow(CtrData)){
        ind = which(substr(colnames(CtrData),1,3)=="Rep")
        CtrData[i,ind] = as.numeric(CtrData[i,ind]) / OD600_ctr[ind]
      }
    }
    
    if(length(conditions) == 4){
      if(useexfactors == T){
        for(i in 1:nrow(ExData)){
          ind = which(substr(colnames(ExData),1,3)=="Rep")
          ExData[i,ind] = as.numeric(ExData[i,ind]) / OD600_ex[ind]
        }
      }
    }
    
  }else{
    
    if(scalebygroup == T){
      #Scale factors - Positive 
      maxpos = max(OD600_pos)
      scalefactor_pos = rep(NA,length(OD600_pos))
      for(i in 1:length(OD600_pos)){
        
        scalefactor_pos[i] = maxpos/OD600_pos[i]
        
      }
      
      #Scale factors - Negative 
      maxneg = max(OD600_neg)
      scalefactor_neg = rep(NA,length(OD600_neg))
      for(i in 1:length(OD600_neg)){
        
        scalefactor_neg[i] = maxneg/OD600_neg[i]
        
      }
      
      #Scale factors - Controls
      
      if(usecontrolfactors == T){
        maxctr = max(OD600_ctr)
        scalefactor_ctr = rep(NA,length(OD600_ctr))
        for(i in 1:length(OD600_ctr)){
          
          scalefactor_ctr[i] = maxctr/OD600_ctr[i]
          
        }
      }
      
      if(length(conditions) == 4){
        if(useexfactors == T){
          
          maxex = max(OD600_ex)
          scalefactor_ex = rep(NA,length(OD600_ex))
          for(i in 1:length(OD600_ex)){
            
            scalefactor_ex[i] = maxex/OD600_ex[i]
            
          }
          
        }
      }
      
      #Positive OD600 normalization
      for(i in 1:nrow(PosData)){
        ind = which(substr(colnames(PosData),1,3)=="Rep")
        PosData[i,ind] = as.numeric(PosData[i,ind]) * scalefactor_pos[ind] #This works correctly
        
      }
      
      #Negative OD600 normalization
      for(i in 1:nrow(NegData)){
        ind = which(substr(colnames(NegData),1,3)=="Rep")
        NegData[i,ind] = as.numeric(NegData[i,ind]) * scalefactor_neg[ind]
        
      }
      
      #Control OD600 normalization
      if(usecontrolfactors == T){
        for(i in 1:nrow(CtrData)){
          ind = which(substr(colnames(CtrData),1,3)=="Rep")
          CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_ctr[ind]
        }
      }
      
      if(length(conditions) == 4){
        if(useexfactors == T){
          for(i in 1:nrow(ExData)){
            ind = which(substr(colnames(ExData),1,3)=="Rep")
            ExData[i,ind] = as.numeric(ExData[i,ind]) * scalefactor_ex[ind]
          }
        }
      }
    }else{
      
      #Scale factors 
      
      if(length(conditions) <= 3){
        if(usecontrolfactors == T){
          maxfactor = max(c(OD600_pos,OD600_neg,OD600_ctr))
        }else{
          maxfactor = max(c(OD600_pos,OD600_neg))
        }
      }else if(length(conditions) == 4){
        
        maxfactor = max(c(OD600_pos,OD600_neg))
        
        if(usecontrolfactors == T){
          if(useexfactors == T){
            maxfactor = max(c(OD600_pos,OD600_neg,OD600_ctr,OD600_ex))
          }else{
            maxfactor = max(c(OD600_pos,OD600_neg,OD600_ctr))
          }
        }else{
          if(useexfactors == T){
            maxfactor = max(c(OD600_pos,OD600_neg,OD600_ex))
          }
        }
        
      }
      
      scalefactor_pos = rep(NA,length(OD600_pos));scalefactor_neg = rep(NA,length(OD600_neg))
      if(usecontrolfactors == T){
        scalefactor_ctr = rep(NA,length(OD600_ctr))
        for(i in 1:length(OD600_ctr)){
          scalefactor_ctr[i] = maxfactor/OD600_ctr[i]
        }
      }
      if(length(conditions) == 4){
        if(useexfactors == T){
          scalefactor_ex = rep(NA,length(OD600_ex))
          for(i in 1:length(OD600_ex)){
            scalefactor_ex[i] = maxfactor/OD600_ex[i]
          }
        }
      }
      for(i in 1:length(OD600_pos)){
        scalefactor_pos[i] = maxfactor/OD600_pos[i]
      }
      for(i in 1:length(OD600_neg)){
        scalefactor_neg[i] = maxfactor/OD600_neg[i]
      }
      
      
      #Positive OD600 normalization
      for(i in 1:nrow(PosData)){
        ind = which(substr(colnames(PosData),1,3)=="Rep")
        PosData[i,ind] = as.numeric(PosData[i,ind]) * scalefactor_pos[ind] #This works correctly
        
      }
      
      #Negative OD600 normalization
      for(i in 1:nrow(NegData)){
        ind = which(substr(colnames(NegData),1,3)=="Rep")
        NegData[i,ind] = as.numeric(NegData[i,ind]) * scalefactor_neg[ind]
        
      }
      
      #Control OD600 normalization
      if(usecontrolfactors == T){
        for(i in 1:nrow(CtrData)){
          ind = which(substr(colnames(CtrData),1,3)=="Rep")
          CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_ctr[ind]
        }
      }
      
      if(length(conditions) == 4){
        if(useexfactors == T){
          for(i in 1:nrow(ExData)){
            ind = which(substr(colnames(ExData),1,3)=="Rep")
            ExData[i,ind] = as.numeric(ExData[i,ind]) * scalefactor_ex[ind]
          }
        }
      }
      
    }
    
    
  }
  
}

########################  Internal Standard  ###########################

if(useIS == T){
  
  if(divideIS == T){
    
    #Positive
    for(i in 1:nrow(PosData)){
      ind = which(substr(colnames(PosData),1,3)=="Rep")
      PosData[i,ind] = as.numeric(PosData[i,ind]) / IS_Abundance_pos[ind]
      
    }
    
    #Negative
    for(i in 1:nrow(NegData)){
      ind = which(substr(colnames(NegData),1,3)=="Rep")
      NegData[i,ind] = as.numeric(NegData[i,ind]) / IS_Abundance_neg[ind]
      
    }
    
    #Control
    if(length(conditions) == 3){
      for(i in 1:nrow(CtrData)){
        ind = which(substr(colnames(CtrData),1,3)=="Rep")
        CtrData[i,ind] = as.numeric(CtrData[i,ind]) / IS_Abundance_ctr[ind]
        
      }
    }
    
    #Extra Group
    if(length(conditions) == 4){
      for(i in 1:nrow(CtrData)){
        ind = which(substr(colnames(CtrData),1,3)=="Rep")
        CtrData[i,ind] = as.numeric(CtrData[i,ind]) / IS_Abundance_ctr[ind]
        
      }
      for(i in 1:nrow(ExData)){
        ind = which(substr(colnames(ExData),1,3)=="Rep")
        ExData[i,ind] = as.numeric(ExData[i,ind]) / IS_Abundance_ex[ind]
        
      }
    }
    
  }else{
    
    if(ISscalebygroup == T){
      #Scale factors - Positive Controls
      maxabund_pos = max(IS_Abundance_pos)
      scalefactor_abund_pos = rep(NA,length(IS_Abundance_pos))
      for(i in 1:length(IS_Abundance_pos)){
        
        scalefactor_abund_pos[i] = maxabund_pos/IS_Abundance_pos[i]
        
      }
      
      #Scale factors - Negative Controls
      maxabund_neg = max(IS_Abundance_neg)
      scalefactor_abund_neg = rep(NA,length(IS_Abundance_neg))
      for(i in 1:length(IS_Abundance_neg)){
        
        scalefactor_abund_neg[i] = maxabund_neg/IS_Abundance_neg[i]
        
      }
      
      #Scale factors - Controls
      if(length(conditions) == 3){
        maxabund_ctr = max(IS_Abundance_ctr)
        scalefactor_abund_ctr = rep(NA,length(IS_Abundance_ctr))
        for(i in 1:length(IS_Abundance_ctr)){
          
          scalefactor_abund_ctr[i] = maxabund_ctr/IS_Abundance_ctr[i]
          
        }
      }
      
      #Scale factors - Extra group
      if(length(conditions) == 4){
        maxabund_ctr = max(IS_Abundance_ctr)
        scalefactor_abund_ctr = rep(NA,length(IS_Abundance_ctr))
        for(i in 1:length(IS_Abundance_ctr)){
          
          scalefactor_abund_ctr[i] = maxabund_ctr/IS_Abundance_ctr[i]
          
        }
        
        maxabund_ex = max(IS_Abundance_ex)
        scalefactor_abund_ex = rep(NA,length(IS_Abundance_ex))
        for(i in 1:length(IS_Abundance_ex)){
          
          scalefactor_abund_ex[i] = maxabund_ex/IS_Abundance_ex[i]
          
        }
      }
      
      #Positive
      for(i in 1:nrow(PosData)){
        ind = which(substr(colnames(PosData),1,3)=="Rep")
        PosData[i,ind] = as.numeric(PosData[i,ind]) * scalefactor_abund_pos[ind]
        
      }
      
      #Negative
      for(i in 1:nrow(NegData)){
        ind = which(substr(colnames(NegData),1,3)=="Rep")
        NegData[i,ind] = as.numeric(NegData[i,ind]) * scalefactor_abund_neg[ind]
        
      }
      
      #Control
      if(length(conditions) == 3){
        for(i in 1:nrow(CtrData)){
          ind = which(substr(colnames(CtrData),1,3)=="Rep")
          CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_abund_ctr[ind]
          
        }
      }
      
      #Extra Group
      if(length(conditions) == 4){
        for(i in 1:nrow(CtrData)){
          ind = which(substr(colnames(CtrData),1,3)=="Rep")
          CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_abund_ctr[ind]
          
        }
        
        for(i in 1:nrow(ExData)){
          ind = which(substr(colnames(ExData),1,3)=="Rep")
          ExData[i,ind] = as.numeric(ExData[i,ind]) * scalefactor_abund_ex[ind]
          
        }
      }
    }else{
      
      if(length(conditions) == 2){
        maxabund = max(c(IS_Abundance_pos,IS_Abundance_neg))
      }else if(length(conditions) == 3){
        maxabund = max(c(IS_Abundance_pos,IS_Abundance_neg,IS_Abundance_ctr))
      }else if(length(conditions) == 4){
        maxabund = max(c(IS_Abundance_pos,IS_Abundance_neg,IS_Abundance_ctr,IS_Abundance_ex))
      }
      
      scalefactor_abund_pos = rep(NA,length(IS_Abundance_pos))
      scalefactor_abund_neg = rep(NA,length(IS_Abundance_neg))
      for(i in 1:length(IS_Abundance_pos)){
        scalefactor_abund_pos[i] = maxabund/IS_Abundance_pos[i]
      }
      for(i in 1:length(IS_Abundance_neg)){
        scalefactor_abund_neg[i] = maxabund/IS_Abundance_neg[i]
      }
      
      #Positive
      for(i in 1:nrow(PosData)){
        ind = which(substr(colnames(PosData),1,3)=="Rep")
        PosData[i,ind] = as.numeric(PosData[i,ind]) * scalefactor_abund_pos[ind]
        
      }
      
      #Negative
      for(i in 1:nrow(NegData)){
        ind = which(substr(colnames(NegData),1,3)=="Rep")
        NegData[i,ind] = as.numeric(NegData[i,ind]) * scalefactor_abund_neg[ind]
        
      }
      
      if(length(conditions) == 3){
        if(ncontrol > 0){
          scalefactor_abund_ctr = rep(NA,length(IS_Abundance_ctr))
          for(i in 1:length(IS_Abundance_ctr)){
            scalefactor_abund_ctr[i] = maxabund/IS_Abundance_ctr[i]
          }
          
          #Control
          for(i in 1:nrow(CtrData)){
            ind = which(substr(colnames(CtrData),1,3)=="Rep")
            CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_abund_ctr[ind]
            
          }
        }
      }
      
      #Extra Group
      if(length(conditions) == 4){
        if(ncontrol > 0){
          scalefactor_abund_ctr = rep(NA,length(IS_Abundance_ctr))
          for(i in 1:length(IS_Abundance_ctr)){
            scalefactor_abund_ctr[i] = maxabund/IS_Abundance_ctr[i]
          }
          
          #Control
          for(i in 1:nrow(CtrData)){
            ind = which(substr(colnames(CtrData),1,3)=="Rep")
            CtrData[i,ind] = as.numeric(CtrData[i,ind]) * scalefactor_abund_ctr[ind]
            
          }
        }
        
        scalefactor_abund_ex = rep(NA,length(IS_Abundance_ex))
        for(i in 1:length(IS_Abundance_ex)){
          
          scalefactor_abund_ex[i] = maxabund/IS_Abundance_ex[i]
          
        }
        
        for(i in 1:nrow(ExData)){
          ind = which(substr(colnames(ExData),1,3)=="Rep")
          ExData[i,ind] = as.numeric(ExData[i,ind]) * scalefactor_abund_ex[ind]
          
        }
        
      }
      
    }
    
  }
  
  
  
}


################################################ Align Experimental Conditions #############################################################

#Build Combined Dataframe
if(length(conditions) == 2){
  DFs = c("PosData","NegData")
  nVOCs = c(nrow(PosData),nrow(NegData))
  maxind = which(nVOCs == max(nVOCs))
  
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  
  tmp = get(DFs[maxind])
  fillDFs = DFs[-maxind]
  rn = rownames(tmp)
}else if(length(conditions) == 3){
  DFs = c("PosData","NegData","CtrData")
  nVOCs = c(nrow(PosData),nrow(NegData),nrow(CtrData))
  maxind = which(nVOCs == max(nVOCs))
  
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  
  tmp = get(DFs[maxind])
  fillDFs = DFs[-maxind]
  rn = rownames(tmp)
}else if(length(conditions) == 4){
  DFs = c("PosData","NegData","CtrData","ExData")
  nVOCs = c(nrow(PosData),nrow(NegData),nrow(CtrData),nrow(ExData))
  maxind = which(nVOCs == max(nVOCs))
  
  if(length(maxind)>1){
    maxind = maxind[1]
  }
  
  tmp = get(DFs[maxind])
  fillDFs = DFs[-maxind]
  rn = rownames(tmp)
}else{
  cat("Error: This code does not support more than 4 experimental groups. \n")
}

if(length(conditions) == 2){
  if(maxind == 1){
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    FinalData = cbind(tmp,NDempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 2){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    FinalData = cbind(PDempty,tmp)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }
}else if(length(conditions) == 3){
  if(maxind == 1){
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    CDempty = data.frame(matrix(NA,nrow(tmp),ncol(CtrData)-1))
    colnames(CDempty) = colnames(CtrData)[1:ncontrol]
    FinalData = cbind(tmp,NDempty,CDempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 2){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    CDempty = data.frame(matrix(NA,nrow(tmp),ncol(CtrData)-1))
    colnames(CDempty) = colnames(CtrData)[1:ncontrol]
    FinalData = cbind(PDempty,tmp,CDempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 3){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    FinalData = cbind(PDempty,NDempty,tmp)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }
}else if(length(conditions) == 4){
  if(maxind == 1){
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    CDempty = data.frame(matrix(NA,nrow(tmp),ncol(CtrData)-1))
    colnames(CDempty) = colnames(CtrData)[1:ncontrol]
    Exempty = data.frame(matrix(NA,nrow(tmp),ncol(ExData)-1))
    colnames(Exempty) = colnames(ExData)[1:nextra]
    FinalData = cbind(tmp,NDempty,CDempty,Exempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 2){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    CDempty = data.frame(matrix(NA,nrow(tmp),ncol(CtrData)-1))
    colnames(CDempty) = colnames(CtrData)[1:ncontrol]
    Exempty = data.frame(matrix(NA,nrow(tmp),ncol(ExData)-1))
    colnames(Exempty) = colnames(ExData)[1:nextra]
    FinalData = cbind(PDempty,tmp,CDempty,Exempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 3){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    Exempty = data.frame(matrix(NA,nrow(tmp),ncol(ExData)-1))
    colnames(Exempty) = colnames(ExData)[1:nextra]
    FinalData = cbind(PDempty,NDempty,tmp,Exempty)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }else if(maxind == 4){
    PDempty = data.frame(matrix(NA,nrow(tmp),ncol(PosData)-1))
    colnames(PDempty) = colnames(PosData)[1:npos]
    NDempty = data.frame(matrix(NA,nrow(tmp),ncol(NegData)-1))
    colnames(NDempty) = colnames(NegData)[1:nneg]
    CDempty = data.frame(matrix(NA,nrow(tmp),ncol(CtrData)-1))
    colnames(CDempty) = colnames(CtrData)[1:ncontrol]
    FinalData = cbind(PDempty,NDempty,CDempty,tmp)
    rownames(FinalData) = rn
    endcol = FinalData$RT
    FinalData = FinalData[,- which(colnames(FinalData) == "RT")]
    FinalData$RT = endcol
  }
  
}


#Fill Combined Data Frame

Pind = seq(1:npos)
Nind = seq(max(Pind)+1,max(Pind)+nneg)
Abindp = seq(1:npos)
Abindn = seq(1:nneg)

if(length(conditions) == 3){
  Cind = seq(max(Nind)+1,max(Nind)+ncontrol)
  Abind2 = seq(1:ncontrol)
}else{
  Cind = NULL
  Abind2 = NULL
}

if(length(conditions) == 4){
  Cind = seq(max(Nind)+1,max(Nind)+ncontrol)
  Abind2 = seq(1:ncontrol)
  Eind = seq(max(Cind)+1,max(Cind)+nextra)
  Abinde = seq(1:nextra)
}

 
for(i in 1:length(fillDFs)){
  
  data = get(fillDFs[i])
  
  if(length(conditions) == 2){
    if(substr(fillDFs[i],1,3) == conditions[1]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Pind] = data[j,Abindp]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions[2]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Nind] = data[j,Abindn]
            
          } 
          
        }
        
      }
      
    }
    
    if((i %% 1) == 0) cat("Condition Aligned:",fillDFs[i],"\n")
    
  }else if(length(conditions) == 3){
    if(substr(fillDFs[i],1,3) == conditions[1]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Pind] = data[j,Abindp]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions[2]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Nind] = data[j,Abindn]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions[3]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Cind] = data[j,Abind2]
            
          } 
          
        }
        
      }
      
    }
    
    if((i %% 1) == 0) cat("Condition Aligned:",fillDFs[i],"\n")
    
  }else if(length(conditions) == 4){
    if(substr(fillDFs[i],1,3) == conditions2[1]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Pind] = data[j,Abindp]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions2[2]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Nind] = data[j,Abindn]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions2[3]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Cind] = data[j,Abind2]
            
          } 
          
        }
        
      }
      
    }else if(substr(fillDFs[i],1,3) == conditions2[4]){
      
      for(j in 1:nrow(data)){
        for(k in 1:nrow(FinalData)){
          
          if(data$RT[j] > FinalData$RT[k] - RTtol2 && data$RT[j] < FinalData$RT[k] + RTtol2 ){
            
            FinalData[k,Eind] = data[j,Abinde]
            
          } 
          
        }
        
      }
      
    }
    
    if((i %% 1) == 0) cat("Condition Aligned:",fillDFs[i],"\n")
    
  }
  
}

#Added 5-21-23
for(i in 1:nrow(FinalData)){
  for(j in 1:(ncol(FinalData)-1)){
    if(!is.na(FinalData[i,j]) && FinalData[i,j] == 0){
      FinalData[i,j] = NA
    }
  }
}

#Converts to numeric matrix
tmp = matrix(as.numeric(unlist(FinalData)),nrow(FinalData))
colnames(tmp) = colnames(FinalData)
rownames(tmp) = rownames(FinalData)

FinalData = tmp

################################################################################################

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
  FinalData = FinalData[-contamindex,]
}


### Optional: Manually remove VOCs in a specific retention time window
if(ClearNoise == T){
  tmpdata = data.frame(FinalData)
  manualrem = rep(NA,nrow(tmpdata))
  count = 1
  
  if(length(Remove_LowRT) == length(Remove_HighRT)){
    for(a in 1:length(Remove_LowRT)){
      for(i in 1:nrow(tmpdata)){
        
        if(tmpdata$RT[i] <= Remove_HighRT[a] && tmpdata$RT[i] >= Remove_LowRT[a]){
            manualrem[count] = i
            count = count+1
          }
        
        }
    }
    
      manualrem = manualrem[!is.na(manualrem)]
      if(length(manualrem)>0){
        tmpdata = tmpdata[-manualrem,]
      }else{
        tmpdata = tmpdata
      }
      
      FinalData = tmpdata
      
  }else{
    cat("ERROR: You have provided an unequal number of retention time pairs for filtering.")
    break
  }
  
  
}

#This is the aligned and normalized data containing all experimental conditions/data
View(FinalData)

############################################  DATAFRAMES FOR STATS  ####################################################

###1) VOCs observed in all conditions and samples
nMissing = rowSums(is.na(FinalData))
AlwaysObservedData = FinalData[which(nMissing == 0),]


###2) VOCs observed in all conditions except the control
if(length(conditions) <= 3){
  nMissing = rowSums(is.na(FinalData[,seq(1,(npos+nneg))]))
  PosNegData = FinalData[which(nMissing == 0),]
}else if(length(conditions) == 4){
  if(originalposnegdata == T){
    nMissing = rowSums(is.na(FinalData[,seq(1,(npos+nneg))]))
    PosNegData = FinalData[which(nMissing == 0),]
  }else{
    tmpvec2 = c(seq(1,(npos+nneg)),seq((npos+nneg+ncontrol)+1,(npos+nneg+ncontrol+nextra)))
    nMissing = rowSums(is.na(FinalData[,tmpvec2]))
    PosNegData = FinalData[which(nMissing == 0),]
  }
  
}


###3) VOCs observed in only the positive, negative, control or extra conditions

if(length(conditions) == 2){
  
    nMissing_p = rowSums(is.na(FinalData[,Pind]))
    nMissing_n = rowSums(is.na(FinalData[,Nind]))
    PosOnlyInd = NegOnlyInd = rep(NA,length(nMissing_p))
    count1=count2=1
    for(i in 1:length(nMissing_p)){
      if(as.numeric(nMissing_p[i]) == 0 && as.numeric(nMissing_n[i]) == 3){
        PosOnlyInd[count1] = i
        count1 = count1+1
      }else if(as.numeric(nMissing_n[i]) == 0 && as.numeric(nMissing_p[i]) == 3){
        NegOnlyInd[count2] = i
        count2 = count2+1
      }
    }
    PosOnlyInd = PosOnlyInd[!is.na(PosOnlyInd)]
    NegOnlyInd = NegOnlyInd[!is.na(NegOnlyInd)]
    
    if(length(PosOnlyInd) != 0){
      PosOnlyData = FinalData[PosOnlyInd,]
    }
    if(length(NegOnlyInd) != 0){
      NegOnlyData = FinalData[NegOnlyInd,]
    }
    
  
}else if(length(conditions) == 3){
  
  if(maintainposnegonly == T){ #Maintains earlier functionality 
    nMissing_p = rowSums(is.na(FinalData[,Pind]))
    nMissing_n = rowSums(is.na(FinalData[,Nind]))
    PosOnlyInd = NegOnlyInd = rep(NA,length(nMissing_p))
    count1=count2=1
    for(i in 1:length(nMissing_p)){
      if(as.numeric(nMissing_p[i]) == 0 && as.numeric(nMissing_n[i]) == 3){
        PosOnlyInd[count1] = i
        count1 = count1+1
      }else if(as.numeric(nMissing_n[i]) == 0 && as.numeric(nMissing_p[i]) == 3){
        NegOnlyInd[count2] = i
        count2 = count2+1
      }
    }
    PosOnlyInd = PosOnlyInd[!is.na(PosOnlyInd)]
    NegOnlyInd = NegOnlyInd[!is.na(NegOnlyInd)]
    
    if(length(PosOnlyInd) != 0){
      PosOnlyData = FinalData[PosOnlyInd,]
    }
    if(length(NegOnlyInd) != 0){
      NegOnlyData = FinalData[NegOnlyInd,]
    }
    
  }else{ #Upgraded functionality to include CtrOnlyData
    nMissing_p = rowSums(is.na(FinalData[,Pind]))
    nMissing_n = rowSums(is.na(FinalData[,Nind]))
    nMissing_c = rowSums(is.na(FinalData[,Cind]))
    PosOnlyInd=rep(NA,length(nMissing_p));NegOnlyInd=rep(NA,length(nMissing_n));CtrOnlyInd=rep(NA,length(nMissing_c))
    count1=count2=count3=1
    for(i in 1:nrow(FinalData)){
      if(as.numeric(nMissing_p[i]) == 0 && as.numeric(nMissing_n[i]) == 3 && as.numeric(nmissing_c[i]) == 3){
        PosOnlyInd[count1] = i
        count1 = count1+1
      }else if(as.numeric(nMissing_n[i]) == 0 && as.numeric(nMissing_p[i]) == 3 && as.numeric(nmissing_c[i]) == 3){
        NegOnlyInd[count2] = i
        count2 = count2+1
      }else if(as.numeric(nmissing_c[i]) == 0 && as.numeric(nMissing_p[i]) == 3 && as.numeric(nmissing_n[i]) == 3){
        CtrOnlyInd[count3] = i
        count3 = count3+1
      }
    }
    PosOnlyInd = PosOnlyInd[!is.na(PosOnlyInd)]
    NegOnlyInd = NegOnlyInd[!is.na(NegOnlyInd)]
    CtrOnlyInd = CtrOnlyInd[!is.na(CtrOnlyInd)]
    
    if(length(PosOnlyInd) != 0){
      PosOnlyData = FinalData[PosOnlyInd,]
    }
    if(length(NegOnlyInd) != 0){
      NegOnlyData = FinalData[NegOnlyInd,]
    }
    if(length(CtrOnlyInd) != 0){
      CtrOnlyInd = FinalData[CtrOnlyInd,]
    }
  }
  
}else if(length(conditions) == 4){
  
  if(maintainposnegonly == T){ #Maintains earlier functionality
    nMissing_p = rowSums(is.na(FinalData[,Pind]))
    nMissing_n = rowSums(is.na(FinalData[,Nind]))
    PosOnlyInd = NegOnlyInd = rep(NA,length(nMissing_p))
    count1=count2=1
    for(i in 1:length(nMissing_p)){
      if(as.numeric(nMissing_p[i]) == 0 && as.numeric(nMissing_n[i]) == 3){
        PosOnlyInd[count1] = i
        count1 = count1+1
      }else if(as.numeric(nMissing_n[i]) == 0 && as.numeric(nMissing_p[i]) == 3){
        NegOnlyInd[count2] = i
        count2 = count2+1
      }
    }
    PosOnlyInd = PosOnlyInd[!is.na(PosOnlyInd)]
    NegOnlyInd = NegOnlyInd[!is.na(NegOnlyInd)]
    
    if(length(PosOnlyInd) != 0){
      PosOnlyData = FinalData[PosOnlyInd,]
    }
    if(length(NegOnlyInd) != 0){
      NegOnlyData = FinalData[NegOnlyInd,]
    }
    
  }else{ #Upgraded functionality to identify control and group 4 only VOCs
    nMissing_p = rowSums(is.na(FinalData[,Pind]))
    nMissing_n = rowSums(is.na(FinalData[,Nind]))
    nMissing_c = rowSums(is.na(FinalData[,Cind]))
    nMissing_e = rowSums(is.na(FinalData[,Eind]))
    PosOnlyInd=rep(NA,length(nMissing_p));NegOnlyInd=rep(NA,length(nMissing_n));CtrOnlyInd=rep(NA,length(nMissing_c));ExOnlyInd==rep(NA,length(nMissing_e)) 
    count1=count2=count3=count4=1
    for(i in 1:nrow(FinalData)){
      if(as.numeric(nMissing_p[i]) == 0 && as.numeric(nMissing_n[i]) == 3 && as.numeric(nmissing_c[i]) == 3 && as.numeric(nmissing_e[i]) == 3){
        PosOnlyInd[count1] = i
        count1 = count1+1
      }else if(as.numeric(nMissing_n[i]) == 0 && as.numeric(nMissing_p[i]) == 3 && as.numeric(nmissing_c[i]) == 3 && as.numeric(nmissing_e[i]) == 3){
        NegOnlyInd[count2] = i
        count2 = count2+1
      }else if(as.numeric(nmissing_c[i]) == 0 && as.numeric(nmissing_p[i]) == 3 && as.numeric(nmissing_n[i]) == 3 && as.numeric(nmissing_e[i]) == 3){
        CtrOnlyInd[count3] = i
        count3 = count3+1
      }else if(as.numeric(nmissing_e[i]) == 0 && as.numeric(nmissing_p[i]) == 3 && as.numeric(nmissing_n[i]) == 3 && as.numeric(nmissing_c[i]) == 3){
        ExOnlyInd[count4] = i
        count4 = count4+1
      }
    }
    PosOnlyInd = PosOnlyInd[!is.na(PosOnlyInd)]
    NegOnlyInd = NegOnlyInd[!is.na(NegOnlyInd)]
    CtrOnlyInd = CtrOnlyInd[!is.na(CtrOnlyInd)]
    ExOnlyInd = ExOnlyInd[!is.na(ExOnlyInd)]
    
    if(length(PosOnlyInd) != 0){
      PosOnlyData = FinalData[PosOnlyInd,]
    }
    if(length(NegOnlyInd) != 0){
      NegOnlyData = FinalData[NegOnlyInd,]
    }
    if(length(CtrOnlyInd) != 0){
      CtrOnlyData = FinalData[CtrOnlyInd,]
    }
    if(length(ExOnlyInd) != 0){
      ExOnlyData = FinalData[ExOnlyInd,]
    }
    
  }
  
}


###4) Keep only the peaks seen in all samples but X
if(filterlowobs == T){
  nMissing = rowSums(is.na(FinalData))
  ThreshObservedData = FinalData[which(nMissing <= Xmissing),]
}


################################################################################################


##################################  Write Files  ###############################################

if(writef == T){
  
  #Positive Data
  finalname = paste("Aligned_Normalized_Positive_",organism,".csv",sep="")
  write.csv(PosData,finalname)
  
  #Negative Data
  finalname = paste("Aligned_Normalized_Negative_",organism,".csv",sep="")
  write.csv(NegData,finalname)
  
  #Control Data
  if(ncontrol > 0){
    finalname = paste("Aligned_Normalized_Control_",organism,".csv",sep="")
    write.csv(CtrData,finalname)
  }
  
  #Group 4 Data
  if(length(conditions) == 4){
    finalname = paste("Aligned_Normalized_Group4_",organism,".csv",sep="")
    write.csv(ExData,finalname)
  }
  
  #Combined Data
  finalname = paste("Aligned_Normalized_AllConditions_",organism,".csv",sep="")
  write.csv(FinalData,finalname)
  
  #VOCs observed in all experimental conditions
  write.csv(AlwaysObservedData,"AlwaysObservedVOCs.csv")
  
  #VOCs observed in all positive and negative experimental conditions
  write.csv(PosNegData,"PosNegVOCs.csv")
  
  #VOCs missing in X experimental conditions
  if(filterlowobs == T){
    write.csv(ThreshObservedData,"ThresholdObservedVOCs.csv")
  }
  
  #VOCs observed in only the positive conditions
  if(length(PosOnlyInd) != 0 ){
    write.csv(PosOnlyData,"PosOnlyVOCs.csv")
  }
  
  #VOCs observed in only the negative conditions
  if(length(NegOnlyInd) != 0){
    write.csv(NegOnlyData,"NegOnlyVOCs.csv")
  }
  
  if(length(CtrOnlyInd) != 0){
    write.csv(CtrOnlyData,"ControlOnlyVOCs.csv")
  }
  
  if(length(conditions) == 4){
    if(length(ExOnlyInd) != 0 ){
      write.csv(ExOnlyData,"Group4OnlyVOCs.csv")
    }
  }
  
}


################################################################################################


####################################################  DATA TRANSFORMATIONS  ##################################################################

if(RunNorm == T){
  
  mydatafornorm = get(normdata)
  
  if(NormMethod == 1){
    
    cat("Data Transformation: Normalization/Transformation not applied. \n")
    TransformedData = mydatafornorm
    ylab = "Raw Abundance"
    if(useIS == T){
      ylab = "Internal Standard Adjusted Abundance"
    }else if(usenormfactors == T){
      ylab = "Scaled/Normalized Abundance"
    }else if(useIS == T && usenormfactors == T){
      ylab = "Internal Standard Adjusted and Normalized Abundance"
    }
    
  }else if(NormMethod == 2){
    
    cat("Data Transformation: PQN applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized Abundance"
    
  }else if(NormMethod == 3){
    
    cat("Data Transformation: PQN + log10 transform applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized and log10 transformed Abundance"
    
  }else if(NormMethod == 4){
    
    cat("Data Transformation: PQN + log10 transform + centering applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = F)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized, log10 transformed, and centered Abundance"
    
  }else if(NormMethod == 5){
    
    cat("Data Transformation: PQN + log10 transform + centering + scaling applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized, log10 transformed, centered, and scaled Abundance"
    
  }else if(NormMethod == 6){
    
    cat("Data Transformation: log10 transformation only. \n")
    TransformedData = log10(mydatafornorm)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed Abundance"
    
  }else if(NormMethod == 7){
    
    cat("Data Transformation: centering only. \n")
    TransformedData = scale(mydatafornorm,center = T,scale = F)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "Centered Abundance"
    
  }else if(NormMethod == 8){
    
    cat("Data Transformation: scaling only. \n")
    TransformedData = scale(mydatafornorm,center = F,scale = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "Scaled Abundance"
    
  }else if(NormMethod == 9){
    
    cat("Data Transformation: zscore only. \n")
    TransformedData = zscore(mydatafornorm,featuresarerows = T,removeNA = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "ZScore Normalized Abundance"
    
  }else if(NormMethod == 10){
    
    cat("Data Transformation: log10 transformation + centering. \n")
    TransformedData = log10(mydatafornorm)
    TransformedData = scale(TransformedData,center = T,scale = F)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed and centered Abundance"
    
  }else if(NormMethod == 11){
    
    cat("Data Transformation: log10 transformation + centering + scaling. \n")
    TransformedData = log10(mydatafornorm)
    TransformedData = scale(TransformedData,center = T,scale = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed, centered, and scaled Abundance"
    
  }else if(NormMethod == 12){
    
    cat("Data Transformation: centering and scaling only. \n")
    TransformedData = scale(mydatafornorm,center = T,scale = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "Centered and Scaled Abundance"
    
  }else if(NormMethod == 13){
    
    cat("Data Transformation: impute missing values only. \n")
    ###  Impute Missing Values (set them equal to something)  ###
    for(i in 1:nrow(mydatafornorm)){
      for(j in 1:ncol(mydatafornorm)){

        if(is.na(mydatafornorm[i,j])){
          mydatafornorm[i,j] = imputevalueto #Arbitrarily low abundance (set according to the middle of the noise)
        }else if(mydatafornorm[i,j] == 0){
          mydatafornorm[i,j] = imputevalueto
        }

      }
    }
    
    TransformedData = mydatafornorm
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    
    ylab = "Raw Abundance with Missing Value Imputation"
    
  }else if(NormMethod == 14){
    
    cat("Data Transformation: impute missing values + log10 transform. \n")
    ###  Impute Missing Values (set them equal to something)  ###
    for(i in 1:nrow(mydatafornorm)){
      for(j in 1:ncol(mydatafornorm)){
        
        if(is.na(mydatafornorm[i,j])){
          mydatafornorm[i,j] = imputevalueto #Arbitrarily low abundance (set according to the middle of the noise)
        }else if(mydatafornorm[i,j] == 0){
          mydatafornorm[i,j] = imputevalueto
        }
        
      }
    }
    
    TransformedData = mydatafornorm
    TransformedData = log10(TransformedData)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed Abundance with Missing Value Imputation"
    
  }else if(NormMethod == 15){
    
    cat("Data Transformation: impute missing values + log10 transform + centering. \n")
    ###  Impute Missing Values (set them equal to something)  ###
    for(i in 1:nrow(mydatafornorm)){
      for(j in 1:ncol(mydatafornorm)){
        
        if(is.na(mydatafornorm[i,j])){
          mydatafornorm[i,j] = imputevalueto #Arbitrarily low abundance (set according to the middle of the noise)
        }else if(mydatafornorm[i,j] == 0){
          mydatafornorm[i,j] = imputevalueto
        }
        
      }
    }
    
    TransformedData = mydatafornorm
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = F)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed and centered Abundance with Missing Value Imputation"
    
  }else if(NormMethod == 16){
    
    cat("Data Transformation: impute missing values + log10 transform + centering + scaling. \n")
    ###  Impute Missing Values (set them equal to something)  ###
    for(i in 1:nrow(mydatafornorm)){
      for(j in 1:ncol(mydatafornorm)){
        
        if(is.na(mydatafornorm[i,j])){
          mydatafornorm[i,j] = imputevalueto #Arbitrarily low abundance (set according to the middle of the noise)
        }else if(mydatafornorm[i,j] == 0){
          mydatafornorm[i,j] = imputevalueto
        }
        
      }
    }
    
    TransformedData = mydatafornorm
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "log10 transformed, centered, and scaled Abundance with Missing Value Imputation"
    
  }else if(NormMethod == 17){
    
    cat("Data Transformation: PQN + zscore normalization applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = zscore(TransformedData,featuresarerows = T,removeNA = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized and Z-Score Transformed Abundance"
    
  }else if(NormMethod == 18){
    
    cat("Data Transformation: PQN + log10 transform + zscore normalization applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    TransformedData = zscore(TransformedData,featuresarerows = T,removeNA = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized, log10 transformed, Z-Score Abundance"
    
  }else if(NormMethod == 19){
    
    cat("Data Transformation: PQN + log10 transform + centering + zscore normalization applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = F)
    TransformedData = zscore(TransformedData,featuresarerows = T,removeNA = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized, log10 transformed, centered, Z-Score Abundance"
    
  }else if(NormMethod == 20){
    
    cat("Data Transformation: PQN + log10 transform + centering + scaling + zscore normalization applied. \n")
    mydatafornorm2 = PREP2.0(mydatafornorm)
    TransformedData = PQN(mydatafornorm2)
    TransformedData = t(TransformedData)
    TransformedData = log10(TransformedData)
    TransformedData = scale(TransformedData,center = T,scale = T)
    TransformedData = zscore(TransformedData,featuresarerows = T,removeNA = T)
    colnames(TransformedData) = colnames(mydatafornorm)
    rownames(TransformedData) = rownames(mydatafornorm)
    TransformedData[,ncol(TransformedData)] = mydatafornorm[,ncol(mydatafornorm)]
    ylab = "PQN Normalized, log10 transformed, centered, scaled, Z-Score Abundance"
    
  }else{
    cat("Warning: You have selected a normalization method that does not exist. Please select an integer between 1 and 20. \n")
  }
  
}

#Convert NAs to 0 (temporarily)
if(RunNorm == T){
  for(i in 1:nrow(TransformedData)){
    for(j in 1:ncol(TransformedData)){
      
      if(is.na(TransformedData[i,j])){
        TransformedData[i,j] = 0
      }
      
    }
  }
  
}


################################################  STATS & PLOTTING  #############################################################

#Possible Dataframe Names:
#1) FinalData - raw data after filtering and normalization; no transformation/PQN
#2) PosNegData - VOCs found in all positive and negative samples; no transformation/PQN
#3) PosOnlyData - VOCs found in ALL positive samples and ZERO negative samples
#4) NegOnlyData - VOCs found in ALL negative samples and ZERO positive samples
#5) AlwaysObservedData - VOCs observed in ALL SAMPLES.
#6) ThreshObservedData - VOCs found in all samples, up to a threshold X. This is total missing observations, not missing observation per experimental condition.
#7) TransformedData - data after performing the desired transformation above
#8) PQNData - PQN normalized data

mydataforstats = get(statdata)

#Or manually after running once
#mydataforstats = PQNData

###########################################   FILTER VOCS FOUND IN REFERENCE LEVEL  ###############################

if(FilterUsingControls == T){
  
  FilterData = data.frame(mydataforstats)
  
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
  
  mydataforstats = mydataforstats[! round(mydataforstats[,ncol(mydataforstats)],2) %in% round(FilterRTs,2),]
  
}


#Converts to numeric matrix
tmp = matrix(as.numeric(unlist(mydataforstats)),nrow(mydataforstats))
colnames(tmp) = colnames(mydataforstats)
rownames(tmp) = rownames(mydataforstats)
mydataforstats = tmp

####Bar plots and statistics 
#Note: the FinalData matrix will always be in the order of Positive, Negative, Control, Group 4.

#Parametric Statistics - Student's t.test - two sided, unpaired, unequal variance
if(uselevels == F){
  
  if(length(conditions) <= 3){
    VOC_pvals_PN = VOC_pvals_PC = VOC_pvals_NC = rep(NA,nrow(mydataforstats))
  }else if(length(conditions) == 4){
    VOC_pvals_PN = VOC_pvals_PC = VOC_pvals_NC = VOC_pvals_PE = VOC_pvals_NE = VOC_pvals_CE = rep(NA,nrow(mydataforstats))
  }
  
  for(i in 1:nrow(mydataforstats)){
    
    if(length(conditions) == 2){
      VOCsummary = t.test(mydataforstats[i,Pind],mydataforstats[i,Nind],alternative = tails,paired = paired,var.equal = eqvar)
      VOC_pvals_PN[i] = VOCsummary$p.value
      
      if(length(Pind) < 3 || length(Nind) < 3){
        cat("Warning: P-value estimate is inaccurate. \n")
      }
    }else if(length(conditions) == 3){
      VOCsummary = t.test(mydataforstats[i,Pind],mydataforstats[i,Nind],alternative = tails,paired = paired,var.equal = eqvar)
      VOC_pvals_PN[i] = VOCsummary$p.value
      
      if(length(Pind) < 3 || length(Nind) < 3){
        cat("Warning: P-value estimate is inaccurate. \n")
      }
      
      if(length(Cind) >= 3){
        VOCsummary2 = t.test(mydataforstats[i,Pind],mydataforstats[i,Cind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PC[i] = VOCsummary2$p.value
        
        VOCsummary3 = t.test(mydataforstats[i,Nind],mydataforstats[i,Cind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NC[i] = VOCsummary3$p.value
      }
    }else if(length(conditions) == 4){
      
      VOCsummary = t.test(mydataforstats[i,Pind],mydataforstats[i,Nind],alternative = tails,paired = paired,var.equal = eqvar)
      VOC_pvals_PN[i] = VOCsummary$p.value
      
      if(length(Pind) < 3 || length(Nind) < 3){
        cat("Warning: P-value estimate is inaccurate. \n")
      }
      
      if(length(Cind) >= 3){
        VOCsummary2 = t.test(mydataforstats[i,Pind],mydataforstats[i,Cind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PC[i] = VOCsummary2$p.value
        
        VOCsummary3 = t.test(mydataforstats[i,Nind],mydataforstats[i,Cind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NC[i] = VOCsummary3$p.value
      }
      
      if(length(Eind) >= 3){
        VOCsummary4 = t.test(mydataforstats[i,Pind],mydataforstats[i,Eind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PE[i] = VOCsummary4$p.value
        
        VOCsummary5 = t.test(mydataforstats[i,Nind],mydataforstats[i,Eind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NE[i] = VOCsummary5$p.value
        
        VOCsummary6 = t.test(mydataforstats[i,Cind],mydataforstats[i,Eind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_CE[i] = VOCsummary6$p.value
      }
    }
    
    
  }
}else{
  if(length(conditions) == 2){
    VOC_pvals_PN = rep(NA,nrow(mydataforstats))
  }else if(length(conditions) == 3){
    VOC_pvals_PN = VOC_pvals_PC = VOC_pvals_NC = rep(NA,nrow(mydataforstats))
  }else if(length(conditions) == 4){
    VOC_pvals_PN = VOC_pvals_PC = VOC_pvals_NC = VOC_pvals_PE = VOC_pvals_NE = VOC_pvals_CE = rep(NA,nrow(mydataforstats))
  }
  
  ntests = length(pairs)/2
  ind = c(1,2)
  
  for(i in 1:ntests){
    currentpair = pairs[ind]
    g1 = currentpair[1]
    g2 = currentpair[2]
    
    if(length(conditions) == 2){
      
      if(g1 == conditions[1]){
        g1ind = Pind
        tind1 = 1
      }else if(g1 == conditions[2]){
        g1ind = Nind
        tind1 = 2
      }
      
      if(g2 == conditions[1]){
        g2ind = Pind
        tind2 = 1
      }else if(g2 == conditions[2]){
        g2ind = Nind
        tind2 = 2
      }
      
    }else if(length(conditions) == 3){
      
      if(g1 == conditions[1]){
        g1ind = Pind
        tind1 = 1
      }else if(g1 == conditions[2]){
        g1ind = Nind
        tind1 = 2
      }else if(g1 == conditions[3]){
        g1ind = Cind
        tind1 = 3
      }
      
      if(g2 == conditions[1]){
        g2ind = Pind
        tind2 = 1
      }else if(g2 == conditions[2]){
        g2ind = Nind
        tind2 = 2
      }else if(g2 == conditions[3]){
        g2ind = Cind
        tind2 = 3
      }
      
    }else if(length(conditions) == 4){
      
      if(g1 == conditions[1]){
        g1ind = Pind
        tind1 = 1
      }else if(g1 == conditions[2]){
        g1ind = Nind
        tind1 = 2
      }else if(g1 == conditions[3]){
        g1ind = Cind
        tind1 = 3
      }else if(g1 == conditions[4]){
        g1ind = Eind
        tind1 = 4
      }
      
      if(g2 == conditions[1]){
        g2ind = Pind
        tind2 = 1
      }else if(g2 == conditions[2]){
        g2ind = Nind
        tind2 = 2
      }else if(g2 == conditions[3]){
        g2ind = Cind
        tind2 = 3
      }else if(g2 == conditions[4]){
        g2ind = Eind
        tind2 = 4
      }
      
    }
    
    if(length(g1ind) < 3 || length(g2ind) < 3){
      cat("Warning: P-value estimate is inaccurate. \n")
    }
    
    if(tind1 == 1 && tind2 == 2){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PN[j] = VOCsummary$p.value
      }
    }else if(tind1 == 2 && tind2 == 1){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PN[j] = VOCsummary$p.value
      }
    }else if(tind1 == 1 && tind2 == 3){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PC[j] = VOCsummary$p.value
      }
    }else if(tind1 == 3 && tind2 == 1){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PC[j] = VOCsummary$p.value
      }
    }else if(tind1 == 1 && tind2 == 4){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PE[j] = VOCsummary$p.value
      }
    }else if(tind1 == 4 && tind2 == 1){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_PE[j] = VOCsummary$p.value
      }
    }else if(tind1 == 2 && tind2 == 3){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NC[j] = VOCsummary$p.value
      }
    }else if(tind1 == 3 && tind2 == 2){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NC[j] = VOCsummary$p.value
      }
    }else if(tind1 == 2 && tind2 == 4){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NE[j] = VOCsummary$p.value
      }
    }else if(tind1 == 4 && tind2 == 2){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_NE[j] = VOCsummary$p.value
      }
    }else if(tind1 == 3 && tind2 == 4){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_CE[j] = VOCsummary$p.value
      }
    }else if(tind1 == 4 && tind2 == 3){
      for(j in 1:nrow(mydataforstats)){
        VOCsummary = t.test(mydataforstats[j,g1ind],mydataforstats[j,g2ind],alternative = tails,paired = paired,var.equal = eqvar)
        VOC_pvals_CE[j] = VOCsummary$p.value
      }
    }
    
    ind = ind + 2
  }
}

if(uselevels == F){
  ###################### Unadjusted
  ## Get the significant VOC names less than the specified alpha value (sigthresh)
  cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and negative \n")
  rownames(mydataforstats)[which(VOC_pvals_PN < sigthresh)]
  
  if(length(Cind) >= 3){
    cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and controls \n")
    rownames(mydataforstats)[which(VOC_pvals_PC < sigthresh)]
    cat("VOCs with unadjusted p-values less than the specified alpha value, comparing negative and controls \n")
    rownames(mydataforstats)[which(VOC_pvals_NC < sigthresh)]
  }
  
  if(length(conditions) == 4){
    if(length(Eind) >= 3){
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_PE < sigthresh)]
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing negative and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_NE < sigthresh)]
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing control and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_CE < sigthresh)]
    }
  }
  
  
  ###################### FDR Adjusted P-values
  ## Get the FDR-adjusted significant VOC names less than the specified alpha value (sigthresh) #Bug fix in v6.8.1
  
  if(length(Cind) < 3){
    VOC_adjp_PN = p.adjust(VOC_pvals_PN,method = MHCmethod)
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and negative \n")
    rownames(mydataforstats)[which(VOC_adjp_PN < sigthresh)]
    
    allptable = data.frame(matrix(NA,nrow=length(VOC_pvals_PN),ncol=2))
    rownames(allptable) = rownames(mydataforstats)
    colnames(allptable) = c("p","adjp")
    allptable$p = VOC_pvals_PN
    allptable$adjp = VOC_adjp_PN
  }
  
  
  
  if(length(Cind) >= 3){
    pvec = c(VOC_pvals_PN,VOC_pvals_PC,VOC_pvals_NC)
    adjpvec = p.adjust(pvec,method=MHCmethod)
    
    allptable = data.frame(matrix(NA,nrow=length(VOC_pvals_PN),ncol=6))
    rownames(allptable) = rownames(mydataforstats)
    colnames(allptable) = c("PNp","PCp","NCp","PNadjp","PCadjp","NCadjp")
    allptable$PNp = VOC_pvals_PN
    allptable$PCp = VOC_pvals_PC
    allptable$NCp = VOC_pvals_NC
    
    tmplen = length(VOC_pvals_PN)
    tmpseq = seq(1,length(pvec),tmplen)
    
    allptable$PNadjp = adjpvec[tmpseq[1]:(tmpseq[2]-1)]
    allptable$PCadjp = adjpvec[tmpseq[2]:(tmpseq[3]-1)]
    allptable$NCadjp = adjpvec[tmpseq[3]:(length(pvec))]
    
    VOC_adjp_PN = allptable$PNadjp
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and negative \n")
    rownames(mydataforstats)[which(VOC_adjp_PN < sigthresh)]
    
    VOC_adjp_PC = allptable$PCadjp
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and controls \n")
    rownames(mydataforstats)[which(VOC_adjp_PC < sigthresh)]
    
    VOC_adjp_NC = allptable$NCadjp
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and controls \n")
    rownames(mydataforstats)[which(VOC_adjp_NC < sigthresh)]
  }
  
  if(length(conditions) == 4){
    if(length(Eind) >= 3){
      
      pvec = c(VOC_pvals_PN,VOC_pvals_PC,VOC_pvals_PE,VOC_pvals_NC,VOC_pvals_NE,VOC_pvals_CE)
      adjpvec = p.adjust(pvec,method=MHCmethod)
      
      allptable = data.frame(matrix(NA,nrow=length(VOC_pvals_PN),ncol=12))
      rownames(allptable) = rownames(mydataforstats)
      colnames(allptable) = c("PNp","PCp","PEp","NCp","NEp","CEp","PNadjp","PCadjp","PEadjp","NCadjp","NEadjp","CEadjp")
      allptable$PNp = VOC_pvals_PN
      allptable$PCp = VOC_pvals_PC
      allptable$PEp = VOC_pvals_PE
      allptable$NCp = VOC_pvals_NC
      allptable$NEp = VOC_pvals_NE
      allptable$CEp = VOC_pvals_CE
      
      tmplen = length(VOC_pvals_PN)
      tmpseq = seq(1,length(pvec),tmplen)
      
      allptable$PNadjp = adjpvec[tmpseq[1]:(tmpseq[2]-1)]
      allptable$PCadjp = adjpvec[tmpseq[2]:(tmpseq[3]-1)]
      allptable$PEadjp = adjpvec[tmpseq[3]:(tmpseq[4]-1)]
      allptable$NCadjp = adjpvec[tmpseq[4]:(tmpseq[5]-1)]
      allptable$NEadjp = adjpvec[tmpseq[5]:(tmpseq[6]-1)]
      allptable$CEadjp = adjpvec[tmpseq[6]:(length(pvec))]
      
      VOC_adjp_PN = allptable$PNadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and negative \n")
      rownames(mydataforstats)[which(VOC_adjp_PN < sigthresh)]
      
      VOC_adjp_PC = allptable$PCadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and controls \n")
      rownames(mydataforstats)[which(VOC_adjp_PC < sigthresh)]
      
      VOC_adjp_NC = allptable$NCadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and controls \n")
      rownames(mydataforstats)[which(VOC_adjp_NC < sigthresh)]
      
      VOC_adjp_PE = allptable$PEadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_PE < sigthresh)]
      
      VOC_adjp_NE = allptable$NEadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_NE < sigthresh)]
      
      VOC_adjp_CE = allptable$CEadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing controls and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_CE < sigthresh)]
      
    }else{
      pvec = c(VOC_pvals_PN,VOC_pvals_PC,VOC_pvals_NC)
      adjpvec = p.adjust(pvec,method=MHCmethod)
      
      allptable = data.frame(matrix(NA,nrow=length(VOC_pvals_PN),ncol=6))
      rownames(allptable) = rownames(mydataforstats)
      colnames(allptable) = c("PNp","PCp","NCp","PNadjp","PCadjp","NCadjp")
      allptable$PNp = VOC_pvals_PN
      allptable$PCp = VOC_pvals_PC
      allptable$NCp = VOC_pvals_NC
      
      tmplen = length(VOC_pvals_PN)
      tmpseq = seq(1,length(pvec),tmplen)
      
      allptable$PNadjp = adjpvec[tmpseq[1]:(tmpseq[2]-1)]
      allptable$PCadjp = adjpvec[tmpseq[2]:(tmpseq[3]-1)]
      allptable$NCadjp = adjpvec[tmpseq[3]:(length(pvec))]
      
      VOC_adjp_PN = allptable$PNadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and negative \n")
      rownames(mydataforstats)[which(VOC_adjp_PN < sigthresh)]
      
      VOC_adjp_PC = allptable$PCadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and controls \n")
      rownames(mydataforstats)[which(VOC_adjp_PC < sigthresh)]
      
      VOC_adjp_NC = allptable$NCadjp
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and controls \n")
      rownames(mydataforstats)[which(VOC_adjp_NC < sigthresh)]
    }
  }
  
}else{
  
  if(sum(is.na(VOC_pvals_PN)) != length(VOC_pvals_PN)){
    cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and negative \n")
    rownames(mydataforstats)[which(VOC_pvals_PN < sigthresh)]
    
    VOC_adjp_PN = p.adjust(VOC_pvals_PN,method = MHCmethod)
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and negative \n")
    rownames(mydataforstats)[which(VOC_adjp_PN < sigthresh)]
  }else{
    VOC_adjp_PN = rep(NA,length(VOC_pvals_PN))
  }
  
  if(sum(is.na(VOC_pvals_PC)) != length(VOC_pvals_PC)){
    cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and controls \n")
    rownames(mydataforstats)[which(VOC_pvals_PC < sigthresh)]
    
    VOC_adjp_PC = p.adjust(VOC_pvals_PC,method = MHCmethod)
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and controls \n")
    rownames(mydataforstats)[which(VOC_adjp_PC < sigthresh)]
  }else{
    VOC_adjp_PC = rep(NA,length(VOC_pvals_PC))
  }
  
  if(sum(is.na(VOC_pvals_NC)) != length(VOC_pvals_NC)){
    cat("VOCs with unadjusted p-values less than the specified alpha value, comparing negative and controls \n")
    rownames(mydataforstats)[which(VOC_pvals_NC < sigthresh)]
    
    VOC_adjp_NC = p.adjust(VOC_pvals_NC,method = MHCmethod)
    cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and controls \n")
    rownames(mydataforstats)[which(VOC_adjp_NC < sigthresh)]
  }else{
    VOC_adjp_NC = rep(NA,length(VOC_pvals_NC))
  }
  
  if(length(conditions) == 4){
    if(sum(is.na(VOC_pvals_PE)) != length(VOC_pvals_PE)){
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing positive and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_PE < sigthresh)]
      
      VOC_adjp_PE = p.adjust(VOC_pvals_PE,method = MHCmethod)
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing positive and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_PE < sigthresh)]
    }else{
      VOC_adjp_PE = rep(NA,length(VOC_pvals_PE))
    }
    
    if(sum(is.na(VOC_pvals_NE)) != length(VOC_pvals_NE)){
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing negative and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_NE < sigthresh)]
      
      VOC_adjp_NE = p.adjust(VOC_pvals_NE,method = MHCmethod)
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing negative and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_NE < sigthresh)]
    }else{
      VOC_adjp_NE = rep(NA,length(VOC_pvals_NE))
    }
    
    if(sum(is.na(VOC_pvals_CE)) != length(VOC_pvals_CE)){
      cat("VOCs with unadjusted p-values less than the specified alpha value, comparing controls and group 4 \n")
      rownames(mydataforstats)[which(VOC_pvals_CE < sigthresh)]
      
      VOC_adjp_CE = p.adjust(VOC_pvals_CE,method = MHCmethod)
      cat("VOCs with FDR-adjusted p-values less than the specified alpha value, comparing controls and group 4 \n")
      rownames(mydataforstats)[which(VOC_adjp_CE < sigthresh)]
    }else{
      VOC_adjp_CE = rep(NA,length(VOC_pvals_CE))
    }
  }
  
}

if(writeadjustedp == T){
  write.csv(allptable,paste(organism,"_pvals_",MHCmethod,"adj.csv",sep=""))
}

if(PlotSigOnly == T){
  require(ggplot2)
  
  if(length(conditions) <= 3){
  if(length(Cind) >= 3){ 
    
    plotind1 = which(VOC_adjp_PN < sigthresh)
    plotind2 = which(VOC_adjp_PC < sigthresh)
    plotind3 = which(VOC_adjp_NC < sigthresh)
    plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3))))
    PNp = VOC_adjp_PN[plotind]
    PCp = VOC_adjp_PC[plotind]
    NCp = VOC_adjp_NC[plotind]
    
    if(AllowUnadjusted == T){
      plotind1 = which(VOC_pvals_PN < sigthresh)
      plotind2 = which(VOC_pvals_PC < sigthresh)
      plotind3 = which(VOC_pvals_NC < sigthresh)
      plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3))))
      PNp = VOC_pvals_PN[plotind]
      PCp = VOC_pvals_PC[plotind]
      NCp = VOC_pvals_NC[plotind]
    }
    
    plotdata = mydataforstats[plotind,]
    
    if(SpecificVOCs == T){
      plotdata = plotdata[rownames(plotdata) %in% VOClist,]
    }
    
    for(i in 1:nrow(plotdata)){
      currentp = c(PNp[i],PCp[i],NCp[i])
      df = data.frame(matrix(NA,length(conditions),3))
      colnames(df) = c("Condition","Average","SE")
      if(length(conditions) == 2){
        df$Condition = c("+DHPPA","-DHPPA")
      }else if(length(conditions) == 3){
        df$Condition = c("+DHPPA","-DHPPA","BHI Control")
      }
      
      for(z in 1:length(currentp)){
        if(is.na(currentp[z])){
          currentp[z] = 1
        }
      }
      
      Pavg = mean(plotdata[i,Pind],na.rm = T)
      Navg = mean(plotdata[i,Nind],na.rm = T)
      Cavg = mean(plotdata[i,Cind],na.rm = T)
      avgs = c(Pavg,Navg,Cavg)
      df$Average = avgs
      
      Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
      Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
      Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
      ses = c(Pse,Nse,Cse)
      df$SE = ses
      
      for(j in 1:nrow(df)){
        for(k in 1:ncol(df)){
          if(is.na(df[j,k])){
            df[j,k] = 0
          }
        }
      }
      
      #Dataframe for plotting
      if(length(conditions) == 2){
        df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA"))
        colors = c(negcolor,poscolor)
      }else if(length(conditions) == 3){
        df$Condition = factor(df$Condition,levels=c("+DHPPA","-DHPPA","BHI Control"))
        colors = c(poscolor,negcolor,ctrcolor)
      }
      
      
      
      ttl = gsub('\\..*','',rownames(plotdata)[i])
      
      #Generate base plot
      p = ggplot(data=df, aes(x=Condition, y=Average)) +
        geom_bar(stat="identity",fill = colors)
      p = p + theme_minimal()
      p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
      p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
      p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
      p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
      p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                          position=position_dodge(.9)) 
      
      #Add p-values
      
      if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
        upperlim = max(avgs + ses)+max(avgs + ses)*0.2
        p = p + ylim(lowerlim,upperlim)
        
        if(currentp[1] < sigthresh){
          
          #Add p-values for non-significant comparisons, if desired
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
          }
          
        }else if(currentp[2] < sigthresh){
          
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
          }
          
        }else if(currentp[3] < sigthresh){
          
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
          }
          
        }
        
      }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
        upperlim = max(avgs + ses)+max(avgs + ses)*0.40
        p = p + ylim(lowerlim,upperlim)
        
        if(currentp[1] < sigthresh && currentp[2] < sigthresh){
          
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
          }
          
        }else if(currentp[1] < sigthresh && currentp[3] < sigthresh){
          
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
          }
          
          
        }else if(currentp[2] < sigthresh && currentp[3] < sigthresh){
          
          if(addpvalues == T){
            upperlim = max(avgs + ses) + max(avgs + ses)*0.6
            p = p + ylim(lowerlim,upperlim)
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
            }
            
          }else{
            
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
            p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
            }
            
            
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
            p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
            if(AllowUnadjusted == T){
              pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }else{
              pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
              p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
            }
            
          }
          
        }
      }else if(length(which(currentp<sigthresh)) == 3){ # if all three " "
        upperlim = max(avgs + ses) + max(avgs + ses)*0.6
        p = p + ylim(lowerlim,upperlim)
        
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        if(AllowUnadjusted == T){
          pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }else{
          pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }
        
        
        p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
        if(AllowUnadjusted == T){
          pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
        }else{
          pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
        }
        
        
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
        p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
        if(AllowUnadjusted == T){
          pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
        }else{
          pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
        }
        
        
      }
      
      print(p)
    }
    
  }else{ #two or less control samples - not recommended for publications.
    
    plotind1 = as.numeric(which(VOC_adjp_PN < sigthresh))
    PNp = VOC_adjp_PN[plotind1]
    
    if(AllowUnadjusted == T){
      plotind1 = as.numeric(which(VOC_pvals_PN < sigthresh))
      PNp = VOC_pvals_PN[plotind1]
    }
    
    plotdata = mydataforstats[plotind1,]
    
    if(SpecificVOCs == T){
      plotdata = plotdata[rownames(plotdata) %in% VOClist,]
    }
    
    for(i in 1:nrow(plotdata)){
      currentp = PNp[i]
      if(length(Cind) >= 1){
        df = data.frame(matrix(NA,length(conditions),3))
        colnames(df) = c("Condition","Average","SE")
        if(length(conditions) == 2){
          df$Condition = c("+DHPPA","-DHPPA")
        }else if(length(conditions) == 3){
          df$Condition = c("+DHPPA","-DHPPA","BHI Control")
        }
        
        Pavg = mean(plotdata[i,Pind],na.rm = T)
        Navg = mean(plotdata[i,Nind],na.rm = T)
        if(length(Cind) == 2){
          Cavg = mean(plotdata[i,Cind],na.rm=T)
        }else{
          Cavg = plotdata[i,Cind]
        }
        
        avgs = c(Pavg,Navg,Cavg)
        df$Average = avgs
        
        Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
        Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
        if(length(Cind) == 2){
          Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
          cat("Warning: This standard error estimate is poor. \n")
        }else{
          Cse = NA
        }
        ses = c(Pse,Nse,Cse)
        df$SE = ses
        
        
        #Dataframe for plotting
        if(length(conditions) == 2){
          df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA"))
          colors = c(negcolor,poscolor)
        }else if(length(conditions) == 3){
          df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA"))
          colors = c(ctrcolor,negcolor,poscolor)
        }
        
      }else{
        df = data.frame(matrix(NA,length(conditions),3))
        colnames(df) = c("Condition","Average","SE")
        df$Condition = c("+DHPPA","-DHPPA")
        
        Pavg = mean(plotdata[i,Pind],na.rm = T)
        Navg = mean(plotdata[i,Nind],na.rm = T)
        avgs = c(Pavg,Navg)
        df$Average = avgs
        
        Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
        Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
        ses = c(Pse,Nse)
        df$SE = ses
        
        
        #Dataframe for plotting
        df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA"))
        colors = c(negcolor,poscolor)
      }
      
      ttl = gsub('\\..*','',rownames(plotdata)[i])
      
      
      
      #Generate base plot
      p = ggplot(data=df, aes(x=Condition, y=Average)) +
        geom_bar(stat="identity",fill = colors)
      p = p + theme_minimal()
      p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
      p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
      p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
      p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
      p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                          position=position_dodge(.9)) 
      
      #Add p-values
      if(currentp < sigthresh){
        upperlim = max(avgs + ses)+max(avgs + ses)*0.2
        p = p + ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        if(AllowUnadjusted == T){
          pstat = paste("p-value:",formatC(currentp,format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }else{
          pstat = paste("FDR p-value:",formatC(currentp,format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }
      }else if(addpvalues == T){
        upperlim = max(avgs + ses)+max(avgs + ses)*0.2
        p = p + ylim(lowerlim,upperlim)
        p = p + ylim(lowerlim,upperlim)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
        p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
        if(AllowUnadjusted == T){
          pstat = paste("p-value:",formatC(currentp,format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }else{
          pstat = paste("FDR p-value:",formatC(currentp,format = "e",digits = 3))
          p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
        }
      }
      
      print(p)
    }
    
  }
  }else if(length(conditions) == 4){
    if(length(Cind) >= 3){ 
      
      plotind1 = which(VOC_adjp_PN < sigthresh)
      plotind2 = which(VOC_adjp_PC < sigthresh)
      plotind3 = which(VOC_adjp_NC < sigthresh)
      plotind4 = which(VOC_adjp_PE < sigthresh)
      plotind5 = which(VOC_adjp_NE < sigthresh)
      plotind6 = which(VOC_adjp_CE < sigthresh)
      plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3,plotind4,plotind5,plotind6))))
      PNp = VOC_adjp_PN[plotind]
      PCp = VOC_adjp_PC[plotind]
      NCp = VOC_adjp_NC[plotind]
      PEp = VOC_adjp_PE[plotind]
      NEp = VOC_adjp_NE[plotind]
      CEp = VOC_adjp_CE[plotind]
      
      if(AllowUnadjusted == T){
        plotind1 = which(VOC_pvals_PN < sigthresh)
        plotind2 = which(VOC_pvals_PC < sigthresh)
        plotind3 = which(VOC_pvals_NC < sigthresh)
        plotind4 = which(VOC_pvals_PE < sigthresh)
        plotind5 = which(VOC_pvals_NE < sigthresh)
        plotind6 = which(VOC_pvals_CE < sigthresh)
        plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3,plotind4,plotind5,plotind6))))
        PNp = VOC_pvals_PN[plotind]
        PCp = VOC_pvals_PC[plotind]
        NCp = VOC_pvals_NC[plotind]
        PEp = VOC_pvals_PE[plotind]
        NEp = VOC_pvals_NE[plotind]
        CEp = VOC_pvals_CE[plotind]
      }
      
      plotdata = mydataforstats[plotind,]
      
      if(SpecificVOCs == T){
        plotdata = plotdata[rownames(plotdata) %in% VOClist,]
      }
      
      for(i in 1:nrow(plotdata)){
        currentp = c(PNp[i],PCp[i],NCp[i],PEp[i],NEp[i],CEp[i])
        ptable = data.frame(matrix(NA,length(currentp),2))
        if(AllowUnadjusted == T){
          colnames(ptable) = c("Comparison","UnadjustedP")
          ptable$UnadjustedP = currentp
          ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Ctr","Neg_vs_Ctr","Pos_vs_Group4","Neg_vs_Group4","Ctr_vs_Group4")
        }else{
          colnames(ptable) = c("Comparison","FDRAdjustedP")
          ptable$FDRAdjustedP = currentp
          ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Ctr","Neg_vs_Ctr","Pos_vs_Group4","Neg_vs_Group4","Ctr_vs_Group4")
        }
        
        df = data.frame(matrix(NA,length(conditions),3))
        colnames(df) = c("Condition","Average","SE")
        df$Condition = c("+DHPPA","-DHPPA","BHI Control","Tungstate Inhibition")
        
        Pavg = mean(plotdata[i,Pind],na.rm = T)
        Navg = mean(plotdata[i,Nind],na.rm = T)
        Cavg = mean(plotdata[i,Cind],na.rm = T)
        Eavg = mean(plotdata[i,Eind],na.rm = T)
        avgs = c(Pavg,Navg,Cavg,Eavg)
        df$Average = avgs
        
        Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
        Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
        Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
        Ese = sd(plotdata[i,Eind])/sqrt(length(plotdata[i,Eind]))
        ses = c(Pse,Nse,Cse,Ese)
        df$SE = ses
        
        for(j in 1:nrow(df)){
          for(k in 1:ncol(df)){
            if(is.na(df[j,k])){
              df[j,k] = 0
            }
          }
        }
        
        #Dataframe for plotting
        df$Condition = factor(df$Condition,levels=c("+DHPPA","-DHPPA","BHI Control","Tungstate Inhibition"))
        colors = c(poscolor,negcolor,ctrcolor,excolor)
        
        #BARPLOTHELP
        #Flip order so that it gives a specific order
        #df$Condition = factor(df$Condition,levels=c("+DHPPA","AFMT Inhibition","-DHPPA","BHI Control"))
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
        
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        if(i == 1){
          cat("Four conditions allow multiple significance bar stylings: Add p-values and significance bars in graphic editors (inkscape or adobe illustrator). \n")
          cat("P-values for each VOC are printed in the Console. \n")
        }
        
        text = paste("P-values for: ",rownames(plotdata)[i],"\n",sep="")
        cat(text)
        print(ptable[,1:2])
        
        #Automatically adjust the plot height
        if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
          upperlim = max(avgs + ses)+max(avgs + ses)*0.40
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 3){ 
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 4){
          upperlim = max(avgs + ses) + max(avgs + ses)*0.8
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 5){
          upperlim = max(avgs + ses) + max(avgs + ses)*1
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 6){
          upperlim = max(avgs + ses) + max(avgs + ses)*1.2
          p = p + ylim(lowerlim,upperlim)
        }
        
        if(addpvalues == T){
          upperlim = max(avgs + ses) + max(avgs + ses)*1.2
          p = p + ylim(lowerlim,upperlim)
        }
        
        print(p)
      }
      
    }else{ #two or less control samples - not recommended for publications.
      
      plotind1 = as.numeric(which(VOC_adjp_PN < sigthresh))
      PNp = VOC_adjp_PN[plotind1]
      plotind2 = as.numeric(which(VOC_adjp_PE < sigthresh))
      PEp = VOC_adjp_PE[plotind2]
      plotind3 = as.numeric(which(VOC_adjp_NE < sigthresh))
      NEp = VOC_adjp_NE[plotind3]
      plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3))))
      
      if(AllowUnadjusted == T){
        plotind1 = as.numeric(which(VOC_pvals_PN < sigthresh))
        PNp = VOC_pvals_PN[plotind1]
        plotind2 = as.numeric(which(VOC_pvals_PE < sigthresh))
        PEp = VOC_pvals_PE[plotind2]
        plotind3 = as.numeric(which(VOC_pvals_NE < sigthresh))
        NEp = VOC_pvals_NE[plotind3]
        plotind = as.numeric(names(table(c(plotind1,plotind2,plotind3))))
        
      }
      
      plotdata = mydataforstats[plotind,]
      
      if(SpecificVOCs == T){
        plotdata = plotdata[rownames(plotdata) %in% VOClist,]
      }
      
      for(i in 1:nrow(plotdata)){
        currentp = c(PNp[i],PEp[i],NEp[i])
        ptable = data.frame(matrix(NA,length(currentp),2))
        ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Group4","Neg_vs_Group4")
        if(AllowUnadjusted == T){
          colnames(ptable) = c("Comparison","UnadjustedP")
          ptable$UnadjustedP = currentp
        }else{
          colnames(ptable) = c("Comparison","FDRAdjustedP")
          ptable$FDRAdjustedP = currentp
        }
        
        if(length(Cind) >= 1){
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","-DHPPA","BHI Control","Tungstate Inhibition")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          if(length(Cind) == 2){
            Cavg = mean(plotdata[i,Cind],na.rm=T)
          }else{
            Cavg = plotdata[i,Cind]
          }
          Eavg = mean(plotdata[i,Eind],na.rm = T)
          
          avgs = c(Pavg,Navg,Cavg,Eavg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          if(length(Cind) == 2){
            Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
            cat("Warning: This standard error estimate is poor.")
          }else{
            Cse = NA
          }
          Ese = sd(plotdata[i,Eind])/sqrt(length(plotdata[i,Eind]))
          ses = c(Pse,Nse,Cse,Ese)
          df$SE = ses
          
          
          #Dataframe for plotting
          df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA","Tungstate Inhibition"))
          colors = c(ctrcolor,negcolor,poscolor,excolor)
          
        }else{
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","-DHPPA","Tungstate Inhibition")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          Eavg = mean(plotdata[i,Eind],na.rm = T)
          avgs = c(Pavg,Navg,Eavg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          ses = c(Pse,Nse,Ese)
          df$SE = ses
          
          
          #Dataframe for plotting
          df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA","Tungstate Inhibition"))
          colors = c(negcolor,poscolor,excolor)
          #BARPLOTHELP
          #Flip order so that it gives a specific order
          #df$Condition = factor(df$Condition,levels=c("+DHPPA","AFMT Inhibition","-DHPPA","BHI Control")
        }
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
      
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        if(i == 1){
          cat("Four conditions allow multiple significance bar stylings: Add p-values and significance bars in graphic editors (inkscape or adobe illustrator). \n")
          cat("P-values for each VOC are printed in the Console. \n")
        }
        
        text = paste("P-values for: ",rownames(plotdata)[i],"\n",sep="")
        cat(text)
        ptable
        
        #Automatically adjust the plot height
        if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
          upperlim = max(avgs + ses)+max(avgs + ses)*0.40
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 3){ 
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }else if(addpvalues == T){
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }
        
        print(p)
      }
      
    }
  }
  
}else{ #Plot all VOCs, regardless of significance
  require(ggplot2)
  
  if(length(conditions) <= 3){
    if(length(Cind) >= 3){ 
      
      plotdata = mydataforstats[startrow:stoprow,]
      
      if(SpecificVOCs == T){
        plotdata = plotdata[rownames(plotdata) %in% VOClist,]
      }
      
      PNp = VOC_adjp_PN[startrow:stoprow]
      PCp = VOC_adjp_PC[startrow:stoprow]
      NCp = VOC_adjp_NC[startrow:stoprow]
      
      if(AllowUnadjusted == T){
        PNp = VOC_pvals_PN[startrow:stoprow]
        PCp = VOC_pvals_PC[startrow:stoprow]
        NCp = VOC_pvals_NC[startrow:stoprow]
      }
      
      for(i in 1:nrow(plotdata)){
        currentp = c(PNp[i],PCp[i],NCp[i])
        df = data.frame(matrix(NA,length(conditions),3))
        colnames(df) = c("Condition","Average","SE")
        df$Condition = c("+DHPPA","-DHPPA","BHI Control")
        
        Pavg = mean(plotdata[i,Pind],na.rm = T)
        Navg = mean(plotdata[i,Nind],na.rm = T)
        Cavg = mean(plotdata[i,Cind],na.rm = T)
        avgs = c(Pavg,Navg,Cavg)
        df$Average = avgs
        
        Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
        Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
        Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
        ses = c(Pse,Nse,Cse)
        df$SE = ses
        
        for(j in 1:nrow(df)){
          for(k in 1:ncol(df)){
            if(is.na(df[j,k])){
              df[j,k] = 0
            }
          }
        }
        
        for(z in 1:length(currentp)){
          if(is.na(currentp[z])){
            currentp[z] = 1
          }
        }
        
        #Dataframe for plotting
        
        df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA"))
        colors = c(ctrcolor,negcolor,poscolor)
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
        
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        #Add p-values
        if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
          
          if(currentp[1] < sigthresh){
            
            #Add p-values for non-significant comparisons, if desired
            if(addpvalues == T){
              upperlim = max(avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
            }
            
          }else if(currentp[2] < sigthresh){
            
            if(addpvalues == T){
              upperlim = max (avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
            }
            
          }else if(currentp[3] < sigthresh){
            
            if(addpvalues == T){
              upperlim = max(avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
            }
            
          }
          
        }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
          upperlim = max(avgs + ses)+max(avgs + ses)*0.40
          p = p + ylim(lowerlim,upperlim)
          
          if(currentp[1] < sigthresh && currentp[2] < sigthresh){
            
            if(addpvalues == T){
              upperlim = max(avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
            }
            
          }else if(currentp[1] < sigthresh && currentp[3] < sigthresh){
            
            if(addpvalues == T){
              upperlim = max(avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
            }
            
            
          }else if(currentp[2] < sigthresh && currentp[3] < sigthresh){
            
            if(addpvalues == T){
              upperlim = max(avgs + ses) + max(avgs + ses)*0.6
              p = p + ylim(lowerlim,upperlim)
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
              }
              
            }else{
              
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
              p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
              }
              
              
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
              p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
              if(AllowUnadjusted == T){
                pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }else{
                pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
                p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.35)
              }
              
            }
            
          }
        }else if(length(which(currentp<sigthresh)) == 3){ # if all three " "
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
          
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }
          
          
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
          p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
          }
          
          
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
          p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
          }
          
          
        }
        
        #In case no comparisons are significant
        if(addpvalues == T){
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
          
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[1],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[1],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }
          
          
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.3),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
          p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.3,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.25),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[3],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[3],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2.5,y=max(avgs + ses)+max(avgs + ses)*0.35)
          }
          
          
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.5),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
          p = p+geom_segment(aes(x=3,y=max(avgs + ses)+max(avgs + ses)*0.5,xend=3,yend=max(avgs + ses)+max(avgs + ses)*0.45),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp[2],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp[2],format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=2,y=max(avgs + ses)+max(avgs + ses)*0.55)
          }
          
        }
        
        print(p)
      }
      
    }else{ #if less than 3 control samples
      
      plotdata = mydataforstats[startrow:stoprow,]
      
      if(SpecificVOCs == T){
        plotdata = plotdata[rownames(plotdata) %in% VOClist,]
      }
      
      PNp = VOC_adjp_PN[startrow:stoprow]
      
      if(AllowUnadjusted == T){
        PNp = VOC_pvals_PN[startrow:stoprow]
      }
      
      for(i in 1:nrow(plotdata)){
        currentp = PNp[i]
        if(length(Cind) >= 1){
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","DHPPA","BHI Control")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          if(length(Cind) == 2){
            Cavg = mean(plotdata[i,Cind],na.rm=T)
          }else{
            Cavg = plotdata[i,Cind]
          }
          
          avgs = c(Pavg,Navg,Cavg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          if(length(Cind) == 2){
            Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
            cat("Warning: This standard error estimate is poor. \n")
          }else{
            Cse = NA
          }
          ses = c(Pse,Nse,Cse)
          df$SE = ses
          
          
          #Dataframe for plotting
          df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA"))
          colors = c(ctrcolor,negcolor,poscolor)
        }else{
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","-DHPPA")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          avgs = c(Pavg,Navg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          ses = c(Pse,Nse)
          df$SE = ses
          
          
          #Dataframe for plotting
          df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA"))
          colors = c(negcolor,poscolor)
        }
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
        
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        #Add p-values
        if(currentp < sigthresh){
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp,format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp,format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }
        }else if(addpvalues == T){
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
          p = p + ylim(lowerlim,upperlim)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.1),size=0.8)
          p = p+geom_segment(aes(x=1,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=1,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          p = p+geom_segment(aes(x=2,y=max(avgs + ses)+max(avgs + ses)*0.1,xend=2,yend=max(avgs + ses)+max(avgs + ses)*0.05),size=0.8)
          if(AllowUnadjusted == T){
            pstat = paste("p-value:",formatC(currentp,format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }else{
            pstat = paste("FDR p-value:",formatC(currentp,format = "e",digits = 3))
            p = p+annotate("text",label=pstat,x=1.5,y=max(avgs + ses)+max(avgs + ses)*0.15)
          }
        }
        
        print(p)
      }
    }
  }else if(length(conditions) == 4){
    if(length(Cind) >= 3){ 
      
      plotdata = mydataforstats[startrow:stoprow,]
      
      if(SpecificVOCs == T){
        plotdata = plotdata[rownames(plotdata) %in% VOClist,]
      }
      
      PNp = VOC_adjp_PN[startrow:stoprow]
      PCp = VOC_adjp_PC[startrow:stoprow]
      NCp = VOC_adjp_NC[startrow:stoprow]
      PEp = VOC_adjp_PE[startrow:stoprow]
      NEp = VOC_adjp_NE[startrow:stoprow]
      CEp = VOC_adjp_CE[startrow:stoprow]
      
      if(AllowUnadjusted == T){
        PNp = VOC_pvals_PN[startrow:stoprow]
        PCp = VOC_pvals_PC[startrow:stoprow]
        NCp = VOC_pvals_NC[startrow:stoprow]
        PEp = VOC_pvals_PE[startrow:stoprow]
        NEp = VOC_pvals_NE[startrow:stoprow]
        CEp = VOC_pvals_CE[startrow:stoprow]
      }
      
      
      for(i in 1:nrow(plotdata)){
        currentp = c(PNp[i],PCp[i],NCp[i],PEp[i],NEp[i],CEp[i])
        ptable = data.frame(matrix(NA,length(currentp),2))
        ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Ctr","Neg_vs_Ctr","Pos_vs_Group4","Neg_vs_Group4","Ctr_vs_Group4")
        if(AllowUnadjusted == T){
          colnames(ptable) = c("Comparison","UnadjustedP")
          ptable$UnadjustedP = currentp
        }else{
          colnames(ptable) = c("Comparison","FDRAdjustedP")
          ptable$FDRAdjustedP = currentp
        }
        
        for(z in 1:length(currentp)){
          if(is.na(currentp[z])){
            currentp[z] = 1
          }
        }
        
        df = data.frame(matrix(NA,length(conditions),3))
        colnames(df) = c("Condition","Average","SE")
        df$Condition = c("+DHPPA","-DHPPA","BHI Control","Tungstate Inhibition")
        
        Pavg = mean(plotdata[i,Pind],na.rm = T)
        Navg = mean(plotdata[i,Nind],na.rm = T)
        Cavg = mean(plotdata[i,Cind],na.rm = T)
        Eavg = mean(plotdata[i,Eind],na.rm = T)
        avgs = c(Pavg,Navg,Cavg,Eavg)
        df$Average = avgs
        
        Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
        Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
        Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
        Ese = sd(plotdata[i,Eind])/sqrt(length(plotdata[i,Eind]))
        ses = c(Pse,Nse,Cse,Ese)
        df$SE = ses
        
        for(j in 1:nrow(df)){
          for(k in 1:ncol(df)){
            if(is.na(df[j,k])){
              df[j,k] = 0
            }
          }
        }
        
        #Dataframe for plotting
        #df$Condition = factor(df$Condition,levels=make.unique(df$Condition))
        
        #BARPLOTHELP
        #Flip order so that it gives a specific order
        df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA","Tungstate Inhibition"))
        colors = c(ctrcolor,negcolor,poscolor,excolor)
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
        
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        if(i == 1){
          cat("Four conditions allow multiple significance bar stylings: Add p-values and significance bars in graphic editors (inkscape or adobe illustrator). \n")
          cat("P-values for each VOC are printed in the Console. \n")
        }
        
        text = paste("P-values for: ",rownames(plotdata)[i],"\n",sep="")
        cat(text)
        ptable
        
        #Automatically adjust the plot height
        if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
          upperlim = max(avgs + ses)+max(avgs + ses)*0.40
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 3){ 
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 4){
          upperlim = max(avgs + ses) + max(avgs + ses)*0.8
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 5){
          upperlim = max(avgs + ses) + max(avgs + ses)*1
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 6){
          upperlim = max(avgs + ses) + max(avgs + ses)*1.2
          p = p + ylim(lowerlim,upperlim)
        }else if(addpvalues == T){
          upperlim = max(avgs + ses) + max(avgs + ses)*1.2
          p = p + ylim(lowerlim,upperlim)
        }
        print(p)
      }
      
    }else{ #two or less control samples - not recommended for publications.
      
      plotdata = mydataforstats[startrow:stoprow,]
      
      PNp = VOC_adjp_PN[startrow:stoprow]
      PEp = VOC_adjp_PE[startrow:stoprow]
      NEp = VOC_adjp_NE[startrow:stoprow]
      
      if(AllowUnadjusted == T){
        PNp = VOC_pvals_PN[startrow:stoprow]
        PEp = VOC_pvals_PE[startrow:stoprow]
        NEp = VOC_pvals_NE[startrow:stoprow]
      }
      
      for(i in 1:nrow(plotdata)){
        currentp = c(PNp[i],PEp[i],NEp[i])
        ptable = data.frame(matrix(NA,length(currentp),2))
        if(AllowUnadjusted == T){
          colnames(ptable) = c("Comparison","UnadjustedP")
          ptable$UnadjustedP = currentp
          ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Group4","Neg_vs_Group4")
        }else{
          colnames(ptable) = c("Comparison","FDRAdjustedP")
          ptable$FDRAdjustedP = currentp
          ptable$Comparison = c("Pos_vs_Neg","Pos_vs_Group4","Neg_vs_Group4")
        }
        
        for(z in 1:length(currentp)){
          if(is.na(currentp[z])){
            currentp[z] = 1
          }
        }
        
        if(length(Cind) >= 1){
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","-DHPPA","BHI Control","Tungstate Inhibition")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          if(length(Cind) == 2){
            Cavg = mean(plotdata[i,Cind],na.rm=T)
          }else{
            Cavg = plotdata[i,Cind]
          }
          Eavg = mean(plotdata[i,Eind],na.rm = T)
          
          avgs = c(Pavg,Navg,Cavg,Eavg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          if(length(Cind) == 2){
            Cse = sd(plotdata[i,Cind])/sqrt(length(plotdata[i,Cind]))
            cat("Warning: This standard error estimate is poor.")
          }else{
            Cse = NA
          }
          Ese = sd(plotdata[i,Eind])/sqrt(length(plotdata[i,Eind]))
          ses = c(Pse,Nse,Cse,Ese)
          df$SE = ses
          
          
          #Dataframe for plotting
          df$Condition = factor(df$Condition,levels=c("BHI Control","-DHPPA","+DHPPA","Tungstate Inhibition"))
          colors = c(ctrcolor,negcolor,poscolor,excolor)
          
        }else{
          df = data.frame(matrix(NA,length(conditions),3))
          colnames(df) = c("Condition","Average","SE")
          df$Condition = c("+DHPPA","-DHPPA","Tungstate Inhibition")
          
          Pavg = mean(plotdata[i,Pind],na.rm = T)
          Navg = mean(plotdata[i,Nind],na.rm = T)
          Eavg = mean(plotdata[i,Eind],na.rm = T)
          avgs = c(Pavg,Navg,Eavg)
          df$Average = avgs
          
          Pse = sd(plotdata[i,Pind])/sqrt(length(plotdata[i,Pind]))
          Nse = sd(plotdata[i,Nind])/sqrt(length(plotdata[i,Nind]))
          ses = c(Pse,Nse,Ese)
          df$SE = ses
          
          
          #Dataframe for plotting
          #df$Condition = factor(df$Condition,levels=make.unique(df$Condition))
          #BARPLOTHELP
          #Flip order so that it gives a specific order
          df$Condition = factor(df$Condition,levels=c("-DHPPA","+DHPPA","Tungstate Inhibition"))
          colors = c(negcolor,poscolor,excolor)
        }
        
        ttl = gsub('\\..*','',rownames(plotdata)[i])
        
        #Generate base plot
        p = ggplot(data=df, aes(x=Condition, y=Average)) +
          geom_bar(stat="identity",fill = colors)
        p = p + theme_minimal()
        p = p + ggtitle(ttl) + ylab(ylab) + xlab("")
        p = p+theme(plot.title = element_text(hjust = 0.5,size = 24))
        p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20))
        p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
        p = p+geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2,
                            position=position_dodge(.9)) 
        
        if(i == 1){
          cat("Four conditions allow multiple significance bar stylings: Add p-values and significance bars in graphic editors (inkscape or adobe illustrator). \n")
          cat("P-values for each VOC are printed in the Console. \n")
        }
        
        text = paste("P-values for: ",rownames(plotdata)[i],"\n",sep="")
        cat(text)
        print(ptable[,1:2])
        
        #Automatically adjust the plot height
        if(length(which(currentp<sigthresh)) == 1){ #if one pairwise test is significant
          upperlim = max(avgs + ses)+max(avgs + ses)*0.2
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 2){ #if two " "
          upperlim = max(avgs + ses)+max(avgs + ses)*0.40
          p = p + ylim(lowerlim,upperlim)
        }else if(length(which(currentp<sigthresh)) == 3){ 
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }else if(addpvalues == T){
          upperlim = max(avgs + ses) + max(avgs + ses)*0.6
          p = p + ylim(lowerlim,upperlim)
        }
        
        print(p)
      }
      
    }
  }
  
  
}
  

####################################################  HEATMAP  #############################################################
#Heatmap (remember to z-score normalize your data before plotting or else the scale will not come out in the correct range (-3,3))

if(PlotHeatmap == T){
  
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
  
  if(length(conditions) <= 3){
    #Set sample colors from user input
    pcacolors = c(rep(pcacolorpos,npos),rep(pcacolorneg,nneg),rep(pcacolorctr,ncontrol))
    
    #Prep PCA Data
    pcadata = mydataforstats[,-which(colnames(mydataforstats) == "RT")]
    
    #Run PCA
    pca = prcomp(t(pcadata))
    PC1 = c(pca$x[,comp1])
    PC2 = c(pca$x[,comp2])
    
    #Get axis bounds
    min1 = min(PC1)
    max1 = max(PC1)
    min2 = min(PC2)
    max2 = max(PC2)
    
    if(pca3d == T){
      require(rgl)
      PC3 = c(pca$x[,comp3])
      min3 = min(PC3)
      max3 = max(PC3)
      rgl::plot3d(PC1,PC2,PC3,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",organism,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5,size=12)
      legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = c("Positive","Negative","Control","Tungstate Inhibition"),col = c(pcacolorpos,pcacolorneg,pcacolorctr,pcacolorex),pch=20,pt.cex = 1.5)
    }else{
      plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",organism,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
      legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = c("Positive","Negative","Control","Tungstate Inhibition"),col = c(pcacolorpos,pcacolorneg,pcacolorctr,pcacolorex),pch=20,pt.cex = 1.5)
    }
    
  }else if(length(conditions) == 4){
    #Set sample colors from user input
    pcacolors = c(rep(pcacolorpos,npos),rep(pcacolorneg,nneg),rep(pcacolorctr,ncontrol),rep(pcacolorex,nextra))
    
    #Prep PCA Data
    pcadata = mydataforstats[,-which(colnames(mydataforstats) == "RT")]
    
    #Run PCA
    pca = prcomp(t(pcadata))
    PC1 = c(pca$x[,comp1])
    PC2 = c(pca$x[,comp2])
    
    #Get axis bounds
    min1 = min(PC1)
    max1 = max(PC1)
    min2 = min(PC2)
    max2 = max(PC2)
    
    if(pca3d == T){
      require(rgl)
      PC3 = c(pca$x[,comp3])
      min3 = min(PC3)
      max3 = max(PC3)
      rgl::plot3d(PC1,PC2,PC3,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",organism,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5,size=12)
      legend(max1 + (0.2*max1),max2 + (0.55*max2),legend = c("Positive","Negative","Control","Tungstate Inhibition"),col = c(pcacolorpos,pcacolorneg,pcacolorctr,pcacolorex),pch=20,pt.cex = 1.5)
    }else{
      plot(PC1,PC2,xlim=c((min1 + (0.2*min1)),(max1 + (0.6*max1))),ylim=c((min2 + (0.2*min2)),(max2 + (0.6*max2))),main=paste("Unsupervised PCA for ",organism,sep=""),col=pcacolors,pch=20,cex=2,cex.lab = 1.5,cex.main = 1.5,cex.axis = 1.5)
      legend(max1 + (0.13*max1),max2 + (0.55*max2),legend = c("Positive","Negative","Control","Tungstate Inhibition"),col = c(pcacolorpos,pcacolorneg,pcacolorctr,pcacolorex),pch=20,pt.cex = 1.5)
    }
    
    #Plot
    
  }
  
  
}


############################  Functional Group Analysis  #############################################################

if(AnalyzeFunctionalGroups == T){
  
  FG = FunctionalGroups(rownames(mydataforstats),flagambiguous)
  
  nFG = length(names(table(FG)))
  
  piecolors = piecolorfnc(nFG)
  
  if(customFGlist == T){
    pie(table(FG),labels = FGlabels,main = paste("Functional Groups for: ",organism,sep=""),col = piecolors)
  }else{
    pie(table(FG),main = paste("Functional Groups for: ",organism,sep=""),col = piecolors)
  }
  
}

#######################################  Prepare Results for Meta Analysis  #############################################################

if(Prep4Meta == T){
  
  #Generate directory, if needed
  str = "/Meta"
  target = paste(MetaDirectory,str,sep="")
  
  if(dir.exists(target)){
    cat("Warning: This directory already exists. A new directory has NOT been created. \n")
  }else{
    dir.create(target)
  }
  
  #Write the processed data matrix to CSV
  setwd(target)
  metafilename = paste(organism,".csv",sep="")
  write.csv(mydataforstats,metafilename)
  
  #Generate reference file with experiment specific parameters / variables
  if(length(conditions) == 2){
    nsamp = c(npos,nneg)
  }else if(length(conditions) == 3){
    nsamp = c(npos,nneg,ncontrol)
  }else if(length(conditions) == 4){
    nsamp = c(npos,nneg,ncontrol,nextra)
  }
  
  reffile = list(conditions,nsamp,organism,MetaDirectory,NormMethod)
  names(reffile) = c("conditions","nsamples","experimentname","dir","NormMethod")
  
  reffilename = paste(organism,"_ref.RData",sep="")
  save(reffile,file=reffilename)
  
}

