#For use by Dr. Barbara Smith Lab and SBHSE at Arizona State University
#Written by: Jarrett Eshima

#Fall 2022

#Most functions have example usage at the end of the section

###################################################################################################
########################################    Z-Score    ############################################
###################################################################################################

#Z-score by feature, simple case. Does not allow for specification of a "reference" level/group

#data is an mxn matrix of numeric values
#featuresarerows: True/False parameter specifying if features/variables/genes/proteins/metabolites are the row names
#removeNA: True/False parameter specifying if NA values should be ignored. Default is true

zscore = function(data,featuresarerows,removeNA){
  
  if(typeof(data) != "double"){
    data2 = matrix(as.numeric(unlist(data)),nrow = nrow(data),ncol=ncol(data))
    rownames(data2) = rownames(data)
    colnames(data2) = colnames(data)
    
    if(typeof(data2) != "double"){
      cat("Warning: values cannot be coerced into numbers.")
    }
    
  }else{
    data2 = data
  }
  
  if(featuresarerows){
    data2 = data2
  }else{
    data2 = t(data2)
  }
  
  if(missing(removeNA)){
    removeNA = T
  }
  
  zscore = matrix(NA,nrow(data2),ncol(data2))
  
  for(i in 1:nrow(data2)){
    for(j in 1:ncol(data2)){
      
      if(length(which(!is.na(data2[i,])))>1){ #This helps deal with sparse observations for plotting purposes, designed to pair with MatrixtoHeatmap()
        zscore[i,j] = (data2[i,j] - mean(data2[i,],na.rm = removeNA))/sd(data2[i,],na.rm = removeNA)
      }else{
        zind = which(!is.na(data2[i,]))
        data2[i,zind] = 0
        zscore[i,j] = data2[i,j]
      }
      
      
    }
  }
  
  return(zscore)
  
}

#Example Usage
# test1 = matrix(rnorm(100,8,2),10,10)
# test2 = matrix(rnorm(100,1,2),10,10)
# test = cbind(test1,test2)
# zscore(test,featuresarerows = T,removeNA = F)
# zscore(test,featuresarerows = T)


###################################################################################################
#######################    Plot a matrix as a heatmap using ggplot2    ############################
###################################################################################################

#Matrix --> heatmap, with reasonable control over graphics (Version 2)

#ADDED: ROW AND COLUMN CLUSTERING

#Dependency: ggplot2
#install.packages("ggplot2")

#Required arguments:
#Data is a matrix containing samples as rows (or columns) and features/genes/proteins/metabolites/variables/etc as columns. Must contain unique row names and column names
#samplesarecols: True/False parameter specifying if the samples/observations are set as the column names
#featureasrows: True/False parameter specifying if the genes/proteins/metabolites/variables should be set as the heatmap rows

#Optional arguments:
#limits: values specifying min/max for heatmap. Accepts string "default" or "Default". If no argument is supplied, the min and max observations are used. 
#colors: values specifying the lower, middle, and upper colors to generate the color scale. Default is blue (low), white (mid), red (high)
#customfeats: list of strings containing all feature names in the original data matrix, in any order (example: customrows = c("Gene4","Gene9","Gene1",etc.))
#customsamps: list of strings containing all sample names in the original data matrix, in any order (example: customrows = c("Patient1","Patient5","Patient7",etc.))
#reverserows: True/False parameter specifying if the rows should be presented in the opposite order (first feature at the top)
#clusterrows: True/False parameter specifying if the rows (features/genes) should be clustered
#clustercols: True/False parameter specifying if the columns (samples/patients/replicates) should be clustered
#clustermethod: Hierarchical clustering method to implement if clusterrows or clustercols is TRUE (Default: ward.D ; see ?hclust)
#clusterdistance: Distance metric to implement if clusterrows or clustercols is TRUE (Default: euclidean ; see ?dist)

#Note: "reverserows" will not work if customfeats or customsamps is provided. If the reverse order is desired, add rev() to customfeat or customsamp, see below:
#customfeats = rev(customfeats)
#customsamps = rev(customsamps)

#Important Note: If customfeats and/or customsamps are provided, hierarchical clustering results are overwritten

#Remember that if you are presenting numeric expression / abundance values, visualization is best on the z-score scale - compare apples to apples!! (Example: limits = c(-4,4))

MatrixtoHeatmap2 = function(data,samplesarecols,featuresasrows,limits,colors,customfeats,customsamps,reverserows,clusterrows,clustercols,clustermethod,distancemethod){
  
  require(ggplot2,quietly = F)
  
  
  if(samplesarecols == T){
    data = data
  }else{
    data = t(data)
  }
  
  samples = colnames(data)
  features = rownames(data)
  
  if(missing(clusterrows)){
    clusterrows = F
  }
  
  if(missing(clustercols)){
    clustercols = F
  }
  
  if(missing(clustermethod)){
    if(clusterrows == T || clustercols == T){
      clustermethod = "ward.D"
      cat("No clustering method selected. Default is euclidean distance. \n")
      cat("Options: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid")
    }
  }
  
  if(missing(distancemethod)){
    if(clusterrows == T || clustercols == T){
      distancemethod = "euclidean"
      cat("No clustering method selected. Default is euclidean distance. \n")
      cat("Options: euclidean, maximum, manhattan, canberra, binary, minkowski")
    }
  }
  
  if(clusterrows == T){
    ord = hclust( dist(data, method = distancemethod), method = clustermethod )$order
    clusterfeats = features[ord]
  }
  
  if(clustercols == T){
    ord2 = hclust( dist(t(data), method = distancemethod), method = clustermethod )$order
    clustersamps = samples[ord2]
  }
  
  if(missing(reverserows)){
    reverserows = F
  }
  
  #Build data.frame for ggplot2
  heatdat = data.frame(matrix(NA,nrow=ncol(data)*nrow(data),ncol=3))
  colnames(heatdat) = c("Module","Response","Value")
  heatdat$Module = rep(samples,nrow(data))
  
  #Get "Response" labels
  resplab = rep(NA,ncol(data)*nrow(data))
  for(i in 1:length(features)){
    
    count2 = i*ncol(data)
    if(i == 1){
      resplab[1:count2] = rep(features[i],ncol(data))
      count1 = count2+1
    }else{
      resplab[count1:count2] = rep(features[i],ncol(data))
      count1 = count2+1
    }
    
  }
  
  heatdat$Response = resplab
  
  
  #Fill ggplot2 matrix
  for(i in 1:nrow(heatdat)){
    
    for(j in 1:nrow(data)){
      
      for(k in 1:ncol(data)){
        
        if(heatdat$Module[i] == colnames(data)[k] && heatdat$Response[i] == rownames(data)[j]){
          
          heatdat$Value[i] = data[j,k]
          
        }
        
      }
      
    }
    
    if((i %% 100) == 0) cat("% Done:",i/nrow(heatdat)*100,"\n")
    
  }
  
  #write.csv(heatdat,"heatmaprows.csv")
  
  #Fix row/column feature order instead of alphabetical (comment this section to leave as alphabetical)
  heatdat$Response = as.character(heatdat$Response)
  heatdat$Module = as.character(heatdat$Module)
  
  if(reverserows == T && ! missing(customfeats)){
    cat("Warning: This function does not allow reversing of rows when a custom feature list is provided. Use rev(customfeats) to accomplish this. \n")
  }
  
  if(reverserows == T && ! missing(customsamps)){
    cat("Warning: This function does not allow reversing of rows when a custom samples list is provided. Use rev(customsamps) to accomplish this. \n")
  }
  
  if(clustercols == T){
    heatdat$Module = factor(heatdat$Module,levels=clustersamps)
  }
  
  if(clusterrows == T){
    heatdat$Response = factor(heatdat$Response,levels=clusterfeats)
  }
  
  if(! missing(customfeats)){
    heatdat$Response = factor(heatdat$Response,levels=customfeats)
    if(clusterrows == T){
      cat("Warning: clustering of rows has been overwritten. \n")
    }
  }
  
  if(! missing(customsamps)){
    heatdat$Module = factor(heatdat$Module,levels=customsamps)
    if(clustercols == T){
      cat("Warning: clustering of columns has been overwritten. \n")
    }
  }
  
  if(missing(limits)){
    limits = c(min(heatdat$Value),max(heatdat$Value))
    cat("Warning: the minimum and maximum observations have been used to set the limits. In the case of some z-score presentations, this may limit interpretation / visualization. \n")
  }
  
  if(missing(colors)){
    colors = c("blue","white","red")
    cat("Default colors selected. \n")
  }
  
  if(reverserows == T){
    heatdat = heatdat[nrow(heatdat):1, ]
  }
  
  #Plotting (add custom graphic options here - title, text size, etc.)
  if(featuresasrows == T){
    p = ggplot(heatdat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
    p = p+scale_fill_gradient2(low=colors[1],mid=colors[2],high=colors[3],limits=c(limits[1],limits[2]))
    p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
    p
  }else{
    p = ggplot(heatdat,aes(Response,Module)) + geom_tile(aes(fill=Value),color="black")
    p = p+scale_fill_gradient2(low=colors[1],mid=colors[2],high=colors[3],limits=c(limits[1],limits[2]))
    p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
    p
  }
  
  
}

#Example Usage
# test = data.frame(matrix(rnorm(200,0,1),ncol = 20,nrow=10)) #this samples from the standard normal (i.e. z-scores)
# rownames(test) = paste("Gene",seq(1,nrow(test)),sep="")
# colnames(test) = paste("Subject",seq(1,ncol(test)),sep="")
# 
# MatrixtoHeatmap2(test,samplesarecols = T,featuresasrows = T)
# MatrixtoHeatmap2(test,samplesarecols = T,featuresasrows = T,clusterrows = T,clustermethod = "ward.D",distancemethod = "euclidean")
# MatrixtoHeatmap2(test,samplesarecols = T,featuresasrows = T,clusterrows = T)
# MatrixtoHeatmap2(test,samplesarecols = T,featuresasrows = T,clusterrows = T,customfeats = c("Gene2","Gene4","Gene6","Gene8","Gene10","Gene1","Gene3","Gene5","Gene7","Gene9") )
# MatrixtoHeatmap2(test,samplesarecols = T,featuresasrows = T,clustercols = T)


###################################################################################################
#######################    Extract Functional Groups from GC-MS Data    ###########################
###################################################################################################

FunctionalGroups = function(list,flagambiguous){
  
  list = tolower(list)
  fg = rep(NA,length(list))
  
  if(flagambiguous == T){
    for(i in 1:length(list)){
      
      nchars = nchar(list[i])
      flag = 0
      
      for(j in 1:nchars){
        
        if(substr(list[i],j,j+2) == "one"){
          fg[i] = "ketone"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+1) == "ol"){
          fg[i] = "alcohol"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+1) == "al"){
          fg[i] = "aldehyde"
          flag = flag+1
        }else if(substr(list[i],j,j+7) == "aldehyde"){
          fg[i] = "aldehyde"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ane"){
          fg[i] = "hydrocarbon"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ene"){
          fg[i] = "double bonded aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+7) == "oic acid"){
          fg[i] = "carboxylic acid"
          flag = flag+1
        }else if(substr(list[i],j,j+6) == "ic acid"){
          fg[i] = "carboxylic acid"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+4) == "ester"){
          fg[i] = "ester"
          flag = flag+1
        }else if(substr(list[i],j,j+9) == "carboxylate"){
          fg[i] = "ester"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+4) == "cyclo"){
          fg[i] = "aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ine"){
          fg[i] = "nitro-aromatic"
          flag = flag+1
        }else if(substr(list[i],j,j+4) == "azole"){
          fg[i] = "nitro-aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ide"){
          fg[i] = "amide"
          flag = flag+1
        }
        
        
        if(flag > 1){
          fg[i] = "Other"
        }else if(flag == 0){
          fg[i] = "Unknown"
        }
        
        
        
      }
      
      flag = 0
      if((i %% 50) == 0) cat("Functional Group Assignment In-Progress, % Done:",i/length(list)*100,"\n")
      
    }
    
    return(fg)
    
  }else{
    
    for(i in 1:length(list)){
      
      nchars = nchar(list[i])
      counter = flag = 0
      
      for(j in 1:nchars){
        
        if(substr(list[i],j,j+2) == "one"){
          fg[i] = "ketone"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+1) == "ol"){
          fg[i] = "alcohol"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+1) == "al"){
          fg[i] = "aldehyde"
          flag = flag+1
        }else if(substr(list[i],j,j+7) == "aldehyde"){
          fg[i] = "aldehyde"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ane"){
          fg[i] = "hydrocarbon"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ene"){
          fg[i] = "double bonded aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+7) == "oic acid"){
          fg[i] = "carboxylic acid"
          flag = flag+1
        }else if(substr(list[i],j,j+6) == "ic acid"){
          fg[i] = "carboxylic acid"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+4) == "ester"){
          fg[i] = "ester"
          flag = flag+1
        }else if(substr(list[i],j,j+9) == "carboxylate"){
          fg[i] = "ester"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+4) == "cyclo"){
          fg[i] = "aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ine"){
          fg[i] = "nitro-aromatic"
          flag = flag+1
        }else if(substr(list[i],j,j+4) == "azole"){
          fg[i] = "nitro-aromatic"
          flag = flag+1
        }
        
        if(substr(list[i],j,j+2) == "ide"){
          fg[i] = "amide"
          flag = flag+1
        }
        
        
        if(flag > 1){
          counter = counter+1
        }
        
        
        
      }
      flag = 0
      if((i %% 50) == 0) cat("Functional Group Assignment In-Progress, % Done:",i/length(list)*100,"\n")
    }
    
    if(counter > 0){
      cat(paste("Warning: There were ",counter," VOCs with an ambiguous functional group assigned. The suffix is used as default.",sep=""),"\n")
    }
    return(fg)
  }
  
  
  
}


###################################################################################################
#########################    Probabilistic Quotient Normalization    ##############################
###################################################################################################

#Reference: https://pubs.acs.org/doi/full/10.1021/ac051632c 

#Steps:
#1) Run PREP()
#2) Run PQN()

PREP2.0 <- function (X) {
  analyte.names <- rownames(X) #Determine the analyte names.
  X <- t(X) #Transpose the data matrix without the first column.
  colnames(X) <- analyte.names #Make the analyte names the names of matrix columns.
  X[is.na(X)] <- 0
  return(as.matrix(X))
}

## Function: Probabilistic Quotient Normalization ##
PQN <- function (X) {
  
  obs <- dim(X)[1]     #Define number of observations.
  dimm <- dim(X)[2]     #Define number of variables (dimensions).
  
  X[0==X] <- 1E-08     #Set zeroes to an arbitrarily small value.
  normRef <- apply(X,2,function(x){median(x[x>1E-08])})     #Define reference spectrum as median for all analytes.
  
  M <- matrix(rep(normRef, each=obs), ncol=length(normRef))     #Convert reference spectrum in matrix equivalent in size to data matrix.
  Q <- X/M     #Divide the concentration of the analyte in each sample by the median value for each analyte.
  Q[0.001 >= Q] <- NA      #Set very small values of "Q" equal to NA for elimination in a subsequent step.   
  for (i in 1:obs) {
    X[i,] <- X[i,]/median(Q[i,], na.rm=TRUE)}     #Divide each analyte in a given sample by the median quotient in that sample.
  
  X[1 >= X] <- 0     #Convert very small normalized values to 0.
  return(matrix(X, nrow=obs, ncol=dimm))
}


#Use
# urinedata=read.csv("C:/Users/Jarrett Eshima/Documents/Research/Spring 2018/Mental Health Clinical Study/Data/Data Processing/Data Files/MentalHealth_MajorityFilterClassifications.csv")
# urinedata=PREP(urinedata)
# nurinedata=PQN(urinedata)

###################################################################################################
#######################    FILTER SILOXANES FROM GC-MS DATA    ####################################
###################################################################################################

#List of VOC names are expected to be the rownames or provided in the FIRST COLUMN
#If userownames = T, then VOC names should be present along the rows. If false, the first column is expected to contain VOC names.

FiltSilox = function(data,userownames){
  contamindex = rep(NA,nrow(data))
  count = 1
  for(i in 1:nrow(data)){
    
    if(userownames == T){
      VOC = rownames(data)[i]
    }else{
      VOC = data[i,1]
    }
    
    
    ncharVOC = nchar(VOC)
    
    for(j in 1:ncharVOC-2){
      
      if(substr(VOC,j,j+5) == "Silane" || substr(VOC,j,j+5) == "silane"){
        contamindex[count] = i
        count = count+1
      }
      
      if(substr(VOC,j,j+7) == "Siloxane" || substr(VOC,j,j+7) == "siloxane"){
        contamindex[count] = i
        count = count+1
      }
      
      if(substr(VOC,j,j+2) == "TMS"){
        contamindex[count] = i
        count = count+1
      }
      
      if(substr(VOC,j,j+4) == "silyl" || substr(VOC,j,j+4) == "Silyl"){
        contamindex[count] = i
        count = count+1
      }
      
      if(substr(VOC,j,j+4) == "TBDMS"){
        contamindex[count] = i
        count = count+1
      }
      
    }
    
  }
  contamindex = contamindex[!is.na(contamindex)]
  cat("Contaminant Indices:",contamindex,"\n")
  
  if(userownames == T){
    cat("Contaminants Removed:",rownames(data)[contamindex],"\n")
  }else{
    cat("Contaminants Removed:",data[contamindex,1],"\n")
  }
  
  
  if(length(contamindex)>0){
    data = data[-contamindex,]
  }
  
  return(data)
}

###################################################################################################
#######################   Quicker working directory   #############################################
###################################################################################################

#Get the working directory from the text copied to your local machines clipboard

#Usage:
#1) Copy path as text to clipboard
#2) setwd(quickdir())

quickdir = function(){
  text=scan("clipboard",what="character")
  if(typeof(text) == "character"){
    tmp = gsub("\\\\","/",text)
    if(length(tmp)>1){
      n = length(tmp)
      str = NA
      for(i in 1:n){
        if(i == 1){
          str = tmp[i]
        }else{
          str = paste(str," ",tmp[i],sep="",collapse = "")
        }
        
      }
      return(str)
    }else{
      return(tmp)
    }
  }else{
    cat("Cannot find characters / string copied to clipboard")
  }
}


###################################################################################################
#######################   Reverse your matrix row order   #########################################
###################################################################################################


bottomsup = function(data){
  
  rnames = rownames(data)
  
  for(i in 1:nrow(data)){
    
    if(i == 1){
      tmp = rbind(data[nrow(data),],data[nrow(data)-i,])
    }else{
      tmp = rbind(tmp,data[nrow(data)-i,])
    }
    
  }
  rownames(tmp) = rev(rnames)
  return(tmp)
  
}

data = matrix(rnorm(100,0,1),10,10)
rownames(data) = paste("Gene",seq(1,10),sep="")
colnames(data) = paste("Sample",seq(1,10),sep="")
newdata = bottomsup(data)
print(newdata)

###################################################################################################
#############   Trick to assist with workspace management during loops ############################
###################################################################################################

#Trick to access atomic vectors and assign variables to workspace in advanced loops

assignatomic = function(x,y,z){
  x[[y]] = z
  return(x)
}

###################################################################################################
##########    PARSE and ORDER functions for data cleanup    #######################################
###################################################################################################

#Function Inputs
#data: m x n matrix containing any numeric values (omic and big data sets) - subjects/samples can be rows or columns, code will automatically adjust
#samplenames: character vector specifying the sample names in the desired order (this can easily be grabbed from the clinical data matrix)
#orderrow: This is a true / false value. If True, sample names are expected to be the row names of the data matrix. If false, sample names are expected to be the column names.

#The order function
ORDER = function(data,samplenames,orderrow){
  
  assignatomic = function(x,y,z){
    x[[y]] = z
    return(x)
  }
  
  nr = nrow(data)
  nc = ncol(data)
  tmp = tmp2 = matrix(NA,nr,nc)
  
  if(orderrow == T){
    for(i in 1:nrow(data)){
      for(j in 1:length(samplenames)){
        if(rownames(data)[i] == samplenames[j]){
          tmp[j,] = data[i,]
        }
      }
    }
    rownames(tmp) = samplenames
    colnames(tmp) = colnames(data)
  }else if(orderrow == F){
    for(i in 1:ncol(data)){
      for(j in 1:length(samplenames)){
        if(colnames(data)[i] == samplenames[j]){
          tmp[,j] = data[,i]
        }
      }
    }
    colnames(tmp) = samplenames
    rownames(tmp) = rownames(data)
  }else{
    cat("Please specify the orderrow arguement")
  }
  
  assignatomic(.GlobalEnv,paste("orderdata"),tmp)
  
}


#Function Inputs
#data: m x n matrix containing any numeric values (omic and big data sets) - subjects/samples can be rows or columns, code will automatically adjust
#groups: Strings specifying the experimental groups used to parse the data. Note that these names must appear somewhere in the sample names or clinical data file
#useclinicaldata: should be set as logical (TRUE/FALSE) - if this is true you will need to provide the clinical data and handle inputs
#clinicaldata: sample/patient meta data containing unique ID codes, instrument type, tissue type, sex, age, disease subtype, etc.
#handle: Name of the category which the groups belong. For example if your groups were "Male,Female", then the handle would be "Sex". Make sure that this matches the column name provided in the clinicaldata object.
#useshortstring: TRUE / FALSE - determines whether to use a string with length 3 for searching your data matrix column/row names. May be useful when abbreviating samples names (i.e. Sample1Pos and Sample1Neg). An exact string match is required, so check capitalization.

#Note: the handle argument is case sensitive. Make sure to double check you use the right capitalization. Please do not use row/column names with multiple instances of the handle... only the last will be selected. 

#Crucial: 
#Data and clinicaldata must be in the same sample order. Use the ORDER function to rearrange data matrix so that data and clinicaldata are presented in the same order. 
#Additionally, if subjects/samples are set as the columns in the data matrix then they will need to be rows in the clinicaldata matrix (vice versa will work too). (you can transpose t() your matrix to simplify this requirement)


#The parse function
PARSE = function(data,groups,useclinicaldata,clinicaldata,handle,useshortstring){
  
  require(stringr)
  assignatomic = function(x,y,z){
    x[[y]] = z
    return(x)
  }
  
  ngroup = length(groups)
  nr = nrow(data)
  nc = ncol(data)
  
  #Create containers - this is required in order to output to global environment
  for(i in 1:ngroup){
    tmpname = paste(groups[i],sep="")
    assign(tmpname,matrix(NA,nr,nc))
  }
  
  
  if(useclinicaldata == T){
    strlen = nchar(handle)
    
    if(nrow(clinicaldata)>0 && ncol(clinicaldata)>0){
      nr_clinical = nrow(clinicaldata)
      nc_clinical = ncol(clinicaldata)
    }else{
      cat("Please provide clinical data. \n")
    }
    
    
    colindex = NA
    for(j in 1:nc_clinical){
      test = colnames(clinicaldata)[j]
      width = nchar(test)
      for(k in 1:(width-strlen)){
        if(str_sub(test,k,k+strlen-1) == handle){
          colindex = j #Column matching the handle string
        }
      }
      
    }
    rowindex = NA
    if(is.na(colindex)){
      for(j in 1:nr_clinical){
        test = colnames(clinicaldata)[j]
        width = nchar(test)
        for(k in 1:(width-strlen)){
          if(str_sub(test,k,k+strlen-1) == handle){
            rowindex = j #row matching the handle string, if no column matches
          }
        }
        
      }
    }
    
    if(is.na(rowindex) && is.na(colindex)){
      cat("Your handle / group name could not be located in the clinical data matrix. Double check that the handle argument has the correct capitalization. \n")
    }
    
    #Parse
    for(i in 1:ngroup){
      tmpname = paste(groups[i],sep="")
      tmpmatrix = get(tmpname)
      
      if(is.na(rowindex)){
        for(j in 1:nr_clinical){
          if(clinicaldata[j,colindex] == groups[i]){
            tmpmatrix[,j] = data[,j] #This only works if data and clinicaldata are in the same sample order.
          }
        }
        colnames(tmpmatrix) = colnames(data)
        rownames(tmpmatrix) = rownames(data)
      }else if(is.na(colindex)){
        for(j in 1:nc_clinical){
          if(clinicaldata[rowindex,j] == groups[i]){
            tmpmatrix[j,] = data[j,] #This only works if data and clinicaldata are in the same sample order.
            
          }
        }
        rownames(tmpmatrix) = rownames(data)
        colnames(tmpmatrix) = colnames(data)
      }
      
      #Clean up NA values in group matrices
      
      if(is.na(rowindex)){
        remove = which(is.na(colSums(tmpmatrix)))
        tmpmatrix = tmpmatrix[,-remove]
      }else{
        remove = which(is.na(colSums(t(tmpmatrix))))
        tmpmatrix = tmpmatrix[-remove,]
      }
      
      
      assignatomic(.GlobalEnv,paste(groups[i]),tmpmatrix) #Function output to global environment (I guess this is atypical but oh well)
      
    }
    
  }else{
    
    #Create containers - this is required in order to output to global environment
    for(i in 1:ngroup){
      tmpname = paste(groups[i],sep="")
      assign(tmpname,matrix(NA,nr,nc))
    }
    
    
    for(i in 1:ngroup){
      
      identifier = groups[i]
      
      if(useshortstring == F){
        idlen = nchar(identifier)
      }
      
      if(useshortstring == T){
        shortid = str_sub(identifier,1,3)
        sidlen = nchar(shortid) #this will always be 3
      }
      
      keeprow = keeprow2 = rep(NA,nr)
      for(j in 1:nrow(data)){
        width = nchar(rownames(data)[j])
        
        if(useshortstring == F){
          for(m in 1:width-idlen+1){
            if(str_sub(rownames(data)[j],m,m+idlen-1) == identifier){
              keeprow[j] = j 
            }
          }
        }
        
        if(useshortstring == T){
          for(n in 1:width-sidlen+1){
            if(str_sub(rownames(data)[j],n,n+sidlen-1) == shortid){
              keeprow2[j] = j
            }
          }
        }
        
        keeprow = keeprow[!is.na(keeprow)]
        keeprow2 = keeprow2[!is.na(keeprow2)]
        if(useshortstring == F){
          parsedata = data[,keeprow]
        }else{
          parsedata = data[,keeprow2]
        }
        
        assignatomic(.GlobalEnv,paste(groups[i]),parsedata)
        
      }
      
      if(useshortstring == F){
        if(length(keeprow) != 0){
          parsedata = data[keeprow,]
          assignatomic(.GlobalEnv,paste(groups[i]),parsedata)
        }
      }
      
      if(useshortstring == T){
        if(length(keeprow2) != 0){
          parsedata = data[keeprow2,]
          assignatomic(.GlobalEnv,paste(groups[i]),parsedata)
        }
      }
      
      if(length(keeprow) == 0 && length(keeprow2) == 0){
        
        keepcol = keepcol2 = rep(NA,nc)
        for(k in 1:ncol(data)){
          width2 = nchar(colnames(data)[k])
          
          if(useshortstring == F){
            for(o in 1:width2-idlen+1){
              if(str_sub(colnames(data)[k],o,o+idlen-1) == identifier){
                keepcol[k] = k
              }
            }
          }
          
          if(useshortstring == T){
            for(p in 1:width2-sidlen+1){
              if(str_sub(colnames(data)[k],p,p+sidlen-1) == shortid){
                keepcol2[k] = k
              }
            }
          }
          
          if(useshortstring == F){
            keepcol = keepcol[!is.na(keepcol)]
            parsedata = data[,keepcol]
          }else{
            keepcol2 = keepcol2[!is.na(keepcol2)]
            parsedata = data[,keepcol2]
          }
          
          assignatomic(.GlobalEnv,paste(groups[i]),parsedata)
        }
        
        if(useshortstring == F){
          if(length(keeprow) == 0 && length(keepcol) == 0 ){
            cat("Double check your group names (and capitalization), column and row names did not contain a matching string. \n")
          }
        }else{
          if(length(keeprow2) == 0 && length(keepcol2) == 0){
            cat("Double check your group names (and capitalization), column and row names did not contain a matching string. \n")
          }
        }
        
      }
      
    }
    
  }
  
}


#Example usage
#Data
nobs = 1000 #proteins/transcripts/metabolites/etc
nsub = 200 #subjects
data = matrix(rnorm(nsub*nobs,0,1),nobs,nsub)
colnames(data) = paste("Subject",seq(1,nsub),sep="")
rownames(data) = paste("Gene",seq(1,nobs),sep="")

#Samplenames
samplenames = paste("Subject",seq(1,nsub),sep="")
samplenames[1] = paste("Subject100")
samplenames[100] = paste("Subject1")
orderrow = F

ORDER(data,samplenames,orderrow = F)

#Clinical Data
clinicaldata = matrix(NA,nsub,10)
colnames(clinicaldata) = c("Sample","Instrument","Accession","Subject","Disease","Duration","Sex","AgeOnset","Symptoms","Genotype")
clinicaldata = data.frame(clinicaldata)
clinicaldata$Sample = colnames(data)
clinicaldata$Instrument = "NovaSeq 2500"
clinicaldata$Accession = "GSEXXXXXXX"
clinicaldata$Subject = "NZO43WQMN2"
clinicaldata$Disease = sample(c("ALS","FTD"),nsub,replace=T)
clinicaldata$Duration = sample(seq(1,100),nsub,replace=T)
clinicaldata$Sex = sample(c("Male","Female"),nsub,replace=T)
clinicaldata$AgeOnset = sample(seq(30,80,1),nsub,replace=T)
clinicaldata$Symptoms = sample(c("Limb Weakness","Upper","Bulbar","Generalized"),nsub,replace = T)
clinicaldata$Genotype = sample(c("C9orf72","SOD1","None"),nsub,replace=T)

#Required arguments
useclinicaldata = T

#Test 1
groups = c("Male","Female")
handle = "Sex"

#Test 2
groups = c("ALS","FTD")
handle = "Disease"

#Test 3
groups = c("C9orf72","SOD1","None")
handle = "Genotype"

#Parse data according to the groups specified
PARSE(data,groups,useclinicaldata = T,clinicaldata,handle)


#Example for PARSE function with no reference file

#Data
nobs = 1000 #proteins/transcripts/metabolites
nsub = 200 #subjects
data = matrix(rnorm(nsub*nobs,0,1),nobs,nsub)
colnames(data) = paste("Subject",seq(1,nsub),sample(c("Positive","Negative"),nsub,replace = T),sep="")
rownames(data) = paste("Gene",seq(1,nobs),sep="")
groups = c("Positive","Negative")
useshortstring = F

PARSE(data,groups,useclinicaldata = F,useshortstring = T)
PARSE(data,groups,useclinicaldata = F,useshortstring = F)


###################################################################################################
#######################    Plot a matrix as a heatmap using ggplot2    ############################
###################################################################################################

#Matrix --> heatmap, with reasonable control over graphics (Version 1 - does not support clustering capabilities)

#Dependency: ggplot2
#install.packages("ggplot2")

#Required arguments:
#Data is a matrix containing samples as rows (or columns) and features/genes/proteins/metabolites/variables/etc as columns. Must contain unique row names and column names
#samplesarecols: True/False parameter specifying if the samples/observations are set as the column names
#featureasrows: True/False parameter specifying if the genes/proteins/metabolites/variables should be set as the heatmap rows

#Optional arguments:
#limits: values specifying min/max for heatmap. Accepts string "default" or "Default". If no argument is supplied, the min and max observations are used. 
#colors: values specifying the lower, middle, and upper colors to generate the color scale. Default is blue (low), white (mid), red (high)
#customfeats: list of strings containing all feature names in the original data matrix, in any order (example: customrows = c("Gene4","Gene9","Gene1",etc.))
#customsamps: list of strings containing all sample names in the original data matrix, in any order (example: customrows = c("Patient1","Patient5","Patient7",etc.))
#reverserows: True/False parameter specifying if the rows should be presented in the opposite order (first feature at the top)
   #Note: "reverserows" will not work if customfeats or customsamps is provided. If the reverse order is desired, add rev() to customfeat or customsamp, see below:
   #customfeats = rev(customfeats)
   #customsamps = rev(customsamps)

#Remember that if you are presenting numeric expression / abundance values, visualization is best on the z-score scale - compare apples to apples!! (Example: limits = c(-4,4))

MatrixtoHeatmap = function(data,samplesarecols,featuresasrows,limits,colors,customfeats,customsamps,reverserows){
  
  require(ggplot2,quietly = F)
  
  
  if(samplesarecols == T){
    data = data
  }else{
    data = t(data)
  }
  
  samples = colnames(data)
  features = rownames(data)
  
  #Build data.frame for ggplot2
  heatdat = data.frame(matrix(NA,nrow=ncol(data)*nrow(data),ncol=3))
  colnames(heatdat) = c("Module","Response","Value")
  heatdat$Module = rep(samples,nrow(data))
  
  #Get "Response" labels
  resplab = rep(NA,ncol(data)*nrow(data))
  for(i in 1:length(features)){
    
    count2 = i*ncol(data)
    if(i == 1){
      resplab[1:count2] = rep(features[i],ncol(data))
      count1 = count2+1
    }else{
      resplab[count1:count2] = rep(features[i],ncol(data))
      count1 = count2+1
    }
    
  }
  
  heatdat$Response = resplab
  
  
  #Fill ggplot2 matrix
  for(i in 1:nrow(heatdat)){
    
    for(j in 1:nrow(data)){
      
      for(k in 1:ncol(data)){
        
        if(heatdat$Module[i] == colnames(data)[k] && heatdat$Response[i] == rownames(data)[j]){
          
          heatdat$Value[i] = data[j,k]
          
        }
        
      }
      
    }
    
    if((i %% 100) == 0) cat("% Done:",i/nrow(heatdat)*100,"\n")
    
  }
  
  #write.csv(heatdat,"heatmaprows.csv")
  
  #Fix row/column feature order instead of alphabetical (comment this section to leave as alphabetical)
  heatdat$Response = as.character(heatdat$Response)
  heatdat$Module = as.character(heatdat$Module)
  
  if(reverserows == T && ! missing(customfeats)){
    cat("Warning: This function does not allow reversing of rows when a custom feature list is provided. Use rev(customfeats) to accomplish this.")
  }
  
  if(reverserows == T && ! missing(customsamps)){
    cat("Warning: This function does not allow reversing of rows when a custom samples list is provided. Use rev(customsamps) to accomplish this.")
  }
  
  
  if(reverserows == T){
    heatdat = heatdat[nrow(heatdat):1, ]
  }
  
  if(missing(customfeats)){
    heatdat$Response = factor(heatdat$Response,levels=make.unique(heatdat$Response))
  }else{
    heatdat$Response = factor(heatdat$Response,levels=customfeats)
  }
  
  if(missing(customsamps)){
    heatdat$Module = factor(heatdat$Module,levels=make.unique(heatdat$Module))
  }else{
    heatdat$Module = factor(heatdat$Module,levels=customsamps)
  }
  
  if(missing(reverserows)){
    reverserows = F
  }
  
  
  if(missing(limits)){
    limits = c(min(heatdat$Value),max(heatdat$Value))
    cat("Warning: the minimum and maximum observations have been used to set the limits. In the case of some z-score presentations, this may limit interpretation / visualization. \n")
  }
  
  
  if(missing(colors)){
    colors = c("blue","white","red")
    cat("Default colors selected. \n")
  }
  
  #Plotting (add custom graphic options here - title, text size, etc.)
  if(featuresasrows == T){
    p = ggplot(heatdat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
    p = p+scale_fill_gradient2(low=colors[1],mid=colors[2],high=colors[3],limits=c(limits[1],limits[2]))
    p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
    p
  }else{
    p = ggplot(heatdat,aes(Response,Module)) + geom_tile(aes(fill=Value),color="black")
    p = p+scale_fill_gradient2(low=colors[1],mid=colors[2],high=colors[3],limits=c(limits[1],limits[2]))
    p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
    p
  }
  
  
}

#Example Usage
test = data.frame(matrix(rnorm(200,0,1),ncol = 20,nrow=10)) #this samples from the standard normal (i.e. z-scores)
rownames(test) = paste("Gene",seq(1,nrow(test)),sep="")
colnames(test) = paste("Subject",seq(1,ncol(test)),sep="")

MatrixtoHeatmap(test,samplesarecols = T,featuresasrows = F,reverserows = T)
MatrixtoHeatmap(test,samplesarecols = T,featuresasrows = T,limits = c(-4,4),colors = c("pink","white","purple"),reverserows = F)


#Custom Order Example
#This is designed to pipe to hierarchical clustering analysis (hclust) or similar methods
customfeats = c("Gene4","Gene8","Gene6","Gene10","Gene3","Gene5","Gene9","Gene2","Gene7","Gene1")
MatrixtoHeatmap(test,samplesarecols = T,featuresasrows = T,limits = c(-3,3),colors = c("green","yellow","blue"),customfeats)
MatrixtoHeatmap(test,samplesarecols = T,featuresasrows = T,limits = c(-3,3),colors = c("green","yellow","blue"),rev(customfeats))

#Axis and legend labels should be cleaned up externally (inkscape or Illustrator)
