########################################################################################
## VAF Data analysis for Alex and James
##
## Purpose: 1) Determine Bonferroni corrected prediction intervals 
##          of background noise from a set control samples
##          (i.e., background samples).
##          
##          
##          2)Compare VAF results from experimental samples
##          (i.e., patient samples) to bacground and identify those
##          that are above background.
##
########################################################################################
## Inputs:
##    User Defned:
##        alphaLevel = the alpha level used to create the 
##        prediction interval.
##
##        sampleTypeControl = vector of IDs for control samples
##
##        sampleTypeExp = vector of IDs for experimental samples
##
##        outFileLoc = location where outputfile should be sent
##
##    Datasets: 
##        df = location of VAF csv dataset with following format,
##        
##
##   Chrom       Loc AO    DP         VAF WT Var Sample ID            Info Expected
##   chr1 115227814 19 45714 0.000415628  G   A      4  1 Orig Background       No
##      
##        sites = location of dataset with site information for 
##        location mapping. Note: no headers are used and the format should be
##        three columns of,
##
##        "site", "chrom#:base pair location start","chrom#:base pair location end"
##
##        example,
##
##        JAK2   chr9:5073733   chr9:5073887
##        TP53-1  chr17:7577504  chr17:7577635
##
########################################################################################
##
## Outputs:
##          CSV file with prediction intervals for backgound noise
##          and how VAF from patient samples compare to background 
##          for a given location.
##
########################################################################################

start_time <- Sys.time()

###################################################################
## Inputs:
##########

## set the alpha level to be used

alphaLevel<-0.001

## create vectors of id numbers that are control or experimental samples
sampleTypeControl<-seq(from=1,to=22,by=1)
sampleTypeExp<-seq(from=23,to=34,by=1)

## set location for where the output dataset should be sent
outFileLoc<-"//ucdenver.pvt/som/CC1/Biostatistics/projects/DegregoriJames/Alex/results/"

## Read Data
df<-read.table("//ucdenver.pvt/som/CC1/Biostatistics/projects/DegregoriJames/Alex/liggettVAFdataWithSampleInfo.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)

## read in a list of the 32 different sites
sites<-read.table("//ucdenver.pvt/som/CC1/Biostatistics/projects/DegregoriJames/Alex/siteSpecificLocations.csv",sep=",",header=FALSE,stringsAsFactors = FALSE)

###################################################################
###################################################################



## create an indicator for background (0) and exerimental samples (1)
df$type<-ifelse(df$ID %in% sampleTypeControl,0, ifelse(df$ID %in% sampleTypeExp,1,NA))


## extract the site information
colnames(sites)<-c("site","start","end")

sites$start<-gsub(".*:","",sites$start)

sites$end<-gsub(".*:","",sites$end)


################################################################
## determine the site for each location (Loc) based on the base pair location
################################################################

df$site<-vapply(1:nrow(df),function(x){
  
  location<-df$Loc[x]
  
  sites$site[which(sites$start <= location & sites$end >= location )]
  
},FUN.VALUE=character(1))





####################################################################################################
## Purpose: summaryVAF function determines the summary information for a set of control 
##          (i.e., background) samples. The function will also fill in locations that are missing
##          in the dataset and set their VAF to 0 , which can be used for control or 
##          experimental samples (e.g., patient) to create a dataset with a complete set of locations
##          
## Inputs:
##
## completeLocAndSite = a vector of the complete set of locations with site 
##                      information as a dataframe with a "location" and "site" column
##
## data = a dataframe to be processed with the format of the following columns
##  Chrom       Loc AO    DP         VAF WT Var Sample ID            Info Expected type  site
##  1  chr1 115227814 19 45714 0.000415628  G   A      4  1 Orig Background       No    0 TIIIA
##
## overallSummary = logical TRUE or FALSE for whether the summary for a set of samples should be provided
##                  at a given locaton.
####################################################################################################





summaryVAF<-function(completeLocSiteAndWT,data,overallSummary=TRUE){
  
  ## list the different var options
  varOpts<-c("A","C","T","G")
  
  
  ## loop over all possible base pair locations
  summaryVAF<-lapply(completeLocSiteAndWT$location,function(location){
    
    
    ## create a subset of the data for that location
    tmp<-data[data$Loc==location,]
    
    
    ## determine the number of records for a given location
    n<-length(unique(tmp$ID[tmp$Loc==location]))
    
    ## check to see if all the records exist. There should be a maximum of 3*# of unique IDs at each site.
    ## for those that don't meet the maximum number create records of 0 for these samples.
    ## This is because the input dataset does not include records for locations where the VAF is 0.
    
    ##number of background
    maxControlSamp<-length(unique(data$ID))*3
    
    if(n!=maxControlSamp){
      
      ## loop over all records for a given id
      dups<-lapply(unique(data$ID), function(x){
        
        ## determine the number of Var for a given location and ID
        numVar<-length(unique(tmp$Var[tmp$ID==x]))
        
        if(!(x %in% unique(tmp$ID)) & nrow(tmp)==0){## if no records exist for a given location and no ids
          
          
          ###################################################################
          ###################################################################
          ## NOTE: for locations that do not exist the WT is not known.
          ## Currently these locations are set to 0 for all 4 base options 
          ## with an unknown "unk" WT indicator in the data
          ###################################################################
          ###################################################################
          
          ## create a copy for a location
          dup<-data[1,]
          
          Chrom<-as.character(completeLocSiteAndWT$Chrom[completeLocSiteAndWT$location==location])
          
          ## update the information to zero
          dup[,c(1,2,3,4,5,6,7,8,9,13)]<-c(Chrom,location,0,0,0,as.character(completeLocSiteAndWT$WT[completeLocSiteAndWT$location==location]),0,-999,x,as.character(completeLocSiteAndWT$site[completeLocSiteAndWT$location==location]))
          
          ## make four copies for the location with the differnt Var options
          dup<-dup[rep(seq_len(nrow(dup)), 4), ]
          
          dup$Var<-varOpts
          
          dup$VAF<-as.numeric(dup$VAF)
          
          
          return(dup)
          
        }else if(!(x %in% unique(tmp$ID)) & nrow(tmp)>=1){## if no records exist for a given location and id but some do a given lcoation so that the WT is known create three copies with zeros
          
          ## create a copy for a location
          dup<-tmp[1,]
          
          ## update the information to zero
          dup[,c(3,4,5,7,8,9)]<-c(0,0,0,0,-999,x)
          
          ## make three copies for the location witht the differnt Var options that are different from the WT
          dup<-dup[rep(seq_len(nrow(dup)), 3), ]
          
          ## set the vars to options other than the WT
          dup$Var<-varOpts[!(varOpts %in% dup$WT)]
          
          
          return(dup)
          
          
          
          
        }else if(numVar == 1){ ## if one record exists for an ID add in the other two var options
          
          ## get the subset of info for a given ID
          dup<-tmp[tmp$ID==x,]
          
          ## create a copy for a location/ID and pdate the information to zero
          tmpRep<-dup[1,]
          tmpRep[,c(3,4,5,7,8,9)]<-c(0,0,0,0,-999,x)
          
          ## make two copies for the location with the differnt Var options that are different from the WT
          tmpRep<-tmpRep[rep(seq_len(nrow(tmpRep)), 2), ]
          
          # set the var what wasn't observed and was the WT  
          tmpRep$Var<-varOpts[!(varOpts %in% c(tmpRep$WT,dup$Var))]
          
          dup<-rbind(dup,tmpRep)
          
          return(dup)
          
        }else if(numVar == 2){ ## if two records exists for an ID add in the other two
          
          ## get the subset of info for a given ID
          dup<-tmp[tmp$ID==x,]
          
          ## create a copy for a location/ID and pdate the information to zero
          tmpRep<-dup[1,]
          tmpRep[,c(3,4,5,7,8,9)]<-c(0,0,0,0,-999,x)
          
          ## make two copies for the location with the differnt Var options that are different from the WT
          tmpRep<-tmpRep[rep(seq_len(nrow(tmpRep)), 1), ]
          
          # set the var what wasn't observed and was the WT  
          tmpRep$Var<-varOpts[!(varOpts %in% c(tmpRep$WT,dup$Var))]
          
          dup<-rbind(dup,tmpRep)
          
          return(dup)
          
        }else if(numVar == 3){ ## if three records exist the return the records
          
          ## get the subset of info for a given ID
          dup<-tmp[tmp$ID==x,]
          
          return(dup)
          
        } else{ ## otherwise throw an error
          
          #stop(print("more than 3 records exist for an ID at a given location"))
          
          stop(print(tmp))
        }
        
        
      })
      
      ## rbind to the dataset and overwrite the existing tmp dataset
      tmp<-do.call("rbind",dups) 
      
      
    }
    
    
    if(overallSummary){  ## check if the summary is wanted for the set of samples
      
      ## find the mean and standard deviation for each VAF at a given location and BP
      finalSum<-lapply(unique(tmp$Var),function(x){
        
        avg<-mean(tmp$VAF[tmp$Var==x])
        
        stdev<-sd(tmp$VAF[tmp$Var==x])
        
        ## deterine the chrom, site and WT information for a given location
        Chrom<-unique(tmp$Chrom[1])
        
        site<-unique(tmp$site[1])
        
        WT<-unique(tmp$WT)
        
        
        
        finalSum<-data.frame("location"=location,"average"=avg,"stdev"=stdev,"Chrom"=Chrom,"site"=site, "WT"=WT,"Var"=x)
        
        return(finalSum)
        
      })
      
      
      finalSum<-do.call("rbind", finalSum)
      
      return(finalSum)
      
    }else{ ## otherwise return the complete set of records with 0 values filled in for WT var combos that didn't have a VAF
      
      return(tmp)
    }
    
  })
  
  summaryVAF<-do.call("rbind",summaryVAF)
  
  return(summaryVAF)
  
}

####################################################################################################
####################################################################################################








#######################################################################################
## Control Group
#######################################################################################



## determine the average VAF for the control group for a given site and BP location
control<-df[df$type==0,]

numControlSamp<-length(unique(control$ID))

## list the different var options
varOpts<-c("A","C","T","G")


## identify all possible locations whether or not they existed in the dataset
completeLoc<-lapply(1:nrow(sites),function(x){
  
  
  ## create a sequnce for all the base pair locatons for a given site
  completeLoc<-seq(from=sites$start[x],to=sites$end[x],by=1)
  
  completeLoc<-data.frame("location"=completeLoc,"site"=rep(sites$site[x],length(completeLoc)))
  
  return(completeLoc)
  
  
})

## combine the lists into a vector
completeLoc<-do.call("rbind",completeLoc)


completeLoc$WT<-vapply(1:nrow(completeLoc),function(x){
  
  location<-completeLoc$location[x]
  
  if(any(df$Loc==location)){
    
  unique(df$WT[df$Loc==location])
    
  }else{
    "ukn"
  }
},FUN.VALUE=character(1))


completeLoc$Chrom<-vapply(1:nrow(completeLoc),function(x){
  
  location<-completeLoc$location[x]
  
  if(any(df$Loc==location)){
    
    unique(df$Chrom[df$Loc==location])
    
  }else{
    "ukn"
  }
},FUN.VALUE=character(1))



summaryControls<-summaryVAF(completeLocSiteAndWT = completeLoc,data=control,overallSummary=TRUE)

summaryControls$locationVAR<-paste(summaryControls$location,summaryControls$Var,sep=".")


###############################################################################
## Check if any locations are missing from the summary of the control dataset
###############################################################################

## create a sequence for each site

missingControls<-lapply(1:nrow(sites),function(x){
  

  ## create a sequnce for all the base pair locatons for a given site
  completeLoc<-seq(from=sites$start[x],to=sites$end[x],by=1)
  
  siteName<-sites$site[x]
  
  ## check to see what locations are missing from the base pair sequence
  missing<-completeLoc[!(completeLoc %in% summaryControls$location[summaryControls$site==siteName])]
  
  ## for each location create a row to be added to the summaryControls
  
  return(data.frame("location"=missing,  "numWithData"=rep(0,times=length(missing)),"average"=rep(0,times=length(missing)),"stdev"=rep(0,times=length(missing)), "Chrom"=rep(summaryControls$Chrom[summaryControls$site==siteName][1],times=length(missing)),  "site"=rep(siteName,times=length(missing))))
  
  
  
})

missingControls<-do.call("rbind",missingControls)


#######################################################################################
#######################################################################################




#######################################################################################
## Experimental Group
#######################################################################################



## determine the average VAF for the experimental  group for a given site and BP location
expGrp<-df[df$type==1,]

numExpSamp<-length(unique(expGrp$ID))

summaryExp<-summaryVAF(completeLocSiteAndWT = completeLoc,data=expGrp,overallSummary = FALSE)

summaryExp$locationVAR<-paste(summaryExp$Loc,summaryExp$Var,sep=".")

#########################################################################################
#########################################################################################



########################################################################################
## Conduct a one-sided t-test for the control group, where H0: mu (controls) = 0
########################################################################################



## the different var options
varOpts<-c("A","C","T","G")


## determine all the different location and Var options
allLocAndVar<-lapply(varOpts,function(x){
  
  paste(completeLoc$location,x,sep=".")
  
})


allLocAndVar<-do.call("c",allLocAndVar)


## keep only the location and VAR combiniations that exist in the experimental ad control groups

allLocAndVar<-allLocAndVar[allLocAndVar %in% c(summaryControls$locationVAR,summaryExp$locationVAR)]


## obtain the p-values for each location using a t test statistic for all possible locations and VARs

testResults<-lapply(allLocAndVar,function(x){

  
  contTmp<-summaryControls[summaryControls$locationVAR==x,]
  
  if(nrow(contTmp)==0){## check to make sure that there are records for a location
    
    return(data.frame("location"=substr(x,start=1,stop=nchar(x)-2),"mean"=0,"SD"=0,"tStat"=NA,"df"=NA,"pval"=1,"site"=NA,"locVar"=x,"Var"=substr(x,start=nchar(x),stop=nchar(x)),"WT"="ukn"))
    
  }else{
  

  ## one sided one sample t-test statistic to determine the probability that the VAF in the control group is greater than 0 
  
  tStat<-(contTmp$average-0)/(contTmp$stdev/sqrt(length(unique(control$ID))))
  

  
  
  degFree<-length(unique(control$ID))-1
  
  
  ## Determine the probability of seeing a t-statistic greater than or equal to that specified
  pval<-pt(q=tStat,df=degFree,lower.tail = FALSE)
  
  ## site info
  site<-unique(contTmp$site)
  
  loc<-unique(contTmp$location)
  
  WT<-unique(contTmp$WT)
  
  Var<-unique(contTmp$Var)
  
  meanVAF<-unique(contTmp$average)
  
  sdVAF<-unique(contTmp$stdev)
  
 
  
  return(data.frame("location"=loc,"mean"=meanVAF,"SD"=sdVAF,"tStat"=tStat,"df"=degFree,"pval"=pval,"site"=site,"locVar"=x,"Var"=Var,"WT"=WT))
    
  }
    
})


testResults<-do.call("rbind",testResults)

#####################################################
## create bonforroni corrected prediction intervals
#####################################################

## calculate the t-value for bonforroni's correction
## number of tests based on the number of different locations assessed
degFree<-length(allLocAndVar)-1

num<-as.numeric(nrow(testResults))

## alpha level
levBon<-1-(alphaLevel/num)



## compute the t-value for a given alpha
tval<-qt(levBon,df=degFree)

## compute the lower bound and set the minimum to 0
testResults$LBonPI<-testResults$mean -(tval*(testResults$SD*sqrt(1+(1/numControlSamp))))

testResults$LBonPI<-ifelse(testResults$LBonPI < 0, 0, testResults$LBonPI)


## compute the upper bound
testResults$UBonPI<-testResults$mean +(tval*(testResults$SD*sqrt(1+(1/numControlSamp))))

########################################################################################
########################################################################################
## merge the experimental samples with the control samples
########################################################################################
########################################################################################

## Combine the prediction interval results with the experimental estimates
## first add the VAR to the location information so there are multiple mutations possible at each location

summaryExp<-summaryExp[,c("Chrom","Loc","VAF","ID","locationVAR")]

colnames(summaryExp)<-c("Chrom","location","VAF","ID","locVar")

results<-merge(testResults,summaryExp,by=c("locVar","location"),all.x=TRUE)

## create an indicator for whether the VAF in the epxerimental group is higher than the upper 99% bon corrected predicton interval
results$sigVAF<-ifelse(results$VAF> results$UBonPI,1,0)


results$location<-as.numeric(results$location)




results$site<-vapply(1:nrow(results),function(x){
  
  location<-results$location[x]
  
  sites$site[which(sites$start <= location & sites$end >= location )]
  
},FUN.VALUE=character(1))

## remove locVar column from dataset
results<-results[,colnames(results)!= "locVar"]

## add what alpha was used  to the dataset
results$alpha<-alphaLevel

## add how many tests were corrected using bonforroni's correction
results$numTest<-num



write.table(results,paste(outFileLoc,"resultsVAF.csv"),sep=",",row.names = FALSE)




end_time <- Sys.time()

end_time - start_time
