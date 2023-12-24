
#(1) Relevant functions for filtering process:######

#Coordinates used:
itm<-"+init=epsg:2039" #ITM coordinations
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#(A) Download data from SQL:

# This function will require a VPN connections to server and
# will give a dataframe that contain all the locations for the tags listed in the ListOfStart dataframe:
#Required impute Parameters:
#1) Start_Time_Str: Starting date in the following format 'YYYY-MM-DD HH:MM:SS' (e.g.,'2020-03-31 06:00:00'), in UTC
#2) End_Time_Str: End date in the following format 'YYYY-MM-DD HH:MM:SS'(e.g.,'2020-11-18 06:00:00'), in UTC
#3) ListOfStart: dataframe with the following variables:
# TAG as int with two or three digits (e.g., 46,47,..,147)
# 4) TagPrefix: The initial value for the atlas tag numbers, the defult is Harod (972006000000)
# 5) quaryVar: The variable required  from the atlas server in a comma seperated chr format, the defult is 'TAG,TIME,X,Y,VARX,VARY,COVXY'.

# Output parameters:
# Dataframe with 17 variables:
# 1) X,Y: Locations as num
# 2) TIME: ATLAS detection time as number format
# 3) TAG: tag name as character with three digits (e.g., "047",..."194")
# 4) date: Date of detection (as Date format)
# 5) VARX,VARY,COVXY
# 6) LON,LAT: Locations as ITM
# 7) dateTime: Real time formate as POSIXct ("2020-07-28 06:00:05")
# 8)  dT: Detection time
# 9) distance
# 10) spd
# 11) stdVarXY
# 12) angl
# 13) traceNorm
# Notes:
# 1) The function will search in the server all the tags specified in ListOfStart
# 2) This function work for the Harod system (in the connection to server)
# This will also give the important atributes for later (speed, angle ext..)

LocationRawData <- function(Start_Time_Str,End_Time_Str,ListOfStart,TagPrefix=972006000000,quaryVar='TAG,TIME,X,Y,VARX,VARY,COVXY') {
  
  # Functions to interact with databases
  # --- Database Connection
  
  #Connection data need to be filled to download from SQL:
  
  dbc <- dbConnect(RMySQL::MySQL(),
                   user = ,            # username 
                   password = ,# password
                   host = ,            # host ip address
                   port=,                        # port Number
                   dbname=)                   # name of data base
  
  
  # --- Set start & end time and convert to ATLAS time
  
  Start_Time_Str_Temp <- as.character.Date(Start_Time_Str)
  ATLAS_Start_Time<-as.numeric(as.POSIXct(Start_Time_Str_Temp,
                                          "%Y-%m-%d %H:%M:%S", tz="UTC"))*1000
  
  End_Time_Str_Temp <- as.character.Date(End_Time_Str)
  ATLAS_End_Time<-as.numeric(as.POSIXct(End_Time_Str_Temp,
                                        "%Y-%m-%d %H:%M:%S", tz="UTC"))*1000
  
  #Make the tag list as full tag name:
  
  FullTag <- as.character(ListOfStart$TAG+TagPrefix)  #Create a list with only the tags
  
  AllTagsLoc <- list() #make an empty list for localization
  
  for (i in 1:length(FullTag)) {
    query = paste('select', quaryVar, 'from LOCALIZATIONS WHERE TAG=',FullTag[i],
                  'AND TIME >=', ATLAS_Start_Time, 'AND TIME <=', ATLAS_End_Time)
    All_Data <- dbGetQuery(dbc,query)
    AllTagsLoc[[i]] <- All_Data
    print(paste(i,"Tags out of", length(FullTag), "are done"))
  }
  
  print(paste("Binding all tags, might take some time"))
  
  AllthetagsLoc <- do.call(rbind.data.frame, AllTagsLoc)
  
  dbDisconnect(dbc)
  
  # Make the locations as ITM
  AllthetagsLoc <-convertSpatial.ITM2WGS84(AllthetagsLoc, xyColNames=c("X","Y"))
  AllthetagsLoc <- as.data.frame(AllthetagsLoc)
  
  # Substruct the tag names
  AllthetagsLoc$TAG <- substr(AllthetagsLoc$TAG, 9, 13)
  
  AllthetagsLoc<-AllthetagsLoc[order(AllthetagsLoc$TAG,AllthetagsLoc$TIME),] #make sure data is sorted chronologically (per tag)
  AddTrubite <- list()
  print(paste("Adding atributes"))
  FullTag2 <- unique(AllthetagsLoc$TAG)
  
  AllthetagsLoc$DateTime<-as.POSIXct((AllthetagsLoc$TIME)/1000, tz="UTC", origin="1970-01-01")
  AllthetagsLoc$date<-as.Date(as.POSIXct((AllthetagsLoc$TIME)/1000, tz="UTC", origin="1970-01-01"))
  
  print(paste("Done!"))
  
  return(AllthetagsLoc)
}

#(B) Filter data by days and locations:

# data - the dataset
# minimumday - Maximum days to remove
# minimumloc - minimum locations to remove

FilterBadDays <- function(data,minimumdays,minimumloc) {
  
  Ring_IDs <- unique(data$Ring_ID)
  onlygooddaysfinal <- list()
  
  for (i in Ring_IDs) {
    d <- subset(data, Ring_ID == i)
    Dat <- unique(d$date)
    Dat <- as.factor(Dat)
    onlygooddays <- list()
    if (length(Dat) <= minimumdays) { #First check if the tag have, if not the tag is removed from analysis
      print(paste("Ring_ID", i, "had less than", minimumdays,"(",length(Dat),")", "days and was REMOVED"))
    }else {
      for (ii in Dat) { #Here the function check each day separately and remove days with less than the threshold locations
        dd <- subset(d, date == ii)
        if (nrow(dd) <= minimumloc) {
          print(paste("Ring_ID", i, "on date", ii, "was REMOVED"))
        } else {
          onlygooddays[[ii]] <- dd
        }
      }
    }
    tt <- do.call(rbind.data.frame, onlygooddays)
    Dat2 <- unique(tt$date)
    Dat2 <- as.factor(Dat2)
    if (length(Dat2) <= minimumdays) {
      print(paste("After filtering Ring_ID", i, "had less than", minimumdays,"(",length(Dat2),")", "days and was REMOVED"))
    }else {
      onlygooddaysfinal[[i]] <- tt
      print(paste("For Ring_ID", i,":",length(unique(tt$date)), "out of",length(unique(d$date)), "days were calculated"))
    }
  }
  
  FilterBadDays <- do.call(rbind.data.frame, onlygooddaysfinal)
  return(FilterBadDays)
  
  
}

# (C) The function is a wrapper to a cpp code (TrackConfidanceVec) saved in file "TrackConf.cpp"
# that calculates the confidence of any point in a track without discarding it
# it runs a loop on all data points per tag and gives higher confidence mark (in the range [0,1,2] )
# to points which have large NBS or are close to confidant points
# Function created by Eitam Arnon.
#http://adv-r.had.co.nz/Rcpp.html
#install.packages('Rcpp')
require('Rcpp')
# knitr:::input_dir()

TrackConfidanceLevelcpp <- function(Data,conectedVel=20,conectedDist=NA,stdlim=80,minNBSforConf2=7,
                                    minNBSforConf1=4,Nconf1forConf2=5)
{
  
  if (!all(c('TAG','X','Y','NBS','TIME') %in% names(Data)))
  {stop("TrackConfidanceLevel needs the variables TAG,X,Y,NBS, and TIME (in miliseconds) to run ")}
  
  if (!('val2' %in% names(Data)))
  {
    if (!all(c('VARX','VARY','COVXY') %in% names(Data)))
    {stop("TrackConfidanceLevel needs either val2 or VARX,VARY,COVXY to run")}
    print("TrackConfidanceLevel calculates stdVarXY as greater eigenvalue")
    Data <- Data %>% mutate(val2= sqrt(((VARX+VARY)+sqrt(VARX^2+VARY^2-2*VARX*VARY+4*COVXY^2))/2)) # greater eigenvalue
  }
  
  Data <- Data %>%  arrange(TIME)
  
  if(!('Conf' %in% names(Data)))
  {Data$Conf=-1}
  Data$aBS <- Data$NBS
  
  listoffilteredData <- list()
  tags <- unique(Data$TAG)
  for (tagInd in 1:length(unique(Data$TAG)))
  {
    tagData <- Data %>% dplyr::filter(TAG==tags[tagInd])
    timediffs=unique(tagData$TIME[1:min(nrow(tagData),1e4)]-lag(tagData$TIME[1:min(nrow(tagData),1e4)]))/1000
    numinalTimeDiff <- min(round(timediffs)[which(round(timediffs)>0.5)])
    if (is.na(conectedDist))
    {conectedDist <- numinalTimeDiff*conectedVel}
    tagData$Conf <-TrackConfidanceVec(as.matrix(tagData %>% dplyr::select(X,Y,aBS,stdVarXY)),
                                      minNBSforConf2,minNBSforConf1,Nconf1forConf2,conectedDist,stdlim)
    listoffilteredData[[tagInd]] <- tagData
  }
  Data <- do.call(rbind.data.frame,listoffilteredData)
  
  Data <- Data %>% dplyr::select(-aBS)
  return(Data)
}

#(D)# A function to assign day number
# Every day begins on DayStartTime (a HH:MM:SS character, in UTC, just like the ATLAS times)
# The input variable data is a data.frame containing a POSIXct time columns (tz="UTC",origin="1970-01-01")
# The data frame time columnmust be named TimeColName
# The day are counted from either as Julian days (Julian=TRUE) or from the first day in the array (default, Julian=FALSE)
# Currently the days are counted from the earliest time in the set without taking into account different Tags
# Returns the same data.frame with an additional column ("DAY")

AssignDayNumber <- function(data,DayStartTime="00:00:00",TimeColName = "TIME",Julian=FALSE){
  datatimes <- as.data.frame(data)
  datatimes <-datatimes[,TimeColName]/1000
  timeshift <- as.numeric(as.POSIXct(paste("1970-01-01",DayStartTime),tz="UTC",origin="1970-01-01"))
  shifteddays <- as.POSIXct(as.numeric(datatimes)-timeshift, tz="UTC", origin="1970-01-01")
  if (Julian)
  {
    data$DAY <- yday(shifteddays)
  }
  else
  {
    data$DAY <- yday(shifteddays)-min(yday(shifteddays))+1
  }
  return(data)
}

# (E) Smoothing ground and fly
#This function take the movement data, separate it to flying and ground locations and then filter each segment seperatly.
#Because the ATALAS system observe better locations when the animal is flying (due to being more visible to more antennas),
#Filtering process need to be different between air and ground locations
#Furthermore, when the bird is on the ground the distance and speed is reduce to minimal by default of the bird biology

#DataForSmoothing: The movement data for further processing
#ind_rts: dataframe used to separate the data to ground and flying segments.
#This dataframe created give the option to found fixed location on the ground to create a radius around
#The locations that exit the radius will be segmented as flying locations, while those inside the radius is
#defined as ground locations.
#This is the data created (can be changed depending on species)

ind_rts<-data.frame("smp_rte"=c(8000), # sample rate, every 8 second
                    "obs_min"=c(5), # minimum locations of within range localizations
                    "p_lim"=c(5), # Tolerance of n localizations outside the buffer
                    "adp_rng"=rep(35)) # adaptive Fixed point buffer radius

#MinFly: minimum points tolerant for flying locations
#MinGround: Minimum points tolerant for ground locations
#MedianSmooth: the K number for smoothing, defined to 3 locations
#Method: smoothing method for ground locations, can be the mean or the median

SmoothGroundFly <- function(DataForSmoothing,ind_rts,MinFly = 2, MinGround = 3, MedianSmooth = 3, Method = "Mean") {
  
  #Arrange the data for further processing:
  DataForSmoothing$SecTime <- DataForSmoothing$TIME/1000 
  DataForSmoothing$Transmitter <- DataForSmoothing$TAG
  DataForSmoothing$TAG <- DataForSmoothing$Ring_ID
  DataForSmoothing <- AddspeedAngleVarxy(as.data.frame(DataForSmoothing)) #Add the speed/angle/varX
  TagNames <- unique(DataForSmoothing$TAG)
  listagging <- list()
  
  for (a in TagNames) {
    Tag <- subset(DataForSmoothing,TAG ==  a)
    Name1<- AssignDayNumber(Tag) #Add days for each tag for further froccessing
    
    days <- unique(Name1$date)
    finisheddays <- list()
    for (b in 1:length(days)) {
      tagday <- subset(Name1,date ==  days[b])
      
      print(paste("Start to work on",a, "on day", days[b]))
      
      
      if (max(dist(as.data.frame(c(tagday[1],tagday[2]), ncol=2, nrow= length(tagday$X), byrow=TRUE))) < 200) {
        #If the data max distance is too small, it is not possible to separate it to segments,
        #Thus if the distances between locations are too small, in advance it defined as 1 ground location (ADP)
        
        print(paste("tag", a, "day", days[b], "have only one ADP"))
        tagday$ADP.ID <- 1
        tagday$ADP.Y <- max(as.data.frame(c(tagday[1],tagday[2]), ncol=2, nrow= length(tagday$X), byrow=TRUE)$Y)
        tagday$ADP.X <- max(as.data.frame(c(tagday[2],tagday[1]), ncol=2, nrow= length(tagday$Y), byrow=TRUE)$X)
        tagday$ADP.start <- head(tagday$TIME,1)
        tagday$ADP.end <- tail(tagday$TIME,1)
        tagday$ADP.duration <- (tagday$ADP.end - tagday$ADP.start)/1000
        tagday$ADP.qlt <- 1
        tagday$seg.ID <- NA
      } else { #If the distance is higher than threshold seperation begin
        
        tagday <- stop_move_analysis(tagday,stopparams = ind_rts) #Seperate to ground and flying locations
        
        print(paste("tag", a, "day", days[b]))
        ADPName <- unique(tagday$ADP.ID)
      }
      ADPName <- unique(tagday$ADP.ID)
      
      tagdayLand <- subset(tagday,ADP.ID  > 0)
      tagdayFly <- subset(tagday,seg.ID  > 0)
      
      ADPlist <- list()
      
      for (i in ADPName) { #Separate all the ground location to filter them
        ADP2 <- subset(tagdayLand, ADP.ID ==i)
        if (nrow(ADP2) <=MinGround){ #If the ground location are below the threshold, remove it
          next
        }
        if (nrow(ADP2) <= 0) {
          next
        }
        
        #atl_filter_covariates, distance_filter, and atl_median_smooth  were used with 'atlastools':
        #atlastools: Pre-processing Tools for High-throughput Animal Tracking Data
        #Pratik Rajan Gupte, 2021
        ADP2T2 <- atl_filter_covariates(data = ADP2 ,
                                        filters = c("spd < 5")) #Filter speed below 5 for ground locations
        
        
        ADP2T2 <- distance_filter(ADP2T2,distThreshold=5*8) #Filter distance below 40 meter for ground locations
        
        if (nrow(ADP2T2) <=MinGround) { #After filtering, remove location that are below threshold
          next
        }
        if (Method == "Median") { #Smooth the ground locations depending on median or mean
          ADP2T2 <- atl_median_smooth(data = ADP2T2,
                                      x = "X", y = "Y",
                                      time = "SecTime",
                                      moving_window = MedianSmooth)
          ADP2T2 <- movMean(data = ADP2T2)
          ADPlist[[i]] <- ADP2T2
        } else if (Method == "Mean") {
          ADP2T2 <- movMean(data = ADP2T2,Weight = c(0.1,0.2,0.4,0.2,0.1))
          ADPlist[[i]] <- ADP2T2
          
        }
        
      }
      
      tagday2 <- do.call(rbind.data.frame, ADPlist)
      
      FlySeg <- unique(tagday$seg.ID)
      Flylistt <- list()
      
      for (ii in FlySeg) { #Work on flying segments after finishing ground segments
        
        FlySegment <- subset(tagdayFly, seg.ID == ii)
        if (nrow(FlySegment) <=MinFly){
          next
        }
        SegFly <- movMean(data = FlySegment) #Only smooth flying segments, no filtering needed
        Flylistt[[ii]] <- SegFly
      }
      tagday3 <- do.call(rbind.data.frame, Flylistt)
      tagday4 <- bind_rows(tagday3,tagday2)
      
      tagday4<-tagday4[order(tagday4$TAG,tagday4$TIME),] #make sure data is sorted chronologically (per tag)
      
      finisheddays[[b]] <- tagday4
    }
    finish1 <- do.call(rbind.data.frame, finisheddays)
    listagging[[a]] <- finish1
  }
  finish2 <- do.call(rbind.data.frame, listagging)
  finish2 <- AddspeedAngleVarxy(as.data.frame(finish2)) #Add atribute after finishing filtering
  
  finish2$TAG <- finish2$Transmitter
  
  return(finish2)
}

# (F) Define stop locations for 'SmoothGroundFly' function:

defineStopPoints<-function(data,
                           ind_rts,
                           time_gap=3600*1000){
  
  
  #Modified from a script of Ingo Schiffner, Emmanuel Lourie & Sivan Margalit
  # A wrapper for AdpFixedPoint & MarkStopLoc functions from toolsForAtlas
  # Main tasks:
  # Loop over several tags with several days each, applying the following functions:
  # AdpFixedPoint - returns a list of track segmentation with stops properties including start/end time, duration, and median location (X/Y)
  # Removing dummy points created by AdpFixedPoint (stops with only a single point)
  # MarkStopLoc - Marking localizations as stop or move segment returning a final data.frame where each point is characterized as either stop (ADP) or movement (seg)
  # Input variables:
  # data - data.frame containing a trajectory with columns $TAG,$TIME,$dT,$DAY,$X,$Y
  # ind_rts - a data.frame with ADP parameters: $smp_rte,$obs_min,$p_lim,$adp_rng
  # output - a data.frame with additional columns presenting the stops-flying segments, and their properties
  
  owl_adp<-NULL
  # 
  #   data<-data%>%
  #     arrange(TAG, TIME)
  
  TagList <- unique(data$TAG)
  
  for (tg in TagList){
    # what is this tag;s sampling frequency?
    freq<-round(min(data$dT[which(data$TAG==tg)],na.rm=TRUE),0)*1000
    if (!any(freq %in% ind_rts$smp_rte))
    { #simpleError(sprintf("your data contains frequency =  %g , which is not included in your ADP parameters ",freq))
      stop(sprintf("your data contains frequency =  %g , which is not included in your ADP parameters ",freq))}
    # Get AdpFixedPoint criterias from the table according to sampling rate:
    smp_rte <- ind_rts$smp_rte[which(ind_rts$smp_rte==freq)] # sample rate
    adp_rng <- ind_rts$adp_rng[which(ind_rts$smp_rte==freq)] # adaptive Fixed point buffer radius
    obs_min <- ind_rts$obs_min[which(ind_rts$smp_rte==freq)] # minimum of within range localizations
    p_lim <- ind_rts$p_lim[which(ind_rts$smp_rte==freq)] #Tolerance of n localizations outside the buffer
    
    
    # created a list of DAYs tracked per tag for looping:
    nDAYs<-unique(data$DAY[which(data$TAG==tg)])
    
    # Start loop per tag-DAY. Notice that data MUST be sorted chrnologically first.
    for (nn in nDAYs){
      DAY_dat<-as.data.frame(data%>%
                               filter((TAG==tg) & (DAY==nn))%>%
                               arrange(TIME))
      
      
      if(nrow(DAY_dat)>20) # only run if DAY data (point for TAG per DAY) is over 100 locations:
      {
        #function to perform track segmentation (from toolsForAtlas)
        AFPList <- AdpFixedPoint (time.vec = DAY_dat$TIME,
                                  x=DAY_dat$X,
                                  y=DAY_dat$Y,
                                  adp_rng=adp_rng,
                                  smp_rte=smp_rte,
                                  obs_min=obs_min,
                                  p_lim=p_lim,
                                  time_gap=time_gap)
        
        #remove dummy points created by AdpFixedPoint
        kdx <- which(AFPList$duration!="1")
        AFPList <- AFPList[kdx,]
        
        
        if (nrow(AFPList)>1){
          AFPList$DAY_number<-nn
          AFPList$TAG<-tg
          owl_adp<-rbind(owl_adp, AFPList)
        }
      }
    }
  }
  # Mark localizations as perch or move segment (from toolsForAtlas)
  filtered_with_ADP<- MarkStopLoc(data, owl_adp)
  return(filtered_with_ADP)
}

stop_move_analysis <- function(data,stopparams="default",DayColName=c("DAY"))
{
  # A wrapper function which defines stops/flights segments and filters the data based on these definations
  # The wrapper activates additional functions including AdpFixedPoint & MarkStopLoc (from toolsForAtlas) which
  # Defines stop points and clears the data based on the segmentation
  # Input data -
  # data - a data.frame containing a trajectory with columns $TAG,$TIME,$dT,$DAY,$X,$Y
  # the "day" variable name can be defined according to the DayColName parameter
  # stopparams - a data.frame defining the stop identification parameters
  # Output - a data.frame with additional columns presenting the stops-flying segments:
  # ADP : a set of variables of each stop:
  # "ID" a set of ID, for each stop starting at the beginning of the data.set
  #  "X","Y", the median for each stop
  #  "start", "end", "duration" times for each stop ( in "TIME" units, usually miliseconds)
  #  "qlt"  some unknown quality parameter
  # seg.ID : a set of ID, for each movement segment starting at the beginning of the data.set
  # indices to match tag frequency:
  # if (stopparams=="default")
  #   ind_rts<-data.frame("smp_rte"=c(8000,4000), # sample rate
  #                       "obs_min"=c(4,8), # minimum of within range localizations
  #                       "p_lim"=c(5,10), #Tolerance of n localizations outside the buffer
  #                       "adp_rng"=rep(35,2)) # adaptive Fixed point buffer radius
  # 
  # else
  ind_rts <- stopparams
  
  
  #Performing track segmentation - adding properties to each of the points
  names(data)[which(colnames(data)==DayColName)] <- "DAY"
  stop_data<-defineStopPoints(data,ind_rts)
  
  #clear the data using RemoveMargins: Abandons previous data for the first exit from the previous
  #day's accommodation point, and after entering the next day's accommodation point
  #  stop_data<-RemoveMargins(stop_data)
  stop_data<-as.data.frame(stop_data)
  
  stop_data<-addLocAttribute(stop_data, locAttributs=c("distanceSpeed", "locQuality")) # function to add attributes for each pair of conseucutive points.
  names(stop_data)[which(colnames(stop_data)=="DAY")] <- DayColName
  return(stop_data)
}

# (G) This function smooths a track by setting any point as weighted average of its neighboring points, according to the Weight variable!
# it is advised to operate the function on separated time-bursts of the data (as wrapped in the function "AvgSmooth"
# The function requires the following input variables:
# dataBurst: a data.frame containing the variables: "X","Y" (coordinates in UTM/ITM format).
# Weight is the weights used for the weighted average , it can have any number of values (recommended to be odd number of values in a symmetric pattern)
# The output is the same data.frame with smoothed coordinates. In the case that replace=False is specified, only the "x", "Y" coloumns will be returned!

movMean <- function(dataBurst,Weight = c(0.25,0.5,0.25),replace=T,rmGeoCoords=T)
{
  Weight <- Weight/sum(Weight)
  X <- as.numeric(stats::filter(dataBurst$X, Weight,method = "convolution",sides = 2))
  Y <- as.numeric(stats::filter(dataBurst$Y, Weight,method = "convolution",sides = 2))
  X[which(is.na(X))] <- dataBurst$X[which(is.na(X))]
  Y[which(is.na(Y))] <- dataBurst$Y[which(is.na(Y))]
  
  if(replace){
    dataBurst$X <- X
    dataBurst$Y <- Y
    if(rmGeoCoords)
    {dataBurst <- subset(dataBurst, select = -c(LAT,LON) )
    print("removing LON and LAT coordinates")}
    return(dataBurst)}
  else
    return(data.frame(X=X,Y=Y))
}

#(H) Subset the data by time
#Data: the dataframe used for subsetting
#SubsetTime: how much to subset

SubsetByTime <- function(Data, SubsetTime) {
  library(dplyr)
  library(data.table)
  library(lubridate)
  Rings <- unique(Data$Ring_ID)
  Listone <- list()
  for (i in Rings) {
    LoopTag <- subset(Data, Ring_ID == i)
    daysfilter <- as.factor(unique(LoopTag$date))
    listtwo <- list()
    for (ii in daysfilter) {
      LoopDay <- subset(LoopTag, date == ii)
      LoopDay <- setDT(LoopDay)[order(LoopDay)]
      
      output <- LoopDay[, .(DateTime = dateTime[1],date = date[1], Ring_ID = Ring_ID[1], X = X[1], Y = Y[1], TIME = TIME[1]) ,
                        by = .(Group = floor_date(dateTime, SubsetTime))]
      listtwo[[ii]] <- output
    }
    dataframetwo <- do.call(rbind.data.frame, listtwo)
    Listone[[i]] <- dataframetwo
    print(paste("Tag", i, "finish"))
  }
  NewLocations <- do.call(rbind.data.frame, Listone)
  NewLocations<-NewLocations[order(NewLocations$Ring_ID,NewLocations$TIME),] #make sure data is sorted chronologically (per tag)
  
  return(NewLocations)
}

#(I) Function to combine the home-range data with a raster layer, and observe which habitats are used in the home-range
# This function give a metadata combining the each individual, and the percentage of habitat types in the home-range
# Data = The movement dataframe
# RasName = The raster layer to create a metadata
# Subset = if the data need subsetting write "TRUE", in our case its false

MetadataHabitatHR <- function(Data,RasName, Subset = FALSE) {
  
  if (Subset == TRUE) {
    Data <- SubsetByTime(Data, SubsetTime = "60 mins")
  }
  #Work on the raster layer given, and combine it with the movement data:
  
  Harod_Habitats2=RasName 
  
  xyt<-subset(Data, select = c(X,Y))
  id<-subset(Data, select = Flag)
  locs1<-id
  coordinates(locs1)<-xyt
  
  #Creating home-range from the movement data:
  
  UD3 <- kernelUD(locs1[,1], h = "href", grid = 500,same4all = FALSE, hlim = c(0.1, 2), kern = "bivnorm", extent = 1,boundary = NULL)
  
  homerange1 <- getverticeshr(UD3, percent = 95)
  homerange1_Core <- getverticeshr(UD3, percent = 50)
  
  homerange2 <- homerange1
  homerange1_Core2 <- homerange1_Core
  
  proj4string(homerange2)<-CRS(itm)
  proj4string(homerange1_Core2)<-CRS(itm)
  
  Raster2 <- projectRaster(Harod_Habitats2, crs=itm) 
  
  Combined2 <- raster::intersect(Raster2,homerange2) #Combining home-range with the raster layer
  
  r3 <- mask(Combined2, homerange2) #Take only the raster from the home-range area
  #plot(r3)
  
  FlagID <- unique(homerange2@data$id)
  BuiltAreas <- list()
  BuiltAreasFull <- list()
  print("Making Metadata")
  # i <- 4
  # Loop for every tag chosen:
  for (i in FlagID) {
    #Taking the data from the home-range:
    FullData <- subset(Data,Flag == i)
    DaysCalculated <- length(unique(FullData$date))
    
    animalnum95 <- homerange2[i,]
    animalnum50 <- homerange1_Core2[i,]
    
    #Adding area for later:
    size95 <- animalnum95$area
    size50 <- animalnum50$area
    #Combining the raster with the data:
    comination1 <- raster::intersect(Raster2,animalnum95)
    comination2 <- raster::intersect(Raster2,animalnum50)
    #Removing unwanted raster locations:
    r4 <- mask(comination1, animalnum95)
    r5 <- mask(comination2, animalnum50)
    # plot(r5)
    
    #Making a new dataframe from 95 and 50 home-range:
    tabtab <- table(round(r4@data@values,0))
    tabtab2 <- table(round(r5@data@values,0))
    
    tab2 <- data.frame(Color = names(tabtab), Count95 = as.integer(tabtab))
    
    tab2 <- tab2[1:9,]
    tab2$Color <- 3:11
    tab2$Count95 <- NA
    
    
    #3 = Main roads
    #4 = dirt road
    #5 = Fields 1
    #6 = Fields 2
    #7 = Fields 3
    #8 = mataim
    #9 = Natural
    #10 = WaterPonds
    #11 = Urban
    
    tab4 <- data.frame(Color = as.integer(names(tabtab)), Count95 = as.integer(tabtab))
    #Adding the habitats. Those lines need to be changed depending what your raster value is:
    tab2$Habitat[tab2$Color == "3"] <- "Main roads"
    tab2$Habitat[tab2$Color == "4"] <- "dirt road"
    tab2$Habitat[tab2$Color == "5"] <- "Fields 1"
    tab2$Habitat[tab2$Color == "6"] <- "Fields 2"
    tab2$Habitat[tab2$Color == "7"] <- "Fields 3"
    tab2$Habitat[tab2$Color == "8"] <- "mataim"
    tab2$Habitat[tab2$Color == "9"] <- "Natural"
    tab2$Habitat[tab2$Color == "10"] <- "WaterPonds"
    tab2$Habitat[tab2$Color == "11"] <- "Urban"
    tab2$Flag <- i
    
    tab2$Count95[match(tab4$Color,tab2$Color)] <- tab4$Count95
    
    tab2$Percent95 <- round(tab2$Count95/sum(tab2$Count95, na.rm=TRUE)*100,digits = 2)
    
    tab3 <- tab2[ -c(1) ]
    
    tab4 <- dcast(tab3, Flag~Habitat, value.var='Percent95')
    tab4$HR <- 95
    tab4$Size <- round(size95,digits = 2)
    
    #habitat 50:
    
    tab2b <- data.frame(Color = names(tabtab2), Count50 = as.integer(tabtab2))
    
    tab2b <- tab2b[1:9,]
    tab2b$Color <- 3:11
    tab2b$Count50 <- NA
    
    tab4b <- data.frame(Color = as.integer(names(tabtab2)), Count50 = as.integer(tabtab2))
    
    
    tab2b$Habitat[tab2b$Color == "3"] <- "Main roads"
    tab2b$Habitat[tab2b$Color == "4"] <- "dirt road"
    tab2b$Habitat[tab2b$Color == "5"] <- "Fields 1"
    tab2b$Habitat[tab2b$Color == "6"] <- "Fields 2"
    tab2b$Habitat[tab2b$Color == "7"] <- "Fields 3"
    tab2b$Habitat[tab2b$Color == "8"] <- "mataim"
    tab2b$Habitat[tab2b$Color == "9"] <- "Natural"
    tab2b$Habitat[tab2b$Color == "10"] <- "WaterPonds"
    tab2b$Habitat[tab2b$Color == "11"] <- "Urban"
    tab2b$Flag <- i
    
    tab2b$Count50[match(tab4b$Color,tab2b$Color)] <- tab4b$Count50
    sum(tab2b$Count50, na.rm=TRUE)
    
    tab2b$Percent50 <- round(tab2b$Count50/sum(tab2b$Count50, na.rm=TRUE)*100,digits = 2)
    
    tab3b <- tab2b[ -c(1) ]
    
    tab4b <- dcast(tab3b, Flag~Habitat, value.var='Percent50')
    tab4b$HR <- 50
    tab4b$Size <- round(size50, digits = 2)
    
    #Combine tab3 with tab3b, and tab4 + tab4b:
    
    tab3_new <- merge(tab3,tab3b,by=c("Habitat","Flag"))
    tab4_new <- rbind(tab4,tab4b)
    tab4_new$Days <- DaysCalculated
    
    BuiltAreas[[i]] <- tab4_new
    BuiltAreasFull[[i]] <- tab3_new
    print(paste("Finished flag number", i))
  }
  
  NumberOfArea <- do.call(rbind.data.frame, BuiltAreas)
  NumberOfAllArea <- do.call(rbind.data.frame, BuiltAreasFull)
  
  # Arranging the data:
  NumberOfArea$Flag <- as.factor(NumberOfArea$Flag)
  NumberOfAllArea$Flag <- as.factor(NumberOfAllArea$Flag)
  NumberOfAllArea$Habitat <- factor(NumberOfAllArea$Habitat, levels = c("Main roads", "dirt road", "Fields 1","Fields 2","Fields 3","mataim","Natural","WaterPonds","Urban"))
  
  return(NumberOfArea) #Single row to 95 and 50
  
}

#(J) Function for day metadata. Create metadata based on ground and flying behaviors
#Specifically, give the max diameter, number of ground locations, time on the ground,time flying, and distance flying

MetadataDays <- function(Data) {
  
  FullList <-list()
  
  for (i in unique(Data$Flag)) {
    Flag1 <- subset(Data, Flag == i)
    DateList <- list()
    for (ii in 1:length(unique(Flag1$date))) { #For each date, calculated separately
      Date1 <- subset(Flag1, date == unique(Flag1$date)[ii])
      MaxDiameter <- max(dist(as.data.frame(c(Date1[1],Date1[2]), ncol=2, nrow= length(Date1$X), byrow=TRUE)))
      
      df <- data.frame(Flag  = unique(Date1$Flag),
                       date = unique(Date1$date),
                       MaxDiameter = MaxDiameter,
                       TAG = unique(Date1$TAG),
                       Ring_ID = unique(Date1$Ring_ID))
      DateList[[ii]] <- df
    }
    FlagList <- do.call(rbind.data.frame, DateList)
    FullList[[i]] <- FlagList
    print(paste("Flag number", i, "Max Dia is done"))
  }
  
  print(paste("Start working on ground behavior"))
  FullList <- do.call(rbind.data.frame, FullList)
  
  #Ground behaviors:
  Ground <- subset(Data,ADP.ID > 0 )
  
  agg1 <- Ground %>% dplyr::select(Flag, date,dateTime, ADP.ID,ADP.duration) %>%
    ddply(Flag~date,summarise,ADPs=length(unique(ADP.ID)),GroundLocations=length(unique(dateTime)),
          MeanGroundDurationMinutes = (mean(ADP.duration, na.rm=TRUE)/1000)/60)
  
  #Flight behaviors:
  rm(Ground)
  
  print(paste("Start working on flying behavior"))
  
  Flying <- subset(Data,seg.ID > 0 )
  agg2 <- Flying %>% dplyr::select(Flag, date,dateTime, distance,spd,seg.ID) %>%
    ddply(Flag~date,summarise,FlyingSegments=length(unique(seg.ID)),TotalDistanceFly=round(sum(distance, na.rm=TRUE),2), 
          MeanSpeedFly = round(mean(spd, na.rm=TRUE),2))
  
  rm(Flying)
  #General behaviors:

  print(paste("Create metadata"))
  agg3 <- merge(agg2, agg1, by= c("Flag","date"))
  FinalMeta <- merge(agg3, FullList, by= c("Flag","date"))
  
  return(FinalMeta)
  
}

#Other functions for behaviors on the ground:

#(I) Functions to create a metadata with ground locations:

#Lap1 - the data to create ADP metadata from:

DurationTimeADPByDate <- function(Lap1) {
  #Create grid of 100 meters:
  #Function that will give every day: (1) Location name (2) Time spend in this location
  
  Name2 <- subset(Lap1, ADP.ID > 0) #Use only the ADP locations:

  Name3 <- NewADPDuration(Name2) 
  
  Name3$ADP.X2 <- round(Name3$ADP.X/100,0)*100
  Name3$ADP.Y2 <- round(Name3$ADP.Y/100,0)*100
  
  Name3$LocName <- paste(Name3$ADP.X2, "_", Name3$ADP.Y2, sep = "")
  
  TotalLocationsDuration <- Name3 %>% dplyr::select(Ring_ID,LocName,date, ADP.ID,ADP.X2, ADP.Y2,ADP.duration) %>%
    ddply(Ring_ID ~ date + LocName, dplyr::summarise, ADP.X2 = unique(ADP.X2),
          ADP.Y2 = unique(ADP.Y2), LocName = unique(LocName),
          ADP.durationMinutes = sum(unique(ADP.duration))) #OK
  
  SumDuration <- TotalLocationsDuration %>% ddply(Ring_ID ~ LocName, dplyr::summarise, SumDurTotalMinutes = sum(ADP.durationMinutes), ADP.X2 = unique(ADP.X2),
                                                  ADP.Y2 = unique(ADP.Y2)
  )
  
  SumDurationMetadata2 <- SumDuration %>% dplyr::group_by(Ring_ID) %>% 
    dplyr::mutate(PercTimeSpend = round(SumDurTotalMinutes/sum(SumDurTotalMinutes)*100,1), GiniIndex = gini.wtd(x=round(SumDurTotalMinutes/sum(SumDurTotalMinutes)*100,1))
    )
  SumDurationMetadata2$VisitsFactor <- NA
  SumDurationMetadata2$VisitsFactor[SumDurationMetadata2$PercTimeSpend >= quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.75))] <- "coreArea"
  SumDurationMetadata2$VisitsFactor[SumDurationMetadata2$PercTimeSpend < quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.75)) &
                                      SumDurationMetadata2$PercTimeSpend > quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.5))] <- "HighArea"
  SumDurationMetadata2$VisitsFactor[SumDurationMetadata2$PercTimeSpend <= quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.5)) &
                                      SumDurationMetadata2$PercTimeSpend > quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.10))] <- "MedianArea"
  SumDurationMetadata2$VisitsFactor[SumDurationMetadata2$PercTimeSpend <= quantile(0:max(SumDurationMetadata2$PercTimeSpend),probs = c(0.10))] <- "LowArea"
  
  TotalLocationsDurationSum <- TotalLocationsDuration %>% dplyr::select(Ring_ID,date,LocName,ADP.durationMinutes) %>%
    ddply(Ring_ID ~ date , dplyr::summarise, NumberOfStops = length(unique(LocName)),
          TimeOnGround = sum(ADP.durationMinutes)) #Still ok
  
  #Combine with dates the core-area and check number of visits to core:
  
  TotalLocationsDuration2 <- merge(SumDurationMetadata2[ ,(!names(SumDurationMetadata2) %in% c("PercTimeSpend", "SumDurTotalMinutes"))],TotalLocationsDuration, by = c("Ring_ID", "ADP.X2","ADP.Y2","LocName"))
  TotalLocationsDuration2 <- TotalLocationsDuration2[order(TotalLocationsDuration2$Ring_ID,TotalLocationsDuration2$date),] #make sure data is sorted chronologically (per tag)
  TotalLocationsDuration2 <- merge(TotalLocationsDurationSum,TotalLocationsDuration2, by = c("Ring_ID", "date")) #Still ok
  
  TotalLocationsDuration3 <- TotalLocationsDuration2 %>% 
    dplyr::select(Ring_ID, date, NumberOfStops, TimeOnGround, LocName, GiniIndex, VisitsFactor, ADP.durationMinutes, ADP.X2, ADP.Y2) %>%
    dplyr::group_by(Ring_ID, date) %>%
    dplyr::mutate(
      NumberOfUniqueStops = unique(NumberOfStops),
      TimeOnGround = unique(TimeOnGround),
      GiniIndex = unique(GiniIndex),
      distinct_Core = sum(VisitsFactor == "coreArea"),
      distinct_HighArea = sum(VisitsFactor == "HighArea"),
      distinct_MedianArea = sum(VisitsFactor == "MedianArea"),
      distinct_LowArea = sum(VisitsFactor == "LowArea"),
      SumCoreArea = sum(ADP.durationMinutes[VisitsFactor == "coreArea"]),
      ADP.X2 = ifelse(any(VisitsFactor == "coreArea"), head(unique(ADP.X2[VisitsFactor == "coreArea"]), 1), NA),
      ADP.Y2 = ifelse(any(VisitsFactor == "coreArea"), head(unique(ADP.Y2[VisitsFactor == "coreArea"]), 1), NA)
    )
  
  TotalLocationsDuration3$month <- month(TotalLocationsDuration3$date)
  TotalLocationsDuration3$Nesting <- NA
  if (nrow(subset(TotalLocationsDuration3, month >=4 & TotalLocationsDuration3$month <= 8)) > 0) {
    TotalLocationsDuration3[TotalLocationsDuration3$month >=4 & TotalLocationsDuration3$month <= 8,]$Nesting <- "Nesting"
  }
  
  if (nrow(subset(TotalLocationsDuration3, month < 4 | month > 8)) > 0) {
    TotalLocationsDuration3[TotalLocationsDuration3$month < 4 | TotalLocationsDuration3$month > 8,]$Nesting <- "NonNesting"
  }
  
  return(TotalLocationsDuration3)
  
} #Function for total duration in each location

#Add time in each ADP:

NewADPDuration <- function(Name2) {
  ADPDur <- Name2 %>% dplyr::select(Flag,date,dateTime,ADP.ID, ADP.start,ADP.end) %>%
    ddply(Flag ~ ADP.ID + date, dplyr::summarise,
          ADP.start = unique(head((dateTime),1)),
          ADP.end = unique(tail((dateTime),1)),
          ADP.duration = round(as.numeric(difftime(tail((dateTime),1), head((dateTime),1), units="mins")),2)
    ) 
  ADPDur <- ADPDur[!duplicated(ADPDur$ADP.start), ]
  Name3 <- merge(ADPDur,Name2[ ,(!names(Name2) %in% c("ADP.duration", "ADP.start","ADP.end"))], by = c("Flag", "date", "ADP.ID"))
  Name3<-Name3[order(Name3$Flag,Name3$date),] #make sure data is sorted chronologically (per tag)
  
  return(Name3)
}

#(2) Filter tags according to given parameters:#########

#LapwingCaptureData: Table with relevant TAGs (a column name TAG is obligate) numbers, together with other relevant informations (such as sex, weight or more)

for (i in IDList) {
  
  #Download the data with SQL (dates as 'YYYY-MM-dd hh:mm:ss':
  Tag1Raw <- LocationRawData('2023-09-11 05:00:00','2023-11-27 05:00:00',subset(LapwingCaptureData, Ring_ID == i),quaryVar='NBS,TAG,TIME,X,Y,VARX,VARY,COVXY') # Download the data I want

  #Add localization, speed, distance and accuracy to the movement data:
  
  DF <-addLocAttribute(Tag1Raw, locAttributs=c( "locQuality")) #This function is part of the 'toolsForAtlas' package
  
  #Remove tags with less than 60 days from analysis, as well as days with less than 1500 locations
  DF2 <-FilterBadDays(DF,60,1500) 

  if (nrow(DF2) > 0) { #If the individual is above the threshold of locations and days the function move to the next step
    
    #Remove locations with wrong localizations (too far from antennas, or not visible enough)
    EitamFilt <- TrackConfidanceLevelcpp(DF2) 
    Listoflapwings <- subset(EitamFilt, Conf == 2) #Take only the data that with confidence of 2
  
    Listoflapwings<- AssignDayNumber(Listoflapwings) #Give day number to the tags for the next filtering
    #The next function help filter locations according to the bird behavior
    
    Listoflapwings <- SmoothGroundFly(Listoflapwings, ind_rts,MinFly = 2, MinGround = 3, MedianSmooth = 3, Method = "Median")
    #After filtering the data, re-assemble the movement data. Use a 100-meter radius to seperate long distance flying:
    
    ind_rts2<-data.frame("smp_rte"=c(8000), # sample rate
                         "obs_min"=c(10), # minimum locations of within range localizations
                         "p_lim"=c(3), # Tolerance of n localizations outside the buffer
                         "adp_rng"=rep(100)) # adaptive Fixed point buffer radius
    
    Listoflapwings <- Listoflapwings[, !(colnames(Listoflapwings) %in% c("ADP.ID", "ADP.X","ADP.Y", "ADP.start","ADP.end", "ADP.duration","ADP.qlt","seg.ID"))]
    Listoflapwings <- AddspeedAngleVarxy(Listoflapwings)
    Listoflapwings <- stop_move_analysis(Listoflapwings,stopparams = ind_rts2)

    #Separate the data to ground segments to create the Home-range calculations:
    
    FilterGroundOnly <- subset(Listoflapwings,ADP.ID > 0 )
    HR60_EveryHour <- SubsetByTime(FilterGroundOnly, SubsetTime = "60 mins") #Subset the ground locations to ~1 hour
    
    ### Change to seasons:
    print(paste("Ring_ID", i, "Metadata habitats"))
    
    HR60_EveryHour$MonthNum <- month(HR60_EveryHour$date) #Add the month to seperate by season
    
    Name1HR60_EveryHourNestingSeason <- subset(HR60_EveryHour,MonthNum >=4 & MonthNum <= 8)
    Name1HR60_EveryHourNonNestingSeason <- subset(HR60_EveryHour,MonthNum < 4 | MonthNum > 8)
    
    #Create metadata with habitats. This is done for each season separately.
    if (nrow(Name1HR60_EveryHourNestingSeason) > 20) {
      MetadataHR_Nesting <- MetadataHabitatHR(Name1HR60_EveryHourNestingSeason,Full_Raster6)
      MetadataHR_Nesting$Season <- "Nesting"
    }
    if (nrow(Name1HR60_EveryHourNonNestingSeason) > 20) {
      MetadataHR_NonNesting <- MetadataHabitatHR(Name1HR60_EveryHourNonNestingSeason,Full_Raster6)
      MetadataHR_NonNesting$Season <- "Non-Nesting"
    }
    if (nrow(Name1HR60_EveryHourNonNestingSeason) > 20 & nrow(Name1HR60_EveryHourNestingSeason) > 20) {
      Full_MetadataHR <- rbind(MetadataHR_Nesting,MetadataHR_NonNesting)
    } 
    if (nrow(Name1HR60_EveryHourNonNestingSeason) > 20 & nrow(Name1HR60_EveryHourNestingSeason) <= 20) {
      Full_MetadataHR <- MetadataHR_NonNesting
    }
    if (nrow(Name1HR60_EveryHourNonNestingSeason) <= 20 & nrow(Name1HR60_EveryHourNestingSeason) > 20) {
      Full_MetadataHR <- MetadataHR_Nesting
    }
    #All seasons:
    #Testing it for all year:
    
    Full_MetadataHRAll <- MetadataHabitatHR(HR60_EveryHour,Full_Raster6)
    Full_MetadataHRAll$Season <- "All"
    Full_MetadataHR <- rbind(Full_MetadataHRAll,Full_MetadataHR)
    
    Full_MetadataHR[is.na(Full_MetadataHR)] <- 0
    Full_MetadataHR$AllFields <- Full_MetadataHR$`Fields 1` + Full_MetadataHR$`Fields 2` + Full_MetadataHR$`Fields 3` + Full_MetadataHR$Natural
    Full_MetadataHR <- Full_MetadataHR[, !(colnames(Full_MetadataHR) %in% c("Fields 1", "Fields 2","Fields 3", "Natural" ))]
    saveRDS(Full_MetadataHR, file = print(paste("FixedMovement/Ring_ID", i, "HR_HabitatMetadata2023.rds")))
   
     print(paste("Ring_ID", i, "Metadata days"))
    MetaDays <- MetadataDays(Listoflapwings) #Create the daily metadata
    
    #Save the daily metadata separately
    saveRDS(MetaDays, file = print(paste("AllMovement/Ring_ID", i, "MetadataMovement2023.rds")))
    
    #Create the ground metadata with ground duration, and COA
    
    print(paste("Starting ring", i, "Ground metadata"))
    
    FlagTwo <- DurationTimeADPByDate(Listoflapwings) #Create the ADP metadata by dates
    FlagTwo<-FlagTwo[order(FlagTwo$Ring_ID,FlagTwo$date),] #make sure data is sorted chronologically (per tag)
    DateNumber <- unique(FlagTwo$date)
    DateData <- list()
    
    for (ii in 1:length(DateNumber)) { #Create the core area
      Date1 <- subset(FlagOne, date == DateNumber[ii])
      Date2 <- subset(FlagTwo, date ==  DateNumber[ii])
      Date1 <- Date1[, (colnames(Date1) %in% c("X", "Y"))]
      coreArea <- subset(Date2, VisitsFactor == "coreArea")
      
      if (nrow(coreArea) <= 0) { #If the lapwing is not in the core area, add NA
        if (!exists('StartLocations')) {
          Date2$MaxDisp <- NA
        } else { #Create max displacement from COA
          LocationsData <- as.data.frame(Date1[, (colnames(Date1) %in% c("X", "Y"))], ncol=2, nrow= length(Date1$X), byrow=TRUE)
          NewMetrix <- rbind(StartLocations[c(9,10)],LocationsData)
          MaxDisplacementCore <- max(as.matrix(dist(NewMetrix))[1,]) #MaxDispFromCore
          Date2$MaxDisp <- NA
          Date2$MaxDisp <- unique(MaxDisplacementCore)
        } } else {
          StartLocations <- as.data.frame(c(coreArea), ncol=2, nrow= length(Date1$X), byrow=TRUE)
          StartLocations <- StartLocations[!duplicated(StartLocations), ]
          colnames(StartLocations)[9] <- "X"
          colnames(StartLocations)[10] <- "Y"
          LocationsData <- as.data.frame(Date1[, (colnames(Date1) %in% c("X", "Y"))], ncol=2, nrow= length(Date1$X), byrow=TRUE)
          
          NewMetrix <- rbind(StartLocations[c(9,10)],LocationsData)
          MaxDisplacementCore <- max(as.matrix(dist(NewMetrix))[1,]) #MaxDispFromCore
          Date2$MaxDisp <- NA
          Date2$MaxDisp <- unique(MaxDisplacementCore)
        }
      DateData[[ii]] <- Date2
    }
    DateData2 <- do.call(rbind.data.frame, DateData)
    
    #Save the new ground locations:
    
    saveRDS(DateData2, file = print(paste("FixedMovement/NewGround Ring", i, "DayGroundMetadata2023WithCoreDis.rds")))
    
  }
  
}
