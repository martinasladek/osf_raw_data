#appropriate for exploding raw data from results for nonquestionnaires
deJSON.data <- function(df,json_name = "data",convert2num=FALSE) {
  orig_col <- ncol(df)
  newdf0 <- as.data.frame(cSplit(df,json_name,sep="},{",direction="long",type.convert=FALSE))
  newdf <- as.data.frame(cSplit(newdf0,json_name,sep=",",direction="wide",type.convert=FALSE))

  #assign variable names to new exploded columns
  one_trial <- newdf[1,orig_col:ncol(newdf)]
  one_trial <- gsub(":.*","",one_trial)
  one_trial <- gsub("\\[","",one_trial)
  one_trial <- gsub("\\{","",one_trial)
  one_trial <- gsub("\"","",one_trial)
  names(newdf)[orig_col:ncol(newdf)] <- one_trial

  for (i in orig_col:ncol(newdf)) {
    newdf[,i] <- gsub("\"","",newdf[,i])
    newdf[,i] <- gsub(".*:","",newdf[,i])
    newdf[,i] <- gsub("}","",newdf[,i])
    newdf[,i] <- gsub("\\{","",newdf[,i])
    newdf[,i] <- gsub("\\[","",newdf[,i])
    newdf[,i] <- gsub("]","",newdf[,i])
    if (convert2num) {
      newdf[,i] <- as.numeric(as.character(newdf[,i]))
    }
  }
  return(newdf)
}

#This function cuts rows from a data frame based on whether any of a specified set of values is present in a specified column. It is intended to cut rows containing practice trials from a series of participant trials.
#Its primary input is a data frame, with each row containing a single participant trial. Because this function does not calculate any participant-specific values (unlike trim), it can be run on a data frame containing multiple participants' trials.
#The field and values parameters together specify the criteria for cutting a line. Field is a string, specifying the column of the dataframe where test vs. practice (or another criterion for cutting a line) is indicated.
#Values is a vector of strings. For each row in the dataframe (i.e. each trial), if any of the strings in values is found in the column specified by field, the row is cut.
#Values is given as a vector rather than a string because some tests use multiple strings to indicate a practice trial (e.g. "prac1", "prac2", "prac3)
#This function returns the dataframe with rows containing any of the strings in values within the field column removed completely.
exclude <- function(df, field, values){ #df is the table of trial-by-trial data, field is the column where practice vs. trial info is stored (e.g. "type"), values are the strings that signify a row to be excluded (e.g. "practice", "prac"){
  for(value in values){
    df <- df[df[,field] != value,]
  }
  return(df)
}

#This function trims outlying values from the reaction time of a participant or series of participants.
#Its primary input is a data frame including a column labeled "rt" that contains trial by trial reaction times.
#The function will return the same data frame with outlying reaction times replaced with NA. What values are marked
#as outlying is determined by the other parameters.
#If a numerical value is given for the nsd parameter, the function will calculate the mean reaction time for each participants' trials and trim all values that are more than nsd standard deviations from that mean for that participant. If no value is given, the function will not trim based on deviation from the individual mean.
#If the zeros parameter is set to TRUE, the function will trim all zero values. This parameter defaults to TRUE, so if zeros is not given a value when calling the function, all zeros will be replaced with NA. If zeros is set to FALSE, zero values will be left unchanged.
#The min_ms parameter gives a minimum value for reaction times; all values below this value will be set to NA. Min_ms defaults to 150, corresponding to the lowest realistic human reaction time (in ms). This value can be set higher (for more time-consuming tasks) or lower (to include more participants). If you do not want to trim below a certain threshold, set min_ms to FALSE.
#The max_ms parameter gives a maximum value for reaction times; all values above this value will be set to NA. Max_ms has a default value of FALSE, so the function will not trim high outlying values unless a numerical value is provided.
#When the showtrimmed parameter is set to TRUE, the function will track how many values are trimmed from each participant's trials and print this number. This parameter defaults to FALSE, so the function will not track the number of values trimmed if showtrimmed = TRUE is not specified.
trim_v2 <- function(df, nsd = FALSE, zeros = TRUE, min_ms = 150, max_ms = FALSE, showtrimmed = FALSE){
  rts <- df$rt
  df_plain <- df
  if(showtrimmed != FALSE){ #tracks the number of trials being trimmed
    ntrimmed = 0
  }
  if(min_ms != FALSE){ #cuts reaction times less than the min_ms cutoff - default is 150 ms
    df[is.na(df$rt) == F & df$rt < min_ms,]$rt <- NA
  }
  #return(df)
  if(max_ms != FALSE & nrow(df[is.na(df$rt) == F & df$rt > max_ms,]) != 0){ #cuts reaction times greater than the max_ms cutoff if one is provided
    df[is.na(df$rt) == F & df$rt > max_ms,]$rt <- NA
  }
  if(zeros != FALSE & nrow(df[is.na(df$rt) == F & df$rt == 0,]) != 0){ #cuts reaction times of 0 (sometimes used to mark timeouts)
    df[is.na(df$rt) == F & df$rt == 0,]$rt <- NA
  }
  ids <- unique(df$id)
  for(i in 1:length(ids)){
    id <- ids[i]
    initntrimmed <- nrow(df_plain[is.na(df_plain$rt) == T & df_plain$id == id,])
    #print(id)
    sd <- sd(df_plain[df_plain$id == id,]$rt, na.rm = T)
    mean <- mean(df_plain[df_plain$id == id,]$rt, na.rm = T)
    if(nsd != F & nrow(df[is.na(df$rt) == F & df$id == id & (df$rt > mean + (nsd * sd) | df$rt < mean - (nsd * sd)),]) != 0){
      df[(is.na(df$rt) == F) & df$id == id & (df$rt > mean + (nsd * sd) | df$rt < mean - (nsd * sd)),]$rt <- NA
    }
    ntrimmed <- nrow(df[is.na(df$rt) == T & df$id == id,])
    #print(paste(ntrimmed - initntrimmed, "trials were trimmed"))
  }
  return(df)
}

#This function drops participants who have more than a threshold number of NA response times. It is meant
#to be used after the trim function has been used to replace outlying response times with NAs. The function
#will print the ID of each participant that is trimmed and the number of outlying and timeout trials
#that this participant had. The df parameter should be a data frame containing one trial in each line - it
#can contain any number of participants. The threshold can be provided as either an absolute number (a value
#1 or above) or a proportion of overall trials (a value from 0 to 1). If you are using an absolute number of
#NAs as the threshold, set the threshold parameter to this number; if you are using a proportion, set the proportion
#parameter.
toomanynas_v2 <- function(df, threshold = FALSE, proportion = FALSE){
  if(threshold == FALSE & proportion == FALSE){
    print("You have not provided exclusion criteria.")
    return(df)
  }
  else{
    count <- 0 #initializes count of trimmed participants
    ids <- unique(df$id)
    for(i in 1:length(ids)){ #checks one participant at a time
      trials <- df[df$id == ids[i],] #selects only trials from a single participant
      narows <- trials[is.na(trials$rt) == TRUE,]
      nas <- nrow(narows) #total number of NA reaction times
      if("state" %in% colnames(narows)){ #only if the test uses timeouts and stores them in the results data
        timeouts <- narows[narows$state == "timeout",] #sorts NA reaction times into those due to timeouts
        fastrts <- narows[narows$state != "timeout",] #and those due to other criteria.
      }
      else{
        fastrts <- narows
        timeouts <- NULL
      }
      if(threshold != FALSE){
        if(nas > threshold){ #trims participant if threshold is exceeded and prints the number of each kind of NA.
          df <- df[df$id != ids[i],]
          count <- count + 1 #counts number of participants trimmed
          cat("Participant ", ids[i], " cut with ", nrow(timeouts), " timeouts and ", nrow(fastrts), " outlying RTs.\n")
        }
      }
      else if(proportion != FALSE){
        if(nas/nrow(trials) > proportion){ #trims participant if proportion is exceeded and prints reason for exclusion
          df <- df[df$id != ids[i],]
          count <- count + 1 #counts number of participants trimmed
          cat("Participant ", ids[i], " cut with ", nrow(timeouts), " timeouts and ", nrow(fastrts), " outlying RTs.\n")
        }
      }
    }
    cat("Number of participants removed: :", count, "\n") #prints total number of participants trimmed
    return(df)#returns data frame containing only trials from participants without NA trials above threshold
  }
}

#This function takes in all of a single participant's trials on simple reaction time and calculates their mean
#reaction time, median reaction time, standard deviation of reaction time, and coefficient of variation of
#reaction time (sd rt/mean rt). It outputs these summary scores as a vector, which also contains participant-level metadata.
#The df parameter should be a data frame containing all trials for a single participant (in order to summarize multiple
#participants, the function must be run in a loop - see below for an example of how to do this).
#The data_begin_column parameter specifies the number of the column where participant-level metadata (e.g. ID, start time)
#ends and trial-level data begins. All data before this column will be saved as metadata and appended to the vector returned
#by the function unaltered.
#The prefix parameter, which takes a string, allows the addition of a prefix string to the column names of the results vector.
#For instance, if using the function for even/odd reliability testing, you could have it label the columns with "odd" or "even."
#This parameter defaults to FALSE, so no prefix will be used if one is not specified.
#The expected_trials parameter takes a numerical value corresponding to the number of trials each participant should have completed.
#If it is set to TRUE, it will add a "complete" value to the results vector that is 1 if the expected number of trials were
#completed and 0 if the true number of trials completed differs from the expected number. This parameter defaults to FALSE, so
#if it is not set to TRUE the "complete" value will be set as NA.
#The report parameter, if set to TRUE, causes the function to print the ID of the participant being run before summarizing. This is
#useful for identifying incorrect or missing data that is causing errors when the function runs. This parameter defaults to FALSE.
summarize_srt <- function(df, data_begin_column, prefix = FALSE, expected_trials = FALSE, report = FALSE){
  if(report == TRUE){ #prints ID of participant whose data is being summarized
    print(df$id[1])
  }
  meta <- df[2,1:(data_begin_column-1)]  #get meta (or person level) data
  if(expected_trials != FALSE){ #Checks true number of trials against expected number (if provided) and fills in 1 if they are the same or 0 if they are not.
    if(nrow(df) == expected_trials){
      complete <- 1
    }
    else{
      complete <- 0
    }
  }
  else{
    complete <- NA
  }
  ntrials <- nrow(df) #Saves the number of trials for this participant
  meanrt <- mean(df$rt, na.rm = TRUE) #Calculates mean reaction time, ignoring NA values
  medianrt <- median(df$rt, na.rm = TRUE) #Calculates median reaction time, ignoring NA values
  sdrt <- sd(df$rt, na.rm = TRUE) #Calculates standart deviation of reaction time, ignoring NA values
  cv <- sdrt/meanrt #Calculates coefficient of variation
  summary <- cbind(meanrt, medianrt, sdrt, cv, ntrials, complete) #Binds all summary scores into one vector with metadata
  if(prefix != FALSE){ #Adds prefix (if specified) to column names
    colnames(summary) <- paste(prefix, colnames(summary), sep = "_")
  }
  summary <- cbind(meta, summary)
  return(summary)
}

#This function takes in all of a participant's trials on digit symbol coding and calculates their
#number of trials, number of correct trials, accuracy (correct trials/total trials), mean reaction time,
#median reaction time, and standard deviation of reaction time. All summary scores other than accuracy are
#calculated on the subset of trials where the participant answered correctly.
#It outputs these summary scores as a vector, which also contains participant-level metadata.
#The df parameter should be a data frame containing all trials for a single participant (in order to summarize
#multiple participants, the function must be run in a loop - see below for an example of how to do this).
#The data_begin_column parameter specifies the number of the column where participant-level metadata
#(e.g. ID, start time) ends and trial-level data begins. All data before this column will be saved as
#metadata and appended unaltered to the vector returned by the function.
#The prefix parameter, which takes a string, allows the addition of a prefix string to the column names
#of the results vector. For instance, if using the function for even/odd reliability testing, you could
#have it label the columns with "odd" or "even." This parameter defaults to FALSE, so no prefix will be
#used if one is not specified.
#The report parameter, if set to TRUE, causes the function to print the ID of the participant being run before summarizing.
#This is useful for identifying incorrect or missing data that is causing errors when the function runs. This parameter defaults to FALSE.
summarize_dsc <- function(df, data_begin_column, prefix = FALSE, report = FALSE){
  if(report == TRUE){ #prints ID of participant whose data is being summarized
    print(df$id[1])
  }
  meta <- df[2,1:(data_begin_column-1)]  #get meta (or person level) data
  ntrials <- nrow(df)
  ncorrect <- nrow(df[df$correct == 1,])
  accuracy <- ncorrect/ntrials #proportion of correct trials
  correctdf <- df[df$correct == 1,] #creates a dataframe of only correct trials - all other summary stats will only include these trials
  meanrt <- mean(correctdf$rt, na.rm = TRUE) #mean reaction time for correct trials
  medianrt <- median(correctdf$rt, na.rm = TRUE) #median reaction time for correct trials
  sdrt <- sd(correctdf$rt, na.rm = TRUE) #sd reaction time for correct trials
  ies <- medianrt/accuracy
  cv <- sdrt/meanrt
  summary <- cbind(ntrials, ncorrect, accuracy, meanrt, medianrt, sdrt, ies, cv) #creates a vector of summary scores
  if(prefix != FALSE){ #adds label to column names if specified
    colnames(summary) <- paste(prefix, colnames(summary), sep = "_")
  }
  summary <- cbind(meta, summary)
  return(summary)
}

#This function takes in all of a single participant's trials on vocabulary and returns a vector
#The df parameter should be a data frame containing data for all trials for a single participant (to run
# this function on more than one particpant, the function must be run in a loop - see below for an example).
#The data_begin_column parameter specifies the number of the column where participant-level metadata
#(e.g. ID, start time) ends and trial-level data begins. All data before this column will be saved as
#metadata and appended unaltered to the vector returned by the function.
#The prefix parameter, which takes a string, allows the addition of a prefix string to the column names
#of the results vector. For instance, if using the function for even/odd reliability testing, you could
#have it label the columns with "odd" or "even." This parameter defaults to FALSE, so no prefix will be
#used if one is not specified.
#The expected_trials parameter takes a numerical value corresponding to the number of trials each participant should have completed.
#If it is set to TRUE, it will add a "complete" value to the results vector that is 1 if the expected number of trials were completed
#and 0 if the true number of trials completed differs from the expected number. This parameter defaults to FALSE, so if it is not set to TRUE the "complete" value will be set as NA.
#The report parameter, if set to TRUE, causes the function to print the ID of the participant being run before summarizing.
#This is useful for identifying incorrect or missing data that is causing errors when the function runs. This parameter defaults to FALSE.
summarize_vocab <- function(df, data_begin_column, prefix = FALSE, expected_trials = FALSE, report = FALSE){
  if(report == TRUE){ #prints ID of participant whose data is being summarized
    print(df$id[1])
  }
  meta <- df[2,1:(data_begin_column-1)]  #get meta (or person level) data
  ncorrect = 0
  if(expected_trials != FALSE){
    if(nrow(df) == expected_trials){
      complete <- 1
    }
    else{
      complete <- 0
    }
  }
  else{
    complete <- NA
  }
  ncorrect <- nrow(df[df$accuracy == 1,])
  ntrials <- nrow(df)
  accuracy <- ncorrect/ntrials #proportion of correct trials
  summary <- cbind(accuracy, ntrials, complete)
  if(prefix != FALSE){
    colnames(summary) <- paste(prefix, colnames(summary), sep = "_")
  }
  summary <- cbind(meta, summary)
  return(summary)
}

#installs packages if not installed, loads package if it has been installed
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}




