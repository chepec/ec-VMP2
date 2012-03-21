source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")
source("/home/taha/chepec/chetex/common/R/common/int2padstr.R")
source("/home/taha/chepec/chetex/common/R/common/hms2seconds.R")
require(lubridate)

##################################################
#################### amp2df ######################
##################################################
amp2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads data from VMP2 potentiostat
   ##   and returns a dataframe with the data and some 
   ##   calculated quantities based on the data.
   ## Usage:
   ##   amp2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeters)
   ## Value: (nneds to be updated)
   ##   Dataframe with the following columns:
   ##   $ sampleid           : chr
   ##   $ sid                : chr
   ##   $ DateTime           : chr
   ##   $ WEarea             : num [square cm]
   ##   $ mode               : num
   ##   $ oxred              : num
   ##   $ error              : num
   ##   $ controlchanges     : num
   ##   $ Nschanges          : num
   ##   $ counterinc         : num
   ##   $ time.second        : num [seconds]
   ##   $ control.volt       : num [volt]
   ##   $ ewe.volt           : num [volt]
   ##   $ current.mA         : num [mA]
   ##   $ charge.mAh         : num [mA hour == 3600 mC]
   ##   $ currentdensity     : num [mA per square cm]
   ##   $ timediff           : num [seconds]
   ##   $ currentdiff        : num [mA]
   ##   $ currentdensitydiff : num [mA per square cm]
   ##   $ dIdt               : num [mA per second]
   ##   $ didt               : num [mA per square cm per second]
   ##   $ charge             : num [mC]
   ##   $ chargedensity      : num [mC per square cm]
   ##   +++ plus two parameter dataframes attached as attributes
   ##   > mstep.param
   ##   > vmp.param
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   # Create a sampleid
   sampleid <- ProvideSampleId(datafilename, implementation = "dirname")
   #
   # Parameter table
   vmp.param <- data.frame(stringsAsFactors = FALSE,
                           matrix(c("sampleid",             "",
                                    "n_HeaderLines",        "^Nb header lines",
                                    "VMPChannel",           "^Run on channel",
                                    "Time",                 "^Acquisition started on",
                                    "ElectrodeMaterial",    "^Electrode material",
                                    "InitialState",         "^Initial state",
                                    "Electrolyte",          "^Electrolyte",
                                    "Comments",             "^Comments",
                                    "ElectrodeSurfaceArea", "", 
                                    "CharacteristicMass",   "^Characteristic mass",
                                    "EquivalentWeight",     "^Equivalent Weight",
                                    "Density",              "^Density"),
                                  ncol = 2, byrow = T))
   names(vmp.param) <- c("parameter", "regexp") 
   vmp.param$value <- as.numeric(NA)
   # Collect params from chifile
   for (i in 1:dim(vmp.param)[1]) {
      # Only grep the value if the regexp is not empty
      if (vmp.param$regexp[i] != "") {
         vmp.param$value[i] <- 
            # ifelse checks length of strsplit() before saving it
            ifelse(length(strsplit(chifile[grep(subset(vmp.param, 
                                                       parameter == vmp.param$parameter[i])$regexp,
                                                chifile,
                                                perl = T,
                                                value = F)], 
                                   " : ")[[1]]) == 1,
                   # if parameter is blank, save "NA"
                   "NA",
                   # if parameter exists, save the part to the right of " : "
                   strsplit(chifile[grep(subset(vmp.param, 
                                                parameter == vmp.param$parameter[i])$regexp, 
                                         chifile, 
                                         perl = T, 
                                         value = F)], 
                            " : ")[[1]][2])
      } else {
         # if the regexp in vmp.param is empty, don't bother grepping
         vmp.param$value[i] <- "NA"
      }
   }
   # Save the current sampleid to vmp.param
   vmp.param$value[which(vmp.param$parameter == "sampleid")] <- sampleid
   # Collect the method from the second line under "Nb header lines"
   vmp.param <- rbind(vmp.param,
                      c("Method", "", 
                        chifile[(grep(subset(vmp.param, parameter == "n_HeaderLines")$regexp,
                                     chifile, perl = T, value = F) + 2)]))
   # Trim extra whitespace from n_HeaderLines value
   vmp.param$value[which(vmp.param$parameter == "n_HeaderLines")] <- 
      sub("\\s+", "", vmp.param$value[which(vmp.param$parameter == "n_HeaderLines")])
   #
   # Create a more resilient id (more unique)
   sid <- paste(paste(int2padstr(ii = month(mdy_hms(vmp.param$value[which(vmp.param$parameter == "Time")])),
                                 pchr="0",
                                 w=2),
                      int2padstr(ii=day(mdy_hms(vmp.param$value[which(vmp.param$parameter == "Time")])),
                                 pchr="0",
                                 w=2),
                      sep = ""),
                sampleid, 
                sep = "-")
   #
   # Multi-step parameter table
   mstep.param <- data.frame(stringsAsFactors = FALSE,
                             matrix(c("Ns",        "^Ns",
                                      "Ei",        "^Ei \\(V\\)",
                                      "Ref",       "^vs\\.",
                                      "ti",        "^ti \\(h:m:s\\)",
                                      "Imax",      "^Imax",
                                      "ImaxUnit",  "^unit Imax",
                                      "Imin",      "^Imin",
                                      "IminUnit",  "^unit Imin",
                                      "dQM",       "^dQM",
                                      "dQMUnit",   "^unit dQM",
                                      "record",    "^record",
                                      "dI",        "^dI",
                                      "dIUnit",    "^unit dI",
                                      "dQ",        "^dQ\\s",
                                      # space after dQ needed to distinguish from dQM
                                      "dQUnit",    "^unit dQ\\s",
                                      # space after dQ needed to distinguish from dQM
                                      "dt",        "^dt \\(s\\)",
                                      "dta",       "^dta \\(s\\)",
                                      "EMin",      "^E range min \\(V\\)",
                                      "EMax",      "^E range max \\(V\\)",
                                      "IRange",    "^I Range", 
                                      "Bandwidth", "^Bandwidth",           
                                      "gotoNs",    "^goto Ns'",
                                      "nccycles",  "^nc cycles"),
                                    ncol = 2, byrow = T))
   names(mstep.param) <- c("parameter", "regexp")
   mstep.param$rowno <- as.numeric(NA)
   # Collect row numbers for multi-step parameters from chifile
   for (i in 1:dim(mstep.param)[1]) {
      mstep.param$rowno[i] <- 
         grep(subset(mstep.param, 
                     parameter == mstep.param$parameter[i])$regexp, 
              chifile, 
              perl = T, 
              value = F)
   }
   # Collect the multi-step parameters
   # Number of steps in experiment (aka columns in multi-step param table)
   n_steps <- length(strsplit(sub("\\s+", " ", 
                                  sub("\\s*$", "",  
                                      sub("^\\s*", "", 
                                          sub(mstep.param$regexp[1], "", 
                                              chifile[mstep.param$rowno][1], 
                                              perl = T)))), 
                              "\\s")[[1]])
   # Prepare names for extra cols
   extra.colnames <- paste("E", 
                           as.character(seq(dim(mstep.param)[2]+1, 
                                            dim(mstep.param)[2]+n_steps)), 
                           sep = "")
   # Expand mstep.param with enough extra columns to fit step-wise params
   mstep.param[, extra.colnames] <- NA
   # Assign step-wise parameters to mstep.param
   for (j in 1:length(extra.colnames)) {
      for (k in 1:dim(mstep.param)[1]) {
         mstep.param[k, which(names(mstep.param) == extra.colnames[j])] <-
            # replace commas with dots
            sub(",", "\\.", 
                # replace multiple whitespace with single whitespace
                strsplit(sub("\\s+", " ", 
                             # trim whitespace at end of string
                             sub("\\s*$", "",  
                                 # trim whitespace at start of string
                                 sub("^\\s*", "", 
                                     # remove the regexp-matching part of the line
                                     sub(mstep.param$regexp[k], "", 
                                         chifile[mstep.param$rowno][k], 
                                         perl = T)))), 
                         "\\s")[[1]][j])
      }
   }
   #
   #
   chidata <- chifile[(as.numeric(vmp.param$value[which(vmp.param$parameter == "n_HeaderLines")]) + 1):length(chifile)]
   # Replace decimal commas with dots
   chidata <- gsub(",", "\\.", chidata)
   zz <- textConnection(chidata, "r")
   data.exp <- data.frame(stringsAsFactors = FALSE,
                          sampleid,
                          sid,
                          vmp.param$value[which(vmp.param$parameter == "Time")],
                          wearea,
                          matrix(scan(zz, 
                                      what = numeric(), 
                                      sep = "\t"),
                                 ncol = 11, 
                                 byrow = T))
   close(zz)
   names(data.exp) <- c("sampleid",
                        "sid",
                        "DateTime",
                        "WEarea",
                        "mode",
                        "oxred",
                        "error",
                        "controlchanges",
                        "Nschanges",
                        "counterinc",
                        "time.second",
                        "control.volt",
                        "ewe.volt",
                        "current.mA",
                        "charge.mAh")
   # Calculate current density
   data.exp$currentdensity <- data.exp$current.mA / wearea
   # Calculate time diff and current diff
   data.exp$timediff <- c(data.exp$time.second[1], diff(data.exp$time.second))
   data.exp$currentdiff <- c(data.exp$current.mA[1], diff(data.exp$current.mA))
   data.exp$currentdensitydiff <- c(data.exp$currentdensity[1], diff(data.exp$currentdensity))
   # Calculate differential of current and current density
   data.exp$dIdt <- data.exp$currentdiff / data.exp$timediff
   data.exp$didt <- data.exp$currentdensitydiff / data.exp$timediff
   # Calculate charge and charge density
   data.exp$charge <- 3600 * data.exp$charge.mAh # convert from mAh to mC
   data.exp$chargedensity <- data.exp$charge / wearea
   #
   # Convert time fields in mstep.param to seconds
   mstep.param[which(mstep.param$parameter == "ti"), extra.colnames] <- 
      hms2seconds(unlist(mstep.param[which(mstep.param$parameter == "ti"), extra.colnames]))
   # Save mstep.param as attrib to data.exp
   attr(data.exp, "mstep.param") <- 
      mstep.param[, which(names(mstep.param) %in% c("parameter", extra.colnames))]
   # Save vmp.param as attrib to data.exp
   attr(data.exp, "vmp.param") <- 
      vmp.param[, which(names(vmp.param) %in% c("parameter", "value"))]
   #
   return(data.exp)
}