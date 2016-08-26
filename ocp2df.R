source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

##################################################
#################### ocp2df ######################
##################################################
ocp2df <- function(datafilename) {
   ## Description:
   ##   Reads time-voltage data (from VMP2 potentiostat)
   ##   and returns a dataframe with the data.
   ## Usage:
   ##   ocp2df(datafilename)
   ## Arguments:
   ##   datafilename: vector with full path to experimental files
   ## Value:
   ##   Dataframe with the following columns:
   ##   $ substrateid     : chr
   ##   $ sampleid        : chr   
   ##   $ time            : num
   ##   $ potential       : num
   #
   
   # Parameter table (add or remove parameters here)
   param.mtx <- 
      matrix(c("n_HeaderLines",          "^Nb header lines",          " : ",
               "VMPChannel",             "^Run on channel",           " : ",
               "AcqStart",               "^Acquisition started on",   " : ",
               "Device",                 "^Device",                   " : ",
               "Computer",               "^Address",                  " : ",
               "ECLabVersion",           "^EC-Lab for windows",       "4",
               "ElectrodeMaterial",      "^Electrode material",       " : ",
               "InitialState",           "^Initial state",            " : ",
               "Electrolyte",            "^Electrolyte",              " : ",
               "Comments",               "^Comments",                 " : ",
               "ElectrodeSurfaceArea",   "",                          " : ",  
               "CharacteristicMass",     "^Characteristic mass",      " : ",
               "EquivalentWeight",       "^Equivalent Weight",        " : ",
               "Density",                "^Density",                  " : ",
               "t",                      "^tR",                       "3", # duration of experiment
               "dEdt",                   "^dER/dt",                   "3",
               "dE",                     "^dER \\(mV\\)",             "3",
               "dtR",                    "^dtR",                      "3",
               "ErangeMin",              "^E range min",              "5",
               "ErangeMax",              "^E range max",              "5"), 
             ncol = 3, byrow = T)
   
   # Construct full parameter dataframe
   param.df <- 
      structure(list(param.mtx[, 1],
                     param.mtx[, 2],
                     param.mtx[, 3],
                     rep("", dim(param.mtx)[1])),
                .Names = c("parameter", 
                           "regexp",
                           # separator: 
                           # either the symbol used by strsplit, 
                           # such as ":", or a number indicating 
                           # the position of the desired value if 
                           # strsplit(sep = " ") is used
                           "separator", 
                           "value"),
                row.names = seq(1, dim(param.mtx)[1]),
                class = "data.frame")    
      
   
   # Read all input datafiles including parameters
   ff <- data.frame(NULL)
   number.ffrow <- 0
   for (f in 1:length(datafilename)) {
      dfile <- file(datafilename[f], "r")
      # read all lines of input file
      chifile <- readLines(dfile, n = -1)
      close(dfile)
      
      substrateid <- 
         ProvideSampleId(datafilename[f], implementation = "dirname")
      sampleid <- 
         paste(strsplit(basename(datafilename[f]), "[_-]")[[1]][c(1, 2, 3)], 
               collapse = "-")
      
      # Collect params from chifile
      for (i in 1:dim(param.df)[1]) {
         # Only grep the value if the regexp is not empty
         if (param.df$regexp[i] != "") {
            # save the matched header row to tmp variable param.row
            param.row <- grep(param.df$regexp[i], chifile, perl = T, value = T)
            if (length(grep("^[0-9]*$", param.df$separator[i])) < 1) {
               # not a number, use the separator as actual separator
               param.df$value[i] <- 
                  strsplit(x = param.row, split = param.df$separator[i])[[1]][2]
               # some param fields are empty in datafile, if so convert to ""
               param.df$value[i] <- 
                  ifelse(is.na(param.df$value[i]), "", param.df$value[i])
            } else {
               # separator is a number, is it as "field-selector"
               param.df$value[i] <- 
                  strsplit(x = param.row, split = "\\s+")[[1]][as.numeric(param.df$separator[i])]
            }
               
         } else {
            # if the regexp in vmp.param is empty, don't bother grepping
            param.df$value[i] <- ""
         }
      }
      # remove extraneous spaces in all param values
      param.df$value <- sub("\\s+$", "", sub("^\\s+", "", param.df$value))
      
      # before including the parameters into the ff-dataframe, we remove all the 
      # no longer required columns and transpose it
      tparam <- as.data.frame(t(param.df[, c("parameter", "value")]))
      colnames(tparam) <- tparam[1, ]
      tparam <- tparam[-1, ]
      
      # the VMP OCP files generally contain a single experiment (range),
      # so this simple approach works
      number.of.header.lines <-
         as.numeric(param.df$value[which(param.df$parameter == "n_HeaderLines")])
      number.of.data.lines <- 
         length(chifile) - number.of.header.lines
      
      ff <- 
         rbind(ff, 
               data.frame(substrateid,
                          sampleid,
                          matrix(scan(textConnection(chifile[(number.of.header.lines + 1):length(chifile)], "r"), 
                                      what = numeric(), 
                                      sep = ""), 
                                 ncol = 4, 
                                 byrow = T),
                          tparam,
                          row.names = seq(number.ffrow + 1, number.ffrow + number.of.data.lines)))
      
      # keep track of total number of datalines (remember, we are inside a for-loop)
      number.ffrow <- number.ffrow + number.of.data.lines
   }

   # rename columns
   names(ff) <- 
      c("substrateid", 
        "sampleid",
        "mode",
        "error",
        "time",
        "potential",
        colnames(tparam))
   
   #
   return(ff)
}
