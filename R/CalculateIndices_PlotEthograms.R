#!/usr/bin/env Rscript


# Author: Aaron M Allen
# Date: 2018.10.09
#     updated: 2021.12.02
#
#
# Description:
# This script takes in the JAABA scores data and computes indices for courtship and plots ethograms.
#


## Input Arguments:
##


## Output
##



## Dependencies:
##     'tidyverse' package
##         install.packages("tidyverse")

## This function was written using R 4.x.x, and tidyverse 1.x.x





## Example usage:

##


















suppressMessages(library("data.table"))
suppressMessages(library("tidyverse"))
suppressMessages(library("cowplot"))




OutputDirectory <- commandArgs(trailingOnly = T)[2]
FileName <- commandArgs(trailingOnly = T)[3]
FPS <- as.numeric(commandArgs(trailingOnly = T)[4])
PropMale <- as.numeric(commandArgs(trailingOnly = T)[5])
OptoLight <- as.logical(str_to_upper(commandArgs(trailingOnly = T)[6]))



message(paste0("OutputDirectory = ",OutputDirectory))
message(paste0("FileName = ",FileName))
message(paste0("FPS = ",FPS))
message(paste0("PropMale = ",PropMale))
message(paste0("OptoLight = ",OptoLight))



LogFile <-file(paste0(OutputDirectory,"/",FileName,"/Logs/CalculateIndicesError.log"))
tryCatch({

    input_dir <- paste0(OutputDirectory,"/",FileName,"/Results/")
    my_data_file <- list.files(paste0(OutputDirectory,"/",FileName,"/Results/")) %>% str_subset("ALLDATA_R.csv.gz")
    message(paste0("Loading ",my_data_file))
    raw_data <- fread(paste0(input_dir,my_data_file),sep = ",",nThread = 8, showProgress = FALSE)

    message(paste0("Wrangling the Data..."))
    raw_data_spread <- raw_data %>%
        select(-Units, -Data_Source) %>%
        spread("Feature","Value") %>%
        mutate(Fly_Id = as.factor(Fly_Id))

    # Calculate Indices Table
    #############################################
    source("../R/calculate_single_indices_table.R")
    message(paste0("Calculating indicies..."))
    calculate_single_indices_table(input = raw_data_spread,
                                    court_init = TRUE,
                                    max_court = TRUE,
                                    predict_sex = TRUE,
                                    frame_rate = FPS,
                                    prop_male = PropMale,
                                    save_path = paste0(OutputDirectory,"/",FileName,"/Results/")
                                    )

    ## Run calculate_indices 2 more times if the optogenetic light was used. One file for
    ## lights-on and one file for light-off.
    # if ( OptoLight == "yes" ) {
    #     calculate_single_indices_table(input = raw_data_spread,
    #                                     court_init = TRUE,
    #                                     max_court = TRUE,
    #                                     predict_sex = TRUE,
    #                                     frame_rate = FPS,
    #                                     prop_male = PropMale,
    #                                     opto_light = "on",
    #                                     save_path = paste0(OutputDirectory,"/",FileName,"/Results/")
    #                                     )
    #     calculate_single_indices_table(input = raw_data_spread,
    #                                     court_init = TRUE,
    #                                     max_court = TRUE,
    #                                     predict_sex = TRUE,
    #                                     frame_rate = FPS,
    #                                     prop_male = PropMale,
    #                                     opto_light = "off",
    #                                     save_path = paste0(OutputDirectory,"/",FileName,"/Results/")
    #                                     )
    #
    # }

    # Ethogram Plots of JAABA Classifiers
    #############################################
    source("../R/plot_jaaba_ethograms.R")
    message(paste0("Plotting ethograms..."))

    FlyId <- as.numeric(unique(raw_data_spread$Fly_Id))
    OddFly <- FlyId[FlyId %% 2 == 1]

    EthoFileName <- paste0(OutputDirectory,"/",FileName,"/Results/",FileName, "_Ethogram.pdf")
    pdf(EthoFileName,width=10,height=7,paper='a4r')
    for (P in OddFly){
        Plot1 <- plot_jaaba_ethograms(raw_data_spread, P)
        Plot2 <- plot_jaaba_ethograms(raw_data_spread, P+1)
        print(plot_grid(Plot1,Plot2, ncol = 1))
    }
    dev.off()

    #############################################

    },
    error=function(e) {
        writeLines(as.character(e), LogFile)
    }
)



sessionInfo()
