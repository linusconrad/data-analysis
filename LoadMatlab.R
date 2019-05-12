# This script is a test to load .mat files exported from abf by matlab
# or .mat files exported by WinWCP

library(tidyverse)
library(magrittr)
library(plyr)
library(purrr)
library(zoo)
library(signal)
library(abf2)
library(R.matlab)

# this is for .mat file written by winWcp

# remove variable names from a list or data frame
# need to do this because otherwise columns are not merged later on
unname = function(df){
  thedataframe =df
  #print(names(df))
  names(thedataframe) = c(NULL)
  thedataframe
}

# rename variables in a two col dataframe I and V
# needs to be done so columns are merged later
name.IV = function(df){
  thedataframe = df
  #print(names(df))
  names(thedataframe) = c("I","V")
  thedataframe
}

addmeta = function(df, metadata, cond){
  # copy the metadata for each row of a dataframe and add them to it
  # "metadata" is a character vector
  # "cond" is a character vector with the names of the metadata variables
  # eg. metadata  = c(120, "TreK mutant 21", "baseline")
  # cond = c("conc potassium", "construct", "condition")
  
  # this creates the rows of metadata
  var =
    ldply(metadata, .fun = rep.int, times = nrow(df)) %>%
    t(.) %>%
    as_tibble(.)
  
  # name the columns correctly
  names(var) = cond
  #  then bind them together
  bind_cols(df, var)
}

# function to load from WinWCP
read.mat.WinWCP = 
  function(file, nsweep){
    raw.mat = 
      readMat(file) %>%
      unname(.)
    
    # this wrangles the data into tidy shape
    tibble.raw =
      # this loads the dataframes in the .mat object and assigns them by the type of trace
      tibble(t = raw.mat[c(seq(1, 2 * nsweep - 1, 2))],
             IV = raw.mat[c(seq(2, 2 * nsweep, 2))],
             sweep = c(1:nsweep))
    tibble.tIV =
      tibble.raw %>%
    # renames the traces and joins them together
      mutate(
        .,
        t = purrr::map(t, ~ data.frame(t = .)),
        t = purrr::map(t, ~ unlist(unname(.))),
        tIV =  purrr::map(IV, ~ data.frame(
          t = t[1],
          I = .[, 1],
          V = .[, 2]
        ))
      ) %>%
      unnest(., tIV)
    
    # as a generic function, just write the filename to the data
    # filenames should be constructed so they contain the metadata and can be extracted accordingly
    # alternatively use a lookuptable and left_join by filename
    
    # join the meta and return data
    addmeta(tibble.tIV, metadata = c(file), cond = c("filename"))
  }

# this is for a .mat file written by Matlab with the script of Manuel
read.mat.Matlab =
  function(file, SR, nsweep){
    raw.mat =
      readMat(file)
    
    # total number of rows of the file at the end
    nrows = raw.mat$File[2] %>%  unlist(.) %>% length(.)
    # number of sample in a sweep
    n.per.sweep = nrows / nsweep
    # time of one sweep (s)
    tend = nrows/ SR / nsweep
    
    tibble.tIV =
      tibble(
        # time in sec
        t = c(rep.int(seq(from = 0, to = tend, length.out = n.per.sweep), nsweep)),
        V = raw.mat$File[2] %>%  unlist(.),
        I = raw.mat$File[3] %>%  unlist(.),
        # assign sweep
        sweep = rep(1:nsweep, each = n.per.sweep)
      )
    # join the meta and return data
    addmeta(tibble.tIV, metadata = c(file), cond = c("filename"))
  }



