# data-analysis
This hosts script to analyse ephys data from raw data to neat summary graphs

This deals with electrophysiological data. The raw data is recorded with WinWCP and written to its native file format.
Then I export them to ACSII and wrangle them in R.
These are thye scripts to do different analysis tasks on multiple files generated like that.

Usually, data files to run analysis on are recorded and exported in a way that they contain all metadata in the filename (concentrations, solutions etc).

A function that loads the data of one file and adds the metadata to each datum is written .
This is fed to ldply and run on all files.
You end up with a big tidyR compatible dataframe to wrangle.

For protocols that contain more samples (long delta-T style protocols, step protocols) a function is written that loads the data and writes some summary for each file to the folder.
Then you again ldply over the summary file to get your data to work on.

