process.SGsmooth.deact =
  
  function(txt, npersweep, sep, cond){
    
    require(signal)       # package for time series data (spectra, filter, etc) 
    require(plyr)         # split apply combine tools
    require(tidyverse)    
    require(magrittr)     # pipe operators
    require(broom)        # used to extract test and model objects in a sane way
    require(readr)        # read csv files
    require(stringr)      # string operation toolbox
    require(minpack.lm)   # Levenberg Markquart algorithm for curve-fitting
    require(nls.multstart)
    require(ggthemes)     # prettier plots
    require(plotrix)      # has sem function
    require(Cairo)        #
    require(purrr)
    require(cowplot)
    
    #print the filename to screen (debugging purposes)
    print(paste0("Start with processing file", txt))
    
    stepdata =
      read_delim(
        txt,
        "\t",
        escape_double = FALSE,
        col_names = c("t", "I", "V"),
        trim_ws = TRUE
      )
    
    print("Finished loading the file")
    #this is the metadata
    txtmetadata = unlist(strsplit(txt, sep))
    
    #we need this numbers to assign sweeps correctly
    nsamples = length(stepdata$t) #number of samples in the file
    #print(nsamples) #this is for debugging only
    nsweeps = nsamples / npersweep # number of sweeps in the file
    #print(nsweeps)
    
    #assign the sweep to the ephysdata
    stepdata$sweep = as.factor(rep(1:nsweeps, each = npersweep))
    
    #transform to physiological convention
    stepdata$I = c(stepdata$I)*-1
    stepdata$V = c(stepdata$V)*-1
    
    print(head(stepdata))
    ########### Assign these before doing all the math
    # change to use this script with another protocol
    t.pre    = 0.1       # start time of the pre-pulse
    t.test   = 0.6       # start time of the test-pulse  
    t.post   = 1.3       # start time of the final tail step (postpulse)
    t.steady = 1.9       # end of the tail step
    
    
    # create a variable for the voltage @ every step
    stepdata %<>%
      group_by(., sweep) %>%
      mutate(., Vpulse = round_any(V[t == 0.65], 1)) %>%  # put a timepoint within the step here to calculate
      ungroup(.)
    
    print(stepdata)
    
    #summarise per sweep values
    stepdata.persweep =
      stepdata %>%
      group_by(sweep, Vpulse) %>%
      summarise(.,
                I.inst              = mean(I[t > (t.test + 0.001) & t < (t.test + 0.002)]),  # Inst current @ pulse
                 I.50              = mean(I[t > (t.test + 0.050) & t < (t.test + 0.051)]),
                 I.100             = mean(I[t > (t.test + 0.100) & t < (t.test + 0.101)]),
                 I.150             = mean(I[t > (t.test + 0.150) & t < (t.test + 0.151)]),
                 I.200             = mean(I[t > (t.test + 0.200) & t < (t.test + 0.201)]),
                 I.400             = mean(I[t > (t.test + 0.400) & t < (t.test + 0.401)]),
                 I.ss               = mean(I[t > (t.post - 0.01)  & t < (t.post - 0.001)]),  #  steady state @ test pulse
                td.test = I.inst - I.ss,
                thresh.test = (I.inst+ I.ss)/2 ,
                I.tailSS            = mean(I[t > (t.steady - 0.01) & t < (t.steady - 0.001)]), #  steady state after pulse
                I.tail              = mean(I[t > (t.post + 0.001) & t < (t.post + 0.002)]),  # inst after pulse
                # measures of the pre-pulse to calculate tapp
                decay.pulse = 1- I.ss / I.inst, # percent decay of the test-pulse current, (deactivation ratio)
                # measures of the pre-pulse to calculate tapp
                I.inst.pre          = mean(I[t > (t.pre + 0.001) & t < (t.pre + 0.002)]),
                I.ss.pre            = mean(I[t > (t.test - 0.01)  & t < (t.test - 0.001)]),
                td.pre              = I.ss.pre - I.inst.pre,
                thresh = (I.inst.pre + I.ss.pre) /2,
                fold.pre = I.ss.pre / I.inst.pre # fold activation on the pre-pulse
      ) %>%
      ungroup(.)
    
    
    # normalise to peak tail current
    stepdata.persweep %<>%
      mutate(., norm.tail = I.tail/I.tail[Vpulse == max(Vpulse)],
             I.tailSSfull = I.tailSS[Vpulse == min(Vpulse)], #fully relaxed current
             td.tail = I.tail - I.tailSSfull,
             thresh.tail = (I.tail + I.tailSSfull)/2)
    
    print(stepdata.persweep)
    
    
    
    plot = list()
    
    plot$timecourse =
      stepdata %>% dplyr::filter(., t > (t.pre-0.05), t < t.steady) %>%
      ggplot(. , aes(t,I, group = Vpulse)) +
      geom_line()
    #
    
    #plot tailcurrents
    plot$tail =
      ggplot(stepdata.persweep, aes(x = Vpulse, y = norm.tail)) +
      geom_point()+
      geom_line() +
      theme(aspect.ratio = 1)
    
    #plot IV
    plot$IV =
      stepdata.persweep %>%
      gather(., key = Measure, value = I, I.ss, I.inst) %>%
      
      ggplot(., aes(x = Vpulse, y = I, shape = Measure))+
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour  = "grey50") +
      geom_vline(xintercept = 0,
                 linetype = 3,
                 colour  = "grey50") +
      scale_shape_manual(values = c(22,25), labels = c("Instantanous", "Steady state")) +#
      geom_point()+
      geom_line(linetype = 2) +
      theme(legend.position = "bottom")
    
    print("Writing sumamry plots to file")
    # 
    # # put them all on one page and write them to a file
    plot_grid(plot$timecourse,
              plot$tail,
              plot$IV,
              align = "h",
              axis = "bt",
              ncol =3) %>%
      ggsave(., filename = paste0(txt, "Plot.pdf"), height = 3.5, width = 10, path = "./tapp deact prot/")
    # 
    
    ##############
    # take some time course elements and apply smoothing to them
    # create a nested data frame separating the different time courses (pre, test, post)
    # into nested columns for each sweep
    
    # create a data frame then call the nesting of the main data set from inside and then column bind them
    # create a data frame then call the nesting of the main data set from inside and then column bind them
    tc.pre = 
      stepdata %>%
      filter(., t > t.pre + 0.001 &
               t < t.test) %>% # this is the pre pulse data (activating pulse)
      group_by(., sweep) %>%
      # set 0 point to instantaneous current
      # normalise to fully relaxed current
      # this works for activating timecourses!
      mutate(
        .,
        offset = mean(I[t > (t.pre + 0.001) &
                          t < (t.pre + 0.002)]),
        Isubs = I - offset,
        norm = mean(Isubs[t > (t.test - 0.01) &
                            t < (t.test - 0.001)]),
        I.normact = Isubs / norm ,
        t.fit = c(t - t.pre + 0.001)
      ) %>%
      nest(.key = "tc.pre")
    
    tc.test =
      stepdata %>%
      filter(., t > t.test + 0.001 &
               t < t.post) %>% # this is the testpulse data
      group_by(., sweep) %>%
      mutate(
        .,
        # set 0 point to steady current
        # normalise to instantanous current
        # this works for deactivating timecourses!
        offset = mean(I[t > t.post - 0.01 &
                          t < t.post - 0.001]),
        Iinst0 = I - offset,
        norm = mean(Iinst0[t > t.test + 0.001 &
                             t < t.test + 0.002]),
        I.normdeact = Iinst0 / norm,
        t.fit = c(t - t.test + 0.001)
      ) %>% #this creates a dataset to do extract the t1/2  from
      group_by(., sweep) %>%
      nest(.key = "tc.test")
    
    
    tc.tail =
      stepdata %>%
      filter(., t > t.post + 0.001 &
               t < t.steady) %>%
      ungroup(.) %>%
      # set 0 point to steady current
      mutate(., offset = mean(I[t > (t.steady - 0.01) &
                                  t < (t.steady - 0.001) & sweep == 24])) %>%
      # normalise to instantanous current
      # this works for deactivating timecourses!
      #filter(., Vpulse > 0) %>%
      group_by(., sweep) %>%
      mutate(
        .,
        Isubs = I - offset,
        norm = mean(Isubs[t > (t.post + 0.001) &
                            t < (t.post + 0.002)]),
        I.tailnorm = Isubs / norm ,
        t.fit = c(t - t.post + 0.001)
      ) %>%
      nest(.key = "tc.tail") 
    
    
    
    #run smoothing on each sweep
    tc.pre %<>%
      mutate(., SGsmooth = purrr::map(tc.pre, ~ sgolayfilt(.$I, p = 3, n = 101, m = 0, ts = 1)))
    
    tc.test %<>%
      mutate(., SGsmooth = purrr::map(tc.test, ~ sgolayfilt(.$I, p = 3, n = 101, m = 0, ts = 1)))
    
    tc.tail %<>%
      mutate(., SGsmooth = purrr::map(tc.tail, ~ sgolayfilt(.$I, p = 3, n = 101, m = 0, ts = 1)))
    
  # get the timepoints
  timepoints =
   unnest(tc.pre, tc.pre, SGsmooth) %>%
     #join the threshold value
     left_join(., stepdata.persweep %>% select(., sweep, thresh)) %>%
    # one row needs to be added because we lose the last point in the lagged time series
     dplyr::mutate(., test.var = c( diff(sign(SGsmooth - thresh)), NA)) %>% 
     #print(.)
     group_by(., sweep) %>%
     dplyr::filter(., test.var  != 0) %>%
     dplyr::summarise(., tapp50.pre = min(t.fit))
  
  timepoints.test =
    unnest(tc.test, tc.test, SGsmooth) %>%
    #join the threshold value
    left_join(., stepdata.persweep %>% select(., sweep, thresh.test)) %>%
    # one row needs to be added because we lose the last point in the lagged time series
    dplyr::mutate(., test.var = c( diff(sign(SGsmooth - thresh.test)), NA)) %>% 
    #print(.)
    group_by(., sweep) %>%
    dplyr::filter(., test.var  != 0) %>%
    dplyr::summarise(., tapp50.test = min(t.fit))
  
  timepoints.tail =
    unnest(tc.tail, tc.tail, SGsmooth) %>%
    #join the threshold value
    left_join(., stepdata.persweep %>% select(., sweep, thresh.tail)) %>%
    dplyr::mutate(., test.var = c( diff(sign(SGsmooth - thresh.tail)), NA)) %>%
    #print(.)
    group_by(., sweep) %>%
    dplyr::filter(., test.var  != 0) %>%
    dplyr::summarise(., tapp50.tail = min(t.fit))
  
  
  
  stepdata.persweep %<>%
    left_join(., timepoints)
  
  stepdata.persweep %<>%
    left_join(., timepoints.test)
  
  stepdata.persweep %<>%
    left_join(., timepoints.tail)
  
  print(stepdata.persweep)

    plot$smooth =
    ggplot(unnest(tc.pre, tc.pre, SGsmooth), aes(x = t.fit, y = I)) +
      geom_line(aes(y = I), colour = "grey50") +
      geom_line(aes(y = SGsmooth), colour = "red") +
      geom_hline(aes(yintercept = thresh), data = stepdata.persweep, colour = "cyan") +
      geom_hline(aes(yintercept = I.inst.pre), data = stepdata.persweep) +
      geom_hline(aes(yintercept = I.ss.pre), data = stepdata.persweep) +
      geom_vline(aes(xintercept = tapp50.pre), data= stepdata.persweep) +
      facet_wrap( ~ Vpulse, ncol = 6,
                  scales = "free",
                  drop = T)
    
    
    plot$smooth.test =
      ggplot(unnest(tc.test, tc.test, SGsmooth), aes(x = t.fit, y = I)) +
      geom_line(aes(y = I), colour = "grey50") +
      geom_line(aes(y = SGsmooth), colour = "red") +
      geom_hline(aes(yintercept = thresh.test), data = stepdata.persweep, color = "cyan") +
      geom_hline(aes(yintercept = I.ss), data = stepdata.persweep) +
      geom_hline(aes(yintercept = I.inst), data = stepdata.persweep) +
      geom_vline(aes(xintercept = tapp50.test), data= stepdata.persweep) +
      facet_wrap( ~ Vpulse, ncol = 6,
                  scales = "free",
                  drop = T)
    
    plot$smooth.tail =
      ggplot(unnest(tc.tail, tc.tail, SGsmooth), aes(x = t.fit, y = I)) +
      geom_line(aes(y = I), colour = "grey50") +
      geom_line(aes(y = SGsmooth), colour = "red") +
      geom_hline(aes(yintercept = thresh.tail), data = stepdata.persweep) +
      geom_hline(aes(yintercept = I.tailSSfull), data = stepdata.persweep) +
      geom_hline(aes(yintercept = I.tail), data = stepdata.persweep) +
      geom_vline(aes(xintercept = tapp50.tail), data= stepdata.persweep) +
      facet_wrap( ~ Vpulse, ncol = 6,
                  scales = "free",
                  drop = T)
    
    
    ggsave(plot$smooth, filename = paste0(txt, "pre-SGsmooth.png"), width = 21, height = 14, path = "./tapp deact prot/")
    ggsave(plot$smooth.test, filename = paste0(txt, "test-SGsmooth.png"), width = 21, height = 14, path = "./tapp deact prot/")
    ggsave(plot$smooth.tail, filename = paste0(txt, "tail-SGsmooth.png"), width = 21, height = 14, path = "./tapp deact prot/")
    
    #write the data to file
    # function to add the meta data 
    addmeta = function(df){
      # copy the meta data for each row
      var =
        ldply(txtmetadata, .fun = rep.int, times = nrow(df)) %>%
        t(.) %>%
        as_tibble(.)
      
      #name the columns correctly
      names(var) = cond
      #  then bind them together
      bind_cols(df, var)
    }
    
    stepdata.persweep %>%
      addmeta(.) %>%
      #select(., id, channel, Pip, sweep, Vpulse, tapp50.pre, I.ss.pre, I.inst.pre, fold.pre) %>%
    write.csv(., file = paste0("./tapp deact prot/" ,str_sub(txt, start = 1, end = -4), "tapp.csv"))
     
  }

# setwd("~/Dphil/data/181214 Isoform selectivity 5K/deact-prot")
# deact.files = list.files(pattern = "txt$")
# 
# dir.create("./tapp deact prot/")
# 
# process.SGsmooth.deact(
#   deact.files[1],
#   npersweep = 39936,
#   sep = "[ ]",
#   cond = c("id", "channel", "Pip", "file")
# )



