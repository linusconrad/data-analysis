# this script fits an bi-exponential function to timecourses of
# a voltage step protocols with a deactivating prepulse

expfit.tail =
  function(txt, npersweep, sep, cond){
    
    require(signal)       # package for time series data (spectra, filter, etc) 
    require(plyr)         # split apply combine tools
    require(pals)
    require(tidyverse)  
    require(magrittr)     # pipe operators
    require(broom)        # used to extract test and model objects in a sane way
    require(readr)        # read csv files
    require(stringr)      # string operation toolbox
    require(minpack.lm)   # Levenberg Markquart algorithm for curve-fitting
    require(nls.multstart)
    require(nlstools)
    require(ggthemes)     # prettier plots
    require(plotrix)      # has sem function
    require(Cairo)        #
    require(purrr)
    
    require(cowplot)
    
    # adhere to my plot style
    theme_linus = theme_few
    theme_set(theme_linus() + theme(legend.position = "bottom",
                                    legend.title = NULL,
                                    strip.placement = "outside"))
    
    
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
    
    #assign the sweep to the ephysdata
    stepdata$sweep = as.factor(rep(1:nsweeps, each = npersweep))
    
    #transform to physiological convention
    stepdata$I = c(stepdata$I)*-1
    stepdata$V = c(stepdata$V)*-1
    
    ########### Assign these before doing all the math
    # change to use this script with another protocol
    t.pre    = 0.1       # start time of the pre-pulse
    t.test   = 0.6       # start time of the test-pulse  
    t.post   = 1.3       # start time of the final tail step (postpulse)
    t.steady = 1.9       # end of the tail step
    
    # create a variable for the voltage @ every step
    stepdata %<>%
      group_by(., sweep) %>%
      mutate(., Vpulse = round_any(V[t == t.test + 0.2], 1)) %>%  # put a timepoint within the step here to calculate
      ungroup(.)
    
    # create a nested data frame separating the different time courses (pre, test, post)
    # into nested columns for each sweep
    
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
               t < t.post,
             Vpulse < -20) %>% # this is the testpulse data
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
    
    # need to substract the current on the last step (-120 Vpulse)
    # to ensure the fully relaxed current is taken for normalisation
    # (not all currents relax in the 0.6 s timeframe, this is slow)
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
      filter(., Vpulse > 0) %>%
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
    
    # fit a 4 parameter polynomial to the data
    # then extract time-points for 5% 50% and 95% completed relaxation
    
    print("Fitting exponentials to the relaxations")
    
    # make time series to calculate the fitted values from
    t.predict = data.frame(t.fit = seq(from = 0,
                                       to = 0.7,
                                       by = 0.0001))
    tc.tail %<>%
      mutate(
        expfit = purrr::map(tc.tail, ~  nls_multstart(
          I.tailnorm ~ A1 * exp(-exp(lrc1) * t.fit) + A2 * exp(-exp(lrc2) * t.fit) + C,
          iter = 300,
          convergence_count = 100,
          start_lower = c(
            A1 = 0.6,
            A2 = 0.3,
            lrc1 = 4,
            lrc2 = 1.5,
            C = -0.1
          ),
          start_upper = c(
            A1 = 0.7,
            A2 = 0.4,
            lrc1 = 6,
            lrc2 = 2.5,
            C = 0.1
          ),
          control = list(maxiter = 100),
          supp_errors = 'Y',
          data = .
        )),
        # take the sumamry output
        fit.summary = purrr::map(expfit, ~ glance(.)),
        # take the tidy output
        estimates = purrr::map(expfit, ~ tidy(.)),
        confints = purrr::map(expfit, ~ as.tibble(confint2(., level = 0.95), rownames = "term")), 
        # make a predicted variable
        exppredict = purrr::map(expfit, ~ augment(., newdata = t.predict))
      )
    
    #print(tc.tail)
    # tc.test %<>%
    #   mutate(
    #     expfit = purrr::map(tc.test, ~ multstart.biexp(I = I.normdeact, t = t.fit, input = .)),
    #     # take the sumamry output
    #     fit.summary = purrr::map(expfit, ~ glance(.)),
    #     # take the tidy output
    #     estimates = purrr::map(expfit, ~ tidy(.)),
    #     # make a predicted variable
    #     exppredict = purrr::map(expfit, ~ augment(., newdata = t.predict))
    #   )
    
    # write the estimates to file, attach metadata
    
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
  
  dir.create("./files tailfit/")  
    
  print("Writing Fit summaries to file")  
  # write it all to file
    tc.tail %>% select(., sweep, fit.summary) %>%
      addmeta(df = .) %>%
      unnest(., fit.summary) %>%
      write.csv(., file =  paste0("./files tailfit/",txt, "tail-expfit-glance.csv"))
    
    
  # merge estimates and confints together
  persweep.tail = 
     tc.tail %>% select(., sweep, estimates, confints) %>%
     addmeta(df = .) %>%
     unnest(., estimates, confints) %>%
     rename(., conf.lo = '2.5 %', conf.hi = `97.5 %`)
  
  #print(persweep.tail)
  # make separate predictions of the two components to plot alongside 
  write.csv(persweep.tail, file =  paste0("./files tailfit/",txt, "tail-expfit.csv"))
  
  # convert to tau instead of rate constant
  convert.lrc =
    function(lnrc){
      1/exp(lnrc)
    }
  
  # manipulate the estiamte data so I can sort it by component and add the component fits to the graph separately
  
  persweep.bycomp =
    persweep.tail %>%
    select(., sweep, term, estimate) %>%
    filter(., term != "C") %>% # this has no "component" and thus is in the way
    # the following reshapes the data as needed
    mutate(., component = parse_number(term),
           parm = str_sub(term, 1, 1)) %>%
    select(.,-term) %>%
    spread(., key = parm, value = estimate) %>%
    group_by(., sweep) %>%
    mutate(.,
           tau = convert.lrc(l),
           speed = as.factor(rank(tau)) %>% revalue(., c('1' = "fast", '2' = "slow"))) %>%
    # add back the constant
    left_join(
      .,
      persweep.tail %>% select(., sweep, term, estimate) %>%
        filter(., term == "C") %>% spread(., key = term, value = estimate)
    ) %>%
    group_by(., sweep, speed) %>%
    # add in the single component prediction
    mutate(., predict.singleexp = purrr::map(A, ~ data.frame(t = t.predict$t.fit,
                                                          predict.singleexp = . * exp(-exp(l) * t.predict$t.fit) + C))) 
  
  print(persweep.bycomp)
  
  # write this to file
  persweep.bycomp %>%
    select(., - predict.singleexp, -component) %>%
    addmeta(.) %>%
    write.csv(., paste0("./files tailfit/",txt, "tail-expfit-bycomp.csv"))
  
  # write the single fit to a frame for plotting
  componentfits = 
  persweep.bycomp %>%
    select(sweep, speed, predict.singleexp) %>%
    unnest(., predict.singleexp)
   
  persweep.tail.tau =
    persweep.tail %>% 
    filter(., term %in% c("lrc1", "lrc2")) %>%
    mutate(., tau = convert.lrc(estimate),
           conf.lo = convert.lrc(conf.lo),
           conf.hi = convert.lrc(conf.hi)) %>%
    group_by(., sweep) %>%
    # in this format it is random whether lerc1 or 2 is the fast one
    mutate(., component = as.factor(rank(tau)) %>% revalue(., c('1' = "fast", '2' = "slow")))

  # plot the fits
    # make an empty list of plots
    plot = list()
    
    #just the tau
    plot$tauestimate =
      persweep.tail.tau %>%
      ggplot(., aes(x = sweep,
                    y = tau)) +
      facet_wrap( ~ component, scales = "free", ncol = 5) +
      geom_hline(yintercept = 0,
                 colour = "grey30",
                 linetype = 3) +
      geom_linerange(aes(y = tau, ymin = conf.lo, ymax = conf.hi)) +
      geom_point() +
      coord_flip()
    
    # all estimates
    plot$allestimate =
      persweep.bycomp %>%
      select(., sweep, speed, A, tau, C) %>%
      gather(., key = "term", value = "estimate", A, tau, C) %>%
      ggplot(., aes(x = sweep,
                    y = estimate,
                    colour = speed)) +
      facet_wrap( ~ term, scales = "free", ncol = 5) +
      #geom_linerange(aes(y = estimate, ymin = conf.lo, ymax = conf.hi)) +
      geom_point() +
      coord_flip()
    
    #timecourse
    plot$exp.tail =
      ggplot(unnest(tc.tail, tc.tail),
             aes(x = t.fit,
                 y = I.tailnorm)) +
      geom_hline(yintercept = 0,
                 colour = "grey30",
                 linetype = 3) +
      geom_line(colour = "grey50") +
      geom_line(data = unnest(tc.tail, exppredict),
                aes(y = .fitted),
                colour = brewer.dark2(4)[4]) +
      geom_point(data = componentfits,
                aes(y = predict.singleexp,
                    x = t,
                    colour = speed),
                size = 0.05
                ) +
      scale_color_brewer(palette = "Dark2") + 
      scale_y_continuous(limits = c(0,1)) +
      scale_x_continuous(limits = c(0,0.6)) +
      facet_wrap( ~ sweep, ncol = 6)
   
    
     print("Writing fit plots to file")
     try(ggsave(plot$exp.tail, filename = paste0(txt, "tail-expfit-plot.png"), width = 21, height = 14, path = "./files tailfit/"))
     try(ggsave(plot$allestimate, filename = paste0(txt, "tail-exppar-plot.pdf"), width = 6, height = 4, path = "./files tailfit/"))
     try(ggsave(plot$tauestimate, filename = paste0(txt, "tail-tau-plot.pdf"), width = 6, height = 4, path = "./files tailfit/"))
    }



setwd("~/Dphil/data/181214 Isoform selectivity 5K/deact-prot/")
# get the list of files to run on
deact.files = list.files(pattern = "txt$")

#test one to debug
# expfit.tail(deact.files[1],
#                   npersweep = 39936,
#                   sep = "[ ]",
#                   cond = c("id", "channel", "Pip", "file"))

## full run
# l_ply(deact.files,
#       expfit.tail,
#       npersweep = 39936,
#       sep = "[ ]",
#       cond = c("id", "channel", "Pip", "file"))
