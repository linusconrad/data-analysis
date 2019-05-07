read.tIV =
  function(txt, npersweep, sep, cond)
    # txt is the filename with a splitable separator, sep, sep must be in regex format
    # npersweep is the number of rows (samples) in a sweep
    # cond is a character vector of column names of the info in filename
  {
    #first read the information from the filename
    txtmetadata = unlist(strsplit(txt, sep))
    #print(txtmetadata) #for debug purposes
    #then read in the actual data
    ephysdata =
      read_delim(
        txt,
        "\t",
        escape_double = FALSE,
        col_names = c("t", "I", "V"),
        trim_ws = TRUE
      )
    #we need this numbers to assign sweeps correctly
    nsamples = length(ephysdata$t) #number of samples in the file
    #print(nsamples) #this is for debugging only
    nsweeps = nsamples / npersweep # number of sweeps in the file
    
    #then create a df with the categorical variables
    variables =
      ldply(as.list(txtmetadata), .fun = rep.int, times = nsamples) %>%
      t(.) %>%
      as_tibble(.)
    #name the columns correctly
    names(variables) = cond
    #add a column to identify the sweep
    variables$sweep = as.factor(rep(1:nsweeps, each = npersweep))
    #bind the two dfs together and create the output
    bind_cols(as_tibble(ephysdata),
              variables)
  }

setwd("~/Dphil/data/181214 Isoform selectivity 5K/ramps")
ramp.files = list.files(pattern = "txt")

# reads in the metadata from the filenames, chacks whether reading the data can go ahead
filecheck = function(filename){
  strsplit(filename, split = " ") %>%
    unlist
  
}

filedf = ldply(.data = ramp.files, .fun = filecheck)


ramp.data = ldply(
  ramp.files,
  read.tIV,
  npersweep = 29696,
  sep = "[ ]",
  cond = c("id", "channel", "Pip", "file"),
  .inform = T
)

# cut out the actual ramp part of the data
ramp.data %<>% filter(., t > 0.4 & t < 1.4)
#reverse to phys convetion
# flip the inside out data to adhere to E-phys conventions
ramp.data$I = c(ramp.data$I)*-1
ramp.data$V = c(ramp.data$V)*-1 

#clean up the cond vector, use only first bit
ramp.data %<>%
  group_by(.,id, file) %>%
  mutate(., cond = str_split(file, pattern = "[.]")[[1]][1]) %>%
  ungroup

#ramp.data %<>% select(., -V1, -file)



# calculate driving force
nernst = # for z = +1
  function(cin, cout) { 
    (8.314 * 293 / 96485) * log((cout / cin)) * 1000
  }

# make a list for objects relating to the ramp data analysis
ramp.obj = list()

# this is the calculated juntion potentials for the rec.
ramp.obj$Pots = data_frame(
  Pip = c("5K-Na", "5K-NMDG", "5K-Li"),
  #Ejunc = c(-5.9,-11.7, -8.3),  # calc
  Ejunc = c(-6.8, -9.8, -2.9),  # empiric
  #junction pots
  EK = c(nernst(140, 5), nernst(140, 5), nernst(140,5))
)

# merge the offsetpotential to the data 
ramp.data %<>%
  left_join(., ramp.obj$Pots)

# correct for liquid junction
ramp.data %<>%
  mutate(., Vcorr = V + Ejunc)

# save this dataset for future use
#save(ramp.data, file = "IsoformsRamp_5K_rawdata-2018.Rda") #data collected before xmas
saveRDS(ramp.data, file = "IsoformsRamp_5K_rawdata-Mar2019.Rds") #all data


# average the sweeps to one Average sweep for each measurement
ramp.obj$byid =
  ramp.data %>%
  group_by(id, t, channel, Pip , cond) %>%
  summarise(., I = mean(I),
            Vcorr = mean (Vcorr),
            EK = mean(EK)) %>%
  ungroup(.)

ramp.obj$byid %<>%
  group_by(., id) %>%
  mutate(., Inorm = I/mean(I[Vcorr > 99.5 & Vcorr < 100.5 & cond == "baseline"]))

# get rid of raw data to save space for now
rm(ramp.data)
# regenerate if needed
#load(file = "IsoformsRamp_5K_rawdata.Rda")
#ramp.data = readRDS("IsoformsRamp_5K_rawdata-Mar2019.Rds")


ramp.obj$bychannel =
  ramp.obj$byid %>%
  group_by(., channel, Pip, cond, t) %>%
  summarise(
    .,
    I.avg         = mean(I),
    I.avg.norm    = mean(Inorm),
    I.merror      = mean(I) - std.error(I),
    I.perror      = mean(I) + std.error(I),
    I.norm.merror = mean(Inorm) - std.error(Inorm),
    I.norm.perror = mean(Inorm) + std.error(Inorm),
    n             = length(I),
    #V             = mean(V),
    Vcorr         = mean (Vcorr),
    EK            = mean(EK)
  )

# calculate some key measures for each measurement (collapse to 1 number/measurement)
# calculate per-patch measures
ramp.obj$avg =
  ramp.obj$byid %>%
  group_by(., id, channel, Pip, cond) %>%
  summarise(.,
            I.0    = mean(I[Vcorr > -0.2 & Vcorr < 0.2]),       # current @ 0 mV
            MP     = mean(Vcorr[I > -1 & I < 1])  # membrane potential
            #out = mean(I[EK > 74.5 & EK < 75.5]), # current with 95 mv driving force
            #inw = mean(I[EK < -74.5 & EK > -75.5]), #-95 driving force
            #rect = abs(out)/abs(inw)
  ) %>%
  ungroup(.)
#calculate fraction of TPA block (~leak)
 
ramp.obj$avg %<>%
  group_by(., id, Pip) %>%
  mutate(., block.TPA = (1 - I.0/I.0[cond == "baseline"]) *100) # percent block

# make an extra column to indicate which patches have a full TPA dataset
ramp.obj$avg %<>%
  left_join(.,
            ramp.obj$avg %>%
              group_by(id) %>%
              summarise(., block.measured = length(MP) == 2)) #if length is 2 there is a TPA measurement, if 1 there is not

# remove the cond. variable and instead merge the TPA block to the MP etc measures determined from the baseline measurement
# this way we can correlate//scatter MP with the TPA block to asses quality of the recording
ramp.obj$avg =
left_join(
  ramp.obj$avg %>% filter(block.TPA == 0) %>% select(.,-cond,-block.TPA),
  ramp.obj$avg %>% filter(block.TPA != 0) %>% select(., id, block.TPA, MP) %>% rename(., MP.TPA = MP)
)
ramp.obj$avg %<>%
  mutate(MP.shift = MP.TPA - MP)

# Plots

#plot IV for each measurent
ramp.obj$IVmeasurement =
  ggplot(ramp.obj$byid,
         aes(
           x = Vcorr,
           y = I,
           colour = as.factor(cond)
         )) +
  facet_wrap(Pip ~ id, scales = "free") +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(xintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(aes(xintercept = MP,
                 colour = as.factor(cond)),
             data = ramp.obj$avg) +
  scale_color_ptol()+
  geom_line()

ramp.obj$IVmeasurement.norm =
  ggplot(ramp.obj$byid,
         aes(
           x = Vcorr,
           y = Inorm,
           colour = as.factor(cond)
         )) +
  facet_wrap(id + channel ~ Pip) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(xintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  scale_color_ptol()+
  geom_line()

# plot comparison in a grid by pipette solution and construct
# add lines for expected reversal
ramp.obj$IVchannel =
  ggplot(ramp.obj$bychannel %>% filter(., channel == "Trek2b-FL"),
         aes(
           x = Vcorr,
           y = (I.avg / 1000),
           colour = as.factor(cond)
         )) +
  labs(x = expression(paste(V + V[{junc}], " (mV)")),
       y = "Average Current, nA") +
  facet_grid(. ~ Pip) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(
    xintercept = c(0, nernst(140, 5)),
    linetype = 3,
    colour  = "grey50"
  ) +
  scale_color_ptol(name = " ", labels = c("Baseline", "TPA block, 1 mM")) +
  geom_ribbon(aes(x = Vcorr,
                   ymin = I.merror /1000,
                   ymax = I.perror /1000,
                   group = cond,
                   colour = NA
                   ),
               fill = "grey90"
               ) +
  geom_line()

# show the current around 0 to gauge reversal
ramp.obj$IVchannel0 =
  ggplot(ramp.obj$bychannel,
         aes(
           x = Vcorr,
           y = I.avg,
           colour = as.factor(cond)
         )) +
  facet_grid(Pip ~ channel) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(
    xintercept = c(0, nernst(140, 5)),
    linetype = 3,
    colour  = "grey50"
  ) +
  scale_color_ptol() +
  scale_y_continuous(limits = c(-200,500)) +
  geom_ribbon(aes(x = Vcorr,
                  ymin = I.merror,
                  ymax = I.perror,
                  group = cond,
                  colour = NA
  ),
  fill = "grey90"
  ) +
  geom_line() 



###
ramp.obj$IVchannel.norm =
  ggplot(ramp.obj$bychannel %>% filter(., channel == "Trek2b-FL"),
         aes(
           x = Vcorr,
           y = I.avg.norm,
           colour = as.factor(cond)
         )) +
  labs(x = expression(paste(V + V[{junc}], " (mV)")),
       y = expression(paste(I/I[{"100 mV"}]))
       ) +
  facet_grid(. ~ Pip) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  #coord_cartesian(xlim = c(-120, -60)) +
  geom_vline(
    xintercept = c(0, nernst(140, 5)),
    linetype = 3,
    colour  = "grey50"
  ) +
  scale_color_ptol(name = " ", labels = c("Baseline", "TPA block, 1 mM")) +
  geom_ribbon(aes(x = Vcorr,
                  ymin = I.norm.merror,
                  ymax = I.norm.perror,
                  group = cond,
                  colour = NA
  ),
  fill = "grey90"
  ) +
  geom_line() 


ramp.obj$plotMP =
  ggplot(ramp.obj$avg %>% filter(., channel == "Trek2b-FL", block.measured == "TRUE"),
         aes(x = Pip,
             #shape = block.measured,
             y = MP,
             colour = Pip)) +
  labs(#x= ,
    y = expression(paste(E[rev], ", (mV)"))) +
  col.pip +
  fill.pip +
  #geom_jitter(width = 0.1) +
   # geom_violin(alpha = 0.5,
   #             fill = "grey50",
   #             colour = NA) +
  
  geom_boxplot(width = 0.2,
              # colour = "grey30",
               aes(fill = Pip),
              alpha = 0.4,
               outlier.size = 1) +
  facet_wrap( ~ channel) +
  geom_hline(
    yintercept = c(0, nernst(140, 5)),
    linetype = 3,
    colour  = "grey50"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        #aspect.ratio = 2,
        legend.position = "none",
        axis.title.x = element_blank()
        )


ramp.obj$I.blocked =
  ggplot(ramp.obj$avg %>% filter(., channel == "Trek2b-FL",block.measured == "TRUE"),
         aes(x = Pip,
             y = block.TPA,
             colour = Pip)) +
  labs(y = expression(paste(I[blocked], ", %"))) +
  #geom_jitter(width = 0.1) +
  #geom_jitter(width = 0.1) +
  # geom_violin(alpha = 0.5,
  #             fill = "grey50",
  #             colour = NA) +
  col.pip +
  fill.pip +
  geom_boxplot(width = 0.2,
               # colour = "grey30",
               aes(fill = Pip),
               alpha = 0.4,
               outlier.size = 1) +
  facet_wrap(~ channel) +
   geom_hline(yintercept = c(100),
              linetype = 3,
             colour  = "grey50") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none"
        #aspect.ratio = 2
        )

# scatter change in MP vs block
ramp.obj$quality.control =
  ggplot(
    ramp.obj$avg %>% filter(., channel == "Trek2b-FL", block.measured == "TRUE"),
    aes(x = MP,
        y = MP.TPA,
        fill = block.TPA)
  ) +
  labs(x = expression(paste(E[rev], ", (mV)")),
       y = expression(paste(E[rev], ", after block (mV)"))) +
  geom_hline(yintercept = c(0, nernst(140, 5)),
             linetype = 3,
             colour  = "grey50") +
  geom_vline(xintercept = c(0),
             linetype = 3,
             colour  = "grey50") +
  scale_x_continuous(limits = c(nernst(140, 5), 0)) +
  scale_y_continuous(limits = c(nernst(140, 5), 0)) +
  #scale_color_viridis(end = 0.9 , option = "D", name = "Blocked, %") +
  scale_fill_gradientn(colours = rev(brewer.greys(300)), guide = "colourbar", name = "Blocked, %") +

geom_abline(
  slope = 1,
  intercept = 0,
  linetype = 3,
  colour  = "grey50"
) +
  points +
  theme(#aspect.ratio = 1,
    legend.position = "right")



ramp.obj$combined = 
plot_grid(ramp.obj$plotMP,
          ramp.obj$I.blocked,
          ramp.obj$quality.control,
          ncol = 3,
          align = "h",
          axis = "bt",
          rel_widths = c(1,1,2.4)
          )

setwd("~/Dphil/data/181214 Isoform selectivity 5K")

ggsave(ramp.obj$combined, width = 7.3, height = 4, filename = "Trek2-Current Rev pot.pdf")
ggsave(ramp.obj$combined, width = 7.3, height = 4, filename = "Trek2-Current Rev pot.png", type = "cairo", dpi = 600)

ggsave(ramp.obj$IVchannel, width = 7.3, height = 4, filename = "Trek2-Current Rev IV.pdf")
ggsave(ramp.obj$IVchannel, width = 7.3, height = 4, filename = "Trek2-Current Rev IV.png", type = "cairo", dpi = 600)


# what is the n for full pair-wise recordings?
ramp.obj$avg %>%
  filter(., block.measured == TRUE, channel == "Trek2b-FL") %>%
  group_by(., channel, Pip) %>%
  summarise(., n = length(unique(id)))

  
ramp.obj$avg %>%
  filter(., block.measured == TRUE) %>%
  group_by(., channel, Pip) %>%
  mutate(., p = unique(id)) %>%
  arrange(., Pip) %>%
  as.data.frame()
  
  
  
  
