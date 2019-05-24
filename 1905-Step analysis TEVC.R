# This script is taking the summaries of the step processing functions and plots the IV curves etc

# make an object for everything of the classic protocol to live in
TEVCstep.obj = list()
TEVCstep.obj$deact_pre = list()

# get the data 
setwd("~/Dphil/data/190512 TECVdata/deact-pre/tapp act prot")

TEVCstep.obj$deact_pre$files =
  list.files(pattern = ".csv$")

TEVCstep.obj$deact_pre$per.sweep =
   ldply(TEVCstep.obj$deact_pre$files, read.csv)

TEVCstep.obj$deact_pre$IVtimepoints =
  TEVCstep.obj$deact_pre$per.sweep %>%
  select(., id, channel, Vpulse, I.0, I.50, I.100, I.150, I.200, I.400, I.600) %>%
  gather(., key = timepoint, value = I, I.0, I.50, I.100, I.150, I.200, I.400, I.600)

# add a timescale, this is how the timepoints map to each other
TEVCstep.obj$deact_pre$timepoints = data.frame(
  timepoint = c("I.0", "I.50", "I.100", "I.150", "I.200", "I.400", "I.600"),
  time = c(3, 50, 100, 150, 200, 400, 600)
)

TEVCstep.obj$deact_pre$IVtimepoints %<>%
  left_join(.,
            TEVCstep.obj$deact_pre$timepoints)

#plot the IVs dep on timepoint
ggplot(TEVCstep.obj$deact_pre$IVtimepoints,
       aes(x = Vpulse, y = I, colour = as.factor(time))) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(xintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  facet_grid(channel ~ time) +
  geom_line(aes(group = id)) 
  #coord_cartesian(ylim = c(-100,100)) + #look at reversal
  #stat_summary(geom = "line")

#fit a high order polynomial to each IV curve
# Do the following:
# 1: fit a polynomial to each IV curves
# 2: solve it for 0, get the crossing point

# make a Voltage vector do do fitting with
V.predict = data.frame(Vpulse = seq(from = -120,
                                   to = 100,
                                   by = 1))

# calculate nernst for 22 C and +1 charge
nernst = # for z = +1
  function(cin, cout) { 
    (8.314 * 293 / 96485) * log((cout / cin)) * 1000
  }

TEVCstep.obj$deact_pre$IVtimepoints %<>%
  left_join(., TEVCstep.obj$deact_pre$IVtimepoints %>%
              group_by(., id, channel, time) %>% 
              dplyr::summarise(., Vcrit = Vpulse[abs(I) == min(abs(I))])
  )

# do the fit
TEVCstep.obj$deact_pre$polyfit =
  TEVCstep.obj$deact_pre$IVtimepoints %>%
  dplyr::mutate(Voi = Vpulse - Vcrit) %>% # this filters for the relevant part of the IV
  dplyr::filter(., abs(Voi) < 30) %>%
  group_by(id, channel, time) %>%
  nest(., .key = IV) %>%
  dplyr::mutate(., IVfit = purrr::map(IV, ~ lm(formula = I ~ stats::poly(Vpulse, 2, raw = TRUE),
                                               data = .)),
                IVfit.summary = purrr::map(IVfit, ~glance(.)),
                # take the tidy output
                IVestimates = purrr::map(IVfit, ~ tidy(.)),
                polypredict = purrr::map(IVfit, ~ augment(x = ., newdata = V.predict)),
                Erev = purrr::map(IVestimates, ~ polyroot(.$estimate))
  )

# extract the correct reversal
TEVCstep.obj$deact_pre$rev =
  TEVCstep.obj$deact_pre$polyfit %>%
  unnest(., Erev) %>% # get rid of the imaginary part of the number
  mutate(., Erev = Re(Erev)) %>%  # left join the critical voltage
  left_join(
    .,
    TEVCstep.obj$deact_pre$IVtimepoints %>%
      group_by(., id, channel, time) %>%
      dplyr::summarise(., Vcrit = Vpulse[abs(I) == min(abs(I))])
  )  %>% # get rid of the solution that are not in the region of the fit
  #dplyr::filter(., Erev )
  dplyr::mutate(., Etest = Erev - Vcrit,
                   bad = abs(Etest) > 10) %>%
  group_by(., id, channel, time) %>%
  dplyr::mutate(., minimum = abs(Etest) == min(abs(Etest))) %>%
  dplyr::filter(., #bad == FALSE,
                minimum == TRUE) %>%
  dplyr::select(., bad, Erev, Vcrit)

# inspect the dubous fits
TEVCstep.obj$deact_pre$rev %>%
  dplyr::filter(., bad == TRUE)

# join to the models etc
TEVCstep.obj$deact_pre$polyfit %<>%
  select(., -Erev) %>%
  left_join(., TEVCstep.obj$deact_pre$rev )

# plot the IVs with fitted reversal to check validity
ggplot(TEVCstep.obj$deact_pre$IVtimepoints  %>% dplyr::filter(., abs(Vpulse-Vcrit) < 30)
       , aes(x = Vpulse,
             y = I,
             colour = time,
             fill = time)) +
  facet_wrap(~ id + channel, scales = "free") +
  scale_color_viridis(begin = 0, end = .8, option  = "B") +
  scale_fill_viridis(begin = 0, end = .8, option  = "B") +
  geom_hline(yintercept = c(0),
             linetype = 3,
             colour  = "grey50") +
  scale_y_continuous(limits = c(-0.5,0.5)) +
  #scale_x_continuous(limits = c(-120,-50)) +
  geom_line(data = unnest(TEVCstep.obj$deact_pre$polyfit, polypredict) #%>% dplyr::filter(., id %in% bad)
            ,
            aes(y = .fitted,
                colour = time,
                group = time)) +
  geom_vline(data = TEVCstep.obj$deact_pre$polyfit #%>% dplyr::filter(., id %in% bad)
             ,
             aes(xintercept = Erev,
                 colour = time))+
  geom_point() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank())

# # make a colour scale to use for the different Pip solutions for consistency
# col.pip = scale_color_brewer(type = "qual", palette = "Set2", name = "Pipette\n Solution")
# fill.pip = scale_fill_brewer(type = "qual", palette = "Set2", name = "Pipette\n Solution")
col.ch = scale_color_brewer(palette = "Accent", name = " ")
fill.ch = scale_fill_brewer(palette = "Accent", name = " ")

# scatter/point specs for consistency
points = geom_point(
  size = 1,
  #alpha = 0.3,
  shape = 21,
  colour = "grey50"
) 


TEVCstep.obj$deact_pre$plot.Erev =
  ggplot(TEVCstep.obj$deact_pre$polyfit #%>% dplyr::filter(., id == "190513cellH")
         ,
         aes(x = time,
                                             y = Erev,
                                             colour = channel)) +
  #facet_grid(id ~ .) +
  labs(title = "Deactivating pre-pulse (-120 mV)",
       x = "Time after jump, (ms)",
       y = expression(paste(E[rev], ", (mV)"))) +
  col.ch +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_hline(yintercept = nernst(90,2),
             linetype = 3,
             colour  = "grey50") +
  scale_y_continuous(limits = c(-109, 0)) +
  geom_line(alpha = 0.2, aes(group = id)) +
  #geom_boxplot(aes(x = as.factor(time))) +
  stat_summary(geom = "pointrange") +
  stat_summary(geom = "line")
#stat_smooth(se = F)

TEVCstep.obj$deact_pre$per.id =
  TEVCstep.obj$deact_pre$per.sweep %>%
  group_by(., id, channel) %>%
  dplyr::summarise(., fold.100 = I.600[sweep == 23] / I.0[sweep == 23],
            prepulse =  "deactivating") # decay of the current on a step from +120 to -120

# add the extrapolated rev pots
TEVCstep.obj$deact_pre$per.id %<>%
  left_join(., 
            TEVCstep.obj$deact_pre$polyfit %>%
              group_by(., id, channel) %>%
             dplyr::summarise(
                .,
                Erev.inst = Erev[time == 3],
                Erev.steady = Erev[time == 600],
                Erev.shift = Erev.inst - Erev.steady, # reversal shift in each patch from inst to ss
                prepulse =  "deactivating"
              )) 
#######################################
# repeat for the other protocol


TEVCstep.obj$act_pre = list()

# get the data 
setwd("~/Dphil/data/190512 TECVdata/act-pre/tapp deact prot/")

TEVCstep.obj$act_pre$files =
  list.files(pattern = ".csv$")

TEVCstep.obj$act_pre$per.sweep =
  ldply(TEVCstep.obj$act_pre$files, read.csv)

TEVCstep.obj$act_pre$IVtimepoints =
  TEVCstep.obj$act_pre$per.sweep %>%
  select(., id, channel, Vpulse, I.0, I.50, I.100, I.150, I.200, I.400, I.600) %>%
  gather(., key = timepoint, value = I, I.0, I.50, I.100, I.150, I.200, I.400, I.600)

# add a timescale, this is how the timepoints map to each other
TEVCstep.obj$act_pre$timepoints = data.frame(
  timepoint = c("I.0", "I.50", "I.100", "I.150", "I.200", "I.400", "I.600"),
  time = c(3, 50, 100, 150, 200, 400, 600)
)

TEVCstep.obj$act_pre$IVtimepoints %<>%
  left_join(.,
            TEVCstep.obj$act_pre$timepoints)

#plot the IVs dep on timepoint
ggplot(TEVCstep.obj$act_pre$IVtimepoints,
       aes(x = Vpulse, y = I, colour = as.factor(time))) +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_vline(xintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  facet_grid(channel ~ time) +
  geom_line(aes(group = id)) 
#coord_cartesian(ylim = c(-100,100)) + #look at reversal
#stat_summary(geom = "line")

#fit a high order polynomial to each IV curve
# Do the following:
# 1: fit a 5 order polynomial to each IV curves
# 2: solve it for 0, get the crossing point

# make a Voltage vector do do fitting with
V.predict = data.frame(Vpulse = seq(from = -120,
                                    to = 100,
                                    by = 1))

TEVCstep.obj$act_pre$IVtimepoints %<>%
  left_join(., TEVCstep.obj$act_pre$IVtimepoints %>%
              group_by(., id, channel, time) %>% 
              dplyr::summarise(., Vcrit = Vpulse[abs(I) == min(abs(I))])
  )

# do the fit
TEVCstep.obj$act_pre$polyfit =
  TEVCstep.obj$act_pre$IVtimepoints %>%
  dplyr::mutate(Voi = Vpulse - Vcrit) %>% # this filters for the relevant part of the IV
  dplyr::filter(., abs(Voi) < 30) %>%
  group_by(id, channel, time) %>%
  nest(., .key = IV) %>%
  dplyr::mutate(., IVfit = purrr::map(IV, ~ lm(formula = I ~ stats::poly(Vpulse, 2, raw = TRUE),
                                               data = .)),
                IVfit.summary = purrr::map(IVfit, ~glance(.)),
                # take the tidy output
                IVestimates = purrr::map(IVfit, ~ tidy(.)),
                polypredict = purrr::map(IVfit, ~ augment(x = ., newdata = V.predict)),
                Erev = purrr::map(IVestimates, ~ polyroot(.$estimate))
  )

# extract the correct reversal
TEVCstep.obj$act_pre$rev =
  TEVCstep.obj$act_pre$polyfit %>%
  unnest(., Erev) %>% # get rid of the imaginary part of the number
  mutate(., Erev = Re(Erev)) %>%  # left join the critical voltage
  left_join(
    .,
    TEVCstep.obj$act_pre$IVtimepoints %>%
      group_by(., id, channel, time) %>%
      dplyr::summarise(., Vcrit = Vpulse[abs(I) == min(abs(I))])
  )  %>% # get rid of the solution that are not in the region of the fit
  #dplyr::filter(., Erev )
  dplyr::mutate(., Etest = Erev - Vcrit,
                bad = abs(Etest) > 10) %>%
  group_by(., id, channel, time) %>%
  dplyr::mutate(., minimum = abs(Etest) == min(abs(Etest))) %>%
  dplyr::filter(., #bad == FALSE,
                minimum == TRUE) %>%
  dplyr::select(., bad, Erev, Vcrit)

# inspect the dubous fits
TEVCstep.obj$act_pre$rev %>%
  dplyr::filter(., bad == TRUE)

# join to the models etc
TEVCstep.obj$act_pre$polyfit %<>%
  select(., -Erev) %>%
  left_join(., TEVCstep.obj$act_pre$rev )

# plot the IVs with fitted reversal to check validity
ggplot(TEVCstep.obj$act_pre$IVtimepoints  %>% dplyr::filter(., abs(Vpulse-Vcrit) < 30)
       , aes(x = Vpulse,
             y = I,
             colour = time,
             fill = time)) +
  facet_wrap(~ id + channel, scales = "free") +
  scale_color_viridis(begin = 0, end = .8, option  = "B") +
  scale_fill_viridis(begin = 0, end = .8, option  = "B") +
  geom_hline(yintercept = c(0),
             linetype = 3,
             colour  = "grey50") +
  scale_y_continuous(limits = c(-0.5,0.5)) +
  #scale_x_continuous(limits = c(-120,-50)) +
  geom_line(data = unnest(TEVCstep.obj$act_pre$polyfit, polypredict) #%>% dplyr::filter(., id %in% bad)
            ,
            aes(y = .fitted,
                colour = time,
                group = time)) +
  geom_vline(data = TEVCstep.obj$act_pre$polyfit #%>% dplyr::filter(., id %in% bad)
             ,
             aes(xintercept = Erev,
                 colour = time))+
  geom_point() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank())

TEVCstep.obj$act_pre$plot.Erev =
  ggplot(TEVCstep.obj$act_pre$polyfit #%>% dplyr::filter(., id == "190513cellH")
         ,
         aes(x = time,
             y = Erev,
             colour = channel)) +
  #facet_grid(id ~ .) +
  labs(title = "Activating pre-pulse (100 mV)",
       x = "Time after jump, (ms)",
       y = expression(paste(E[rev], ", (mV)"))) +
  col.ch +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour  = "grey50") +
  geom_hline(yintercept = nernst(90,2),
             linetype = 3,
             colour  = "grey50") +
  #scale_y_continuous(limits = c(-109, 0)) +
  geom_line(alpha = 0.2, aes(group = id)) +
  #geom_boxplot(aes(x = as.factor(time))) +
  stat_summary(geom = "pointrange") +
  stat_summary(geom = "line")
#stat_smooth(se = F)

TEVCstep.obj$act_pre$per.id =
  TEVCstep.obj$act_pre$per.sweep %>%
  group_by(., id, channel) %>%
  dplyr::summarise(.,
                   steady.m120 = I.600[sweep == 1] / I.0[sweep == 1] * 100, # persistent inward current at -120mV
                   tapp.deact.m120 = tapp50.test[sweep == 1],
                   tapp50.tail.100 = tapp50.tail[sweep == 23],
                   prepulse = "activating") # decay of the current on a step from +120 to -120

# add the extrapolated rev pots
TEVCstep.obj$act_pre$per.id %<>%
  left_join(., 
            TEVCstep.obj$act_pre$polyfit %>%
              group_by(., id, channel) %>%
              dplyr::summarise(
                .,
                Erev.inst = Erev[time == 3],
                Erev.steady = Erev[time == 600],
                Erev.shift = Erev.inst - Erev.steady # reversal shift in each patch from inst to ss
              )) 

# Prepare data sets for plotting
# join the per .id data
# Join the perpatch data together

TEVCstep.obj$measures.per.patch =
  bind_rows(TEVCstep.obj$deact_pre$per.id,
            TEVCstep.obj$act_pre$per.id)

# join the per sweep data together 

TEVCstep.obj$deact_pre$per.sweep %<>%
  dplyr::mutate(., protocol = rep.int("deact. pre p.", length(id)))

TEVCstep.obj$act_pre$per.sweep %<>%
  dplyr::mutate(., protocol = rep.int("act. pre p.", length(id)))

# make a big list for the plots and associated objects
TEVCstep.obj$plots = list()

TEVCstep.obj$plots$gridIVdata =
  bind_rows(TEVCstep.obj$deact_pre$per.sweep,
            TEVCstep.obj$act_pre$per.sweep)

# Continue to assign all the plots in the separate script

# check if some records//files got lost on the way
# all data ids, if its not in here it is irrelevant
records = list()
records$act_pre = tibble(id = TEVCstep.obj$act_pre$per.sweep$id %>% as.character(.) %>% unique,
                         act_pre = rep.int(TRUE, length(id)))
records$deact_pre = tibble(id = TEVCstep.obj$deact_pre$per.sweep$id %>% as.character(.) %>% unique,
                           deact_pre = rep.int(TRUE, length(id)))

records$ramp = tibble(id = ramp.obj$data$id %>% unique(),
                           ramp = rep.int(TRUE, length(id)))
records$dt = tibble(id = deltat.obj$summaries.gathered$id %>% as.character(.) %>% unique(),
                      dt = rep.int(TRUE, length(id)))
  
records$all =
  data.frame(
    id =
      c(
        TEVCstep.obj$act_pre$per.sweep$id %>% as.character(.) %>% unique,
        TEVCstep.obj$deact_pre$per.sweep$id %>% as.character(.) %>% unique,
        ramp.obj$data$id %>% unique(),
        deltat.obj$summaries.gathered$id %>% as.character(.) %>% unique()
      ) %>%
      unique(.)
  )

records$all %<>%
  left_join(., records$act_pre)

records$all %<>%
  left_join(., records$deact_pre)
records$all %<>%
  left_join(., records$ramp)
records$all %<>%
  left_join(., records$dt)


