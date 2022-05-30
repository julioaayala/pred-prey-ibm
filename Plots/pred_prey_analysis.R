# Julio Ayala
# February 2022
# Script to plot results from the predator-prey matlab program
# Folders with csv files need to be located in
# Pred_prey_1morph, Pred_prey_1morph_sexual, Prey
# Outputs pdf files with plots for radiations, final number of species, and population sizes

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpattern)

# Initial variables
tstep <- 500         # Timestep for calculation of clusters/species by intermediate generations
bandwidth <- 0.4     # Kernel bandwidth for kernel density estimator

get_turns <- function(trait_dist, turns, mut_step = 0.0625){
  ## Function to get filtered clusters based on a KDE
  ## Receives a distribution, inflection points, and mutational step.
  ## Returns a new vector of filtered inflection points (i.e. clusters)
  old <- turns
  new <- c()
  maxpeak <- max(trait_dist$y)
  mindensity <- 0.2
  # Filter on density, minimum 0.2*max_density
  while (length(old)!=length(new)){
    if (length(new)>0){
      old <- new
    }
    new <- c()
    for (i in 1:length(old)){
      if (i%%2==1){ # Only eval for maxima
        new <- c(new,old[i]) # Add the peak
        if (trait_dist$y[old[i]]/maxpeak > mindensity){
          new <- c(new,old[i-1])
        }
      }
    }
    new <- sort(new)
  }

  ## Filtering criteria on distance between clusters (Trait value), with 5 mutational
  ## steps of distance as a minimum
  next_iter <- F
  new <- c()
  while (length(old)!=length(new)){
    if (length(new)>0){
      old <- new
    }
    new <- c()    # Filtered indices list
      next_iter <- F
      for (i in 1:length(old)){
        if (!next_iter){
          if (i%%2==1){ # Only eval for maxima and points not set to discard
            if (i+2 <= length(old)){ #Eval with next point
              diff <- trait_dist$x[old[i+2]] - trait_dist$x[old[i]]
              if (diff > (mut_step*5)){
                # Add the peak if the distance is >3*mut
                keep <- c()
                if(i!=1){
                  keep <- old[(i-1):(i+1)]
                }
                else{
                  keep <- old[i:(i+1)]
                }
              }
              else {
                # Get the mean position between the two points, weighted by its density
                meanpos <- as.integer(weighted.mean(c(old[i],old[(i+2)]),c(trait_dist$y[old[i]],trait_dist$y[old[(i+2)]])))
                keep <- c(meanpos)
                if(i!=1){
                  keep <- c(keep,old[(i-1)])
                }
                if(i+3!=length(old)){
                  allnext <- c(old[(i+3):length(old)])
                  keep  <- c(keep, allnext)
                }
                next_iter <- T
              }
            }
            else { # Add last point
              if(i!=1){
                keep <- old[(i-1):i]
              }
              else{
                keep <- old[i]
              }
            }
            new <- c(new,keep)
          }
        }
      }
      new <- sort(unique(new))
  }
  return(new)
}

assigncluster <- function(trait_value, t, distributions, timestep=500) {
  # Function to assign a cluster and return the cluster number, as well as range of the cluster
  # Receives the trait value, generation (t), and array of distributions.
  # Returns a cluster number.
  cluster <- 1
  distindex <- round(t/timestep) +1 # Get the index of the distribution based on

  dist_t <- distributions[[distindex]]
  # mut_step <- 0.0625
  if (length(dist_t)>0){
    for (i in 1:length(dist_t)){ # Find cluster in distribution based on limits
      # if (i%%2==1 & between(trait_value,
      #             (dist_t[i] - 7*mut_step),
      #             (dist_t[i] + 7*mut_step))){ # Assign cluster based on proximity to maxima
      #   cluster <- i
      # }
      if (i>1 & i%%2==0 &  trait_value>dist_t[i]){
        cluster <- cluster+1
      }
    }

    if (length(dist_t)>1){
      cluster <- cluster+length(dist_t)-2
    }
  }


  return(cluster)
}


## Plotting
## Predator-prey 1 morph ---------------------------------
files_pred1 <- list.files('Pred_prey_1morph/', pattern = "\\.csv$")
pdf("Pred_1morph_asexual.pdf", width = 8, height = 6) # Create output file
finished_1morph <- c()
# Dataframe to contain each simulation parameters
filesinfo <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                        attack=double(),g=double(),pmutpred=double(),capred=double(),
                        sexual_reproduction = integer(), fname=character(),
                        populationspred=integer(), morphs = integer(),
                        lastngen=integer(), popsize=integer(), populationsprey=integer())

# Dataframe to contain each simulation parameters (Every 10k gen)
filesinfoextended <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                                attack=double(),g=double(),pmutpred=double(),capred=double(),
                                sexual_reproduction = integer(), fname=character(),
                                populationspred=integer(), morphs = integer(), t=integer(),type=character())


for (i in 1:length(files_pred1)){ # For all files
  filetitle <- files_pred1[i];
  filecontents <- read.csv(paste('Pred_prey_1morph/',files_pred1[i], sep = ""), header = F,
                           sep = '\t', col.names = c('t','alpha','n', 'fitness', 'type'))
  filevalues <- strsplit(files_pred1[i],'_')[[1]]
  filevalues <- filevalues[c(4,6,8,10,12,14,18,22)]
  last <- filecontents %>% filter(type=='pred') %>% filter(t==max(t))
  only_pred <- filecontents %>% filter(type=='pred')
  only_prey <- filecontents %>% filter(type=='prey')
  dists_t_pred <- list() # List of trait distributions at time t
  dists_t_prey <- list() # List of trait distributions at time t
  oldturns_t <- c()
  for (j in c(1,seq(tstep,10000,tstep))){
    if (j<=max(only_pred$t)){
      if ((j%%tstep)==0){
        temp_pred <- only_pred %>% filter(t==j)
        densitykernel_t <- density(temp_pred$alpha, width = bandwidth, weights = temp_pred$n/sum(temp_pred$n))
        deltaY_t <- diff(densitykernel_t$y)
        turns_t <- which(deltaY_t[-1] * deltaY_t[-length(deltaY_t)] < 0) + 1
        turns_t <- get_turns(densitykernel_t,turns_t)
        if (length(turns_t)>5){
          turns_t <- oldturns_t
        }
        npopulations_t <- (length(turns_t)+1)/2

        dists_t_pred[[length(dists_t_pred)+1]] <- densitykernel_t$x[turns_t]
        filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,
                                                           ifelse(as.integer(npopulations_t)==0,1,as.integer(npopulations_t)), 1,
                                                           j,"pred")
        oldturns_t <- turns_t
      }
      else{
        dists_t_pred[[length(dists_t_pred)+1]] <- c(2.0) # Initial populations
      }
    }
    else { # Add extinctions
      filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,0, 1,j,"pred")
    }
  }

  ## Repeat for prey
  for (j in c(1,seq(tstep,10000,tstep))){
    if (j<=max(only_prey$t)){
      if ((j%%tstep)==0){
        temp_prey <- only_prey %>% filter(t==j)
        densitykernel_t <- density(temp_prey$alpha, width = 0.5, weights = temp_prey$n/sum(temp_prey$n))
        deltaY_t <- diff(densitykernel_t$y)
        turns_t <- which(deltaY_t[-1] * deltaY_t[-length(deltaY_t)] < 0) + 1
        turns_t_corrected <- get_turns(densitykernel_t,turns_t,mut_step = 0.125)
        # TODO: temporary fix in split points
        if(length(turns_t_corrected)==0){
          npopulations_t <- (length(turns_t)+1)/2
          dists_t_prey[[length(dists_t_prey)+1]] <- densitykernel_t$x[turns_t]
        }
        else{
          npopulations_t <- (length(turns_t_corrected)+1)/2
          dists_t_prey[[length(dists_t_prey)+1]] <- densitykernel_t$x[turns_t_corrected]
        }

        filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,
                                                           ifelse(as.integer(npopulations_t)==0,1,as.integer(npopulations_t)), 1,
                                                           j,"prey")
      }
      else{
        dists_t_prey[[length(dists_t_prey)+1]] <- c(2.0) # Initial populations
      }
    }
    else { # Add extinctions
      filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,0, 1,j,"prey")
    }
  }

  #Test
  npopulationsprey <- filesinfoextended[nrow(filesinfoextended),]$populationspred
  ##
  allsamples <- rep(last$alpha, last$n)


  densitykernel <- density(last$alpha, width = bandwidth, weights = last$n/sum(last$n))
  deltaY <- diff(densitykernel$y)
  turns <- which(deltaY[-1] * deltaY[-length(deltaY)] < 0) + 1
  turns <- get_turns(densitykernel,turns)
  npopulations <- ceiling((length(turns)+1)/2)
  filesinfo[nrow(filesinfo)+1,] <- c(filevalues, filetitle,
                                     as.integer(npopulations), 1,
                                     filecontents[nrow(filecontents),]$t,
                                     sum(last$n),npopulationsprey)
  plottitle <- bquote(sigma[gamma]~"="~.(filevalues[3])~";"~mu[pred]~"= "~.(formatC(filevalues[6], format = "e", digits = 2))~";"~g~"= "~.(filevalues[5]))
  if (10000 %in% filecontents$t) {
    only_pred$tgroup <- round(only_pred$t/50)*50
    only_prey$tgroup <- round(only_prey$t/50)*50

    # A
    only_pred_filtered <- only_pred %>%
      rowwise() %>%
      mutate(cluster = assigncluster(alpha, t, distributions=dists_t_pred, timestep = tstep))

    # B
    only_prey_filtered <- only_prey %>%
      rowwise() %>%
      mutate(cluster = assigncluster(alpha, t, distributions=dists_t_prey, timestep = tstep))


    only_pred_filtered_1 <- only_pred_filtered %>%
      group_by(tgroup,cluster,type) %>%
      summarise(popsize = sum(n),
                min_trait = min(alpha),
                max_trait = max(alpha),
                mean_alpha = weighted.mean(alpha,n),
                mean_fitness = weighted.mean(fitness,n), .groups = 'drop')

    only_prey_filtered_1 <- only_prey_filtered %>%
      group_by(tgroup,cluster,type) %>%
      summarise(popsize = sum(n),
                min_trait = min(alpha),
                max_trait = max(alpha),
                mean_alpha = weighted.mean(alpha,n),
                mean_fitness = weighted.mean(fitness,n), .groups = 'drop')

    only_pred_filtered_1$cluster <- as.factor(only_pred_filtered_1$cluster)
    only_pred_filtered_1$type <- "pred"
    only_prey_filtered_1$cluster <- as.factor(only_prey_filtered_1$cluster + length(levels(only_pred_filtered_1$cluster)))
    only_prey_filtered_1$type <- "prey"

    grouped_sim <- rbind(only_pred_filtered_1,only_prey_filtered_1)

    pl <- grouped_sim %>%
      ggplot(aes(x=tgroup, y=mean_alpha, color=type, fill=cluster, ymin=min_trait, ymax=max_trait)) +
      geom_path(size=0.7, alpha=1) +
      geom_ribbon(data=only_pred_filtered_1, aes(x=tgroup, y=mean_alpha, fill=cluster, ymin=min_trait, ymax=max_trait),linetype = 0, alpha=0.5) +
      geom_ribbon(data=only_prey_filtered_1, aes(x=tgroup, y=mean_alpha, fill=cluster, ymin=min_trait, ymax=max_trait),linetype = 0, alpha=0.2, color="#00BCF5") +
      scale_fill_manual(values=alpha(c(rep("#F57A00", length(levels(only_pred_filtered_1$cluster))),rep("#00BCF5", 10)),1)) +
      scale_color_manual(values=c("#F57A00","#00BCF5")) +
      labs(title = plottitle, colour="") +
      scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = seq(0,10000,1000)) +
      ylab("Trait value") +
      xlab("Time") +
      theme_bw() +
      theme(legend.position="bottom", legend.text = element_text(size=11)) +
      guides(fill=FALSE, color = guide_legend(override.aes = list(size = 1, alpha=1, fill=c("#F57A00","#00BCF5"))))

    dk <- data.frame(x = densitykernel$x, y = densitykernel$y)

    pl2 <- ggplot(dk, aes(x=x, y=y)) + # KDE plot
      geom_density(stat = "identity", fill=alpha("#F57A00",alpha = 0.3), color="#F57A00") +
      theme_bw() +
      theme(axis.title.x=element_blank(), # Remove all legends and gridlines
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(t = 28, r = 10, b = 68,l = 0)) +
      coord_flip()


    if (ceiling(npopulations)==3){
      lims_1 <- layer_scales(pl)$y$range$range # Extract limits for reescaling
      lims_2 <- layer_scales(pl2)$x$range$range
      pl <- pl + ylim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2])) #Reescale limits
      pl2 <- pl2 + xlim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2]))
      print(filetitle)
      print(grid.arrange(pl,pl2,ncol=2,  widths=c(2.0, 0.25))) # Plot radiation and density
      #print(npopulations)
    }
  }
}
dev.off()

pdf("Pred_1morph_asexual_summary.pdf", width = 8, height = 6) # Create output file
# Summary of number of morphs in the last generation
predpreysummary_asexual <- filesinfo %>%
  mutate(sigmaalpha = as.numeric(sigmaalpha),
         pmutprey = as.numeric(pmutprey),
         pmutpred = as.numeric(pmutpred),
         sigmagamma = as.numeric(sigmagamma),
         attack = as.numeric(attack),
         g = as.numeric(g),
         capred = as.numeric(capred),
         sexual_reproduction = as.integer(sexual_reproduction),
         populationspred = as.integer(populationspred),
         morphs = as.integer(morphs),
         lastngen = as.integer(lastngen),
         popsize = as.integer(popsize),
         populationsprey = as.integer(populationsprey)) %>%
  mutate(populationspred = replace(populationspred, populationspred==0, 1)) %>%
  mutate(populationspred = replace(populationspred, lastngen<10000, 0))

predpreysummary_asexual$extinction <- ifelse(predpreysummary_asexual$populationspred==0,T,F)

# Fixing format of extended info df
filesinfoextended$t <- as.integer(filesinfoextended$t)
filesinfoextended$populationspred <- as.integer(filesinfoextended$populationspred)

finalmorphs <- predpreysummary_asexual %>%
  filter(pmutpred %in% c(1e-5,1e-4,1e-3)) %>%
  ggplot(aes(x=as.factor(g), y=sigmagamma, fill=as.factor(populationspred))) +
  geom_tile() +
  facet_wrap(~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "Pmut: 5e-05",
    "1e-05" = "Pmut: 1e-05",
    "5e-04" = "Pmut: 5e-04",
    "1e-04" = "Pmut: 1e-04",
    "0.005" = "Pmut: 5e-03",
    "0.001" = "Pmut: 1e-03")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  theme_bw() +
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("#FFFFFF","#F57A00"))(max(predpreysummary_asexual$populationspred)+1), na.value="#EEEEEE", name="Morphs")+
  theme(strip.background = element_rect(fill="white"))

finalmorphs5 <- predpreysummary_asexual %>%
  filter(pmutpred %in% c(5e-5,5e-4,5e-3)) %>%
  ggplot(aes(x=as.factor(g), y=sigmagamma, fill=as.factor(populationspred))) +
  geom_tile() +
  facet_wrap(~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "Pmut: 5e-05",
    "1e-05" = "Pmut: 1e-05",
    "5e-04" = "Pmut: 5e-04",
    "1e-04" = "Pmut: 1e-04",
    "0.005" = "Pmut: 5e-03",
    "0.001" = "Pmut: 1e-03")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  theme_bw() +
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("#FFFFFF","#F57A00"))(max(predpreysummary_asexual$populationspred)+1), na.value="#EEEEEE", name="Morphs")+
  theme(strip.background = element_rect(fill="white"))


finalpopsize <- predpreysummary_asexual %>%
  filter(pmutpred %in% c(1e-5,1e-4,1e-3)) %>%
  ggplot(aes(x=as.factor(g), y=sigmagamma, fill=popsize)) +
  geom_tile() +
  facet_wrap(~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "Pmut: 5e-05",
    "1e-05" = "Pmut: 1e-05",
    "5e-04" = "Pmut: 5e-04",
    "1e-04" = "Pmut: 1e-04",
    "0.005" = "Pmut: 5e-03",
    "0.001" = "Pmut: 1e-03")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"))


print(finalmorphs)
print(finalpopsize)
filesinfotmp <- filesinfoextended
filesinfoextended <- filesinfotmp
filesinfoextendedpred <- filesinfoextended %>%
  filter(type=="pred")%>%
  group_by(fname) %>%
  mutate(speciation = populationspred - lag(populationspred, default = 0))

filesinfoextendedpred <- merge(x = filesinfoextendedpred, y = predpreysummary_asexual[ ,c("fname", "extinction", "populationspred")], by = "fname", all.x=TRUE)
names(filesinfoextendedpred)[names(filesinfoextendedpred)=="populationspred.x"] <- "populationspred"


filesinfoextendedpred$pmutpred <- as.numeric(filesinfoextendedpred$pmutpred)
filesinfoextendedpred <- filesinfoextendedpred %>%
  group_by(fname,populationspred) %>%
  mutate(mint=min(t)) %>%
  ungroup()

firstsplit <- filesinfoextendedpred %>% filter(speciation==1 & populationspred==2 & populationspred.y>1)
secondsplit <- filesinfoextendedpred %>% filter(speciation==1 & populationspred==3 & populationspred.y>2)
firstsplit$event <- 'First split'
secondsplit$event <- 'Second split'

splits <- union(firstsplit, secondsplit)
splits$event <- as.factor(splits$event)

splits %>%
  filter(as.numeric(pmutpred) %in% c(1e-04, 0.001, 1e-05) & type=="pred" & extinction==F) %>%
  ggplot(aes(x=as.factor(g), y=sigmagamma, fill=mint)) +
  geom_tile() +
  facet_grid(event ~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "Pmut: 5e-05",
    "1e-05" = "Pmut: 1e-05",
    "5e-04" = "Pmut: 5e-04",
    "1e-04" = "Pmut: 1e-04",
    "0.005" = "Pmut: 5e-03",
    "0.001" = "Pmut: 1e-03",
    "First split" = "First split",
    "Second split" = "Second split")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  scale_fill_stepsn(n.breaks = 10, colours = rev(viridis::viridis(10)), name = "Split generation") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))
# MODIFY COLORS
dev.off()


# Pred 1 morph sexual reproduction ---------------------------------

files_pred1 <- list.files('Pred_prey_1morph_sexual/', pattern = "\\.csv$")
pdf("Pred_1morph_sexual_1.pdf", width = 8, height = 6) # Create output file
finished_1morph <- c()
# Dataframe to contain each simulation parameters
filesinfo <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                        attack=double(),g=double(),pmutpred=double(),capred=double(),
                        sexual_reproduction = integer(), fname=character(),
                        populationspred=integer(), morphs = integer(),
                        lastngen=integer(), popsize=integer(), populationsprey=integer())

# Dataframe to contain each simulation parameters (Every 10k gen)
filesinfoextended <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                                attack=double(),g=double(),pmutpred=double(),capred=double(),
                                sexual_reproduction = integer(), fname=character(),
                                populationspred=integer(), morphs = integer(), t=integer(), type=character())


for (i in 1:length(files_pred1)){
  filetitle <- files_pred1[i];
  filecontents <- read.csv(paste('Pred_prey_1morph_sexual/',files_pred1[i], sep = ""), header = F,
                           sep = '\t', col.names = c('t','alpha','n', 'fitness', 'type'))
  filevalues <- strsplit(files_pred1[i],'_')[[1]]
  filevalues <- filevalues[c(4,6,8,10,12,14,18,22)]
  last <- filecontents %>% filter(type=='pred') %>% filter(t==max(t))
  only_pred <- filecontents %>% filter(type=='pred')
  only_prey <- filecontents %>% filter(type=='prey')
  dists_t_pred <- list() # List of trait distributions at time t
  dists_t_prey <- list() # List of trait distributions at time t
  oldturns_t <- c()
  for (j in c(2,seq(tstep,10000,tstep))){
    if (j<=max(only_pred$t)){
      if ((j%%tstep)==0){
        temp_pred <- only_pred %>% filter(t==j)
        densitykernel_t <- density(temp_pred$alpha, width = bandwidth, weights = temp_pred$n/sum(temp_pred$n))
        deltaY_t <- diff(densitykernel_t$y)
        turns_t <- which(deltaY_t[-1] * deltaY_t[-length(deltaY_t)] < 0) + 1
        turns_t <- get_turns(densitykernel_t,turns_t)
        if (length(turns_t)>5){ # Limit to 3 morphs
          turns_t <- oldturns_t
        }
        npopulations_t <- (length(turns_t)+1)/2
        dists_t_pred[[length(dists_t_pred)+1]] <- densitykernel_t$x[turns_t]
        filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,
                                                           ifelse(as.integer(npopulations_t)==0,1,as.integer(npopulations_t)), 1,
                                                           j,"pred")
        oldturns_t <- turns_t
      }
      else{
        dists_t_pred[[length(dists_t_pred)+1]] <- c(2.0) # Initial populations
      }
    }
    else { # Add extinctions
      filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,0, 1,j,"pred")
    }
  }

  ## Repeat for prey
  for (j in c(1,seq(tstep,10000,tstep))){
    if (j<=max(only_prey$t)){
      if ((j%%tstep)==0){
        temp_prey <- only_prey %>% filter(t==j)
        densitykernel_t <- density(temp_prey$alpha, width = 0.5, weights = temp_prey$n/sum(temp_prey$n))
        deltaY_t <- diff(densitykernel_t$y)
        turns_t <- which(deltaY_t[-1] * deltaY_t[-length(deltaY_t)] < 0) + 1
        turns_t_corrected <- get_turns(densitykernel_t,turns_t,mut_step = 0.125)
        # TODO: temporary fix in split points, meanwhile, add the previous value
        if(length(turns_t_corrected)==0){
          prev <- dists_t_prey[[length(dists_t_prey)]]
          dists_t_prey[[length(dists_t_prey)+1]] <- prev
        }
        else{
          npopulations_t <- (length(turns_t)+1)/2
          dists_t_prey[[length(dists_t_prey)+1]] <- densitykernel_t$x[turns_t]
        }

        filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,
                                                           ifelse(as.integer(npopulations_t)==0,1,as.integer(npopulations_t)), 1,
                                                           j,"prey")
      }
      else{
        dists_t_prey[[length(dists_t_prey)+1]] <- c(2.0) # Initial populations
      }
    }
    else { # Add extinctions
      filesinfoextended[nrow(filesinfoextended)+1,] <- c(filevalues, filetitle,0, 1,j,"prey")
    }
  }


  allsamples <- rep(last$alpha, last$n)


  densitykernel <- density(last$alpha, width = bandwidth, weights = last$n/sum(last$n))
  deltaY <- diff(densitykernel$y)
  turns <- which(deltaY[-1] * deltaY[-length(deltaY)] < 0) + 1
  turns <- get_turns(densitykernel,turns)
  npopulations <- ceiling((length(turns)+1)/2)
  npopulationsprey <- filesinfoextended[nrow(filesinfoextended),]$populationspred
  filesinfo[nrow(filesinfo)+1,] <- c(filevalues, filetitle,
                                     as.integer(npopulations), 1,
                                     filecontents[nrow(filecontents),]$t,
                                     sum(last$n), npopulationsprey)
  plottitle <- bquote(sigma[gamma]~"="~.(filevalues[3])~";"~mu[pred]~"= "~.(formatC(filevalues[6], format = "e", digits = 2))~";"~g~"= "~.(filevalues[5])~choosiness~"= "~.(filevalues[7]))
  if (10000 %in% filecontents$t) {
    only_pred$tgroup <- round(only_pred$t/50)*50
    only_prey$tgroup <- round(only_prey$t/50)*50

    # A
    only_pred_filtered <- only_pred %>%
      rowwise() %>%
      mutate(cluster = assigncluster(alpha, t, distributions=dists_t_pred, timestep = tstep))

    # B
    only_prey_filtered <- only_prey %>%
      rowwise() %>%
      mutate(cluster = assigncluster(alpha, t, distributions=dists_t_prey, timestep = tstep))


    only_pred_filtered_1 <- only_pred_filtered %>%
      group_by(tgroup,cluster,type) %>%
      summarise(popsize = sum(n),
                min_trait = min(alpha),
                max_trait = max(alpha),
                mean_alpha = weighted.mean(alpha,n),
                mean_fitness = weighted.mean(fitness,n), .groups = 'drop')

    only_prey_filtered_1 <- only_prey_filtered %>%
      group_by(tgroup,cluster,type) %>%
      summarise(popsize = sum(n),
                min_trait = min(alpha),
                max_trait = max(alpha),
                mean_alpha = weighted.mean(alpha,n),
                mean_fitness = weighted.mean(fitness,n), .groups = 'drop')

    only_pred_filtered_1$cluster <- as.factor(only_pred_filtered_1$cluster)
    only_pred_filtered_1$type <- "pred"
    only_prey_filtered_1$cluster <- as.factor(only_prey_filtered_1$cluster + length(levels(only_pred_filtered_1$cluster)))
    only_prey_filtered_1$type <- "prey"

    grouped_sim <- rbind(only_pred_filtered_1,only_prey_filtered_1)

    pl <- grouped_sim %>%
      ggplot(aes(x=tgroup, y=mean_alpha, color=type, fill=cluster, ymin=min_trait, ymax=max_trait)) +
      geom_path(size=0.7, alpha=1) +
      geom_ribbon(data=only_pred_filtered_1, aes(x=tgroup, y=mean_alpha, fill=cluster, ymin=min_trait, ymax=max_trait),linetype = 0, alpha=0.5) +
      geom_ribbon(data=only_prey_filtered_1, aes(x=tgroup, y=mean_alpha, fill=cluster, ymin=min_trait, ymax=max_trait),linetype = 0, alpha=0.2, color="#00BCF5") +
      scale_fill_manual(values=alpha(c(rep("#F57A00", length(levels(only_pred_filtered_1$cluster))),rep("#00BCF5", 10)),1)) +
      scale_color_manual(values=c("#F57A00","#00BCF5")) +
      labs(title = plottitle, colour="") +
      scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = seq(0,10000,1000)) +
      ylab("Trait value") +
      xlab("Time") +
      theme_bw() +
      theme(legend.position="bottom", legend.text = element_text(size=11)) +
      guides(fill=FALSE, color = guide_legend(override.aes = list(size = 1, alpha=1, fill=c("#F57A00","#00BCF5"))))

    dk <- data.frame(x = densitykernel$x, y = densitykernel$y)

    pl2 <- ggplot(dk, aes(x=x, y=y)) + # KDE plot
      geom_density(stat = "identity", fill=alpha("#F57A00",alpha = 0.3), color="#F57A00") +
      theme_bw() +
      theme(axis.title.x=element_blank(), # Remove all legends and gridlines
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(t = 28, r = 10, b = 68,l = 0)) +
      coord_flip()


    if (npopulations>3){
      lims_1 <- layer_scales(pl)$y$range$range # Extract limits for reescaling
      lims_2 <- layer_scales(pl2)$x$range$range
      pl <- pl + ylim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2])) #Reescale limits
      pl2 <- pl2 + xlim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2]))
      print(filetitle)
      print(grid.arrange(pl,pl2,ncol=2,  widths=c(2.0, 0.25))) # Plot radiation and density
      #print(npopulations)
    }
  }
}

dev.off()
pdf("Pred_1morph_sexual_summary.pdf", width = 8, height = 6) # Create output file

predpreysummary_sexual <- filesinfo %>%
  mutate(sigmaalpha = as.numeric(sigmaalpha),
         pmutprey = as.numeric(pmutprey),
         pmutpred = as.numeric(pmutpred),
         sigmagamma = as.numeric(sigmagamma),
         attack = as.numeric(attack),
         g = as.numeric(g),
         capred = as.numeric(capred),
         sexual_reproduction = as.integer(sexual_reproduction),
         populationspred = as.integer(populationspred),
         morphs = as.integer(morphs),
         lastngen = as.integer(lastngen),
         popsize = as.integer(popsize)) %>%
  mutate(populationspred = replace(populationspred, populationspred==0, 1)) %>%
  mutate(populationspred = replace(populationspred, lastngen<10000, 0))

filesinfoextended$t <- as.integer(filesinfoextended$t)
filesinfoextended$populationspred <- as.integer(filesinfoextended$populationspred)
predpreysummary_sexual$populationspred[predpreysummary_sexual$populationspred>3] <- 3
predpreysummary_sexual$extinction <- ifelse(predpreysummary_sexual$populationspred==0,T,F)

predpreysummary_asexual$capred[predpreysummary_asexual$capred=="0"] <- "-1"

predpreysummary <- rbind(predpreysummary_sexual,predpreysummary_asexual)

# predpreysummary$populationspred[predpreysummary$populationspred==0] <- NaN

finalmorphs <- predpreysummary %>%
  # filter(pmutpred %in% c(1e-5,1e-4,1e-3)) %>%
  ggplot(aes(x=g, y=sigmagamma, fill=as.factor(populationspred))) +
  geom_tile() +
  geom_tile_pattern(aes(pattern=extinction), colour="black", alpha=0.5) +
  facet_grid(as.numeric(capred) ~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "5e-05",
    "1e-05" = "Low",
    "5e-04" = "5e-04",
    "1e-04" = "Medium",
    "0.005" = "5e-03",
    "0.001" = "High",
    "-1" = "Clonal reproduction",
    "0" = "Random mating",
    "6.9078" = " Low choosiness",
    "44.3614" = "High choosiness")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  theme_bw() +
  # scale_fill_manual(drop=FALSE, values=colorRampPalette(c("#ffffff","#F57A00"))(max(predpreysummary$populationspred)+1), na.value="#FFFFFF", name="Morphs")+
  scale_fill_manual(drop=FALSE, values=c("#FFFFFF","#698CC7","#4D5DA3","#373385"), na.value="#FFFFFF", name="Final\nmorphs")+
  theme(strip.background = element_rect(fill="white")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "", breaks = NULL, labels = NULL)) +
  scale_x_continuous(breaks=c(0.55,0.60, 0.65),sec.axis = sec_axis(~ . , name = "Predator mutation rate", breaks = NULL, labels = NULL)) +
  scale_pattern_manual(values = c('none','none')) +
  guides(pattern=FALSE)

finalmorphs

finalpopsize <- predpreysummary %>%
  # filter(pmutpred %in% c(1e-5,1e-4,1e-3)) %>%
  ggplot(aes(x=g, y=sigmagamma, fill=popsize)) +
  geom_tile() +
  geom_tile_pattern(aes(pattern=extinction), colour="black", alpha=0.5) +
  facet_grid(as.numeric(capred) ~ pmutpred, labeller = as_labeller(c(c(
    "5e-05" = "5e-05",
    "1e-05" = "Low",
    "5e-04" = "5e-04",
    "1e-04" = "Medium",
    "0.005" = "5e-03",
    "0.001" = "High",
    "-1" = "Clonal reproduction",
    "0" = "Random mating",
    "6.9078" = " Low choosiness",
    "44.3614" = "High choosiness")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  theme_bw() +
  # scale_fill_manual(drop=FALSE, values=colorRampPalette(c("#ffffff","#F57A00"))(max(predpreysummary$populationspred)+1), na.value="#FFFFFF", name="Morphs")+
  scale_fill_stepsn(n.breaks = 10, colours = rev(viridis::viridis(10)), name = "Population \nsize") +
  theme(strip.background = element_rect(fill="white")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "", breaks = NULL, labels = NULL)) +
  scale_x_continuous(breaks=c(0.55,0.60, 0.65),sec.axis = sec_axis(~ . , name = "Predator mutation rate", breaks = NULL, labels = NULL)) +
  scale_pattern_manual(values = c('none','none')) +
  guides(pattern=FALSE)

print(finalmorphs)
print(finalpopsize)
dev.off()

# Speciation events
# Filter only predators and find speciation events by the difference of t and t-1
filesinfoextendedpred <- filesinfoextended %>%
  filter(type=="pred")%>%
  group_by(fname) %>%
  mutate(speciation = populationspred - lag(populationspred, default = 0))

filesinfoextendedpred <- merge(x = filesinfoextendedpred, y = predpreysummary_sexual[ ,c("fname", "extinction", "populationspred")], by = "fname", all.x=TRUE)
names(filesinfoextendedpred)[names(filesinfoextendedpred)=="populationspred.x"] <- "populationspred"

# Find the first occurence of the speciation event
filesinfoextendedpred$pmutpred <- as.numeric(filesinfoextendedpred$pmutpred)
filesinfoextendedpred <- filesinfoextendedpred %>%
  group_by(fname,populationspred) %>%
  mutate(mint=min(t)) %>%
  ungroup()

# Get first and second speciation event
firstsplit <- filesinfoextendedpred %>% filter(speciation==1 & populationspred==2 & populationspred.y>1)
secondsplit <- filesinfoextendedpred %>% filter(speciation==1 & populationspred==3 & populationspred.y>2)
firstsplit$event <- 'First split'
secondsplit$event <- 'Second split'

splits <- union(firstsplit, secondsplit)
splits$event <- as.factor(splits$event)

firstsplit %>%
  filter(as.numeric(pmutpred) %in% c(1e-04, 0.001, 1e-05) & type=="pred" & extinction==F) %>%
  ggplot(aes(x=as.factor(g), y=sigmagamma, fill=mint)) +
  geom_tile() +
  facet_grid(capred ~ pmutpred , labeller = as_labeller(c(c(
    "5e-05" = "Pmut: 5e-05",
    "1e-05" = "Pmut: 1e-05",
    "5e-04" = "Pmut: 5e-04",
    "1e-04" = "Pmut: 1e-04",
    "0.005" = "Pmut: 5e-03",
    "0.001" = "Pmut: 1e-03",
    "First split" = "First split",
    "Second split" = "Second split")))) +
  ylab("Predator niche width") +
  xlab("Predator efficiency") +
  scale_fill_stepsn(n.breaks = 10, colours = rev(viridis::viridis(10)), name = "Split generation") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Prey only ---------------------------------
files_pred1 <- list.files('Prey/', pattern = "\\.csv$")
pdf("Prey.pdf", width = 8, height = 6) # Create output file
finished_1morph <- c()
# Dataframe to contain each simulation parameters
filesinfo <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                        attack=double(),g=double(),pmutpred=double(),capred=double(),
                        sexual_reproduction = integer(), fname=character(),
                        populationspred=integer(), morphs = integer(),
                        lastngen=integer(), popsize=integer())

# Dataframe to contain each simulation parameters (Every 10k gen)
filesinfoextended <- data.frame(sigmaalpha=double(),pmutprey=double(),sigmagamma=double(),
                                attack=double(),g=double(),pmutpred=double(),capred=double(),
                                sexual_reproduction = integer(), fname=character(),
                                populationspred=integer(), morphs = integer(), t=integer(), type=character())


for (i in 1:length(files_pred1)){
  filetitle <- files_pred1[i];
  filecontents <- read.csv(paste('Prey/',files_pred1[i], sep = ""), header = F,
                           sep = '\t', col.names = c('t','alpha','n', 'fitness', 'type'))
  filevalues <- strsplit(files_pred1[i],'_')[[1]]
  filevalues <- filevalues[c(3,5)]
  last <- filecontents %>% filter(type=='prey') %>% filter(t==max(t))
  only_prey <- filecontents %>% filter(type=='prey')
  dists_t_prey <- list() # List of trait distributions at time t

  ## Get clusters every tstep generations
  for (j in c(1,seq(tstep,10000,tstep))){
    if (j<=max(only_prey$t)){
      if ((j%%tstep)==0){
        temp_prey <- only_prey %>% filter(t==j)
        densitykernel_t <- density(temp_prey$alpha, width = 0.5, weights = temp_prey$n/sum(temp_prey$n))
        deltaY_t <- diff(densitykernel_t$y)
        turns_t <- which(deltaY_t[-1] * deltaY_t[-length(deltaY_t)] < 0) + 1
        turns_t_corrected <- get_turns(densitykernel_t,turns_t,mut_step = 0.125)
        # TODO: temporary fix in split points, meanwhile, add the previous value
        if(length(turns_t_corrected)==0){
          prev <- dists_t_prey[[length(dists_t_prey)]]
          dists_t_prey[[length(dists_t_prey)+1]] <- prev
        }
        else{
          npopulations_t <- (length(turns_t)+1)/2
          dists_t_prey[[length(dists_t_prey)+1]] <- densitykernel_t$x[turns_t]
        }

      }
      else{
        dists_t_prey[[length(dists_t_prey)+1]] <- c(2.0) # Initial populations
      }
    }
  }

  allsamples <- rep(last$alpha, last$n)

  # Calculate clusters
  densitykernel <- density(last$alpha, width = bandwidth, weights = last$n/sum(last$n))
  deltaY <- diff(densitykernel$y)
  turns <- which(deltaY[-1] * deltaY[-length(deltaY)] < 0) + 1
  turns <- get_turns(densitykernel,turns,mut_step = 0.125)
  npopulations <- ceiling((length(turns)+1)/2)

  plottitle <- bquote(sigma[alpha]~"="~.(filevalues[1])~";"~mu[prey]~"= "~.(formatC(filevalues[2], format = "e", digits = 2)))
  if (10000 %in% filecontents$t) {
    only_prey$tgroup <- round(only_prey$t/50)*50

    only_prey_filtered <- only_prey %>%
      rowwise() %>%
      mutate(cluster = assigncluster(alpha, t, distributions=dists_t_prey, timestep = tstep))

    only_prey_filtered_1 <- only_prey_filtered %>%
      group_by(tgroup,cluster,type) %>%
      summarise(popsize = sum(n),
                min_trait = min(alpha),
                max_trait = max(alpha),
                mean_alpha = weighted.mean(alpha,n),
                mean_fitness = weighted.mean(fitness,n))

    only_prey_filtered_1$cluster <- as.factor(only_prey_filtered_1$cluster)
    only_prey_filtered_1$type <- "prey"

    grouped_sim <- only_prey_filtered_1

    pl <- grouped_sim %>%
      ggplot(aes(x=tgroup, y=mean_alpha, color=type, fill=cluster, ymin=min_trait, ymax=max_trait)) +
      geom_path(size=0.7, alpha=1) +
      geom_ribbon(data=only_prey_filtered_1, aes(x=tgroup, y=mean_alpha, fill=cluster, ymin=min_trait, ymax=max_trait),linetype = 0, alpha=0.5, color="#00BCF5") +
      scale_fill_manual(values=alpha(c(rep("#00BCF5", 10)),1)) +
      scale_color_manual(values=c("#00BCF5")) +
      labs(title = plottitle, colour="") +
      scale_x_continuous(labels = function(x) format(x, scientific = TRUE), breaks = seq(0,10000,1000)) +
      ylab("Trait value") +
      xlab("Time") +
      theme_bw() +
      theme(legend.position="none") +
      guides(fill=FALSE, color = guide_legend(override.aes = list(size = 1, alpha=1, fill=c("#00BCF5"))))

    dk <- data.frame(x = densitykernel$x, y = densitykernel$y)

    pl2 <- ggplot(dk, aes(x=x, y=y)) + # KDE plot
      geom_density(stat = "identity", fill=alpha("#00BCF5",alpha = 0.3), color="#00BCF5") +
      theme_bw() +
      theme(axis.title.x=element_blank(), # Remove all legends and gridlines
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(t = 28, r = 10, b = 68,l = 0)) +
      coord_flip()


    if (npopulations>0){
      lims_1 <- layer_scales(pl)$y$range$range # Extract limits for reescaling
      lims_2 <- layer_scales(pl2)$x$range$range
      pl <- pl + ylim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2])) #Reescale limits
      pl2 <- pl2 + xlim(min(lims_1[1],lims_2[1]),max(lims_1[2],lims_2[2]))
      print(filetitle)
      print(grid.arrange(pl,pl2,ncol=2,  widths=c(2.0, 0.25))) # Plot radiation and density
    }
  }
}
dev.off()
