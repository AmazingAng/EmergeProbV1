# useful functions

## 1. Generate Spike trains
# generate input and output spike train
poisson_gen <- function(fr, dt = 1,t_tot = 10000){
  # Poisson Spike Generator
  spike_train <- rep(0, t_tot)
  # Producing the spike train.
  rand_number <- runif(t_tot)
  # firing rate should adjust from s to ms
  spike_train[rand_number <= fr/1000*dt] <- 1
  return(spike_train)
}

# generate poisson spike train given input firing rates
gen_poisson_mat <- function(pat, ttot){
  spikel = matrix(0,ttot, pre_n)
  for(t in 1:ttot){
    f=pat+5
    # Poisson Spike Generator
    lgn.spike <- as.numeric(runif(pre_n) <= f/1000*dt)
    #image(matrix(lgn.spike, 30,30), useRaster = T)
    #Sys.sleep(0.1)
    spikel[t, ] = lgn.spike
  }
  return(spikel)
}

# generate spikes from gaussian firing rate given input place
gen_ff_spike <- function(theta, ttot){
  # generate input spike given theta and ttot
  ff.mu = theta /pi*pre_n +1
  spikel = matrix(0,ttot, pre_n)
  g.x = 1:pre_n
  #	for drifting grating eq.5
  for(t in 1:ttot){
    f=ff.bg + ff.A * exp( -pmin((g.x-ff.mu)^2, (g.x-pre_n-ff.mu)^2, (g.x+pre_n-ff.mu)^2)/(2*ff.sigma^2) )
    f[f<0] = 0
    # Poisson Spike Generator
    lgn.spike <- as.numeric(runif(pre_n) <= f/1000*dt)
    #image(matrix(lgn.spike, 30,30))
    #Sys.sleep(0.1)
    spikel[t, ] = lgn.spike
  }
  return(spikel)
}

## 2. generate connection function between neurons
RandomConnFunc <- function(n,  n2 = NULL, k = 1, self = F){
  # Calculates connectivity matrix with uniform probability from population n2 to (->) population n
  # ARGS:
  # k: expected degree
  # RETURNS:
  # conmat =  connectivity matrix
  # width of the square
  if(is.null(n2)){
    n2 = n
  }
  nsub = n*n2
  # 2-d gaussian
  p=k/n2;#generate gaussian
  if(self == T){
    diag(p) <- 0
  }
  # generate r.v.
  pV= matrix(runif(nsub), nrow = n, ncol = n2)
  conmat = (p > pV)
  return(conmat)
}
GaussianWeightFunc <- function(n, n2 = NULL,pref.e = NULL,sigma.e = NULL, noise = 0.1, g.max = 1){
  # Calculates weight matrix with gaussian probability from population n2 to (->) population n
  # ARGS:
  # k: expected degree
  # RETURNS:
  # conmat =  weight matrix
  
  # width of the square
  if(is.null(n2)){
    n2 = n
  }
  if(is.null(pref.e)){
    # initialize pref randomly if it is undefined
    pref.e = runif(ne, 0,pi)
  }
  if(is.null(sigma.e)){
    # initialize pref randomly if it is undefined
    sigma.e = rep(ne/5, ne)
  }
  g.x = 1:n2
  g.noise = noise
  WeightMat = matrix(0, nrow = n, ncol = n2)
  for(i in 1:n){
    # 1-d gaussian
    g.center = pref.e[i]/pi*nl
    g.sigma = sigma.e[i]
    g.weight <- g.noise + g.max * exp( -pmin((g.x-g.center)^2, (g.x-nl-g.center)^2, (g.x+nl-g.center)^2)/(2*g.sigma^2))
    WeightMat[i,] = g.weight
  }
  
  # renormalize the mean to 1
  #g.mean = apply(MARGIN = 1,WeightMat, mean)
  #WeightMat = sweep(WeightMat, 1, g.mean, FUN = "/")
  
  return(WeightMat)
}

## Function to test tuning curve
network_test <- function(WeightMat){
  ## input variable:
  # WeightMat: Weight Matrix
  # tepoch   : Represenation time for each stimuli
  # ttot     : total test time
  
  ## 1. Generate Test Input
  spikel = matrix(0,tepoch*step, nl)
  # test input 8 orientations from 0 to 180 (each 22.5 degree)
  for(s in 1:step){
    spikel[(tepoch*(s-1)+1):(tepoch*s), ] = gen_ff_spike(theta.set[s], tepoch)
    
  }
  ye = matrix(0,ttot, ne) # membrane voltage trace
  yi =  matrix(0,ttot, ni) # membrane voltage trace
  spikee = matrix(0,ttot, ne) # spikings
  spikei = matrix(0,ttot, ni) # spikings
  B = exp(-1/taul)
  
  # recordings
  IEe.record =  matrix(0,ttot, ne)
  IEi.record =  matrix(0,ttot, ne)
  IEb.record =  matrix(0,ttot, ne)
  
  IIe.record =  matrix(0,ttot, ni)
  IIi.record =  matrix(0,ttot, ni)
  IIb.record =  matrix(0,ttot, ni)
  
  # start simulation
  
  for(t in 1:ttot){
    if(t %% tepoch == 1){
      # reset membrane potential for different input stimuli
      ye0 = runif(ne, min = vleak, max = vth)
      yi0 = runif(ni, min = vleak, max = vth)
      #browser()
      IEe = gle*(conle*WeightMat)%*%spikel[t,] /tausyn 
      IEi = 0
      ye[t,] = t((vleak + B*(ye0-vleak)) + IEe*(E_glu - ye0) + IEi*(E_gaba - ye0))
      IIe = gli*conli%*%spikel[t,]/tausyn 
      IIi = 0
      yi[t,] = t((vleak + B*(yi0-vleak)) + IIe * (E_glu - yi0) + IIi* (E_gaba - yi0))
    }else{
      IEe = Bsyn*IEe+(gle*(conle*WeightMat)%*%spikel[t,] + gee*conee%*%spikee[t-1,])/tausyn 
      IEi = Bsyn*IEi + gie*conie%*%spikei[t-1,]/tausyn 
      ye[t,] = (vleak + B*(ye[t-1,]-vleak)) + t(IEe* (E_glu - ye[t-1,]) + IEi* (E_gaba - ye[t-1,]) )
      IIe = Bsyn*IIe + (gli*conli%*%spikel[t,] + gei*conei%*%spikee[t-1,])/tausyn
      IIi = Bsyn*IIi + gii*conii%*%spikei[t-1,]/tausyn
      yi[t,] = (vleak + B*(yi[t-1,]-vleak)) + t(IIe* (E_glu - yi[t-1,]) + IIi* (E_gaba - yi[t-1,]))
    }
    spikee[t,] = (ye[t,] > vth)
    spikei[t,] = (yi[t,] > vth)
    ye[t,] = ye[t,]*(ye[t,] <= vth) + vreset*(ye[t,] > vth)
    yi[t,] = yi[t,]*(yi[t,] <= vth) + vreset*(yi[t,] > vth)
    # record
    if(t != 1){
      IEe.record[t,] = t(IEe* (E_glu - ye[t-1,]))
      IEi.record[t,] = t(IEi* (E_gaba - ye[t-1,]) )
      IIe.record[t,] = t(IIe* (E_glu - yi[t-1,]))
      IIi.record[t,] = t(IIi* (E_gaba - yi[t-1,]))
    }
  }
  return(list(spikee = spikee, spikei = spikei, ye=ye, yi=yi, IEe.record=IEe.record, IEi.record= IEi.record,
              IIe.record=IIe.record, IIi.record=IIi.record))
}

network_test_noise <- function(WeightMat){
  ## add noise to the input when measure the tuning curve
  ## input variable:
  # WeightMat: Weight Matrix
  # tepoch   : Represenation time for each stimuli
  # ttot     : total test time
  
  ## 1. Generate Test Input
  spikel = matrix(0,tepoch*step, nl)
  
  tepoch_noise = tepoch/10 # resample the orientation noise for each time interval
  # test input 8 orientations from 0 to 180 (each 22.5 degree)
  for(s in 1:step){
    for(s_noise in 1:10){
      spikel[(tepoch*(s-1)+1+tepoch_noise*(s_noise-1)):(tepoch*(s-1) +tepoch_noise*s_noise), ] = gen_ff_spike(theta.set[s]+rnorm(1, sd = 15/180*pi), tepoch_noise)
    }
  }
  ye = matrix(0,ttot, ne) # membrane voltage trace
  yi =  matrix(0,ttot, ni) # membrane voltage trace
  spikee = matrix(0,ttot, ne) # spikings
  spikei = matrix(0,ttot, ni) # spikings
  B = exp(-1/taul)
  
  # recordings
  IEe.record =  matrix(0,ttot, ne)
  IEi.record =  matrix(0,ttot, ne)
  IEb.record =  matrix(0,ttot, ne)
  
  IIe.record =  matrix(0,ttot, ni)
  IIi.record =  matrix(0,ttot, ni)
  IIb.record =  matrix(0,ttot, ni)
  
  # start simulation
  
  for(t in 1:ttot){
    if(t %% tepoch == 1){
      # reset membrane potential for different input stimuli
      ye0 = runif(ne, min = vleak, max = vth)
      yi0 = runif(ni, min = vleak, max = vth)
      #browser()
      IEe = gle*(conle*WeightMat)%*%spikel[t,] /tausyn 
      IEi = 0
      ye[t,] = t((vleak + B*(ye0-vleak)) + IEe*(E_glu - ye0) + IEi*(E_gaba - ye0))
      IIe = gli*conli%*%spikel[t,]/tausyn 
      IIi = 0
      yi[t,] = t((vleak + B*(yi0-vleak)) + IIe * (E_glu - yi0) + IIi* (E_gaba - yi0))
    }else{
      IEe = Bsyn*IEe+(gle*(conle*WeightMat)%*%spikel[t,] + gee*conee%*%spikee[t-1,])/tausyn 
      IEi = Bsyn*IEi + gie*conie%*%spikei[t-1,]/tausyn 
      ye[t,] = (vleak + B*(ye[t-1,]-vleak)) + t(IEe* (E_glu - ye[t-1,]) + IEi* (E_gaba - ye[t-1,]) )
      IIe = Bsyn*IIe + (gli*conli%*%spikel[t,] + gei*conei%*%spikee[t-1,])/tausyn
      IIi = Bsyn*IIi + gii*conii%*%spikei[t-1,]/tausyn
      yi[t,] = (vleak + B*(yi[t-1,]-vleak)) + t(IIe* (E_glu - yi[t-1,]) + IIi* (E_gaba - yi[t-1,]))
    }
    spikee[t,] = (ye[t,] > vth)
    spikei[t,] = (yi[t,] > vth)
    ye[t,] = ye[t,]*(ye[t,] <= vth) + vreset*(ye[t,] > vth)
    yi[t,] = yi[t,]*(yi[t,] <= vth) + vreset*(yi[t,] > vth)
    # record
    if(t != 1){
      IEe.record[t,] = t(IEe* (E_glu - ye[t-1,]))
      IEi.record[t,] = t(IEi* (E_gaba - ye[t-1,]) )
      IIe.record[t,] = t(IIe* (E_glu - yi[t-1,]))
      IIi.record[t,] = t(IIi* (E_gaba - yi[t-1,]))
    }
  }
  return(list(spikee = spikee, spikei = spikei, ye=ye, yi=yi, IEe.record=IEe.record, IEi.record= IEi.record,
              IIe.record=IIe.record, IIi.record=IIi.record))
}
## function for analyze the tuning curve
calc_fr <- function(spikes, test.t){
  # calclate firing rate from spike train
  # spikes is a column vector of spike trains from a neuron, test.t is the time for a single trial (epoch, in ms)
  unname(tapply(spikes, (seq_along(spikes)-1) %/% test.t, sum))/(test.t/1000)
}

ori_tuning_loss <- function(par, observed,stim.ang){
  pref <- par[1]
  sd <- par[2]
  B <- par[3]
  A <- par[4]
  rhat <- A + B * exp(-0.5 * (ang_ori_diff(stim.ang - pref)/sd)^2) 
  return(sum((observed - rhat)^2))
}

ang_ori_diff <-function(x){
  ang <- pmin(abs(x), abs(x-180), abs(x+180))
  return(ang)
}

ang_ori <- function(x){
  ang <- x%%180
  return(ang)
}

get_tuning_single <- function(fr.single, theta.set, step, type = "ori"){
  # get tuning function by fitting a gaussian function for firing rates of a gaussian function
  stim.ang <- theta.set/pi*180
  pref.guess <- stim.ang[which.max(fr.single)]
  optim.result <- optim(c(pref.guess, 90, 5,5), ori_tuning_loss, method="L-BFGS-B", lower = c(-180, 1, 0, 0), upper = c(360, 90, 100,100), control=list(factr=1e7), observed = fr.single, stim.ang = stim.ang)
  par = optim.result$par
  
  par[1] <- ang_ori(par[1])
  return(par)
}

tuning_estimate <- function(par, stim.ang){
  # calculate estimated tuning curve giving tuning parameters
  pref <- par[1]
  sd <- par[2]
  B <- par[3]
  A <- par[4]
  rhat <- A + B * exp(-0.5 * (ang_ori_diff(stim.ang - pref)/sd)^2) 
}

calc_osi <- function(par){
  # calculate orientation selectivity index
  rbest = tuning_estimate(par, par[1])
  rnull = tuning_estimate(par,ang_ori(par[1]+90))
  osi = (rbest-rnull)/(rbest+0.0001)
  return(osi)
}

# circular variance (CVAR)
calc_cvar <- function(rates, theta.set){
  osi  = sqrt((rates%*%sin(2*(theta.set)))^2+(rates%*%cos(2*(theta.set)))^2) / (sum(rates)+0.0001);
}

# distribution of CV of ISI
calc_cvofisi <- function(spiketrain){
  # get coefficient of variation 
  isi = diff(which(spiketrain == 1))
  cv = sd(isi)/mean(isi)
}

# Plot tuning of different input
calc_max <- function(spikes, test.t){
  # calclate maximum conductance for each orientation
  # spikes is a column vector of spike trains from a neuron, test.t is the time for a single trial (epoch, in ms)
  unname(tapply(spikes, (seq_along(spikes)-1) %/% test.t, max))
}
calc_mean <- function(spikes, test.t){
  # calclate mean input for each orientation
  # spikes is a column vector of spike trains from a neuron, test.t is the time for a single trial (epoch, in ms)
  unname(tapply(spikes, (seq_along(spikes)-1) %/% test.t, mean))
}

plot_input_tuning <- function(bin.width = 10){
  IEeFroml <-rowsum(spikel%*%t(conle*WeightMat), (1:ttot-1) %/% bin.width)
  IEeFrome <- rowsum(spikee%*%t(conee), (1:ttot-1) %/% bin.width)
  IIeFroml <-rowsum(spikel%*%t(conli), (1:ttot-1) %/% bin.width)
  IIeFrome <-rowsum(spikee%*%t(conei), (1:ttot-1) %/% bin.width)
  
  conducEel <- apply(X = IEeFroml, MARGIN = 2, calc_mean, test.t = tepoch/bin.width) * gle
  conducEee <- apply(X = IEeFrome, MARGIN = 2, calc_mean, test.t = tepoch/bin.width) * gee
  conducIel <- apply(X = IIeFroml, MARGIN = 2, calc_mean, test.t = tepoch/bin.width) * gli
  conducIee <- apply(X = IIeFrome, MARGIN = 2, calc_mean, test.t = tepoch/bin.width) * gei
  conducEe <- apply(X = IEe.record, MARGIN = 2, calc_mean, test.t = tepoch)
  conducEi <- apply(X = IEi.record, MARGIN = 2, calc_mean, test.t = tepoch)
  conducIe <- apply(X = IIe.record, MARGIN = 2, calc_mean, test.t = tepoch)
  conducIi <- apply(X = IIi.record, MARGIN = 2, calc_mean, test.t = tepoch)
  
  par(mfrow = c(5,5),mai = c(0.3,0.3,0.1,0.1))
  for(i in 1:min(ne, 50)){
    ts.plot(cbind(conducEe[,i], conducEi[,i], conducEel[,i], conducEe[,i]+conducEi[,i]), col = c("black", "red", "green", "blue"))
  }
  for(i in 1:min(ni,50)){
    
    ts.plot(cbind(conducIe[,i], conducIi[,i], conducIel[,i], conducIe[,i]+conducIi[,i]), col = c("black", "red", "green", "blue"))
    
  }
}

# distribution of CV of ISI
calc_cvofisi <- function(spiketrain){
  # get coefficient of variation 
  isi = diff(which(spiketrain == 1))
  cv = sd(isi)/mean(isi)
}
## Analyze
analyze_tuning <- function(test){
  spikee = test$spikee
  spikei = test$spikei
  ye= test$ye
  yi = test$yi
  IEe.record = test$IEe.record
  IEi.record = test$IEi.record
  IIe.record = test$IIe.record
  IIi.record = test$IIi.record
  
  fre <- apply(X = spikee, MARGIN = 2, calc_fr, test.t = tepoch)
  fri <- apply(X = spikei, MARGIN = 2, calc_fr, test.t = tepoch)
  #test = get_tuning_single(fre[,1], theta.set, step)
  tuninge <- apply(fre, MARGIN = 2, FUN = get_tuning_single, theta.set = theta.set, step= step)
  tuningi <- apply(fri, MARGIN = 2, FUN = get_tuning_single, theta.set = theta.set, step= step)
  # osi (R(θbest)-R(θbest+90°))/R(θbest)
  osie <- apply(tuninge, MARGIN = 2, FUN = calc_osi)
  osii <- apply(tuningi, MARGIN = 2, FUN = calc_osi)
  par(mfrow = c(2,1))
  hist(osie, main = "Histogram of OSI for E neurons")
  hist(osii, main = "Histogram of OSI for I neurons")
  dev.off()
  # sigma (width parameter)
  par(mfrow = c(2,1))
  hist(tuninge[2,], main = "Histogram of sigma for E neurons")
  hist(tuningi[2,], main = "Histogram of sigma for I neurons")
  # cvar
  cvare <- apply(fre, MARGIN = 2, FUN = calc_cvar, theta.set = theta.set)
  cvari <- apply(fri, MARGIN = 2, FUN = calc_cvar, theta.set = theta.set)
  if(sum(!is.na(cvare)!=0)&sum(!is.na(cvari)!=0)){
    par(mfrow = c(2,1))
    hist(cvare, main = "Histogram of CVAR for E neurons")
    hist(cvari, main = "Histogram of CVAR for I neurons")
    dev.off()
  }
  xaxis.ang = seq(0, 180, by = 22.5)
  # report to pdf
  #setwd("G:/Research/Orientation Development/Code/Rewrite Spiking")
  if(exists("pdf.name")){
    pdf(paste(pdf.name,"tunings",".pdf", sep = ""))
  }else{
    pdf(paste("tunings",format(Sys.time(), "%y%m%d%H%M%S"),".pdf", sep = ""))
  }
  fre.max = apply(fre,MARGIN = 2, max)
  fri.max = apply(fri,MARGIN = 2, max)
  fre.mean = apply(fre,MARGIN = 2, mean)
  fri.mean = apply(fri,MARGIN = 2, mean)
  # preferred distribution of the selective neurons
  par(mfrow = c(2,2))
  hist(tuninge[1, ], main = "E neuron", breaks = xaxis.ang)
  if(sum(cvare>0.1) >0){
    hist(tuninge[1, cvare >0.1 & fre.max>1],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
    legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
  }
  if(sum(cvare>0.15) >0){
    hist(tuninge[1, cvare >0.15 & fre.max>1],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
    legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
  }
  if(sum(cvare>0.2) >0){
    hist(tuninge[1, cvare >0.2 & fre.max>1],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
    legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
  }
  mtext('Preferred Distribution', side = 3, outer = F, line = 2.5, cex = 1.3)
  hist(tuningi[1, ], main = "I neuron")
  
  # Senior and freshman
  if(nfresh > 0){
    hist(tuninge[1, fresh.ind], main = "Fresh neuron", breaks = xaxis.ang)
    if(sum(cvare>0.1) >0){
      hist(tuninge[1, cvare >0.1 & fre.max>1 & 1:ne %in% fresh.ind],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    if(sum(cvare >0.15 & fre.max>1& 1:ne %in% fresh.ind) >0){
      hist(tuninge[1, cvare >0.15 & fre.max>1& 1:ne %in% fresh.ind],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    if(sum(cvare >0.2 & fre.max>1& 1:ne %in% fresh.ind) >0){
      hist(tuninge[1, cvare >0.2 & fre.max>1& 1:ne %in% fresh.ind],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    mtext('Preferred Distribution', side = 3, outer = F, line = 2.5, cex = 1.3)
    
    hist(tuninge[1, -fresh.ind], main = "Senior neuron", breaks = xaxis.ang)
    if(sum(cvare>0.1) >0){
      hist(tuninge[1, cvare >0.1 & fre.max>1 & !(1:ne %in% fresh.ind)],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    if(sum(cvare>0.15) >0){
      hist(tuninge[1, cvare >0.15 & fre.max>1& !(1:ne %in% fresh.ind)],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    if(sum(cvare>0.2) >0){
      hist(tuninge[1, cvare >0.2 & fre.max>1& !(1:ne %in% fresh.ind)],col=rgb(0.5,0.5,0.5,0.5), add= T, breaks = xaxis.ang)
      legend("topright", legend = c("Feshman Neuron", "Senior Neuron"),fill=c(rgb(0.8,0.8,0.8,0.2), rgb(0.5,0.5,0.5,0.7)))
    }
    mtext('Preferred Distribution', side = 3, outer = F, line = 2.5, cex = 1.3)
    
  }

  
  # bias index
  par(mfrow = c(1,1))
  # horizontal bias index
  sel_thres = 0.1
  data = tuninge[1,cvare > sel_thres]
  hbi = sum( data>75 &data <105)/(sum(data <15 | data>165)+0.0001)
  # cardinal bias index
  cbi = (sum( data>75 &data <105)+sum(data <15 | data>165)) / 
    (sum( data>30 &data <60)+sum(data >120 & data<150) + 0.0001)
  barplot(c(hbi, cbi), names.arg =c("hbi", "cbi"), main = "bias index")
  
  # maximum firing rate
  par(mfrow = c(2,2))
  hist(fre.max)
  hist(fri.max)
  hist(fre.mean)
  hist(fri.mean)
  
  if(nfresh >0){
    hist(fre.max[fresh.ind], main =  "fresh")
    hist(fre.max[-fresh.ind], main =  "senior")
  }


  par(mfrow = c(2,1))
  CVe = apply(spikee[6001:8000,], MARGIN = 2, FUN = calc_cvofisi)
  CVi = apply(spikei[6001:8000,], MARGIN = 2, FUN = calc_cvofisi)
  if(sum(!is.na(CVe)!=0)&sum(!is.na(CVi)!=0)){
    hist(CVe)
    hist(CVi)
  }
  
  
  # tuning distributions
  par(mfrow = c(2,1))
  hist(osie, main = "Histogram of OSI for E neurons")
  hist(osii, main = "Histogram of OSI for I neurons")
  # sigma (width parameter)
  par(mfrow = c(2,1))
  hist(tuninge[2,], main = "Histogram of sigma for E neurons")
  hist(tuningi[2,], main = "Histogram of sigma for I neurons")
  # Circular variance
  par(mfrow = c(2,1))
  hist(cvare, main = "Histogram of CVAR for E neurons")
  hist(cvari, main = "Histogram of CVAR for I neurons")
  
  # Plot tuning curve
  par(mfrow = c(5,5),mai = c(0.2,0.2,0.2,0.2))
  xaxis.ang = seq(0, 180, by = 10)
  for(i in 1:ne){
    plot(theta.set/pi*180, fre[,i],xlim = c(0,180))
    lines(xaxis.ang,tuning_estimate(tuninge[,i],xaxis.ang), col = "red")
    
  }
  for(i in 1:ni){
    plot(theta.set/pi*180, fri[,i], xlim = c(0,180))
    lines(xaxis.ang,tuning_estimate(tuningi[,i],xaxis.ang ), col = "red")
    
  }
  dev.off()
  return(list(fre = fre, fri =fri, tuninge = tuninge, tuningi = tuningi, cvare = cvare, cvari = cvari, osie = osie, osii = osii))
}

