# Code for Main Model: Development of orientation preference after eye-opening
# 2020/02/11
## !Update: input noise is applied to the input orientation

rm(list = ls())
setwd(choose.dir()) # set working directory to the folder of the project (~/EmergeProbV1)

source("utils.R")
library(circular)
library(distr)

dt = 1 # ms
set.seed(1234)

# Parallel
library(parallel)
library(foreach)
library(doParallel)

print(paste("# of cores detected", detectCores()))# detect available cores automatically
numCores  = 25 # set number of cores manually, I used 25 cores on a server, you may use less cores 
registerDoParallel(numCores)  # use multicore, set to the number of cores

############################################################
## 1. Input ##
# prior distribution, adapted from Xue-Xin & Stocker 2015
c0 = 1/(2*pi-2) 
# prior density distribution
den<-function(x){
  de<-2-abs(sin(2*x));
  de*c0
}
dist <-AbscontDistribution(d=den, low1 = 0, up1 = pi)
rdist <- r(dist)
# test = rdist(500)
# hist(test, breaks = 12)

# input parameters
repeat_n = 25 # times of repetition
pre_n = 100 # input neurons
nl = pre_n
t_tot = 1500000 # ms
trial_t = 100 # ms time for each input presentation
trial_n = t_tot/trial_t
ff.bg = 10 # back ground firing rate (gaussian function)
ff.A = 40 # peak firing rate (gaussian function)
ff.sigma = 10 # sigma (gaussian function)

############################################################
## 2. Output ##
# output parameters
ne = 80
ni = floor(ne * 0.25)
tausyn = 3        	# tausyn 3
Bsyn = exp((-1./tausyn) *dt) 

taul = 20	 	# taul   25
B = exp(-1/taul)
vth = -50 # threshold voltage (mV)
vleak = -70 # mV, leaky potential
vreset = -70 # mV, reset potential
E_glu =  0.0  # mV AMPA reversal potential
E_gaba = -70.0   # mVGABA reversal potential

kee = ne*0		# degree of connection from E to E: total number * probability
kei = ne*0.5      # kei
kie = ni*0.5      # kie
kii = ni*0.75     # kii
kle = nl       	  # degree of connection from input to E	
kli = nl*0.4    	# kil	
L1_sum = nl
gee = 0           # conductance from E to E
gie = 0.06        # gie  
gei = 0.02        # gei   
gii = 0.02        # gii  
gle = 0.012       # conductance from input to E
gli = 0.006       # gli   

# params for Initial condition
pref.theta = pi # 90 degree
kappa = 0.8 # spread of the von-mises distribution
sigma.min = 30 # min sigma of the initial FF connection
sigma.max = 45 # max sigma of the initial FF connection
L1_sum = nl # sum of synaptic weights, param for weight rescaling 

# params for test orientation selectivity
step = 8
theta.set = (0:(step-1))*pi/step
tepoch = 1000
ttot = tepoch * step

# plot rate for 90 degree input
par(mfrow = c(2,2))
g.x=1:pre_n
ff.mu = pre_n/2
f=ff.bg + ff.A * exp( -pmin((g.x-ff.mu)^2, (g.x-pre_n-ff.mu)^2, (g.x+pre_n-ff.mu)^2)/(2*ff.sigma^2) )
# plot(g.x, f, type = "l", xlab = "input rate")

WeightMat.list = vector("list", repeat_n) # initialize recorder

results_list = foreach(repeat_i = 1:repeat_n, .packages = c("circular"))%dopar%{
  # generate input spikes
  spikel = matrix(0,t_tot, pre_n)
  # generate variable 
  input_pat = rdist(trial_n) 
  # add stimulus noise
  input_pat = input_pat+ rnorm(trial_n,sd = 15/180*pi)
  input_pat[input_pat<0] = input_pat[input_pat<0] +pi
  input_pat[input_pat>pi] = input_pat[input_pat>pi] -pi
  
  for(i in 1:trial_n){
    spikel[(trial_t*(i-1)+1):(trial_t*i), ] = gen_ff_spike(input_pat[i], trial_t)
  }
  image(spikel[1:1000,])
  
  # output spikes
  spikee = matrix(0, t_tot, ne)
  spikei = matrix(0, t_tot, ni)
  ye = matrix(vleak,t_tot, ne) # membrane voltage trace
  yi =  matrix(vleak,t_tot, ni) # membrane voltage trac
  ye0 = runif(ne, min = vleak, max = vth)
  yi0 = runif(ni, min = vleak, max = vth)
  
  # connection matrix
  conle =  RandomConnFunc(ne, nl, k = kle)
  conee = RandomConnFunc(ne, ne, k = kee)*matrix(runif(ne*ne),ne,ne)
  if(ni>0){
    conli = RandomConnFunc(ni, nl, k = kli)
    conei = RandomConnFunc(ni, ne, k = kei)
    conie = RandomConnFunc(ne, ni, k = kie)
    conii = RandomConnFunc(ni, ni, k = kii)
  }
  image(conei)
  
  ############################################################
  ## 4. Plasticity
  # weight matrix: only for plastic weight:
  w_min = 0 # lower bound on connected weights
  w_max = 2 # upper bound on connected weights
  pref.e =as.numeric(rvonmises(ne, circular(pref.theta), kappa))/2
  hist(pref.e, breaks = 12)
  
  sigma.e = runif(ne, min=sigma.min, max = sigma.max)/180*nl
  WeightMat = GaussianWeightFunc(ne, nl, pref.e,sigma.e)

  #WeightMat = matrix(1, nrow = ne, ncol = nl) # row: post neuron, col pre neuron
  WeightMat = WeightMat * conle
  qnfactor = L1_sum / apply(WeightMat, 1, sum) 
  WeightMat= WeightMat * qnfactor
  
  WeightMat.init = WeightMat
  pdf.name = "init dual"
  test_result.init = network_test(WeightMat.init)
  analyze_results.init = analyze_tuning(test_result.init)
  ## recording the preference over time: each column is the preference at time point
  
  pref_record = analyze_results.init$tuninge[1,]
  # start simulation
  # -- plst param
  tm_plst = 34 # tau - (ms)
  tp_plst = 14 # tau + (ms)
  mp_ratio = 1.05 # ratio between LTD and LTP
  
  Bm_plst =  exp((-1./tm_plst) *dt)
  Bp_plst =  exp((-1./tp_plst) *dt)
  
  A_ltp = 5.e-3 #
  A_ltd = mp_ratio*tp_plst*A_ltp/tm_plst
  
  ym_plst<- rep(0,ne)
  yp_plst<- rep(0,pre_n)
  
  par(mfrow = c(3,3))
  pal.1=colorRampPalette(c("blue", "green"), space="rgb")
  
  for(t in 1:(t_tot)){
    ############################################################
    ## generate Output spikes
    if(t == 1){
      
      #browser()
      IEe = gle*(conle*WeightMat)%*%spikel[t,] /tausyn
      IEi = 0
      ye[1,] = t((vleak + B*(ye0-vleak)) + IEe*(E_glu - ye0) + IEi*(E_gaba - ye0))
      IIe = gli*conli%*%spikel[t,]/tausyn
      IIi = 0
      yi[1,] = t((vleak + B*(yi0-vleak)) + IIe*(E_glu - yi0)+ IIi*(E_gaba - yi0))
      
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
    
    # plasticity
    ## Hidden variables
    # When a post-synaptic spike occurs
    ym_plst = Bm_plst*ym_plst*(1-spikee[t,]) + A_ltd*spikee[t,]
    # When a pre-synaptic spike occurs
    yp_plst = Bp_plst*yp_plst*(1-spikel[t,]) + A_ltp*spikel[t,]
    ## Update Weight
    WeightMat = WeightMat + (spikee[t,]%*%t(yp_plst) - ym_plst%*%t(spikel[t,]))
    
    # ## bounds
    WeightMat = WeightMat*conle
    # ## Normalize  multiplicatively, keeping the quadratic sum constant
    qnfactor = L1_sum / apply(WeightMat, 1, sum)
    WeightMat= WeightMat * qnfactor
    # WeightMat= (WeightMat <= w_min) * w_min + (WeightMat > w_min & WeightMat < w_max)* WeightMat + (WeightMat >= w_max)* w_max
    
    # if(t %% 50000 == 0 ){
    #   plot.ts(t(WeightMat),  plot.type = "single", col = pal.1(ne), xlab = "pre_rate", ylab =  "final weight",  xaxt = 'n')
    #   axis(1, at = seq.int(0, pre_n, by = 10), labels = seq.int(0, pre_n, by = 10))
    #   legend("topleft",legend = 1:ne, col = pal.1(ne), lty = 1,cex= 0.5)
    #   # pref recoding
    #   analyze_results.tmp = analyze_tuning(network_test(WeightMat))
    #   pref_record = cbind(pref_record, analyze_results.tmp$tuninge[1,])
    #   
    #   print(t)
    # }
    if( t == t_tot/2){
      # recording during learning
      WeightMat.mid = WeightMat
    }
  }
  # test results
  # mid
  test_result.mid = network_test(WeightMat.mid)
  pdf.name = "mid dual"
  analyze_results.mid = analyze_tuning(test_result.mid)
  # final
  test_result = network_test(WeightMat)
  pdf.name = "final dual"
  analyze_results = analyze_tuning(test_result)
  print(repeat_i)
  # save data for i-th repeat
  result = list(analyze_results.init = analyze_results.init, analyze_results.mid=analyze_results.mid, analyze_results = analyze_results,
                WeightMat.init = WeightMat.init, WeightMat.mid= WeightMat.mid, WeightMat = WeightMat, 
                test_result.init=test_result.init, test_result.mid=test_result.mid, test_result = test_result,
                conle = conle, conli = conli, conee = conee, conei = conei, conie = conie, conii = conii)
  
  #WeightMat.list[[repeat_i]] = WeightMat- WeightMat.init
}
stopImplicitCluster()

save(list  = c( "results_list"), 
     file = paste("data/results_fig24.R", sep = ""))


########## Paper Figures of the model
###################################################################### 
### Figure 3: Emergence of cardinal bias on distribution orientation selectivity in the cortical network 	
###################################################################### 
load("data/results_fig24.R")
par(lwd = 2)
layout.matrix <- matrix(c(1,2, 3, 4,3,5,6,7), nrow = 2, ncol = 4)
layout(mat = layout.matrix) # Widths of the two columns
#layout.show(7)
# prior distribution, adapted from Xue-Xin & Stocker 2015
c0 = 1/(2*pi-2) 
# prior density distribution
den<-function(x){
  de<-2-abs(sin(2*x));
  de*c0
}
x_den = seq(0,pi,by = 0.01)
y_den = den(x_den)
plot(seq(0,180, length.out = length(x_den)), y_den*pi/180, main ="a. Prior", type = "l", xaxt="n",lwd = 2,
     xlab = "Orientation", ylab = "Density",frame.plot = F, ylim = c(0,0.01))
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))


# Distribution of preference, E
hist_n = 9
hist_bin_width = 20
hist_record_init = matrix(0, nrow = hist_n, ncol = repeat_n)
hist_record_mid = matrix(0, nrow = hist_n, ncol = repeat_n)
hist_record_final = matrix(0, nrow = hist_n, ncol = repeat_n)
# init
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # mid
  hist_record_mid[,repeat_i] = hist((results$analyze_results.mid$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}

hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_mid = apply(hist_record_mid, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_init), ylim = c(0,0.3), col = "red", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation", main = "d. Distrbution of PO")
arrows(barCenters, rowMeans(hist_record_init) , barCenters,
       rowMeans(hist_record_init) + hist_sd_init, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

# Sample trace
sample_neuron_id = 2 # 2 15 20 21 34 41 48
test_sep = 100
test_indicator = (seq_along(results_list[[1]]$test_result.init$spikee[,sample_neuron_id])-1) %% 1000 <= (test_sep-1)
plot.ts(1:sum(test_indicator),results_list[[1]]$test_result.init$spikee[test_indicator,sample_neuron_id]*0.5 + 2, type = "l", 
        ylim = c(0,3), col = "red", xlab = "Time (ms)", ylab = "Spike", main = "b. Sample Traces")
lines(1:sum(test_indicator),results_list[[1]]$test_result.mid$spikee[test_indicator,sample_neuron_id]*0.5 + 1, col = "dark green")
lines(1:sum(test_indicator),results_list[[1]]$test_result$spikee[test_indicator,sample_neuron_id]*0.5, col = "blue")
abline(v = test_sep * 1:7, lty = 2)


barplot(rowMeans(hist_record_mid), ylim = c(0,0.3), col = "dark green",width = 1, space  = 0)
arrows(barCenters, rowMeans(hist_record_mid) , barCenters,
       rowMeans(hist_record_mid) + hist_sd_mid, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

barplot(rowMeans(hist_record_final), ylim = c(0,0.3), col = "blue",width = 1, space  = 0)
arrows(barCenters, rowMeans(hist_record_final) , barCenters,
       rowMeans(hist_record_final) + hist_sd_final, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

# Sample tuning curve
xaxis.ang = seq(1,180,1)
plot(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results.init$tuninge[,sample_neuron_id],xaxis.ang), 
     type = "l", col = "red",xlim = c(0,180),lwd = 2, main = "c. Sample Tuning Curve",
     xlab = "Orientation", ylab = "Spikes Rates (Hz)",frame.plot = F, xaxt="n")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
lines(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results.mid$tuninge[,sample_neuron_id],xaxis.ang),lwd = 2, col = "dark green")
lines(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results$tuninge[,sample_neuron_id],xaxis.ang), lwd = 2,col = "blue")

# cumulativ plot of osi
cvare_init = c()
cvare_mid = c()
cvare_final = c()
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  cvare_init = c(cvare_init, results$analyze_results.init$cvare)
  cvare_mid = c(cvare_mid, results$analyze_results.mid$cvare)
  cvare_final = c(cvare_final, results$analyze_results$cvare)
}
plot(ecdf(cvare_init), col = "red", main = "e. Orientation Selectivity (1-CV)", ylab = "Cumulative", xlim = c(0,1))
lines(ecdf(cvare_mid), col = "dark green")
lines(ecdf(cvare_final), col = "blue")
# KS test
ks.test(cvari_init, cvari_mid)
ks.test(cvari_mid, cvari_final)
###################################################################### 
# Figure 4: Emergence of cardinal bias on distribution orientation selectivity in the cortical network 	
###################################################################### 
par(lwd = 2)
layout.matrix <- matrix(c(1,2, 3, 1,4,5,6,7,8), nrow = 3, ncol = 3)
layout(mat = layout.matrix) # Widths of the two columns
layout.show(8)
## Sample trace
sample_neuron_id = 1
test_sep = 100
test_indicator = (seq_along(results_list[[1]]$test_result.init$spikei[,sample_neuron_id])-1) %% 1000 <= (test_sep-1)
plot.ts(1:sum(test_indicator),results_list[[1]]$test_result.init$spikei[test_indicator,sample_neuron_id]*0.5 + 2, 
        type = "l", ylim = c(0,3), col = "red", xlab = "Time (ms)", ylab = "Spike", main = "a. Sample Traces of I")
lines(1:sum(test_indicator),results_list[[1]]$test_result.mid$spikei[test_indicator,sample_neuron_id]*0.5 + 1, col = "dark green")
lines(1:sum(test_indicator),results_list[[1]]$test_result$spikei[test_indicator,sample_neuron_id]*0.5, col = "blue")
abline(v = test_sep * 1:7, lty = 2)

tuning_estimate_normalize <- function(par, stim.ang){
  # calculate estimated tuning curve giving tuning parameters
  pref <- par[1]
  sd <- par[2]
  B <- par[3]
  A <- par[4]
  rhat <- A + B * exp(-0.5 * (ang_ori_diff(stim.ang - pref)/sd)^2) 
  rhat = rhat/max(rhat)
}
## Population Tuning Curve for E and I
# E initialization
pop_response_e.init = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
pop_response_e.mid = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
pop_response_e.final = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
# I initialzation
pop_response_i.init = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
pop_response_i.mid = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
pop_response_i.final = matrix(0, nrow = length(xaxis.ang), ncol = repeat_n)
for(repeat_i in 1:repeat_n){
  pop_response_e.init[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results.init$tuninge, 2, tuning_estimate_normalize, xaxis.ang))
  pop_response_e.mid[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results.mid$tuninge, 2, tuning_estimate_normalize, xaxis.ang))
  pop_response_e.final[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results$tuninge, 2, tuning_estimate_normalize, xaxis.ang))

  pop_response_i.init[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results.init$tuningi, 2, tuning_estimate_normalize, xaxis.ang))
  pop_response_i.mid[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results.mid$tuningi, 2, tuning_estimate_normalize, xaxis.ang))
  pop_response_i.final[,repeat_i] = rowMeans(apply(results_list[[repeat_i]]$analyze_results$tuningi, 2, tuning_estimate_normalize, xaxis.ang))

}
# E pop tuning curve
plot(xaxis.ang, rowMeans(pop_response_e.init), xaxt= "n",type = "l" ,xlim = c(0,180), ylim = c(0,1),
     xlab = "Orientation", ylab = "Normalized Response",frame.plot = F, col = "red", main = "c. E Population Tuning Curve")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
lines(xaxis.ang,  rowMeans(pop_response_e.mid), col = "dark green")
lines(xaxis.ang,  rowMeans(pop_response_e.final), col = "blue")

# Distribution of preference, I
hist_n = 9
hist_bin_width = 20
hist_record_init = matrix(0, nrow = hist_n, ncol = repeat_n)
hist_record_mid = matrix(0, nrow = hist_n, ncol = repeat_n)
hist_record_final = matrix(0, nrow = hist_n, ncol = repeat_n)
# init
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuningi[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ni
  # mid
  hist_record_mid[,repeat_i] = hist((results$analyze_results.mid$tuningi[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ni
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuningi[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ni
}
hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_mid = apply(hist_record_mid, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
# I histogram of preferred orientation init
barCenters = barplot(rowMeans(hist_record_init), ylim = c(0,0.55), col = "red", xaxt = "n",
                     width = 1, space  = 0, xlab = "Preferred Orientation", main = "f. Distrbution of PO")
arrows(barCenters, rowMeans(hist_record_init) , barCenters,
       rowMeans(hist_record_init) + hist_sd_init, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

# I pop tuning curve
plot(xaxis.ang, rowMeans(pop_response_i.init), xaxt= "n",type = "l" ,xlim = c(0,180), ylim = c(0,1),frame.plot = F, col = "red", 
     xlab = "Orientation", ylab = "Normalized Response",main = "d. I Population Tuning Curve")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
lines(xaxis.ang,  rowMeans(pop_response_i.mid), col = "dark green")
lines(xaxis.ang,  rowMeans(pop_response_i.final), col = "blue")

# I histogram of preferred orientation mid
barplot(rowMeans(hist_record_mid), ylim = c(0,0.4), col = "dark green",width = 1, space  = 0)
arrows(barCenters, rowMeans(hist_record_mid) , barCenters,
       rowMeans(hist_record_mid) + hist_sd_mid, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

# Sample tuning curve
xaxis.ang = seq(1,180,1)
plot(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results.init$tuningi[,sample_neuron_id],xaxis.ang), 
     type = "l", col = "red",xlim = c(0,180), ylim = c(0,80),lwd = 2, main = "b. Sample Tuning Curve",
     xlab = "Orientation", ylab = "Spikes Rates (Hz)",frame.plot = F, xaxt="n")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
lines(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results.mid$tuningi[,sample_neuron_id],xaxis.ang),lwd = 2, col = "dark green")
lines(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results$tuningi[,sample_neuron_id],xaxis.ang), lwd = 2,col = "blue")

# cumulativ plot of osi
cvari_init = c()
cvari_mid = c()
cvari_final = c()
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  cvari_init = c(cvari_init, results$analyze_results.init$cvari)
  cvari_mid = c(cvari_mid, results$analyze_results.mid$cvari)
  cvari_final = c(cvari_final, results$analyze_results$cvari)
}
plot(ecdf(cvari_init), col = "red", main = "f. Orientation Selectivity (1-CV)", ylab = "Cumulative", xlab = "gOSI", xlim = c(0,1), do.points = F)
lines(ecdf(cvari_mid), col = "dark green",do.points=F)
lines(ecdf(cvari_final), col = "blue",do.points=F)
# KS test
ks.test(cvari_init, cvari_mid)
ks.test(cvari_mid, cvari_final)

# I histogram of preferred orientation final
barplot(rowMeans(hist_record_final), ylim = c(0,0.4), col = "blue",width = 1, space  = 0)
arrows(barCenters, rowMeans(hist_record_final) , barCenters,
       rowMeans(hist_record_final) + hist_sd_final, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))
