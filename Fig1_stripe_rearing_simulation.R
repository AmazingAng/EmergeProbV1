## Code for Figure 1 stripe rearing
## 2020/02/16

rm(list = ls())
setwd(choose.dir()) # set working directory to the folder of the project (~/EmergeProbV1)
source("utils.R")
library(circular)

set.seed(1234)
dt = 1 # ms

# Parallel
library(parallel)
library(foreach)
library(doParallel)

print(paste("# of cores detected", detectCores()))# detect available cores automatically
numCores  = 25 # set number of cores manually, I used 25 cores on a server, you may use less cores 
registerDoParallel(numCores)  # use multicore, set to the number of cores

############################################################
## 1. Input ##
# prior density distribution (von mises distribution, center at 90 degree)
for(lens_name in c(0,45,90,135)){
  #lens_name = 90# in the code 90 is the horizontal orientation
  lens_ori = pi * lens_name/180 *2 #pi*2/4 *2 # 90
  lens_kappa = 0.5
  test = dvonmises(seq(0,2*pi,by = 0.01),  circular(lens_ori), lens_kappa)
  # plot(seq(0,pi,by = 0.005), test, main ="prior")
  
  # input parameters
  repeat_n = 25 # number of repetition
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
  ne = 80 # number of exictatory neurons
  ni = floor(ne * 0.25) # number of inhibitory neurons
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
  kei = ne*0.5        	# kei
  kie = ni*0.5        	# kie
  kii = ni*0.75        	# kii
  kle = nl       	# kel	25
  kli = nl*0.4    	        # kil	25
  L1_sum = nl
  gee = 0             # conductance from E to E   0.25	0.70    0.2
  gie = 0.06            # gie   1.3     2.6
  gei = 0.02             # gei   1.5     0.75
  gii = 0.02            # gii   2.2     5.5
  gle = 0.012          # gle  0.012
  gli = 0.006           # gli   0.09    0.12
  
  # params for Initial condition
  pref.theta = pi # 90 degree, peak of the von-mises distribution
  kappa = 0.5 # spread of the von-mises distribution
  sigma.min = 30 # min sigma of the initial FF connection
  sigma.max = 45 # max sigma of the initial FF connection
  L1_sum = nl # sum of synaptic weights, param for weight rescaling 
  
  # params for test orientation selectivity
  step = 8
  theta.set = (0:(step-1))*pi/step
  tepoch = 1000
  ttot = tepoch * step
  
  
  # plot rate for 90 degree input
  # par(mfrow = c(2,2))
  # g.x=1:pre_n
  # ff.mu = pre_n/2
  # f=ff.bg + ff.A * exp( -pmin((g.x-ff.mu)^2, (g.x-pre_n-ff.mu)^2, (g.x+pre_n-ff.mu)^2)/(2*ff.sigma^2) )
  # plot(g.x, f, type = "l", xlab = "input rate")
  # # plot of prior
  # test = dvonmises(seq(0,2*pi,by = 0.01),  circular(lens_ori), lens_kappa)
  # plot(seq(0,pi,by = 0.005), test, type ="l",main = "prior", xlab = "orientation", ylab = "density")
  # 
  
  # WeightMat.list = vector("list", repeat_n) # initialize recorder
  
  
  results_list = foreach(repeat_i = 1:repeat_n, .packages = c("circular", "distr"))%dopar%{
    # generate input spikes
    spikel = matrix(0,t_tot, pre_n)
    input_pat = as.numeric(rvonmises(trial_n, circular(lens_ori), lens_kappa))/2
    # add stimulus noise
    input_pat = input_pat+ rnorm(trial_n,sd = 15/180*pi)
    input_pat[input_pat<0] = input_pat[input_pat<0] +pi
    input_pat[input_pat>pi] = input_pat[input_pat>pi] -pi
    for(i in 1:trial_n){
      spikel[(trial_t*(i-1)+1):(trial_t*i), ] = gen_ff_spike(input_pat[i], trial_t)
    }
    
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
    # image(conei)
    
    ############################################################
    ## 4. Plasticity
    # weight matrix: only for plastic weight:
    w_min = 0 # lower bound on connected weights
    w_max = 2 # upper bound on connected weights
    
    # prior distribution, adapted from Xue-Xin & Stocker 2015
    c0 = 1/(2*pi-2) 
    # prior density distribution
    den_lens_init<-function(x){
      de<-2-abs(sin(2*x));
      de*c0
    }
    
    # initial condition
    pref.e = runif(ne, max = pi)
    sigma.e = runif(ne, min=sigma.min, max = sigma.max)/180*nl
    WeightMat = GaussianWeightFunc(ne, nl, pref.e,sigma.e)

    #WeightMat = matrix(1, nrow = ne, ncol = nl) # row: post neuron, col pre neuron
    WeightMat = WeightMat * conle
    qnfactor = L1_sum / apply(WeightMat, 1, sum) 
    WeightMat= WeightMat * qnfactor
    
    WeightMat.init = WeightMat
    test_result.init = network_test(WeightMat.init)
    analyze_results.init = analyze_tuning(test_result.init)
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
      
      # if(t %% 10000 == 0 ){
      #   plot.ts(t(WeightMat),  plot.type = "single", col = pal.1(ne), xlab = "pre_rate", ylab =  "final weight",  xaxt = 'n')
      #   axis(1, at = seq.int(0, pre_n, by = 10), labels = seq.int(0, pre_n, by = 10))
      #   legend("topleft",legend = 1:ne, col = pal.1(ne), lty = 1,cex= 0.5)
      #   print(t)
      # }
      
    }
    # print(repeat_i)
    test_result = network_test(WeightMat)
    analyze_results = analyze_tuning(test_result)
    # save data for i-th repeat
    result = list(analyze_results.init = analyze_results.init, analyze_results = analyze_results, WeightMat.init = WeightMat.init, WeightMat = WeightMat, 
                  test_result.init=test_result.init, test_result = test_result,
                  conle = conle, conli = conli, conee = conee, conei = conei, conie = conie, conii = conii)
    
  }
  print(lens_name)
  save(list  = c( "results_list"),
       file = paste("data/lens_results_uniforminit_", lens_name,".RData", sep = ""))
  
}
stopImplicitCluster()

###################################################################### 
## Some sample plots
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
hist_record_final = matrix(0, nrow = hist_n, ncol = repeat_n)
# init
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  # init
  hist_record_init[,repeat_i] = hist((results$analyze_results.init$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
  # final
  hist_record_final[,repeat_i] = hist((results$analyze_results$tuninge[1,]+10) %% 180 , breaks = seq(0,180,length.out = 10), plot = F)$counts/ne
}

hist_sd_init = apply(hist_record_init, MARGIN = 1, FUN = sd)
hist_sd_final = apply(hist_record_final, MARGIN = 1, FUN = sd)
barCenters = barplot(rowMeans(hist_record_init), ylim = c(0,0.3), col = "red", xaxt = "n",width = 1, space  = 0,
                     xlab = "Preferred Orientation", main = "d. Distrbution of PO")
arrows(barCenters, rowMeans(hist_record_init) , barCenters,
       rowMeans(hist_record_init) + hist_sd_init, lwd = 2, angle = 90,
       code = 3, length = 0.05)
axis(side = 1, at = c(0.5,2,3.5, 5, 6.5, 8), labels = c(-90,-60,-30,0,30,60))

# Sample trace
sample_neuron_id = 4
test_sep = 100
test_indicator = (seq_along(results_list[[1]]$test_result.init$spikee[,sample_neuron_id])-1) %% 1000 <= (test_sep-1)
plot.ts(1:sum(test_indicator),results_list[[1]]$test_result.init$spikee[test_indicator,sample_neuron_id]*0.5 + 2, type = "l", 
        ylim = c(0,3), col = "red", xlab = "Time (ms)", ylab = "Spike", main = "b. Sample Traces")
lines(1:sum(test_indicator),results_list[[1]]$test_result$spikee[test_indicator,sample_neuron_id]*0.5, col = "blue")
abline(v = test_sep * 1:7, lty = 2)


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
lines(xaxis.ang,tuning_estimate(results_list[[1]]$analyze_results$tuninge[,sample_neuron_id],xaxis.ang), lwd = 2,col = "blue")

# cumulative plot of osi
cvare_init = c()
cvare_final = c()
for(repeat_i in 1:repeat_n){
  results = results_list[[repeat_i]]
  cvare_init = c(cvare_init, results$analyze_results.init$cvare)
  cvare_final = c(cvare_final, results$analyze_results$cvare)
}
plot(ecdf(cvare_init), col = "red", main = "e. Orientation Selectivity (1-CV)", ylab = "Cumulative")
lines(ecdf(cvare_final), col = "blue")
