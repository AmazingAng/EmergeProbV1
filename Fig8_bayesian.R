# Code for figure 8
# 2021/01/28

rm(list = ls())
setwd(choose.dir()) # set working directory to the folder of the project (~/EmergeProbV1)

load("data/results_fig24.R")
source("utils.R")
library(circular)
library(distr)
library(radial)
library(plotrix)
gaussian_example_curve <- function(ff.mu){
  f=exp( -pmin((tuning_x-ff.mu)^2, (tuning_x-tuning_n-ff.mu)^2, (tuning_x+tuning_n-ff.mu)^2)/(2*tuning_sd^2) )
}

## Plot
par(mfcol = c(2,4), lwd = 2)
plot.new()
text(0.5,0.5,"A Schematic plot of decoding model", cex = 1)

## B population vector
tuning_x = as.numeric(qvonmises(seq(1e-5,1-1e-5,length.out = 50), kappa = 0.8)+pi)/2/pi*180
tuning_n = length(tuning_x)
tuning_sd = 20
tuning_center = 45
tuning_rate =  exp( -pmin((tuning_x-tuning_center)^2, (tuning_x-180-tuning_center)^2, (tuning_x+180-tuning_center)^2)/(2*tuning_sd^2) )
oldpar<-radial.plot(tuning_rate,tuning_x/180*pi*2,line.col="black", main = "B Population vector decoder",
                    lwd=2,rad.col="lightblue", labels = NULL, xlab = "", ylab = "", show.grid.labels=F,show.radial.grid=F)
text(0.5,0, "Firing rate", cex = 1.2)
text(0,1.1, "Orientation", cex = 1.2)

## C Tuning curvse of encoding population
tuning_x = 0:179
tuning_n = 180
tuning_sd = 25

tuning_pref_init = as.numeric(qvonmises(seq(1e-5,1-1e-5,length.out = 19), kappa = 0.8)+pi)/2/pi*180
for(pref_i in tuning_pref_init){
  if(pref_i == tuning_pref_init[1]){
    plot(tuning_x, gaussian_example_curve(pref_i), type = "l", xlab = "Orientation", ylab = "Normalized response",
         ylim = c(0,1), main = "C Tuning curvse of encoding population", lwd = 1, col = "red",xaxt = "n")
    axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
    
  }else{
    lines(tuning_x, gaussian_example_curve(pref_i), lwd = 2, col = "red")
  }
}

# dual peak distribution
c0 = 1/(2*pi-2) 
# prior density distribution
den_dual<-function(x){
  de<-2-abs(sin(2*x));
  de*c0
}
dist_dual <-AbscontDistribution(d=den_dual, low1 = 0, up1 = pi)
qdist_dual <- q(dist_dual)
tuning_pref_dual = qdist_dual(seq(1e-5,1-1e-5,length.out = 19))/pi*180

for(pref_i in tuning_pref_dual){
  if(pref_i == tuning_pref_dual[1]){
    plot(tuning_x, gaussian_example_curve(pref_i), type = "l", xlab = "Orientation", ylab = "Normalized response",
         ylim = c(0,1), lwd = 1, col = "blue",xaxt = "n")
    axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
    
  }else{
    lines(tuning_x, gaussian_example_curve(pref_i), lwd = 2, col = "blue")
  }
}




# D and E: Estimation variability and bias of the decoding model
dt = 1
pre_n = 100 # input neurons
nl = pre_n
ff.bg = 10 # back ground firing rate (gaussian function)
ff.A = 40 # peak firing rate (gaussian function)
ff.sigma = 10 # sigma (gaussian function)

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

step = 36 # number of stimuli
tepoch = 5000 # time for each stimuli representationa
theta.set = (0:(step-1))*pi/step
ttot = tepoch * step
t_integrate = 200 # ms integration time

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

# Parallel
library(circular)
library(parallel)
library(foreach)
library(doParallel)

numCores <- detectCores()
numCores  = 25
registerDoParallel(numCores)  # use multicore, set to the number of our cores
repeat_n = length(results_list)

get_ori_from_sincos = function(sinx, cosx){
  if(cosx >= 0&sinx>=0){
    asin(sinx)/pi*180
  }else if(cosx >= 0&sinx<0){
    asin(sinx)/pi*180+ 2*180
  }else if(cosx<0 & sinx>=0){
    180 - asin(sinx)/pi*180
  }else if(cosx<0 & sinx<0){
    180 -  asin(sinx)/pi*180
  }
}

calc_popvec = function(x, dir_theta){
  output = (x%*%dir_theta) / sum(x)
  cos_output = output[1]/sqrt(sum(output^2))
  sin_output = output[2]/sqrt(sum(output^2))
  output_dir =get_ori_from_sincos(sin_output, cos_output)/2
}

popvec_list = foreach(repeat_i = 1:repeat_n, .packages = c("circular"))%dopar%{
  conle =  results_list[[repeat_i]]$conle
  conee =  results_list[[repeat_i]]$conee
  conli =  results_list[[repeat_i]]$conli
  conei =  results_list[[repeat_i]]$conei
  conie =  results_list[[repeat_i]]$conie
  conii =  results_list[[repeat_i]]$conii
  
  WeightMat = results_list[[repeat_i]]$WeightMat
  WeightMat.init = results_list[[repeat_i]]$WeightMat.init
  # ground truth
  ground_truth  = rep(theta.set/pi*180,each = tepoch/t_integrate )
  
  # init
  test_result_init = network_test(WeightMat.init)
  analyze_result_init = analyze_tuning(test_result_init)
  spikee_init = test_result_init$spikee
  pref_init = analyze_result_init$tuninge[1,]
  fre_init <- apply(X = spikee_init, MARGIN = 2, calc_fr, test.t = t_integrate)/(tepoch/t_integrate)
  # Orientation Preference
  pref_theta_init = pref_init/180*pi*2
  dir_theta_init = cbind(cos(pref_theta_init), sin(pref_theta_init))
  # Calculate estimated orientation
  estimated_popvec_init = apply(X = fre_init, MARGIN = 1, FUN = calc_popvec, dir_theta = dir_theta_init)
  estimated_mean_init =  tapply(estimated_popvec_init/180*pi*2, INDEX = rep(1:step, each = tepoch/t_integrate), mean.circular)/2
  estimated_mean_init[estimated_mean_init<0] = estimated_mean_init[estimated_mean_init<0] + pi
  # Calculate estimation error
  error_init = estimated_popvec_init-ground_truth
  # make the error as error of orientation variable: distance 178 to 0 should be -2, not -178
  error_init[error_init >= 90]= error_init[error_init >= 90]-180
  error_init[error_init <= -90]= error_init[error_init <= -90]+180
  # calculate estimation bias and standard deviation
  sd_error_init = tapply(error_init,  INDEX = rep(1:step, each = tepoch/t_integrate), sd)
  bias_error_init = tapply(error_init, INDEX = rep(1:step, each = tepoch/t_integrate), mean)
  
  # final
  test_result_final = network_test(WeightMat)
  analyze_result_final = analyze_tuning(test_result_final)
  spikee_final = test_result_final$spikee
  pref_final = analyze_result_final$tuninge[1,]
  fre_final <- apply(X = spikee_final, MARGIN = 2, calc_fr, test.t = t_integrate)/(tepoch/t_integrate)
  # Orientation Preference
  pref_theta_final = pref_final/180*pi*2
  dir_theta_final = cbind(cos(pref_theta_final), sin(pref_theta_final))
  # Calculate estimated orientation
  estimated_popvec_final = apply(X = fre_final, MARGIN = 1, FUN = calc_popvec, dir_theta = dir_theta_final)
  estimated_mean_final =  tapply(estimated_popvec_final/180*pi*2, INDEX = rep(1:step, each = tepoch/t_integrate), mean.circular)/2
  estimated_mean_final[estimated_mean_final<0] = estimated_mean_final[estimated_mean_final<0] + pi
  # Calculate estimation error
  error_final = estimated_popvec_final-ground_truth
  # make the error as error of orientation variable: distance 178 to 0 should be -2, not -178
  error_final[error_final >= 90]= error_final[error_final >= 90]-180
  error_final[error_final <= -90]= error_final[error_final <= -90]+180
  # calculate estimation bias and standard deviation
  sd_error_final = tapply(error_final,  INDEX = rep(1:step, each = tepoch/t_integrate), sd)
  bias_error_final = tapply(error_final, INDEX = rep(1:step, each = tepoch/t_integrate), mean)
  
  popvec_result = list(estimated_mean_init = estimated_mean_init,error_init = error_init,  bias_error_init = bias_error_init, sd_error_init = sd_error_init,
                       estimated_mean_final = estimated_mean_final,error_final = error_final,  bias_error_final = bias_error_final, sd_error_final = sd_error_final)
}
  
stopImplicitCluster()

mean_record_init = matrix(0, nrow = step, ncol = repeat_n)
bias_record_init = matrix(0, nrow = step, ncol = repeat_n)
sd_record_init = matrix(0, nrow = step, ncol = repeat_n)
mean_record_final = matrix(0, nrow = step, ncol = repeat_n)
bias_record_final = matrix(0, nrow = step, ncol = repeat_n)
sd_record_final = matrix(0, nrow = step, ncol = repeat_n)

for(repeat_i in 1:repeat_n){
  popvec_result = popvec_list[[repeat_i]]
  # init
  mean_record_init[,repeat_i] = popvec_result$estimated_mean_init
  bias_record_init[,repeat_i] = popvec_result$bias_error_init
  sd_record_init[,repeat_i] = popvec_result$sd_error_init
  
  # final
  mean_record_final[,repeat_i] = popvec_result$estimated_mean_final
  bias_record_final[,repeat_i] = popvec_result$bias_error_final
  sd_record_final[,repeat_i] = popvec_result$sd_error_final
}

plot(theta.set/pi*180, rowMeans(sd_record_init), col = "red", type="l", xlab = "Orientation", ylab = "Variability (deg)",
     xaxt = "n", main = "D Estimation variability", ylim = c(0, 5.5))
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))

plot(theta.set/pi*180, rowMeans(sd_record_final), col = "blue", type="l", xlab = "Orientation", ylab = "Variability (deg)",
     xaxt = "n", ylim = c(0, 5.5))
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))

plot(theta.set/pi*180, rowMeans(bias_record_init), col = "red", type="l", xlab = "Orientation", ylab = "Bias (deg)", 
     xaxt = "n", main = "E Estimation bias")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))
plot(theta.set/pi*180, rowMeans(bias_record_final), col = "blue", type="l", xlab = "Orientation", ylab = "Bias (deg)",
     xaxt = "n")
axis(side = 1, at = c(-90,-60,-30,0,30,60)+90, labels = c(-90,-60,-30,0,30,60))

