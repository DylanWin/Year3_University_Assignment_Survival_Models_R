#
#     R-code to generate data for Assignment 1.  This code MUST be
#     placed at the start of your own R-script.  You must edit
#     the argument to the set.seed( ) function to fit your own
#     registration number
#
set.seed(0633)
#
#     Generate sample A times with censorng
#
LenA <- 100
TimeA <- sort(round(rexp(LenA, 0.5), digits = 10))
CensorA <- rbinom(LenA,1,0.5)
cbind(TimeA, CensorA)
#
#     Generate sample B times with censorng
#
LenB <- 100
TimeB <- sort(round(rexp(LenB, 0.5), digits = 10))
CensorB <- rbinom(LenB,1,0.6)
cbind(TimeB, CensorB)
#
#     Combine data
#
Time <- c(TimeA, TimeB)
Censor <- c(CensorA, CensorB)
Group <- factor(c(rep(1,LenA), rep(2,LenB)))
cbind(Time, Censor, Group)

####################################################
# Please insert your R code after this line
####################################################
# F79SU Survival Models
# Assignment 1

# Package reuired to create Kaplan Meier model
library(survival)
# # install package coin##
# library(coin)
# library(plyr)

# Displaying data extracted from given R code
# Data of patience taking Drug A(Group 1) and Drug B(Group 2)
# (ascending time of survival time)
# spacing<- extra space added betwwen data from Drug A and B, better display
spacing<-rep('  ', 50)
info<-data.frame(TimeA[1:50], CensorA[1:50],  spacing, TimeA[51:100], CensorA[51:100], spacing, TimeB[1:50], CensorB[1:50], spacing, TimeB[51:100], CensorB[51:100])
colnames(info)<-c('Time A', 'Censor A', '','Time A', 'Censor A', '', 'Time B', 'Censor B', '','Time B', 'Censor B')
info

# Kapan Meier model for Drug A
Surv(TimeA, CensorA)
K_M_model_A<-survfit(Surv(TimeA, CensorA)~1, conf.type='plain')
names(K_M_model_A)
summary(K_M_model_A)
plot(K_M_model_A, xlab='Time (years) to death of patient', ylab='Probability of survival', main='Plot of Kaplan Meier model of test on Drug A')
K_M_model_A_unbound<-survfit(Surv(TimeA, CensorA)~1, conf.type='none')

# Kaplan Meier model for Drug B
Surv(TimeB, CensorB)
K_M_model_B<-survfit(Surv(TimeB, CensorB)~1, conf.type='plain')
names(K_M_model_B)
summary(K_M_model_B)
plot(K_M_model_B, xlab='Time (years) to death of patient', ylab='Probability of survival', main='Plot of Kaplan Meier model of test on Drug B')

# Plot both Kaplan Meier model of Drug A and Drug B together
K_M_model_A_unbound<-survfit(Surv(TimeA, CensorA)~1, conf.type='none')
K_M_model_B_unbound<-survfit(Surv(TimeB, CensorB)~1, conf.type='none')
plot(K_M_model_A_unbound, main='Plot of Kaplan Meier model of test on Drug A and Drug B',
     xlab='Time (years) to death of patient', ylab='Probability of survival',
     col='blue')
lines(K_M_model_B_unbound, col='red')
legend('topright', legend=c('Drug A', 'Drug B'), col=c('blue', 'red'), lty=c(1,1), cex=0.8)

K_M_model_A
K_M_model_B

# Prob of patient taking Drug A surviving for more than 2.3 years
Ans_1<-summary(survfit(Surv(TimeA, CensorA)~1), times=2.3)
Ans_1
#Answer is 0.591

# Prob of patient taking Drug B surviving for less than/equal to 3.1 years
Ans_2<-summary(survfit(Surv(TimeB, CensorB)~1), times=3.1 )
Ans_2
#Answer is 1-0.362<- 0.638

# Cox Model on both drug type (A and B)
# Group<-c(rep(1,100), rep(2,100)), a drug group indicator, 100 for Drug A(Group 1) then 100 for Drug B(Group 2)
# Time<-Survival time of all patience (all Drug A then B, ascending time)
# Censor<-Observation censor of all patience (all Drug A then B, ascending time)
# corresponding to patient survival times
Cox_Model_A_to_B<-coxph(Surv(Time, Censor)~Group)
summary(Cox_Model_A_to_B)
Cox_Model_A_to_B
survdiff(Surv(Time, Censor)~Group)
beta_hat<-Cox_Model_A_to_B$coefficients

# Create dataframe to store survival time(ascending order), observation censor and drug group indicator
# Note that here data no longer displayed by all Drug A then Drug B, but by ascending time
# In short, it is a mix of rug A and Drug B depending on survival time of patient
data_df <- as.data.frame(matrix(0, ncol = 3, nrow = 200))
data_df[1]<-Time
data_df[2]<-Censor
data_df[3]<-as.vector(as.integer(Group))
colnames(data_df)<-c('Time', 'Censor', 'Group')
data_df<-data_df[order(data_df$Time),]
data_df

# Create dataframe to store nominator and denominator values of of both Drug A and B patients and finally drug group indicator
# A_nom being nominator of Drug A (integer only part)
# B_nom being nominator of Drug B (exponential(integer) only part)
# A_denom being denominator of Drug A (integer only part)
# B_denom being denominator of Drug B (exponential(integer) part)
info_df <- as.data.frame(matrix(0, ncol = 5, nrow = 200))
info_df[5]<-data_df$Group
colnames(info_df)<-c('A_nom', 'B_nom', 'A_denom', 'B_denom', 'Group')

# Initial patients alive for Drug A and B is 100 each
A_denom<-100
B_denom<-100
# Deaths for start period would be zero
censor_A_nom<-0
censor_B_nom<-0

# Compute nominator and denominator part of Drug A and B separately to be merged later
for (i in 1:200){
  # If at that time point, the survival time of the patient is censored
  if (data_df[i, 2]==0){
    # If the patient at that time point is from Group 1 (Drug A)
    if(data_df[i,3]==1){
      censor_A_nom<-censor_A_nom+1
      # If the patient at that time point is from Group 2 (Drug B)
    }else if(data_df[i,3]==2){
      censor_B_nom<-censor_B_nom+1
      
    }
    # If at that time point, the survival time of the patient is censored
  }else if(data_df[i,2]==1){
    # If the patient at that time point is from Group 1 (Drug A)
    if(data_df[i,3]==1){
      
      info_df[i,1]<-1
      
      info_df[i,2]<-0
      
      A_denom<-A_denom-(censor_A_nom)
      info_df[i,3]<-A_denom
      
      B_denom<-B_denom-(censor_B_nom)
      info_df[i,4]<-B_denom
      
      censor_A_nom<-1
      censor_B_nom<-0
      # If the patient at that time point is from Group 2 (Drug B) 
    }else if(data_df[i,3]==2){
      
      info_df[i,1]<-0
      
      info_df[i,2]<-1
      
      A_denom<-A_denom-(censor_A_nom)
      info_df[i,3]<-A_denom
      
      B_denom<-B_denom-(censor_B_nom)
      info_df[i,4]<-B_denom
      
      censor_A_nom<-0
      censor_B_nom<-1
    }
  }
}
# Display data
info_df
# Remove rows where both nominator of Drug A and Drug B are 0
# Menas that at that data time point the patient survival time is censored
info_df<-info_df[!((info_df$A_nom==0)&(info_df$B_nom==0)), ]
# Display data
info_df
# Remove Group column
info_df<-info_df[-5]

# To remove extra constants from denominator Exp: 5exp(x)+5 proportionate to exp(x)+1
for (i in 1:length(info_df[,3])){
  if(info_df[i,3]==info_df[i,4]){
    info_df[i,3]<-1
    info_df[i,4]<-1
  }
  
}


# To plot partial log likelihood against beta values range(0.01 to 1 by 0.01 step)
partial_log_lik_plot<-function(n, info, expected_beta){
  beta<-numeric(n)
  values<-numeric(n)
  
  for (i in 1:n){
    
    beta[i]<-i/100
    
    # Attaining log likelihood values by substituting each remaining rows values in the info dataframe
    # into the formula below to get the final value for the particular beta value
    for (j in 1:length(info[,3])){
      # Formula for each component
      # log(B): (B_nom*beta[i])-log(A_denom+(B_denom*beta[i]))
      x<-info[j,3]+info[j,4]*exp(beta[i])
      term<-((info[j,2]*beta[i]))-log(x)
      values[i]<-values[i]+term
    }
  }
  plot(beta, values, ylab='partial_loglikelihood',main='Partial Log Likelihood for beta between 0.1 and 1', type='l')
  # Plot vertical line of beta_hat value obtained from Cox Regression summary data (0.1911)
  abline(v=expected_beta, col='red', lty='solid')
}
partial_log_lik_plot(100, info_df, beta_hat)

# To plot score function value against beta values range(0.01 to 1 by 0.01 step)
score_func_plot<-function(n, info, expected_beta){
  beta<-numeric(n)
  score_value<-numeric(n)
  
  # Attaining score function values by substituting each remaining rows values in the info dataframe
  # into the formula below to get the final value for the particular beta value
  for (i in 1:n){
    # Formula for each component
    # dlog(B)/dB : B_nom-(B_denom*beta[i])/(A_denom+(B_denom*beta[i]))
    beta[i]<-i/100
    
    for (j in 1:length(info[,3])){
      
      term<-info[j,2]-(info_df[j,4]*exp(beta[i]))/(info[j,3]+info[j,4]*exp(beta[i]))
      score_value[i]<-score_value[i]+term
    }
  }
  plot(beta, score_value, main='Score function for beta between 0.1 and 1', type='l')
  # Plot vertical line of beta_hat value obtained from Cox Regression summary data (0.1911)
  abline(v=expected_beta, col='red', lty='solid')
}
score_func_plot(100, info_df, beta_hat)

# To check score function values with particular beta value
score_func<-function(info, expected_beta){
  score_value<-0
  for (j in 1:length(info[,3])){
    term<-info[j,2]-(info[j,4]*exp(expected_beta))/(info[j,3]+info[j,4]*exp(expected_beta))
    score_value<-score_value+term
    
  }
  score_value
}
# Using beta_hat value
score_func(info_df, beta_hat)





