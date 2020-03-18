################################
## State-transition matrices ##
################################
# WEEKLY transitions
Sgerm = 0.5
Sseed = 0.2
Sadult = 0.05
GermMort = 0.01
SeedMort = 0.01
AdultMort = 0.001

# initial pop matrix (i.e. sown seeds)
N_0 <- matrix(c(100, 0, 0, 0, 0), ncol=1)

# state-transition matrix at time 1
ps <- matrix(c(1-Sgerm, Sgerm, 0, 0, 0,#1st column- seed
               0, 1-Sseed-GermMort, Sseed, 0, GermMort,#2nd- germ
               0, 0, 1-Sadult-SeedMort, Sadult, SeedMort,#3rd- sdl
               0, 0, 0, 1-AdultMort, AdultMort,#4th- adult
               0, 0, 0, 0, 1), #5th- mortality
             nrow=5, ncol=); ps

# projecting 12 time-steps (weeks?) into the future
N_1 <- ps %*% N_0; N_1
N_2 <- ps %*% N_1; N_2
N_3 <- ps %*% N_2; N_3
N_4 <- ps %*% N_3; N_4
N_5 <- ps %*% N_4; N_5
N_6 <- ps %*% N_5; N_6
N_7 <- ps %*% N_6; N_7
N_8 <- ps %*% N_7; N_8
N_9 <- ps %*% N_8; N_9
N_10 <- ps %*% N_9; N_10
N_11 <- ps %*% N_10; N_11
N_12 <- ps %*% N_11; N_12

# using exponentiation for projection instead
t=12
#install.packages("expm")
library(expm)
N_12_exp <- (ps %^% 12) %*% N_0; N_12_exp

# pull out abundances over time
N_seed_bank <- c(N_0[1], N_1[1], N_2[1], N_3[1], N_4[1], N_5[1], N_6[1], N_7[1], N_8[1], N_9[1], N_10[1], N_11[1], N_12[1])
N_germ <- c(N_0[2], N_1[2], N_2[2], N_3[2], N_4[2], N_5[2], N_6[2], N_7[2], N_8[2], N_9[2], N_10[2], N_11[2], N_12[2])
N_seedlings <- c(N_0[3], N_1[3], N_2[3], N_3[3], N_4[3], N_5[3], N_6[3], N_7[3], N_8[3], N_9[3], N_10[3], N_11[3], N_12[3])
N_adults <- c(N_0[4], N_1[4], N_2[4], N_3[4], N_4[4], N_5[4], N_6[4], N_7[4], N_8[4], N_9[4], N_10[4], N_11[4], N_12[4])
N_dead <- c(N_0[5], N_1[5], N_2[5], N_3[5], N_4[5], N_5[5], N_6[5], N_7[5], N_8[5], N_9[5], N_10[5], N_11[5], N_12[5])

# plot abundances over time for each stage
plot(N_seed_bank, col="purple", main="Abundance per stage class over time", xlab="Week", ylab="Abundance", sub="purple=SB, red=germ, blue=seedl, green=adult, black=dead")
lines(N_seed_bank, col="purple")
points(N_germ, col="blue")
lines(N_germ, col="blue")
points(N_seedlings, col="darkgreen")
lines(N_seedlings, col="darkgreen")
points(N_adults, col="red")
lines(N_adults, col="red")
points(N_dead, col="black")
lines(N_dead, col="black")
