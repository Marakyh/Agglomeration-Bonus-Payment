## Topologies of the networks ##

load("./Data/NBP_A.Rda")

#Subset
NBP_A <- NBP_A[NBP_A$District_ID %in% c(242:246), ]

#Control
NBP_Com0 = NBP_A[NBP_A$Rule == 0,]
net0 = Graph(NBP_Com0)
E <- ecount(net0)
N <- vcount(net0)
DD0 <- degree_distribution(net0, cumulative = FALSE) #https://igraph.org/r/doc/degree.html

#Treatment
NBP_Com1 = NBP_A[NBP_A$Rule == 1,]
net1 = Graph(NBP_Com1)
E1 <- ecount(net1)
N1 <- vcount(net1)
DD1 <- degree_distribution(net1, cumulative = FALSE)


# Generate a sequence of degrees from 0 up to the length of DD0 minus 1
degrees0 = 0:(length(DD0) - 1)
degrees1 = 0:(length(DD1) - 1)

plot(degrees0, DD0, type = "l", col = "blue", xlab = "Degree", ylab = "Degree Distribution", main = "Degree Distribution", xlim = c(0, max(degrees0)))
lines(degrees1, DD1, col = "red") # Add DD1 to the same plot
legend("topright", legend = c("Preservation area", "Networking area"), col = c("blue", "red"), lty = 1, bty = "n")

#------------
#More precise:  
NBP_A0 = NBP_A[NBP_A$Rule == 0,]
Municipality_ID_sorted = sort(unique(NBP_A0$Municipality_ID))
E0 = numeric(length(Municipality_ID_sorted))
N0 = numeric(length(Municipality_ID_sorted))

for (i in 1: length(Municipality_ID_sorted)) {
  NBP_Com = NBP_A0[NBP_A0$Municipality_ID == Municipality_ID_sorted[i],]
  net = Graph(NBP_Com)
  if(is(net,"igraph")) {
    E0[i] <- ecount(net)
    N0[i] <- vcount(net)
  } else {
    E0[i] <- 0
    N0[i] <- 0
  }
}

E0_sum <- sum(E0)
N0_sum <- sum(N0)

#Compute Gamma Treatment
NBP_A1 = NBP_A[NBP_A$Rule == 1,]
Municipality_ID_sorted = sort(unique(NBP_A1$Municipality_ID))
E1 = numeric(length(Municipality_ID_sorted))
N1 = numeric(length(Municipality_ID_sorted))

for (i in 1: length(Municipality_ID_sorted)) {
  NBP_Com = NBP_A1[NBP_A1$Municipality_ID == Municipality_ID_sorted[i],]
  net = Graph(NBP_Com)
  if(is(net,"igraph")) {
    E1[i] <- ecount(net)
    N1[i] <- vcount(net)
  } else {
    E1[i] <- 0
    N1[i] <- 0
    
  }
}

E1_sum <- sum(L1)
N1_sum <- sum(N1)




