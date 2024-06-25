#This function takes a dataset with sites and computes an adjacency matrix of 100 m 
#connections. It then calculates how many of the 100 m connections are shared with 
#a neighbor. 
#The input "Con" is a dataset with sites

NeighborCon = function (Con) {
  
  Con = Con[order(Con$Farm_ID),]
  adj = st_is_within_distance(Con$geometry, dist = 100, sparse = FALSE)
  
  Farm = Con$Farm_ID
  
  if (length(unique(Farm)) == 1) {   
    return(0)  #Checks whether there is only one farm in this area, then returns 0
  } else {
  
  a = as.data.frame(table(Farm))[,2] #Counts how many sites belong to a specific farm
  M = adj*1 #transforms the adjacency matrix to a 0/1 matrix
  
  #Connections with others
  current = 0                   #A counter
  my_list_Horizontal = list()   #Initializes list
  my_list_Vertical = list()     #Initializes list
  my_list_Symmetric = list()    #Initializes list
  
  z = dim(M)[2]                 #Dimension of adj.matrix (number of sites), length(Con) would be the same
  
  for (j in 1:length(a)){       #lenght(a) is the number of individual farmers
    if (j == 1){                #First Farm Special Case
      x = 1 + current
      y = a[j] + current
      my_list_Vertical[[j]] = 0
      my_list_Horizontal[[j]] = M[x:y,(y+1):z]
      my_list_Symmetric[[j]] = M[x:y,x:y]
      my_list_Symmetric[[j]][lower.tri(my_list_Symmetric[[j]], diag=TRUE)] = 0
      current = y
      
    } else if (j == length(a)){   #Last Farm Special Case
      x = 1 + current
      y = a[j] + current
      my_list_Vertical[[j]] = M[1:(x-1),x:y]
      my_list_Horizontal[[j]] = 0
      my_list_Symmetric[[j]] = M[x:y,x:y]
      my_list_Symmetric[[j]][lower.tri(my_list_Symmetric[[j]], diag=TRUE)] = 0
      current = y
      
    } else {                      #Intermediate Farm Special Case
      x = 1 + current
      y = a[j] + current
      my_list_Vertical[[j]] = M[1:(x-1),x:y]
      my_list_Horizontal[[j]] = M[x:y,(y+1):z]
      my_list_Symmetric[[j]] = M[x:y,x:y]
      my_list_Symmetric[[j]][lower.tri(my_list_Symmetric[[j]], diag=TRUE)] = 0 
      
      current = y
    }
  }
  
  Count_own = list()
  Count_others = list()
  
  for (i in 1:length(a)) {    
    Count_own[[i]] = sum(unlist(my_list_Symmetric[[i]]))
    Count_others[[i]] = sum(unlist(my_list_Vertical[[i]], my_list_Horizontal[[i]]))
  }
  
  Count = list(Count_own, Count_others)
  
  Con_areas = do.call(rbind.data.frame, Count)
  Con_areas = t(Con_areas)
  rownames(Con_areas) = as.data.frame(table(Farm))[,1]
  colnames(Con_areas) = c("internal","external")
  
  return(Con_areas)
  }
}