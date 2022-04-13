# set global theme for all plots

th<-theme(plot.title = element_text(size = 12, face = "bold"),axis.title=element_text(size=12,color = "black"),
          axis.text.x = element_text(size=10, color = "black"),axis.text.y = element_text(size=10, color = "black"))

#melt dataframe of time and ODs (used in follow_ups.R)
myfun <- function(x){
  melt(x,id.vars = "Time_h", variable.name = "wells", value.name= "OD")
}

# function for creating annotation from file name (used in follow_ups.R)

create_annot_single <- function(input_file){
  
  elements<-strsplit(input_file, "_", perl = TRUE)
  elements
  #sp<-strsplit(ID, "/", perl = TRUE)
  Bug_ID<-elements[[1]][2]
  Replicate_no<-elements[[1]][3]
  Plate_no<-elements[[1]][1]
  
  annotation<-rbind(c(Bug_ID, Replicate_no, Plate_no))
  annotation
  
  return(annotation)
  
}

create_annot_CB <- function(input_file){
  
  elements<-strsplit(input_file, "_", perl = TRUE)
  elements
  #sp<-strsplit(ID, "/", perl = TRUE)
  Bug_ID<-elements[[1]][3]
  Replicate_no<-elements[[1]][4]
  Plate_no<-elements[[1]][2]
  
  annotation<-rbind(c(Bug_ID, Replicate_no,Plate_no))
  #annotation
  
  return(annotation)
  
}

#background (media) correction

bg_sub <- function (X) {
  
  #subtract 1st time point or any minimum if lower than 1st (to avoid negative values)
  X$OD<-as.double(X$OD)
  min_tp <- X %>% dplyr::group_by(wells) %>%
    dplyr::summarise(mintp = min(OD))
  
  Y<-dplyr::inner_join (X, min_tp)
  head(Y)
  firsttp <- subset (Y,  Time_h == 0) %>% dplyr::group_by (wells) %>%
    dplyr::summarise(zero_tp = OD)
  
  test <- dplyr::inner_join (Y, firsttp) %>%
    mutate(OD_to_subtract = pmin(zero_tp, mintp), adj_OD = round(OD - OD_to_subtract, digits = 3))
  
  class(test$adj_OD)
  return(test)
  
  #final$adj_OD = round((final$OD - final$OD_to_subtract), digits = 3)
  
  #final$adj_OD[final$adj_OD<0] <- 0
  
  #return(final)
  
}

# Calculate AUC by trapezoidal rule (this is correct for for unequal spaced data)

calcAUC <- function(adj_OD, Time_h) {
  N <- length(adj_OD)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (adj_OD[s1] + adj_OD[s2])/2 * (Time_h[s2]-Time_h[s1]) )
}


