library("R.matlab")

rsfMRI<- readMat("./data/HCP_Desikan_rfMRI_REST1_BOLD.mat")

time_series<- list()
time_series_retest<- list()

subject_id<- list()

for(s in 1:length(rsfMRI$HCP.All.BOLD)){
  flat<- rsfMRI$HCP.All.BOLD[,s][[1]][[1]][,,1]
  time_series[[s]]<- (flat$dtseriesLR)
  time_series_retest[[s]]<- (flat$dtseriesRL)
  subject_id[[s]]<- flat$subjectid[1]
}


for(s in 1:length(rsfMRI$HCP.All.BOLD)){
  time_series[[s]] <- t(apply(time_series[[s]], 1, function(x)(x-mean(x))/sd(x)))
  time_series_retest[[s]] <- t(apply(time_series_retest[[s]], 1, function(x)(x-mean(x))/sd(x)))
}


HCP_desikan_data<- list("subject_id"=subject_id, "time_series"=time_series, "time_series_retest"=time_series_retest)

save( HCP_desikan_data,file = "HCP_desikan_data.RDa")
