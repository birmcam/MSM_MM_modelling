### fit transitions from CVD: 'cvd->cvd+xxx'
# load package
library(flexsurv)
library(rlist)
# load data, converted to data frame named 'x' (detailed code omitted)
# 
N<-5
x$ethnicity<-replace(x$ethnicity, x$ethnicity=="Mixed", "Other")

# data preparation for flexsurv (clock reset, time unit in years)
cp_summary_2<-function(a) {
  transtime<-vector("list")
  gender<-vector("list")
  age<-vector("list")
  imd<-vector("list")
  ethnicity<-vector("list")
  for (i in 1:(N-1)) {
    e_ind<-which((((x[,a]<x[,setdiff((1:N),a)[i]]) & (pmax(x[,a],x[,setdiff((1:N),a)[i]])<pmin(x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[1]],x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[2]],x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[3]],na.rm = T)))|((x[,a]<x[,setdiff((1:N),a)[i]]) & (is.na(pmin(x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[1]],x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[2]],x[,setdiff((1:N),c(a,setdiff((1:N),a)[i]))[3]],na.rm = T))))) & (pmin(x[,a],x[,setdiff((1:N),a)[i]])>x$cohort_entry_date))
    transtime[[i]]<-as.numeric(x[e_ind,setdiff((1:N),a)[i]]-x[e_ind,a])/365.25
    gender[[i]]<-x$gender[e_ind]-1
    age[[i]]<-(as.numeric(x[e_ind,a]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10 ## in 10 years; treat age as pc
    imd[[i]]<-x$imd[e_ind] #
    ethnicity[[i]]<-x$ethnicity[e_ind]
  }
  e_ind<-which((x[,a]<x$dod) & (is.na(pmin(x[,setdiff((1:N),a)[1]],x[,setdiff((1:N),a)[2]],x[,setdiff((1:N),a)[3]],x[,setdiff((1:N),a)[4]],na.rm = T))) & (is.na(x$dod)==F) & (pmin(x[,a],x$dod)>x$cohort_entry_date))
  transtime[[i+1]]<-as.numeric(x$dod[e_ind]-x[e_ind,a])/365.25
  gender[[i+1]]<-x$gender[e_ind]-1
  age[[i+1]]<-(as.numeric(x[e_ind,a]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+1]]<-x$imd[e_ind] #
  ethnicity[[i+1]]<-x$ethnicity[e_ind]
  e_ind<-which((x[,a]<x$cohort_exit_date) & (is.na(pmin(x[,setdiff((1:N),a)[1]],x[,setdiff((1:N),a)[2]],x[,setdiff((1:N),a)[3]],x[,setdiff((1:N),a)[4]],x$dod,na.rm = T))) & (pmin(x[,a],x$cohort_exit_date)>x$cohort_entry_date))
  transtime[[i+2]]<-as.numeric(x$cohort_exit_date[e_ind]-x[e_ind,a])/365.25
  gender[[i+2]]<-x$gender[e_ind]-1
  age[[i+2]]<-(as.numeric(x[e_ind,a]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+2]]<-x$imd[e_ind] #
  ethnicity[[i+2]]<-x$ethnicity[e_ind]
  sumy<-list(transtime=transtime,gender=gender,age=age,imd=imd,ethnicity=ethnicity,n=length(unlist(transtime)))
  return(sumy)
}

# fitting using flexsurv (knot point selection guided by aic)
smy<-cp_summary_2(1)
trans<-vector("list")
for (i in 1:(N-1)){
  mm_trans<-data.frame(id=1:smy$n,years=c(smy$transtime[[i]],smy$transtime[[setdiff((1:N),i)[1]]],smy$transtime[[setdiff((1:N),i)[2]]],smy$transtime[[setdiff((1:N),i)[3]]],smy$transtime[[setdiff((1:N),i)[4]]],smy$transtime[[6]]),gender=relevel(factor(c(smy$gender[[i]],smy$gender[[setdiff((1:N),i)[1]]],smy$gender[[setdiff((1:N),i)[2]]],smy$gender[[setdiff((1:N),i)[3]]],smy$gender[[setdiff((1:N),i)[4]]],smy$gender[[6]])),ref="1"),age=c(smy$age[[i]],smy$age[[setdiff((1:N),i)[1]]],smy$age[[setdiff((1:N),i)[2]]],smy$age[[setdiff((1:N),i)[3]]],smy$age[[setdiff((1:N),i)[4]]],smy$age[[6]]),imd=factor(c(smy$imd[[i]],smy$imd[[setdiff((1:N),i)[1]]],smy$imd[[setdiff((1:N),i)[2]]],smy$imd[[setdiff((1:N),i)[3]]],smy$imd[[setdiff((1:N),i)[4]]],smy$imd[[6]])),ethnicity=relevel(factor(c(smy$ethnicity[[i]],smy$ethnicity[[setdiff((1:N),i)[1]]],smy$ethnicity[[setdiff((1:N),i)[2]]],smy$ethnicity[[setdiff((1:N),i)[3]]],smy$ethnicity[[setdiff((1:N),i)[4]]],smy$ethnicity[[6]])),ref = "White"),status=c(rep(1,length(smy$transtime[[i]])),rep(0,smy$n-length(smy$transtime[[i]]))),trans=rep(1,smy$n))
  aic<-NA
  trans0<-vector("list")
  for (k in 0:5) {
    trans0[[k+1]]<-flexsurvspline(Surv(years, status) ~ age+gender+imd+ethnicity, inits=c(gamma0=1,gamma1=0.1), subset = (trans == 1), data=mm_trans, k=k, scale="hazard")
    aic[k+1]<-trans0[[k+1]]$AIC
  }
  trans[[i]]<-trans0[[which.min(aic)]]
}

# store results
list.save(trans, 'trans_1.rdata')