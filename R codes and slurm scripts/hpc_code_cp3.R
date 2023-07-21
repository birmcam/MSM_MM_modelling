### fit transitions with 2 pre-existing conditions (one-at-a-time)
# load package
library(flexsurv)
library(rlist)
# load data, converted to data frame named 'x' (detailed code omitted)
# 
N<-5 # total number of conditions
x$ethnicity<-replace(x$ethnicity, x$ethnicity=="Mixed", "Other") # merge "Mixed" and "Other" eth groups

# data preparation for flexsurv (clock reset, time unit in years)
# we use number to index the conditions 1: CVD, 2: T2D, 3: CKD: 4, MH, 5:HF
cp_summary_3<-function(a,b) {
  transtime<-vector("list")
  gender<-vector("list")
  age<-vector("list")
  imd<-vector("list")
  ethnicity<-vector("list")
  for (i in 1:(N-2)) {
    e_ind<-which((((pmax(x[,a],x[,b])<x[,setdiff((1:N),c(a,b))[i]]) & (pmax(x[,a],x[,b],x[,setdiff((1:N),c(a,b))[i]])<pmin(x[,setdiff((1:N),c(a,b,setdiff((1:N),c(a,b))[i]))[1]],x[,setdiff((1:N),c(a,b,setdiff((1:N),c(a,b))[i]))[2]],na.rm = T)))|((pmax(x[,a],x[,b])<x[,setdiff((1:N),c(a,b))[i]]) & (is.na(pmin(x[,setdiff((1:N),c(a,b,setdiff((1:N),c(a,b))[i]))[1]],x[,setdiff((1:N),c(a,b,setdiff((1:N),c(a,b))[i]))[2]],na.rm = T))))) & (pmax(x[,a],x[,b])>x$cohort_entry_date))
    transtime[[i]]<-as.numeric(x[e_ind,setdiff((1:N),c(a,b))[i]]-pmax(x[e_ind,a],x[e_ind,b]))/365.25
    gender[[i]]<-x$gender[e_ind]-1
    age[[i]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
    imd[[i]]<-x$imd[e_ind]
    ethnicity[[i]]<-x$ethnicity[e_ind]
  }
  e_ind<-which((pmax(x[,a],x[,b])<x$dod) & (is.na(pmin(x[,setdiff((1:N),c(a,b))[1]],x[,setdiff((1:N),c(a,b))[2]],x[,setdiff((1:N),c(a,b))[3]],na.rm = T))) & (is.na(x$dod)==F) & (pmax(x[,a],x[,b])>x$cohort_entry_date))
  transtime[[i+1]]<-as.numeric(x$dod[e_ind]-pmax(x[e_ind,a],x[e_ind,b]))/365.25
  gender[[i+1]]<-x$gender[e_ind]-1
  age[[i+1]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+1]]<-x$imd[e_ind]
  ethnicity[[i+1]]<-x$ethnicity[e_ind]
  e_ind<-which((pmax(x[,a],x[,b])<x$cohort_exit_date) & (is.na(pmin(x[,setdiff((1:N),c(a,b))[1]],x[,setdiff((1:N),c(a,b))[2]],x[,setdiff((1:N),c(a,b))[3]],x$dod,na.rm = T))) & (pmax(x[,a],x[,b])>x$cohort_entry_date))
  transtime[[i+2]]<-as.numeric(x$cohort_exit_date[e_ind]-pmax(x[e_ind,a],x[e_ind,b]))/365.25
  gender[[i+2]]<-x$gender[e_ind]-1
  age[[i+2]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+2]]<-x$imd[e_ind]
  ethnicity[[i+2]]<-x$ethnicity[e_ind]
  sumy<-list(transtime=transtime,gender=gender,age=age,imd=imd,ethnicity=ethnicity,n=length(unlist(transtime)))
  return(sumy)
}

a<-c(1,1,1,1,2,2,2,3,3,4)
b<-c(2,3,4,5,3,4,5,4,5,5)
trans_smy<-vector("list")
# fitting using flexsurv (knot point selection guided by aic)
for (l in 1:10) {
  smy<-cp_summary_3(a[l],b[l])
  trans<-vector("list")
  for (i in 1:(N-2)){
    mm_trans<-data.frame(id=1:smy$n,years=c(smy$transtime[[i]],smy$transtime[[setdiff((1:N),i)[1]]],smy$transtime[[setdiff((1:N),i)[2]]],smy$transtime[[setdiff((1:N),i)[3]]],smy$transtime[[5]]),gender=relevel(factor(c(smy$gender[[i]],smy$gender[[setdiff((1:N),i)[1]]],smy$gender[[setdiff((1:N),i)[2]]],smy$gender[[setdiff((1:N),i)[3]]],smy$gender[[5]])),ref="1"),age=c(smy$age[[i]],smy$age[[setdiff((1:N),i)[1]]],smy$age[[setdiff((1:N),i)[2]]],smy$age[[setdiff((1:N),i)[3]]],smy$age[[5]]),imd=factor(c(smy$imd[[i]],smy$imd[[setdiff((1:N),i)[1]]],smy$imd[[setdiff((1:N),i)[2]]],smy$imd[[setdiff((1:N),i)[3]]],smy$imd[[5]])),ethnicity=relevel(factor(c(smy$ethnicity[[i]],smy$ethnicity[[setdiff((1:N),i)[1]]],smy$ethnicity[[setdiff((1:N),i)[2]]],smy$ethnicity[[setdiff((1:N),i)[3]]],smy$ethnicity[[5]])),ref="White"),status=c(rep(1,length(smy$transtime[[i]])),rep(0,smy$n-length(smy$transtime[[i]]))),trans=rep(1,smy$n))
    aic<-NA
    trans0<-vector("list")
    for (k in 0:5) {
      trans0[[k+1]]<-flexsurvspline(Surv(years, status) ~ age+gender+imd+ethnicity,inits=c(gamma0=0.5,gamma1=0.1), subset = (trans == 1), data=mm_trans, k=k, scale="hazard")
      aic[k+1]<-trans0[[k+1]]$AIC
    }
    trans[[i]]<-trans0[[which.min(aic)]]
  }
  trans_smy[[l]]<-trans
}

# store results
list.save(trans_smy, 'trans_cp3.rdata')