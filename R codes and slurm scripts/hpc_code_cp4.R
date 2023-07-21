### fit transitions with 3 pre-existing conditions (one-at-a-time)
# load package
library(flexsurv)
library(rlist)
# load data, converted to data frame named 'x' (detailed code omitted)
#
N<-5
x$ethnicity<-replace(x$ethnicity, x$ethnicity=="Mixed", "Other")

# data preparation for flexsurv (clock reset, time unit in years)
# we use number to index the conditions 1: CVD, 2: T2D, 3: CKD: 4, MH, 5:HF
cp_summary_4<-function(a,b,c) {
  transtime<-vector("list")
  gender<-vector("list")
  age<-vector("list")
  imd<-vector("list")
  ethnicity<-vector("list")
  for (i in 1:(N-3)) {
    e_ind<-which((((pmax(x[,a],x[,b],x[,c])<x[,setdiff((1:N),c(a,b,c))[i]]) & (pmax(x[,a],x[,b],x[,c],x[,setdiff((1:N),c(a,b,c))[i]])<x[,setdiff((1:N),c(a,b,c,setdiff((1:N),c(a,b,c))[i]))]))|((pmax(x[,a],x[,b],x[,c])<x[,setdiff((1:N),c(a,b,c))[i]]) & (is.na(x[,setdiff((1:N),c(a,b,c,setdiff((1:N),c(a,b,c))[i]))])))) & (pmax(x[,a],x[,b],x[,c])>x$cohort_entry_date))
    transtime[[i]]<-as.numeric(x[e_ind,setdiff((1:N),c(a,b,c))[i]]-pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c]))/365.25
    gender[[i]]<-x$gender[e_ind]-1
    age[[i]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10  
    imd[[i]]<-x$imd[e_ind]
    ethnicity[[i]]<-x$ethnicity[e_ind]
  }
  e_ind<-which((pmax(x[,a],x[,b],x[,c])<x$dod) & is.na(pmin(x[,setdiff((1:N),c(a,b,c))[1]],x[,setdiff((1:N),c(a,b,c))[2]],na.rm = T)) & (pmax(x[,a],x[,b],x[,c])>x$cohort_entry_date))
  transtime[[i+1]]<-as.numeric(x$dod[e_ind]-pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c]))/365.25
  gender[[i+1]]<-x$gender[e_ind]-1
  age[[i+1]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+1]]<-x$imd[e_ind]
  ethnicity[[i+1]]<-x$ethnicity[e_ind]
  e_ind<-which((pmax(x[,a],x[,b],x[,c])<x$cohort_exit_date) & is.na(pmin(x[,setdiff((1:N),c(a,b,c))[1]],x[,setdiff((1:N),c(a,b,c))[2]],x$dod,na.rm = T)) & (pmax(x[,a],x[,b],x[,c])>x$cohort_entry_date))
  transtime[[i+2]]<-as.numeric(x$cohort_exit_date[e_ind]-pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c]))/365.25
  gender[[i+2]]<-x$gender[e_ind]-1
  age[[i+2]]<-(as.numeric(pmax(x[e_ind,a],x[e_ind,b],x[e_ind,c])-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind])/10
  imd[[i+2]]<-x$imd[e_ind]
  ethnicity[[i+2]]<-x$ethnicity[e_ind]
  sumy<-list(transtime=transtime,gender=gender,age=age,imd=imd,ethnicity=ethnicity,n=length(unlist(transtime)))
  return(sumy)
}

# excluding transitions emit from d3d4d5 (fitted separately using hpc_code_345.R)
a<-c(1,1,1,1,1,1,2,2,2)
b<-c(2,2,2,3,3,4,3,3,4)
c<-c(3,4,5,4,5,5,4,5,5)
trans_smy<-vector("list")
# fitting using flexsurv (knot point selection guided by aic)
for (l in 1:9) {
  smy<-cp_summary_4(a[l],b[l],c[l])
  trans<-vector("list")
  for (i in 1:(N-3)){
    mm_trans<-data.frame(id=1:smy$n,years=c(smy$transtime[[i]],smy$transtime[[setdiff((1:N),i)[1]]],smy$transtime[[setdiff((1:N),i)[2]]],smy$transtime[[4]]),gender=relevel(factor(c(smy$gender[[i]],smy$gender[[setdiff((1:N),i)[1]]],smy$gender[[setdiff((1:N),i)[2]]],smy$gender[[4]])),ref="1"),age=c(smy$age[[i]],smy$age[[setdiff((1:N),i)[1]]],smy$age[[setdiff((1:N),i)[2]]],smy$age[[4]]),imd=factor(c(smy$imd[[i]],smy$imd[[setdiff((1:N),i)[1]]],smy$imd[[setdiff((1:N),i)[2]]],smy$imd[[4]])),ethnicity=relevel(factor(c(smy$ethnicity[[i]],smy$ethnicity[[setdiff((1:N),i)[1]]],smy$ethnicity[[setdiff((1:N),i)[2]]],smy$ethnicity[[4]])),ref = "White"),status=c(rep(1,length(smy$transtime[[i]])),rep(0,smy$n-length(smy$transtime[[i]]))),trans=rep(1,smy$n))
    aic<-NA
    trans0<-vector("list")
    for (k in 0:5) {
      trans0[[k+1]]<-flexsurvspline(Surv(years, status) ~ age+gender+imd+ethnicity,inits=c(gamma0=0.1,gamma1=0.1), subset = (trans == 1), data=mm_trans, k=k, scale="hazard")
      aic[k+1]<-trans0[[k+1]]$AIC
    }
    trans[[i]]<-trans0[[which.min(aic)]]
  }
  trans_smy[[l]]<-trans
}

# store results
list.save(trans_smy, 'trans_cp4.rdata')