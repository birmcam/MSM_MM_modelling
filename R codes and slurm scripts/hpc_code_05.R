### fit transition 'none->hf'
# load package
library(flexsurv)
library(rlist)
# load data, converted to data frame named 'x' (detailed code omitted)
# 
N<-5
x$ethnicity<-replace(x$ethnicity, x$ethnicity=="Mixed", "Other")

# data preparation for flexsurv (clock reset, time unit in years)
transtime<-vector("list")
gender<-vector("list")
imd<-vector("list")
ethnicity<-vector("list")
for (i in 1:N) {
  e_ind<-which((((x[,i]==pmin(x[,1],x[,2],x[,3],x[,4],x[,5],na.rm = T)) & (x[,i]<pmin(x[,setdiff((1:N),i)[1]],x[,setdiff((1:N),i)[2]],x[,setdiff((1:N),i)[3]],x[,setdiff((1:N),i)[4]],na.rm = T)))|((x[,i]==pmin(x[,1],x[,2],x[,3],x[,4],x[,5],na.rm = T)) & is.na(pmin(x[,setdiff((1:N),i)[1]],x[,setdiff((1:N),i)[2]],x[,setdiff((1:N),i)[3]],x[,setdiff((1:N),i)[4]],na.rm = T)))) & (x[,i]>x$cohort_entry_date) & (x[,i]<"2020-5-11"))
  transtime[[i]]<-as.numeric(x[e_ind,i]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind]
  gender[[i]]<-x$gender[e_ind]-1
  imd[[i]]<-x$imd[e_ind] #
  ethnicity[[i]]<-x$ethnicity[e_ind]
}
e_ind<-which((is.na(pmin(x[,1],x[,2],x[,3],x[,4],x[,5],na.rm = T))) & (is.na(x$dod)==F))
transtime[[i+1]]<-as.numeric(x$dod[e_ind]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind]
gender[[i+1]]<-x$gender[e_ind]-1
imd[[i+1]]<-x$imd[e_ind] #
ethnicity[[i+1]]<-x$ethnicity[e_ind]
e_ind<-which((is.na(pmin(x[,1],x[,2],x[,3],x[,4],x[,5],x$dod,na.rm = T))))
transtime[[i+2]]<-as.numeric(x$cohort_exit_date[e_ind]-x$cohort_entry_date[e_ind])/365.25+x$age_entry[e_ind]
gender[[i+2]]<-x$gender[e_ind]-1
imd[[i+2]]<-x$imd[e_ind] #
ethnicity[[i+2]]<-x$ethnicity[e_ind]
sumy<-list(transtime=transtime,gender=gender,imd=imd,ethnicity=ethnicity,n=length(unlist(transtime)))

# fitting using flexsurv (knot point selection guided by aic)
smy<-sumy
i<-5
mm_trans<-data.frame(id=1:smy$n,years=c(smy$transtime[[i]],smy$transtime[[setdiff((1:N),i)[1]]],smy$transtime[[setdiff((1:N),i)[2]]],smy$transtime[[setdiff((1:N),i)[3]]],smy$transtime[[setdiff((1:N),i)[4]]],smy$transtime[[6]],smy$transtime[[7]]),gender=relevel(factor(c(smy$gender[[i]],smy$gender[[setdiff((1:N),i)[1]]],smy$gender[[setdiff((1:N),i)[2]]],smy$gender[[setdiff((1:N),i)[3]]],smy$gender[[setdiff((1:N),i)[4]]],smy$gender[[6]],smy$gender[[7]])),ref="1"),imd=factor(c(smy$imd[[i]],smy$imd[[setdiff((1:N),i)[1]]],smy$imd[[setdiff((1:N),i)[2]]],smy$imd[[setdiff((1:N),i)[3]]],smy$imd[[setdiff((1:N),i)[4]]],smy$imd[[6]],smy$imd[[7]])),ethnicity=relevel(factor(c(smy$ethnicity[[i]],smy$ethnicity[[setdiff((1:N),i)[1]]],smy$ethnicity[[setdiff((1:N),i)[2]]],smy$ethnicity[[setdiff((1:N),i)[3]]],smy$ethnicity[[setdiff((1:N),i)[4]]],smy$ethnicity[[6]],smy$ethnicity[[7]])),ref = "White"),status=c(rep(1,length(smy$transtime[[i]])),rep(0,smy$n-length(smy$transtime[[i]]))),trans=rep(1,smy$n))
aic<-NA
trans0<-vector("list")
for (k in 0:5) {
  trans0[[k+1]]<-flexsurvspline(Surv(years, status) ~ gender+imd+ethnicity, inits=c(gamma0=1,gamma1=0.1), subset = (trans == 1), data=mm_trans, k=k, scale="hazard")
  aic[k+1]<-trans0[[k+1]]$AIC
}
trans<-trans0[[which.min(aic)]]

# store results
list.save(trans, 'trans_05.rdata')