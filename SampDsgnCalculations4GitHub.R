## Determination of Sampling Probabilities, extraction of sample, and calculations to show that sampling scheme is beneficial ##

#np = number of physicians
#npcp = number of PCPs
#nach = number of acute care hospitals
#nmg = number of medical groups
#nos = number of owner subsidiaries
#drs_in_multiple_cps = Proportion of physicians in multiple corporate parents
#drs_in_multiple_achs_or_mgs = Proportion of physicians in multiple hospitals or medical groups within corporate parent

#Set-up paths to input data and output data and files
library(foreign)
rsource='/Volumes/OMalleyJ/Dartmouth/Biostatistics/U19_DataCore/Code'
setwd(rsource)
datdir='../Data/' #Make these subdirectories if you want to use the same directory structure as me.
outdir='../Output/' 

#read hospital/practice input data and determine number of independent ach's and mg's
data=read.csv(paste(datdir,'jachmg-2017-02-27.csv',sep="")) #ach_or_mg: 0=hospital and 1=medical group
data=data[is.na(data$id),]
n1ach=sum(data$ach_or_mg==0) #Number of independent hospitals
n1mg=sum(data$ach_or_mg==1) #Number of independent medical groups

#Read in system-level data, clean data, and specify parameters
data=read.csv(paste(datdir,'jcps-joss-2017-02-26.csv',sep=""))
data$ska=as.numeric(data$ska)
N=5700 #Total number of surveys that can afford
mx1ach=100 #Number of indep hospitals to sample
mx1mg=800 #Number of independent medical groups to sample

#Scenario: 
scen=1 #1 = log1.9
       #2 = linear
       #3 = square root

#Useful computations on data to know quantities in advance
ncp=sum(!is.na(data$nos))
table(data$nos[1:ncp])
ind=is.na(data$nos)
data$nos[ind]=0

#Align owner subsidiary (OS) and corporate parent (CP) information on relevant records
oscp=matrix(0,nrow=nrow(data),ncol=3)
for (i in (ncp+1):nrow(data)) {
 ind=(data$id==data$parent_id[i])
 oscp[i,]=unlist(data[ind,c("nach","nmg","npcp")]) #Addition of CP info to OS
 oscp[ind,]=oscp[ind,]+unlist(data[i,c("nach","nmg","npcp")]) #Addition of OS info to CP
}
pertot=(data$nach+data$nmg)/(data$nach+data$nmg+oscp[,1]+oscp[,2])
hetfact=2*(data$nos>0)*(pertot^2+(1-pertot)^2)
data$othnach=oscp[,1]
data$othnmg=oscp[,2]
data$othnpcp=oscp[,3]
data$pertot=pertot
data$hetfact=hetfact

mxnach=max(data$nach)
mxnmg=max(data$nmg)
ind=(data$nach>1)
qnach=quantile(data$nach[ind],c(0.2,0.4,0.6,0.8))
ind=(data$nmg>1)
qnmg=quantile(data$nmg[ind],c(0.2,0.4,0.6,0.8))

#Utilization score evaluation: Constrained above by mxach and mxmg, constrained to have at least 1 ach and 1 mg if system or os contains either
ind=(data$nmg<data$othnmg)

mxach=10 #Maximum number of hospitals per system: 5, 10
mxmg=10 #Maximum number of medical groups per system: 5, 10
mxos=3 #Maximum number of owner subsidaries per system
data$unach=data$nach-data$othnach*(data$nos>0)
data$unmg=data$nmg-data$othnmg*(data$nos>0)
data$unpcp=data$npcp-data$othnpcp*(data$nos>0)
ind=(data$unpcp<0) #Fox for two CP's with insufficient number of PCPs
print(data[ind,])
data$unpcp[ind]=0
#In each scenario we require 10 surveys when the unit count is 177
if (scen==1) {
 ind=(data$unach>4)
 sach=data$unach
 sach[ind]=round(4+(log(data$unach[ind],1.9)-log(4,1.9)))
 ind=(data$unmg>4)
 smg=data$unmg
 smg[ind]=round(4+(log(data$unmg[ind],1.9)-log(4,1.9)))	
} else if (scen==2) {
 ind=(data$unach>4)
 sach=data$unach
 sc=(log(177,1.9)-log(4,1.9))/(177-4)
 sach[ind]=round(4+sc*(data$unach[ind]-4))
 ind=(data$unmg>4)
 smg=data$unmg
 smg[ind]=round(4+sc*(data$unmg[ind]-4))	
} else {
 ind=(data$unach>4)
 sach=data$unach
 sc=(log(177,1.9)-log(4,1.9))/(sqrt(177)-sqrt(4))
 sach[ind]=round(4+sc*(sqrt(data$unach[ind])-sqrt(4)))
 ind=(data$unmg>4)
 smg=data$unmg
 smg[ind]=round(4+sc*(sqrt(data$unmg[ind])-sqrt(4)))	 
}
ind=(data$nos>0)
sos=data$nos
sos[ind]=round(1+log(data$nos[ind],6))
score=1+sach+smg #Add 1 for corporate parent to reflect number of surveys used

#Report values of nach and nmg at which sach and smg increase
val=4:max(c(data$unach,data$umg))
lval=round(4+log(val,1.9)-log(4,1.9))
pos=match(4:10,lval)
complextrans=c(0:3,val[pos])
if (scen==2) {
 lval=round(4+sc*(val-4)); pos=match(4:10,lval); c(0:3,val[pos]) #check linear
} else if (scen==3) {
 lval=round(4+sc*(sqrt(val)-sqrt(4))); pos=match(4:10,lval); c(0:3,val[pos]) #check square root
}
  
#Desirability (sampling) score evaluation
ind=!is.na(data$parent_id) #Compute penalty for having an OS
cpospen=data$nos
cpospen[ind]=data$nos[data$parent_id[ind]]
ncomp=(1+data$unach+data$unmg)/(1+ifelse(cpospen>0,1,0)) #Sum organizations but halve prob if system includes an OS
denom=sum(ncomp)
prob=ncomp/denom #Just for systems with cps: Prob of inclusion on a given draw

#Compute expected number of surveys sampled per system
Nr=N-mx1ach-mx1mg #Remaining surveys

#Simulation to make sampling plots reported in manuscript
p=seq(0,1,0.1) #Response rate of individual unit (range from 0 to 1); assume the same for all units
np=length(p)
DsgnEffProbs=matrix(0,nrow=np,ncol=18)
prand=Nr/(nrow(data)+sum(data$unach)+sum(data$unmg)) #Probability under completely random sampling
for (j in 1:length(p)) {
 #Evaluate sample inclusion indicators and probabilities
 probsamp=rep(0,length(prob))
 totnsampOS=rep(0,length(prob))
 ge1nsampOS=rep(0,length(prob))
 prsampOS=rep(0,length(prob))
 survall=rep(0,length(prob))
 mxsurvall=rep(0,length(prob))
 prsurvall=rep(0,length(prob))
 nsim=10000 #Number of simulations (10000 is plenty to overcome Monte Carlo error)
 for (sim in 1:nsim) {
  indsamp=rep(0,length(prob))
	
  #Draw sample
  sys=data$id
  samp=sample(sys,replace=FALSE,prob=prob) #Order is random
 
  #Order CP's so that they always occur before OS's
  nsamp=samp
  n=length(samp)
  pos=0
  for (i in 1:n) {
   pos=pos+1 #position in nsamp
   if (samp[i]>ncp) { #Insert CP immediately before OS
    nsamp[pos:(pos+1)]=c(data$parent_id[samp[i]],samp[i])
    pos=pos+1
   } else if (data$nos[samp[i]]>0){ #Insert a randomly selected OS immediately after CP
    ind=(data$parent_id==samp[i] & data$id>ncp) #Selects on OS's
    if (sum(ind)>1) {
  	 sampos=samp[(i+1):n] #Remaining organizations
 	 potmatch=match(data$id[ind],sampos)
 	 if (sum(!is.na(potmatch))>0) {
 	  rsamp=sampos[min(match(data$id[ind],sampos),na.rm=TRUE)] #Take first occuring OS of system in list or organizations remaining to sample
 	  nsamp[pos:(pos+1)]=c(samp[i],rsamp)
 	  pos=pos+1
 	 } else {
 	  nsamp[pos]=samp[i]
 	 }
    } else {
 	 rsamp=data$id[ind]
 	 nsamp[pos:(pos+1)]=c(samp[i],rsamp)
     pos=pos+1
    }
   } else {
    nsamp[pos]=samp[i]
   }
  }
  nsamp=unique(nsamp) #Remove duplicates
 
  #Bisection method to find sample size
  low=1; high=nrow(data);
  f0=sum(score[nsamp[1:low]])-Nr
  f1=sum(score[nsamp[1:high]])-Nr
  cont=ifelse(f0<0 & f1>0,1,0)
  while (cont==1) {
   mid=round((low+high)/2)
   fmid=sum(score[nsamp[1:mid]])-Nr
   cont=(high-low>1)
   if (fmid<0) {
    low=mid
    f0=mid
   } else {
    high=mid
    f1=fmid
   }
  }

  #Generate sampling probability as long run overage over indicator of whether unit is sampled
  probsamp[nsamp[1:mid]]=probsamp[nsamp[1:mid]]+1
 
  #Sampling design properties: Completeness of CP over OS (these and the following quantities are also evaluated as long-run averages)
  # Number of OS's sampled by CP's
  indsamp[nsamp[1:mid]]=1
  sampandOS=p[j]*indsamp*(1-is.na(data$parent_id)) #Don't want to count CPs themselves
  nsampOS=tapply(sampandOS,data$parent_id,sum)
  Notreturn=1-sampandOS
  nsamp1OS=1-tapply(Notreturn,data$parent_id,prod) 
  idCP=as.numeric(names(nsampOS)) #ID's of corporate parent
  totnsampOS[idCP]=totnsampOS[idCP]+nsampOS #Number of OS per complex system
  ge1nsampOS[idCP]=ge1nsampOS[idCP]+ifelse(nsampOS>0,1,0) #Any OS per complex system
  prsampOS[idCP]=prsampOS[idCP]+nsampOS/data$nos[idCP] #Proportion of OS per complex system
  avesampOS=sum(nsampOS)/sum(data$nos[data$nos>0])
  aveprsampOS=sum(nsamp1OS)/length(data$nos[data$nos>0]) #Mean proportion OS sampled and returned

  #Sampling design properties: Completeness of CP and OS over Hosp, Pract (response rate is 100p%)
  indachmg=(sach>0 & smg>0) #Indicates systems or OS with both ach's and mg's
  survall=indachmg*indsamp*p[j]*(1-(1-p[j])^sach)*(1-(1-p[j])^smg)
  prsurvall=sum(survall)/sum(indachmg)
 
  #The above stratified by CP and OS
  indCP=is.na(data$parent_id)
  prsurvstrat=c(sum(survall[indCP])/sum(indachmg[indCP]),sum(survall[!indCP])/sum(indachmg[!indCP]))
  
  #Occurrence of complete CP and at least one complete OS 
  indall=(data$nos>0 & sach>0 & smg>0) #Just indicates complex systems
  u=unique(data$parent_id[!indCP])
  probsOS=rep(0,sum(indCP))
  survOS=survall[!indCP]
  for (uid in 1:length(u)) {
   ind=(data$parent_id[!indCP]==u[uid])
   probsOS[u[uid]]=1-prod(1-survOS[ind])
  }
  probsTot=survall[indCP]*probsOS
  prsurvCP1OS=sum(probsTot[indall])/sum(indall)
  
  #Calculation of probabilities under completely random sampling
  pind=prand*p[j]
  avesampOSr=pind
  aveprsampOSr=sum(1-(1-pind*(data$nos>0))^data$nos)/length(data$nos[data$nos>0]) #Mean proportion OS sampled
  survallr=indachmg*pind*(1-(1-pind)^data$unach)*(1-(1-pind)^data$unmg)
  prsurvallr=sum(survallr)/sum(indachmg)
  prsurvstratr=c(sum(survallr[indCP])/sum(indachmg[indCP]),sum(survallr[!indCP])/sum(indachmg[!indCP]))
  probsOSr=rep(0,sum(indCP))
  survOSr=survallr[!indCP]
  for (uid in 1:length(u)) {
   ind=(data$parent_id[!indCP]==u[uid])
   probsOSr[u[uid]]=1-prod(1-survOSr[ind])
  }
  probsTotr=survallr[indCP]*probsOSr
  prsurvCP1OSr=sum(probsTotr[indall])/sum(indall)
  
  #Calculation of probabilities holding number of systems sampled and allocation per system fixed (treat
  # whether sample ach and mg at the system as independent events so that if sample system don't necessarily
  # sample ach's or mg's)
  pindc=mean(indsamp)*p[j]
  avesampOSc=pindc
  aveprsampOSc=sum(1-(1-pindc*(data$nos>0))^data$nos)/length(data$nos[data$nos>0]) #Mean proportion OS sampled
  prsurvallc=pindc^3
  prsurvstratc=c(pindc^3,pindc^3)
  probsTotc=prsurvallc*(1-(1-prsurvallc)^data$nos)
  prsurvCP1OSc=sum(probsTotc[indall])/sum(indall)
  DsgnEffProbs[j,]=DsgnEffProbs[j,]+c(avesampOS,aveprsampOS,prsurvall,prsurvstrat,prsurvCP1OS,avesampOSr,aveprsampOSr,prsurvallr,prsurvstratr,prsurvCP1OSr,avesampOSc,aveprsampOSc,prsurvallc,prsurvstratc,prsurvCP1OSc)  
 }
 probsamp=probsamp/nsim
 totnsampOS=totnsampOS/nsim
 ge1nsampOS=ge1nsampOS/nsim
 prsampOS=prsampOS/nsim
}
DsgnEffProbs=DsgnEffProbs/nsim #Dividing by nsim determines the long-run averages
DsgnEffProbs[,1:6]/DsgnEffProbs[,7:12]
DsgnEffProbs[,1:6]/DsgnEffProbs[,13:18]

## Make and save plots showing comparative performance of coupled sampling ##

titles=c("Overall proportion OS surveys returned","Probability of receiving any OS survey","Completeness: corporate parent or owner subsidiary","Completeness: corporate parent","Completeness: owner subsidiary","Completeness: corporate parent and owner subsidiary")
for (i in 1:6) {
 par(mfrow=c(1,1), srt=0, mai=c(0.4, 0.4, 0.3, 0.1), mgp=c(2,1,0))
 plot(p,DsgnEffProbs[,i],main=titles[i],xlab='Response Probability',ylab='Probability of sampling event',col=1,type='l')
 lines(p,DsgnEffProbs[,i+6],type='l',lty=1,col=2)
 lines(p,DsgnEffProbs[,i+2*6],type='l',lty=1,col=3)
 legend(0,max(DsgnEffProbs[,i]),legend=c("Cooperative sampling","Completely random sampling","Uncooperative surveys with same allocation"),lty=c(1,1,1),col=c(1,2,3),bty="n")
 dev.copy2pdf(file=paste(outdir,paste(paste('SampleEventProb',i,sep=''),'pdf',sep='.'),sep=''), width=6, height=6)
}
par(mfrow=c(2,1), srt=0, mai=c(0.7, 0.7, 0.3, 0.1), mgp=c(2,1,0))
wantplots=c(2,1,3,6,4,5); filenames=c('Owner-subsidaryCapture.pdf','HospitalPracticeCapture.pdf','CPvsOSCompleteness.pdf')
for (i in 1:length(wantplots)) {
 j=wantplots[i]
 plot(p,DsgnEffProbs[,j],main=titles[j],xlab='Response Probability',ylab='Probability of sampling event',col=1,type='l')
 lines(p,DsgnEffProbs[,j+6],type='l',lty=1,col=2)
 lines(p,DsgnEffProbs[,j+2*6],type='l',lty=1,col=3)
 legend(0,max(DsgnEffProbs[,j]),legend=c("Cooperative sampling","Completely random sampling","Uncooperative surveys with same allocation"),lty=c(1,1,1),col=c(1,2,3),bty="n")
 if (i==2 || i==4 || i==6) {
  dev.copy2pdf(file=paste(outdir,filenames[i/2],sep=''), width=6, height=6)
 }
}

## Finalize calculations and write output to files ##

#Output sampling probabilities: Check that the orders of data and probsamp are synchronized
data$sampscore=prob
data$smg=smg
data$sach=sach
data$score=score
data$probsamp=probsamp
data$totnsampOS=totnsampOS
data$ge1nsampOS=ge1nsampOS #Prob that at least one OS is sampled at CP
data$prsampOS=prsampOS #Mean of proportion of OS sampled
tfile=paste(datdir,paste(paste('SampDsgnPropprob',scen,sep=""),'.csv',sep=""),sep="")
write.csv(data,file=tfile,row.names=FALSE)

#Report on sampling design from CP perspective
datacp=data[1:ncp,]
scorecp=score[1:ncp]
probsampcp=probsamp[1:ncp]
hasos=(datacp$nos>0)
dataos=data[(ncp+1):nrow(data),]
scoreos=score[(ncp+1):nrow(data)]
probsampos=probsamp[(ncp+1):nrow(data)]

#Summary statistics for Desirability score
orgtype=c(rep(1,1121)+hasos,rep(0,390))
tapply(ncomp,orgtype,'mean')
 
#OS utilization
data$parent_id[is.na(data$parent_id)]=data$id[is.na(data$parent_id)] #To allow summaries by CP
score4os=tapply(score*(data$id>ncp),data$parent_id,sum) #Total score of owner subsidiaries per system
n4os=tapply(probsamp*(data$id>ncp),data$parent_id,sum) #Expected number of owner subsidiaries sampled per system
surv4os=tapply(probsamp*score*(data$id>ncp),data$parent_id,sum) #Expected total surveys used over owner subsidiaries within system

#SKA utilization for system
data$ska=as.numeric(data$ska)
n4ska=tapply(data$ska,data$parent_id,sum) #Total number of SKAs per system
nunitsska=tapply(data$ska*(data$nach+data$nmg),data$parent_id,sum) #Total number of SKA units per system
expnska=tapply(probsamp*data$ska,data$parent_id,sum) #Expected number of SKAs sampled per system
score4ska=tapply(score*data$ska,data$parent_id,sum)
expsska=tapply(score*probsamp*data$ska,data$parent_id,sum) #Expected number of surveys per system

#Total utilization of systems
sachtot=tapply(sach,data$parent_id,sum) #Summed ach score for whole system
smgtot=tapply(smg,data$parent_id,sum) #Summed mg score for whole system
scoretot=tapply(score,data$parent_id,sum) #Summed score for whole system
ntot=tapply(probsamp,data$parent_id,sum) #Expected number of cp's or owner subsidiaries sampled per system
survtot=tapply(probsamp*(score-1),data$parent_id,sum) #Expected total surveys used for units within system

#How many systems of each type are sampled by the big categories
nhosp=c(1,0,mean(datacp$nach[!hasos]),mean(datacp$nach[hasos]))
nmg=c(0,1,mean(datacp$nmg[!hasos]),mean(datacp$nmg[hasos]))
nos=c(0,0,mean(datacp$nos[!hasos]),mean(datacp$nos[hasos]))
totunits=c(n1ach,n1mg,sum(datacp$nach[!hasos]+datacp$nmg[!hasos]),sum(datacp$nach[hasos]+datacp$nmg[hasos]))
npcps=c(NA,NA,sum(datacp$npcp[!hasos]),sum(datacp$npcp[hasos]))/totunits #PCPs per unit
mnscoretot=c(1,1,mean(scoretot[!hasos]),mean(scoretot[hasos])) #Total number of surveys used
freq=c(n1ach,n1mg,length(probsampcp[!hasos]),length(probsampcp[hasos])) #Number of systems
nsys=c(mx1ach,mx1mg,sum(probsampcp[!hasos]),sum(probsampcp[hasos])) #Number of systems expected to be sampled
nsysos=c(0,0,sum(n4os[!hasos]),sum(n4os[hasos])) #Number of systems expected to be sampled
persys=100*nsys/freq #Proportion of systems expected to be sampled
freqos=c(0,0,0,sum(data$nos[hasos]))
peros=100*nsysos/c(1,1,1,freqos[4])
numska=c(0,0,sum(n4ska[!hasos]),sum(n4ska[hasos]))
expska=c(0,0,sum(expnska[!hasos]),sum(expnska[hasos]))
perska=100*expska/c(1,1,numska[3:4])
nunits=c(0,0,sum(nunitsska[!hasos]),sum(nunitsska[hasos]))
expsurvska=c(0,0,sum(expsska[!hasos]),sum(expsska[hasos]))
persurvska=100*expsurvska/c(1,1,nunits[3:4])
expsurvtot=c(mx1ach,mx1mg,sum(survtot[!hasos]),sum(survtot[hasos])) #Proportion of units expected to be sampled (including cp)
tab_all=cbind(systype=seq(1,4),nhosp,nmg,nos,npcps,freq,nsys,persys,freqos,nsysos,peros,totunits,mnscoretot,expsurvtot,perunits=100*expsurvtot/totunits)
tab_all[,c(2:5,7:11,13:15)]=round(tab_all[,c(2:5,7:11,13:15)],digits=1)
tfile=paste(outdir,paste(paste('DesignTableAlln',scen,sep=""),'.csv',sep=""),sep="")
write.csv(tab_all,file=tfile,row.names=FALSE)

#Separate proportion of units by hospitals and practices
totunitsi=c(sum(datacp$nach[!hasos]),sum(datacp$nmg[!hasos]),sum(datacp$nach[hasos]),sum(datacp$nmg[hasos]))
survtotach=tapply(probsamp*sach,data$parent_id,sum) #Expected total surveys used for hospitals
survtotmg=tapply(probsamp*smg,data$parent_id,sum) #Expected total surveys used for practices
expsurvtoti=c(sum(survtotach[!hasos]),sum(survtotmg[!hasos]),sum(survtotach[hasos]),sum(survtotmg[hasos]))
perunitsi=100*expsurvtoti/totunitsi
100*expsurvtot/totunits

#Systems with corporate parents stratified by deciles of complexity score
type=10*round(sach[1:ncp])+round(smg[1:ncp]) #Type of cp system (same as for OS below if nos=0, otherwise quintiles)
type[hasos]=100
intwidth=0.1
q=quantile(scoretot[hasos],seq(0,1,intwidth))
for (i in 1:(length(q)-1)) {
 type[hasos]=type[hasos]+ifelse(scoretot[hasos]>=q[i],1,0)
}
nhosp=tapply(datacp$nach,type,mean)
nmg=tapply(datacp$nmg,type,mean)
nos=tapply(datacp$nos,type,mean)
totunits=tapply(datacp$nach+datacp$nmg,type,sum)
npcps=tapply(datacp$npcp,type,sum)/totunits
freq=tapply(probsampcp,type,length)
nsys=tapply(probsampcp,type,sum) #Number of systems expected to be sampled
persys=100*nsys/freq #Proportion of systems expected to be sampled
freqos=tapply(datacp$nos,type,sum)
nsysos=tapply(n4os,type,sum)
persysos=100*nsysos/(1*(nsysos==0)+freqos) #Proportion of os expected to be sampled
nsystot=tapply(ntot,type,sum)
mnscoretot=tapply(scoretot,type,mean)
expsurvtot=tapply(survtot,type,sum)
numska=tapply(datacp$ska,type,sum)
expska=tapply(probsampcp*datacp$ska,type,sum)
perska=100*expska/numska
tab_type_cp=cbind(cptype=sort(unique(type)),nhosp,nmg,nos,npcps,freq,nsys,persys,freqos,nsysos,persysos,totunits,mnscoretot,expsurvtot,perunits=100*expsurvtot/totunits)
tab_type_cp[,c(2:5,7:8,10:15)]=round(tab_type_cp[,c(2:5,7:8,10:15)],digits=1)
ind=(tab_type_cp[,1]>=100)
tab_type_cp=tab_type_cp[ind,] #Restriction to just complex systems
tfile=paste(outdir,paste(paste('DesignTableCPn',scen,sep=""),'.csv',sep=""),sep="")
write.csv(tab_type_cp,file=tfile,row.names=FALSE)

#Report on pooled CPs and OS
#Stratify by score for number of ach
type=sach
nhosp=tapply(data$unach,type,mean)
nmg=tapply(data$unmg,type,mean)
nos=tapply(data$nos,type,mean)
totunits=tapply(data$unach+data$unmg,type,sum)
npcps=tapply(data$unpcp,type,sum)/totunits
mnscore=tapply(score,type,mean) #Average number of surveys used
freq=tapply(probsamp,type,length)
nsys=tapply(probsamp,type,sum) #Number of systems expected to be sampled
persys=100*nsys/freq #Proportion of systems expected to be sampled
expsurv=tapply(probsamp*(score-1),type,sum) #Proportion of hospitals and practices expected to be sampled
numska=tapply(data$ska,type,sum)
expska=tapply(data$probsamp*data$ska,type,sum)
perska=100*expska/numska
tab_type_cpos=cbind(type=sort(unique(type)),nhosp,nmg,nos,npcps,freq,nsys,persys,totunits,mnscore,expsurv,perunits=100*expsurv/totunits)
tab_type_cpos[,c(2:5,7:8,10:12)]=round(tab_type_cpos[,c(2:5,7:8,10:12)],digits=1)
tfile=paste(outdir,paste(paste('DesignTableCPOS_ACH',scen,sep=""),'.csv',sep=""),sep="")
write.csv(tab_type_cpos,file=tfile,row.names=FALSE)

#Distribution of number of hospitals and practices per corporate owner
par(mfrow=c(2,1), srt=0, mai=c(0.4, 0.6, 0.3, 0.1), mgp=c(2,1,0))
hist(data$nach,breaks=c(seq(0,50),60,300),main='Number of hospitals per corporate owner',xlab='Number',ylab='Count',xlim=c(0,30),right=FALSE,freq=TRUE,plot=TRUE)
hist(data$nmg,breaks=c(seq(0,50),60,100,300),main='Number of practices per corporate owner',xlab='Number',ylab='Count',right=FALSE,freq=TRUE,xlim=c(0,30),plot=TRUE)
dev.copy2pdf(file=paste(outdir,'NHospPractDist.pdf',sep=''), width=6, height=6)

#Stratify by score for number of mg
type=smg
nhosp=tapply(data$unach,type,mean)
nmg=tapply(data$unmg,type,mean)
nos=tapply(data$nos,type,mean)
totunits=tapply(data$unach+data$unmg,type,sum)
npcps=tapply(data$unpcp,type,sum)/totunits
mnscore=tapply(score,type,mean) #Average number of surveys used
freq=tapply(probsamp,type,length)
nsys=tapply(probsamp,type,sum) #Number of systems expected to be sampled
persys=100*nsys/freq #Proportion of systems expected to be sampled
expsurv=tapply(probsamp*score,type,sum) #Proportion of units expected to be sampled
numska=tapply(data$ska,type,sum)
expska=tapply(data$probsamp*data$ska,type,sum)
perska=100*expska/numska
tab_type_cpos=cbind(type=sort(unique(type)),nhosp,nmg,nos,npcps,freq,nsys,persys,totunits,mnscore,expsurv,perunits=100*expsurv/totunits)
tab_type_cpos[,c(2:5,7:8,10:12)]=round(tab_type_cpos[,c(2:5,7:8,10:12)],digits=1)
tfile=paste(outdir,paste(paste('DesignTableCPOS_MG',scen,sep=""),'.csv',sep=""),sep="")
write.csv(tab_type_cpos,file=tfile,row.names=FALSE)

#Stratify by both ach and mg
type=0*c(rep(1,ncp),rep(0,(nrow(data)-ncp)))+100*round(sach)+round(smg) #Type of cp or os (if make first value = 10000 get separate tables)
nhosp=tapply(data$unach,type,mean)
nmg=tapply(data$unmg,type,mean)
nos=tapply(data$nos,type,mean)
totunits=tapply(data$unach+data$unmg,type,sum)
npcps=tapply(data$unpcp,type,sum)/totunits
mnscore=tapply(score,type,mean) #Average number of surveys used
freq=tapply(probsamp,type,length)
nsys=tapply(probsamp,type,sum) #Number of systems expected to be sampled
persys=100*nsys/freq #Proportion of systems expected to be sampled
expsurv=tapply(probsamp*score,type,sum) #Proportion of units expected to be sampled
numska=tapply(data$ska,type,sum)
expska=tapply(data$probsamp*data$ska,type,sum)
perska=100*expska/numska
tab_type_cpos=cbind(type=sort(unique(type)),nhosp,nmg,nos,npcps,freq,nsys,persys,totunits,mnscore,expsurv,perunits=100*expsurv/totunits)
tab_type_cpos[,c(2:5,7:8,10:12)]=round(tab_type_cpos[,c(2:5,7:8,10:12)],digits=1)
tfile=paste(outdir,paste(paste('DesignTableCPOSn',scen,sep=""),'.csv',sep=""),sep="")
write.csv(tab_type_cpos,file=tfile,row.names=FALSE)

#Output expected number of surveys spent
tfile=paste(outdir,'ExpSurvn.csv',sep="")
write.csv(survtot,file=tfile,row.names=FALSE)


