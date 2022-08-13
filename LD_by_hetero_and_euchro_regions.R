install.packages('NAM',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(NAM)
install.packages('tidyr',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(tidyr)
install.packages('parallel',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(parallel)
install.packages('MASS',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(MASS)

load("Geno_ok.Rdata")
SNP_name<- colnames(gen)
a<- data.frame(colnames(gen))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))
colnames(info)[2]<-"chromosome"

info<- as.data.frame(cbind(SNP_name, info$chr, info$Position))
colnames(info)<- c("SNP","Chr","Position")
info$Chr<-as.numeric(as.character(info$Chr))

library(NAM)
range<-read.csv("HetEuregionChr.csv")
#Chr1<-c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")
b_numbers<-seq(1,20)

info$Position<-as.numeric(as.character(info$Position))
library(NAM)

funct<- function(b){
  dat<-subset(info,Chr==b)
  labs<-cut(dat$Position , breaks=c(range[b,2],range[b,3],range[b,4],range[b,5],  range[b,6]),
            labels<-c("Het","Eu","Het","Eu"))
  dat$Position<-labs
  if (any(dat$Position== "Het")){
    dat1<-subset(dat,Position=="Het")
    SNP_list<-dat1$SNP
    col.num <- which(colnames(gen) %in% SNP_list)
    gentemp<-gen[,c(col.num)]
    LD<- LD(gentemp)
    final<-as.data.frame(as.table(LD))
    c<-rep("Het",nrow(final))
    final<-cbind(final, c)
  }else {
    print("none of them")}
}

results<- mclapply(b_numbers, funct, mc.cores =20)
save(results, file="resultsLDfinalHet.RData")

##########euchromatin
load("Geno_ok.Rdata")
SNP_name<- colnames(gen)
a<- data.frame(colnames(gen))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))
colnames(info)[2]<-"chromosome"

info<- as.data.frame(cbind(SNP_name, info$chr, info$Position))
colnames(info)<- c("SNP","Chr","Position")
info$Chr<-as.numeric(as.character(info$Chr))

library(NAM)
range<-read.csv("HetEuregionChr.csv")
#Chr1<-c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")
b_numbers<-seq(1,20)

info$Position<-as.numeric(as.character(info$Position))
library(NAM)

funct<- function(b){
  dat<-subset(info,Chr==b)
  labs<-cut(dat$Position , breaks=c(range[b,2],range[b,3],range[b,4],range[b,5],  range[b,6]),
            labels<-c("Het","Eu","Het","Eu"))
  dat$Position<-labs
  if (any(dat$Position== "Eu")){
    dat1<-subset(dat,Position=="Eu")
    SNP_list<-dat1$SNP
    col.num <- which(colnames(gen) %in% SNP_list)
    gentemp<-gen[,c(col.num)]
    LD<- LD(gentemp)
    final1<-as.data.frame(as.table(LD))
    d<-rep("Eu",nrow(final1))
    final1<-cbind(final1, d)
  }else {
    "none of them"}
}

results<- mclapply(b_numbers, funct, mc.cores =20)
save(results, file="resultsLDfinaEu.RData")



#########working with the results##################################################################################
load("resultsLDfinalEu.RData")

Eu<-data.frame(matrix(ncol=7,nrow=0))
for (i in 1:20){
  chr<-i
  a<-as.data.frame(results[[i]])
  a$Chr<- rep(i,nrow(a))
  SNP_name1<- as.data.frame(a$Var1)
  SNP_name2<- as.data.frame(a$Var2)
  colnames(SNP_name1)<-"name"
  colnames(SNP_name2)<-"name"
  
  info<- SNP_name1%>% separate(name,c("x","chr","Position", "y","z"))
  info<- as.data.frame(cbind(info$chr, info$Position))
  colnames(info)<- c("Chr","PositionSNP1")
  
  info1<- SNP_name2%>% separate(name,c("x","chr","Position", "y","z"))
  info1<- as.data.frame(cbind(info1$chr, info1$Position))
  colnames(info1)<- c("Chr","PositionSNP2")
  
  
  a$PositionSNP1<- info$PositionSNP1
  a$PositionSNP2<- info1$PositionSNP2
  
  Eu<-rbind(Eu,a)
  
}
rm(results)

##eucromatine
load("resultsLDfinalHet.RData")

Het<-data.frame(matrix(ncol=7,nrow=0))
for (i in 1:20){
  chr<-i
  a<-as.data.frame(results[[i]])
  a$Chr<- rep(i,nrow(a))
  SNP_name1<- as.data.frame(a$Var1)
  SNP_name2<- as.data.frame(a$Var2)
  colnames(SNP_name1)<-"name"
  colnames(SNP_name2)<-"name"
  
  info<- SNP_name1%>% separate(name,c("x","chr","Position", "y","z"))
  info<- as.data.frame(cbind(info$chr, info$Position))
  colnames(info)<- c("Chr","PositionSNP1")
  
  info1<- SNP_name2%>% separate(name,c("x","chr","Position", "y","z"))
  info1<- as.data.frame(cbind(info1$chr, info1$Position))
  colnames(info1)<- c("Chr","PositionSNP2")
  
  
  a$PositionSNP1<- info$PositionSNP1
  a$PositionSNP2<- info1$PositionSNP2
  
  Het<-rbind(Het,a)
  
}
rm(results)




Het$PositionSNP1<- as.numeric(as.character(Het$PositionSNP1))
Het$PositionSNP2<- as.numeric(as.character(Het$PositionSNP2))
Het$diff<-abs(Het$PositionSNP1- Het$PositionSNP2)
#Het<- subset(Het, diff!=0)
Het<- subset(Het, diff< 10000000) ##subsetting SNPs pairs that are < 10000000 bp
write.csv(Het, file="heterochromatine.csv")

Eu$PositionSNP1<- as.numeric(as.character(Eu$PositionSNP1))
Eu$PositionSNP2<- as.numeric(as.character(Eu$PositionSNP2))
Eu$diff<-abs(Eu$PositionSNP1- Eu$PositionSNP2)
#Eu<- subset(Eu, diff!=0)
Eu<- subset(Eu, diff< 10000000)
write.csv(Eu, file="Eurochromatine.csv")

###############obtaining averages 
Het<-read.csv("heterochromatine.csv")
Het<- Het[,-1]
dfr<- Het[,c(8,3)]
colnames(dfr)<-c("dist","rsq")
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=100000))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
write.csv(dfr1, file="average_LDhet.csv")

LDHetbyChr<- list()
for (i in 1:20){
  dat<-subset(Het,Chr== i)
  dfr<- dat[,c(8,3)]
  colnames(dfr)<-c("dist","rsq")
  dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=100000))
  dfr <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
  dfr$Chr<- rep(i,nrow(dfr))
  LDHetbyChr[[i]]<- dfr
  
}

save(LDHetbyChr, file="LDbyHetCHR.RData")

####################EU

#Eu<- Eu[,-1]
dfr<- Eu[,c(7,2)]
colnames(dfr)<-c("dist","rsq")
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=100000))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
write.csv(dfr1, file="average_LDEU.csv")

LDEubyChr<- list()
for (i in 1:20){
  dat<-subset(Eu,Chr== i)
  dfr<- dat[,c(7,2)]
  colnames(dfr)<-c("dist","rsq")
  dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=100000))
  dfr <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
  dfr$Chr<- rep(i,nrow(dfr))
  LDEubyChr[[i]]<- dfr
  
}

save(LDEubyChr, file="LDbyCHREU.RData")



#################LD decay for heterochomatin
##########################
load("LDbyCHREU.RData")

EuDecay<-data.frame(matrix(ncol=2,nrow=0)) 
for (i in 1:20){
  dfr1<- LDEubyChr[[i]]
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  dfr1$mid<- dfr1$mid/1000
  maxld<-max(dfr1$mean)
  #You could eleucubrate if it's better to use the maximum ESTIMATED value of LD
  #In that case you just set: maxld<-max(fpoints) 
  h.decay<-maxld/2
  half.decay.distance<-dfr1$mid[which.min(abs(dfr1$mean-h.decay))]
  decay<-cbind(i,half.decay.distance)
  decay<-as.data.frame(decay)
  colnames(decay)<-c("chr","Decay rate")
  EuDecay<-rbind(EuDecay,decay)
}


load("LDbyHetCHR.RData")

HetDecay<-data.frame(matrix(ncol=2,nrow=0)) 
for (i in 1:20){
  dfr1<- LDHetbyChr[[i]]
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  dfr1$mid<- dfr1$mid/1000
  maxld<-max(dfr1$mean)
  #You could eleucubrate if it's better to use the maximum ESTIMATED value of LD
  #In that case you just set: maxld<-max(fpoints) 
  h.decay<-maxld/2
  half.decay.distance<-dfr1$mid[which.min(abs(dfr1$mean-h.decay))]
  decay<-cbind(i,half.decay.distance)
  decay<-as.data.frame(decay)
  colnames(decay)<-c("chr","Decay rate")
  HetDecay<-rbind(HetDecay,decay)
}



final_decay<-merge(EuDecay,HetDecay,by="chr")
write.csv(final_decay, file="LDdecay_chr.csv")  

dfr1Eu$LD.data<-as.numeric(as.character(dfr1Eu$LD.data))
dfr1Eu$distance<-as.numeric(as.character(dfr1Eu$distance))
maxld<-max(dfr1Eu$LD.data)

#You could eleucubrate if it's better to use the maximum ESTIMATED value of LD
#In that case you just set: maxld<-max(fpoints) 
h.decay<-maxld/2

half.decay.distance<-dfr1Eu$distance[which.min(abs(dfr1Eu$LD.data-h.decay))]
mean(dfr1Eu$LD.data)
#####
#################LD decay for hetero
distance<-dfr1$mid
LD.data<-dfr1$mean
n<-323
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=200))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))

###########################

df1<-data.frame(distance,fpoints)
maxld1<-max(LD.data)
#You could eleucubrate if it's better to use the maximum ESTIMATED value of LD
#In that case you just set: maxld<-max(fpoints) 
h.decay1<-maxld1/2
half.decay.distance1<-df$distance[which.min(abs(df1$fpoints-h.dec))]


#Chro<-c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20") ##chromosomes per region
  info$Chr<-as.numeric(as.character(info$Chr))
sum_regions<- data.frame(matrix(ncol=4,nrow=0))
for (i in 1:20) {
  dat<-subset(info,Chr==i)
  labs<-cut(dat$Position , breaks=c(range[i,2],range[i,3],range[i,4],range[i,5],  range[i,6]),
            labels<-c("Het","Eu","Het","Eu"))
  dat$region<-labs
  sum_regions<-rbind(sum_regions,dat)
}

write.csv(sum_regions, file="SNPsheteucreg.csv")
library(dplyr)
stats<-aggregate(sum_regions$SNP, by=list(sum_regions$Chr,sum_regions$region),FUN=length)
write.csv(stats, file="heteuc%byChr.csv")


######marker density ###########
SNPsChr<-data.frame(matrix(ncol=3,nrow=0))
for (i in 1:20) {
  dat<-subset(info,Chr==i)
  labs<-cut(dat$Position , breaks=c(range[i,2],range[i,3],range[i,4],range[i,5],  range[i,6]),
            labels<-c("Het","Eu","Het","Eu"))
  dat$region<-labs
  Euc<-subset(dat, region=="Eu")
  Het<-subset(dat, region=="Het")
  len<-max(Euc$Position)
  len1<-max(Het$Position)
  mark<-nrow(Euc)
  mark1<-nrow(Het)
  dat1<-cbind(len,mark,len1,mark1,i)
  colnames(dat1)<-c("ChrLenEu","NumberMarkersEU","ChrLenHet","NumberMarkersHet","Chr")
  SNPsChr<-rbind(SNPsChr,dat1)
}
write.csv(SNPsChr, file="SNPs_density.csv")



###########

SNPsChr<-data.frame(matrix(ncol=3,nrow=0))
for (i in 1:20) {
  dat<-subset(info,Chr==i)
  len<-max(dat$Position)
  mark<-nrow(dat)
  dat1<-cbind(len,mark,i)
  colnames(dat1)<-c("ChrLen","NumberMarkers","Chr")
  SNPsChr<-rbind(SNPsChr,dat1)
}
write.csv(SNPsChr, file="SNPs_density.csv")




#split between heterochromatin and euchromatin################# from PLINK data

data3<-read.table('LD_dataset1.ld', header=T)
range<-read.csv("HetEuregionChr.csv")


info$Position<-as.numeric(as.character(info$Position))
library(NAM)

final<-as.data.frame(matrix(ncol=9,nrow = 0))
for (i in 1:20){
  dat<-subset(data3,CHR_A==i & CHR_B==i)
  labs<-cut(dat$BP_A , breaks=c(range[i,2],range[i,3],range[i,4],range[i,5],  range[i,6]),
            labels<-c("Het","Eu","Het","Eu"))
  dat$region_A<-labs
  labs1<-cut(dat$BP_B , breaks=c(range[i,2],range[i,3],range[i,4],range[i,5],  range[i,6]),
             labels<-c("Het","Eu","Het","Eu"))
  dat$region_B<-labs1
  final<-rbind(final,dat)
}


Het<-subset(final,final$region_A== "Het" & final$region_B=="Het")
Het$dist<- abs(Het$BP_A-Het$BP_B)
dfr<- Het[,c(10,7)]
colnames(dfr)<-c("dist","rsq")
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=10000))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

dfr1$region<- rep("Het",nrow(dfr1))
write.csv(dfr1,file="hetLD_plink.csv")
####################EU

Eu<-subset(final,final$region_A== "Eu" & final$region_B=="Eu")
Eu$dist<- abs(Eu$BP_A-Eu$BP_B)
dfr2<- Eu[,c(10,7)]
colnames(dfr2)<-c("dist","rsq")
dfr2$distc <- cut(dfr2$dist,breaks=seq(from=min(dfr2$dist),to=max(dfr2$dist)+1,by=10000))
dfr2 <- dfr2 %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))

dfr2 <- dfr2 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


dfr2$region<-rep("Eu", nrow(dfr2))
write.csv(dfr2,file="EuLD_plink.csv")


distance<-(dfr2$mid)/1000
LD.data<- dfr2$mean
plot(distance, LD.data, col="red", pch=20)
lines(distance, LD.data,col="Red", lwd=2)
region<- rep("Heterochromatin",nrow(dfr1))
dfr1Het<-cbind(distance,LD.data,region)
dfr1Het<-as.data.frame(dfr1Het)
dfr1Het$LD.data<-as.numeric(as.character(dfr1Het$LD.data))
plot(dfr1Het$distance, dfr1Het$LD.data, col="red", pch=20)


###the same but with eucromatine
dfr1<- read.csv("average_LDEU.csv")
library(stringr)
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


distance<-(dfr1$mid)/1000
LD.data<- dfr1$mean
plot(distance, LD.data, col="red", pch=20)
lines(distance, LD.data,col="Red", lwd=2)
region<- rep("Euchromatin",nrow(dfr1))
dfr1Eu<-cbind(distance,LD.data,region)
dfr1Eu<-as.data.frame(dfr1Eu)

final_dat<-rbind(dfr1Het,dfr1Eu)
mean(final_dat$LD.data)
final_dat$distance<-as.numeric(as.character(final_dat$distance))
final_dat$LD.data<-as.numeric(as.character(final_dat$LD.data))
#####plot

##########makeing  a line plot with two variables 

####
png("LD_plots_wholegenome.png", width = 900, height = 700)
ggplot(final_dat, aes(x=distance, y=LD.data, group=region, color=region)) +
  geom_point()+
  scale_color_brewer(palette = "Set1")+
  #geom_line( col="red", alpha=0.8)+
  labs(x="Distance (kb)",y=expression(LD~(r^{2}))) +
  theme_bw(base_size = 18)
dev.off()


#average by chromosome
LDHetbyChr<- list()
for (i in 1:20){
  dat<-subset(Het,Chr== i)
  dfr<- dat[,c(8,3)]
  colnames(dfr)<-c("dist","rsq")
  dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist),to=max(dfr$dist)+1,by=100000))
  dfr <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq, na.rm = TRUE),median=median(rsq,na.rm = TRUE))
  dfr$Chr<- rep(i,nrow(dfr))
  LDHetbyChr[[i]]<- dfr
  
}

save(LDHetbyChr, file="LDbyHetCHR.RData")

