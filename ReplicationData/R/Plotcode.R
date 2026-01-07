#Code for ploting

medsea=1;
library(ggplot2)
library(fanplot)
########### fan plot###################;
  pdf('EUplot.pdf',width = 8, height = 6);plot(NULL, xlim = c(1, Tmax), ylim = c(0,1),xlab="Period",ylab="Herfindahl Index");fan(herfEU,roffset=0.5);dev.off()
  pdf('CNplot.pdf',width = 8, height = 6);plot(NULL, xlim = c(1, Tmax), ylim = c(0,1),xlab="Period",ylab="Herfindahl Index");fan(herfCN,roffset=0.5);dev.off()

  
########### median plot##################;
  EuHI.mean <- data.frame(x=1:Tmax,HI=apply(herfEU,2,mean),region="Europe")
  CnHI.mean <- data.frame(x=1:Tmax,HI=apply(herfCN,2,mean),region="China")
  MideastHI.mean <- data.frame(x=1:Tmax,HI=apply(herfMideast,2,mean),region="Middle East")
  IndiaHI.mean <- data.frame(x=1:Tmax,HI=apply(herfIndia,2,mean),region="India")
  SEAsiaHI.mean <- data.frame(x=1:Tmax,HI=apply(herfSEAsia,2,mean),region="Southeast Asia")
  EuropaXHI.mean <- data.frame(x=1:Tmax,HI=apply(herfEuropaX,2,mean),region="EuropaX")
  
  if(medsea==1){
    MedHI.mean <- data.frame(x=1:Tmax,HI=apply(herfMed,2,mean),region="Mediterranean Region")
    
    Herf_MedEach <- rbind(CnHI.mean,EuHI.mean,MedHI.mean)
  }else{
    Herf_MedEach <- rbind(CnHI.mean,EuHI.mean)
  }

  
  Herf_MedEach_Other <- rbind(MideastHI.mean,IndiaHI.mean,SEAsiaHI.mean,EuropaXHI.mean)
  
  plot2<-ggplot(data = Herf_MedEach, aes(x=x, y=HI)) + geom_path(aes(colour=factor(region)),size=1)+ylim(ymin=0,ymax=1)+xlab("Period")+ylab("Herfindahl Index")+theme_bw()+scale_x_continuous(breaks=seq(0,Tmax,Tmax/5))
  ggsave(plot2, file="meanEach_EUCN.pdf", width=8, height=4,scale=1,dpi=300)
  
  plot2<-ggplot(data = Herf_MedEach_Other, aes(x=x, y=HI)) + geom_path(aes(colour=factor(region)),size=1)+ylim(ymin=0,ymax=1)+xlab("Period")+ylab("Herfindahl Index")+theme_bw()+scale_x_continuous(breaks=seq(0,Tmax,Tmax/5))
  ggsave(plot2, file="meanEach_Others.pdf", width=8, height=4,scale=1,dpi=300)
  

  ############original pixel of Chinese empire##################;
  CNorg=as.data.frame(table(CNmax.pixel))
  colnames(CNorg)=c("id","Freq")
  CNorg$fid=as.numeric(as.character(CNorg$id))-1
  CNorg$id=NULL
  write.csv(CNorg,"ChinaOrigin.csv")



