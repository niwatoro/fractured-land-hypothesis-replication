#Main Code
#You should run the Start.R instead of this

medsea=1;

dir.create(folder)
outfolder=paste0(folder,roundI,"/");
dir.create(outfolder)

pixelno=nrow(EAdta)
EAdta$coastal[is.na(EAdta$coastal)] <- 0
pixel.coast=as.array(EAdta$coastal)

EAdta$medcoast[is.na(EAdta$medcoast)] <- 0
pixel.medcoast=as.array(EAdta$medcoast)


pixel.EU=which(as.array(EAdta$Europe)==1)
pixel.CN=which(as.array(EAdta$China)==1)
pixel.India=which(as.array(EAdta$indiancont)==1)
pixel.Mideast=which(as.array(EAdta$mideast)==1)
pixel.SoutheastAsia=which(as.array(EAdta$seasia)==1)
pixel.EuropaX=which(as.array(EAdta$Europe)==1)



pixel.AfricaET=which(as.array(EAdta$AfricaET)==1)
pixel.AfricaWT=which(as.array(EAdta$AfricaWT)==1)
pixel.AfricaME=which(as.array(EAdta$AfricaME)==1)
pixel.AfricaNH=which(as.array(EAdta$AfricaNH)==1)
pixel.AfricaSH=which(as.array(EAdta$AfricaSH)==1)

pixel.AmericaNH=which(as.array(EAdta$AmericaNH)==1)
pixel.AmericaSH=which(as.array(EAdta$AmericaSH)==1)
pixel.AmericaCL=which(as.array(EAdta$AmericaCL)==1)




if(medsea!=1){EAdta$medsea=0};
EAdta$medsea[is.na(EAdta$medsea)] <- 0
medsea.region=which(as.array(EAdta$medsea)==1)



pixel.ruggedness=as.array(EAdta$ElevationS);



EAdta$tmin[is.na(EAdta$tmin)] <- 0
pixel.temp=as.array(EAdta$tmin)*1;
pixel.cold=log(pmax(9-pixel.temp,1));
Theta_cold=Theta_cold/log(9-quantile(pixel.temp,0.1))


EAdta$tmax[is.na(EAdta$tmax)] <- 0
pixel.temp2=as.array(EAdta$tmax)*1;
pixel.hot=log(pmax(pixel.temp2-21,1));
Theta_hot=Theta_hot/log(quantile(pixel.temp2,0.9)-21)




############Theta######################;
Theta_rugged=1/quantile(pixel.ruggedness,0.9)*Theta_rugged_90th;
Theta_sea=0;

Theta0=1
if(Theta_rugged_90th==0 & Theta_cold==0 & Theta_hot==0){Theta0=0}


##############Y settings#################;

if(yielddta=="YKK10"){pixel.Y=as.array(EAdta$YKK10)}
if(yielddta=="YCSI"){pixel.Y=as.array(EAdta$YCSI)}
if(yielddta=="YGAEZ"){pixel.Y=as.array(EAdta$YGAEZ)}
if(yielddta=="YGAEZ4"){pixel.Y=as.array(EAdta$YGAEZ4)}


if(UniformY==1){
  pixel.Y=as.array(rep(0.5,pixelno))
}

pixel.Y[is.na(pixel.Y)]<-0

EAdta$startKK10[is.na(EAdta$startKK10)]<-0
EAdta$startDSMW[is.na(EAdta$startDSMW)]<-0
EAdta$startloess[is.na(EAdta$startloess)]<-0
EAdta$startHSWD[is.na(EAdta$startHSWD)]<-0



pixel.steppeeast=as.array(EAdta$steppeeast) 
id.steppeeast=which(pixel.steppeeast==1)

pixel.steppe=as.array(EAdta$steppe_all)
id.steppe=which(pixel.steppe==1)






if(river==1){
  EAdta$river[is.na(EAdta$river)] <- 0
  pixel.rivers=cbind(EAdta$river,EAdta$river)
  pixel.riverdummy=(apply(pixel.rivers,1,sum)>=1)*1
}

################Self-defined Functions########################
pixelno=nrow(EAdta)
nation.territory=as.list(seq(1:pixelno))

if(start_with_regime==1){
  regimeid=as.array(EAdta$regimeid)
  nation.territory=unstack(data.frame(seq(1:pixelno),regimeid))
}

nation.orgpixel=seq(1:pixelno)
nationno=length(nation.territory)

theta<-function(pixel){
  return(pmax(Theta_rugged*pixel.ruggedness[pixel]+Theta_cold*pixel.cold[pixel]+Theta_hot*pixel.hot[pixel]+Theta_steppe*pixel.steppe[pixel],0))
}




fid=seq(1:pixelno)-1


riverconnected<-function(pixel1,pixel2){
  if(length(pixel1)>1){
    connected=apply(pixel.rivers[pixel1,]*pixel.rivers[pixel2,],1,sum)>0
  }else{
    connected=sum(pixel.rivers[pixel1,]*pixel.rivers[pixel2,])>0
  }
  return(connected)
}


intersect0<-function(x,y) x[x %in% y]
setdiff0<-function(x,y) x[! x %in% y]


RdSel<-function(x) ifelse(length(x)>1,sample(x,1),x)

RdSelWar<-function(x) {
  if(length(x)>6){return(sample(x,6))}
  else{return(x)}
}

RdSort<-function(x){
  if(length(x)>1){return(sample(x))}
  else{return(x)}
}

duplicatedR<-function(x){
  rorder=RdSort(1:length(x))
  duplicated(x[rorder])[order(rorder)]
}

lengthNA<-function(x){
  if(is.na(x[1])){
    return(0)
  }else{
    return(length(x))
  }
}

sep.prob<-function(pixelvec){
  return(pmin(beta*theta(pixelvec)*nation.allborderL[(pixel.belong$nation.fid[unlist(pixelvec)])],1))
}



SplitDisconnect<-function(x){
  pixels=unlist(x)
  setsize=length(pixels)
  if(setsize==1){return(x)}else{
    y=pixels[1]
    yset=c(y)
    splitlist=NA
    addno=1
    renameset=pixels
    while(addno<=setsize){
      newmember=setdiff(intersect0(unlist(pixel.borderM[y]),pixels),yset)
      addno=addno+length(newmember)
      if(length(newmember)>0){
        y=newmember
        yset=c(yset,newmember)
      }else{
        if(is.na(splitlist[1])){
          splitlist=list(yset)
        }else{
          splitlist=c(splitlist,list(yset))
        }
        renameset=setdiff(renameset,yset)
        y=renameset[1]
        yset=c(y)
        addno=addno+length(y)
      }
    }
    return(splitlist)
  }
}


PixelLandConnected<-function(target,territory){
  territory=unlist(territory)
  connectedlist=intersect(target,territory)
  y=target
  while(length(y)>0){
    y_new=setdiff(intersect(unlist(pixel.borderM_land[y]),territory),connectedlist)
    connectedlist=c(connectedlist,y_new)
    y=y_new
  }
  return(connectedlist)
}



costY<-function(nationid,topixel){
  connectpixel=mapply(function(x,y) PixelLandConnected(x,y),topixel,nation.territory[nationid], SIMPLIFY=FALSE)
  connectY=sapply(connectpixel,function(x) sum(pixel.Y[x]))
  unconnectY=nation.resource[nationid]-connectY
  return(connectY+sigma*unconnectY)
}

costYs<-function(nationids,topixels){
  if(length(nationids)==1){
    nationids=rep(nationids,length(topixels))
  }
  return(mapply(costY,nationids,topixels))
}


drawing<-function(){
  map<-readOGR(dsn=paste0(mapfile,".shp"))
  
  pixel.belong=stack(setNames(nation.territory, seq_along(nation.territory)))
  colnames(pixel.belong)[1]="pixelid"
  colnames(pixel.belong)[2]="nation.fid"
  pixel.belong=pixel.belong[order(pixel.belong$pixelid),]
  
  mapU <- unionSpatialPolygons(map,pixel.belong$nation.fid)
  map.df <- fortify(mapU)
  map.df$id=as.numeric(map.df$id)
  
  nationsize=sapply(nation.territory,length)
  nationsize.rank=length(nationsize)-rank(nationsize,ties.method = "random")
  map.df$colorid=nationsize.rank[map.df$id] %% 9
  
  myColors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#636363')
  colScale <- scale_fill_manual(values = myColors)
  colScaleline<-scale_colour_manual(values = myColors)
  
  ggp <- ggplot()
  if(drawEurasia==1){
    ggp <- ggp + as_reference(geom_polygon(data=map.df, aes(x=long, y=lat,group=group,fill=as.factor(colorid)),show.legend=FALSE),id="polygon") + colScale
    ggp <- ggp + as_reference(geom_path(data=map.df, aes(x=long, y=lat,group=group)),id="path")+colScaleline
    ggp <- ggp + with_blend(geom_polygon(data=map_eurasia,aes(x = long, y = lat,group=group)),bg_layer="polygon",blend_type="in",flip_order=TRUE)
    ggp <- ggp + with_blend(geom_polygon(data=map_eurasia,aes(x = long, y = lat,group=group)),bg_layer="path",blend_type="in",flip_order=TRUE)
    ggp <- ggp + coord_sf(xlim=c(-1528000,13353993),ylim=c(138340.8,8043417),expand=FALSE)
  }else{
    ggp <- ggp + geom_polygon(data=map.df, aes(x=long, y=lat,group=group,fill=as.factor(colorid)),show.legend=FALSE) + colScale
    ggp <- ggp + geom_path(data=map.df, aes(x=long, y=lat,group=group))+colScaleline
  }
  ggp <- ggp + labs(title="")+theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  ggsave(ggp, file=paste(outfolder,"Map",timer,".jpeg",sep=""), width=8, height=4,scale=2.5,dpi=300)
  
}

############################################





Chi_cold=(1/(1-Ydiscount_cold)-1)/log(pmax(9+4,1))
Chi_hot=(1/(1-Ydiscount_hot)-1)/log(pmax(29.5-21,1))
pixel.Y=pixel.Y/(1+Chi_cold*pixel.cold+Chi_hot*pixel.hot)



pixel.Y.unrest=pixel.Y

alpha=ifelse(UniformY,1,1/quantile(pixel.Y.unrest[union(pixel.EU,pixel.CN)],0.95))#0.95

pixel.Y=pixel.Y+ymin


##################Main######################

timer=0;
Hihis=c();Hihis.EU=c();Hihis.CN=c();Hihis.Med=c();Hihis.Mideast=c();Hihis.India=c();Hihis.SEAsia=c();Hihis.EuropaX=c();
Empireno.CN=c();Empireno.EU=c();Empireno.India=c();Empireno.Mideast=c();Empireno.SoutheastAsia=c();Empireno.EuropaX=c()
Empireno.AfricaET=c();Empireno.AfricaWT=c();Empireno.AfricaME=c();Empireno.AfricaNH=c();Empireno.AfricaSH=c();
Empireno.AmericaNH=c();Empireno.AmericaSH=c();Empireno.AmericaCL=c();
regimeno.EU=c();
regimesize.EU=matrix(,nrow=10,ncol=Tmax);
regimeno.CN=c();
regimesize.CN=matrix(,nrow=10,ncol=Tmax);
regimeno.Med=c();
regimesize.Med=matrix(,nrow=10,ncol=Tmax);
regimeno.Mideast=c();
regimesize.Mideast=matrix(,nrow=10,ncol=Tmax);
regimeno.India=c();
regimesize.India=matrix(,nrow=10,ncol=Tmax);
regimeno.SEAsia=c();
regimesize.SEAsia=matrix(,nrow=10,ncol=Tmax);
regimeno.EuropaX=c();
regimesize.EuropaX=matrix(,nrow=10,ncol=Tmax);

regimeno.AfricaET=c();
regimesize.AfricaET=matrix(,nrow=10,ncol=Tmax);
regimeno.AfricaWT=c();
regimesize.AfricaWT=matrix(,nrow=10,ncol=Tmax);
regimeno.AfricaME=c();
regimesize.AfricaME=matrix(,nrow=10,ncol=Tmax);
regimeno.AfricaNH=c();
regimesize.AfricaNH=matrix(,nrow=10,ncol=Tmax);
regimeno.AfricaSH=c();
regimesize.AfricaSH=matrix(,nrow=10,ncol=Tmax);
regimeno.AmericaNH=c();
regimesize.AmericaNH=matrix(,nrow=10,ncol=Tmax);
regimeno.AmericaSH=c();
regimesize.AmericaSH=matrix(,nrow=10,ncol=Tmax);
regimeno.AmericaCL=c();
regimesize.AmericaCL=matrix(,nrow=10,ncol=Tmax);




IntersectCN=c();IntersectEU=c();UnionEU=c();UnionCN=c();

Maize1=which(EAdta$Maize1==1)
Maize2=which(EAdta$Maize2==1 & EAdta$Maize1!=1)
Maize3=which(EAdta$Maize3==1 & EAdta$Maize1!=1 & EAdta$Maize2!=1)

america=which(EAdta$Americas==1 & EAdta$Maize1!=1 & EAdta$Maize2!=1 & EAdta$Maize3!=1)
africa=which(EAdta$Africa==1)
australia=which(EAdta$Australia==1)
japan=which(EAdta$Japan==1)

startyr=rep(0,pixelno)
startyr[Maize1]=Maize1_Start
startyr[Maize2]=Maize2_Start
startyr[Maize3]=Maize3_Start

startyr[america]=America_Start
startyr[africa]=Africa_Start
startyr[australia]=Australia_Start
startyr[japan]=Japan_Start



while (TRUE){
  if(timer>=Tmax){
    break;
  }else if(nationno<=1){
    Hihis[timer:Tmax]=1;
    Hihis.EU[timer:Tmax]=1;
    Hihis.CN[timer:Tmax]=1;
    Hihis.Med[timer:Tmax]=1;
    Hihis.Mideast[timer:Tmax]=1;
    Hihis.India[timer:Tmax]=1;
    Hihis.SEAsia[timer:Tmax]=1;
    Hihis.EuropaX[timer:Tmax]=1;
    break;
  }
  timer=timer+1;
  print(paste0("run:",roundI," T:",timer));
  if(file.exists(paste0(mainfolder,"Output/","sim-",roundI,"-T-",timer-1))){
    file.rename(paste0(mainfolder,"Output/","sim-",roundI,"-T-",timer-1),paste0(outfolder,"sim-",roundI," T-",timer))
  }else{
    file.create(paste0(mainfolder,"Output/","sim-",roundI,"-T-",timer))
  }

  RatioY0=as.array(rep(1,pixelno));
  if(BC1000Ratio=="startKK10"){RatioY0=as.array(EAdta$startKK10)}
  if(BC1000Ratio=="startDSMW"){RatioY0=as.array(EAdta$startDSMW)}
  if(BC1000Ratio=="startloess"){RatioY0=as.array(EAdta$startloess)}
  if(BC1000Ratio=="startHSWD"){RatioY0=as.array(EAdta$startHSWD)}
  

  if(timer==1){
    pixel.Y.500=pixel.Y;
    pixel.Y.unrest.500=pixel.Y.unrest;
    increment=pixel.Y*(1-RatioY0)/500

    
    pixel.Y=pixel.Y*RatioY0
    pixel.Y=pmax(pixel.Y,0.000001)
    increment2=pixel.Y.unrest*(1-RatioY0)/500
    
    pixel.Y.unrest=pixel.Y.unrest*RatioY0
    pixel.Y.unrest=pmax(pixel.Y.unrest,0.000001)
    
    
  }else{
    pixel.Y=pmin(pixel.Y+increment*(timer>=startyr),pixel.Y.500)
    pixel.Y.unrest=pmin(pixel.Y.unrest+increment2*(timer>=startyr),pixel.Y.unrest.500)
  }


  if(!is.na(y_steppe)){
    pixel.Y[id.steppe]=y_steppe
  }
  if(!is.na(y_steppe)){
    pixel.Y.unrest[id.steppe]=y_steppe
  }
  
  
  if(!is.na(y_steppeeast)){
    pixel.Y[id.steppeeast]=y_steppeeast
  }
  if(!is.na(y_steppeeast)){
    pixel.Y.unrest[id.steppeeast]=y_steppeeast
  }
  
  
  if (drawmap & timer%%Tplot==0){
    drawing()
  }

  ############### Make War ###################;
  
  nation.resource=sapply(nation.territory,function(x) sum(pixel.Y[x]))
  pixel.belong=stack(setNames(nation.territory, seq_along(nation.territory)))
  colnames(pixel.belong)[1]="pixelid"
  colnames(pixel.belong)[2]="nation.fid"
  pixel.belong=pixel.belong[order(pixel.belong$pixelid),]
  

  pixel.unrestprob=pmax(pmin(alpha*pixel.Y.unrest,1),unrest_min);
  pixel.unrestprob[africa]=pixel.unrestprob[africa];
  pixel.unrestprob[EAdta$Americas==1]=pixel.unrestprob[EAdta$Americas==1];
  
  pixel.unrestprob[timer<startyr]=unrest_min
  
  
  pixel.unrest=which(apply(pixel.unrestprob,1,function(x) rbinom(1, 1, x))==1)
  
  
  pixel.unrest=setdiff(pixel.unrest,pixel.single)
  pixel.unrest.neighbors=lapply(pixel.unrest,function(x) unique(unlist(pixel.borderM[x])))
  
  
  if(conflictmech=="random"){
    pixel.unrest.neighbors.land=lapply(pixel.unrest,function(x) unique(unlist(pixel.borderM_land[x])))
    pixel.unrest.neighbors.sea=mapply(setdiff,pixel.unrest.neighbors,pixel.unrest.neighbors.land,SIMPLIFY=FALSE)
    pixel.unrest.landborder=lapply(pixel.unrest.neighbors.land,length)
    
    pixel.unrest.neighbors.land.Y=lapply(pixel.unrest.neighbors.land,function(x) rep(1,length(x)))
    pixel.unrest.neighbors.land.Ysum=lapply(pixel.unrest.neighbors.land.Y,sum)
    pixel.unrest.neighbors.land.prob=lapply(pixel.unrest.neighbors.land.Y,function(x) x/sum(x))
    pixel.unrest.neighbors.land.prob=mapply(function(x,y) x*y/6,pixel.unrest.neighbors.land.prob,pixel.unrest.landborder,SIMPLIFY=FALSE)
    
        pixel.unrest.neighbors.sea.Y=lapply(pixel.unrest.neighbors.sea,function(x) rep(1,length(x)))
    pixel.unrest.neighbors.sea.Ysum=lapply(pixel.unrest.neighbors.sea.Y,sum)
    pixel.unrest.neighbors.sea.prob=lapply(pixel.unrest.neighbors.sea.Y,function(x) x/sum(x))
    pixel.unrest.neighbors.sea.prob=mapply(function(x,y) (alpha_sea)*x*(6-y)/6,pixel.unrest.neighbors.sea.prob,pixel.unrest.landborder,SIMPLIFY=FALSE)
    
    
    pixel.unrest.neighbors=mapply(function(x,y) c(x,y),pixel.unrest.neighbors.land,pixel.unrest.neighbors.sea,SIMPLIFY=FALSE)
    pixel.unrest.neighbors.prob=mapply(function(x,y) c(x,y),pixel.unrest.neighbors.land.prob,pixel.unrest.neighbors.sea.prob,SIMPLIFY=FALSE)
    
    pixel.unrest.neighbors.prob=lapply(pixel.unrest.neighbors.prob,function(x) x*min(length(x),6)/6)
    pixel.unrest.target=mapply(function(x,y) ifelse(length(x)>1,sample(x,1,prob=y),x),pixel.unrest.neighbors,pixel.unrest.neighbors.prob,SIMPLIFY = TRUE)
    
    unrest.belongnation=as.numeric(pixel.belong$nation.fid[pixel.unrest])
    
    noseawar.prob=(6-unlist(pixel.unrest.landborder))/6*(1-alpha_sea)
    noseawar.id=which(apply(as.array(noseawar.prob),1,function(x) rbinom(1, 1, x))==1)
    
  }
  if(conflictmech=="expgain"){

    unrest.neighbors.belongnation=lapply(pixel.unrest.neighbors,function(x) pixel.belong$nation.fid[x])
    unrest.belongnation=as.numeric(pixel.belong$nation.fid[pixel.unrest])
    unrest.enemyneighbors0=mapply(function(x,y,z) x[y!=z],pixel.unrest.neighbors,unrest.neighbors.belongnation,unrest.belongnation)
    unrest.noenemy=sapply(unrest.enemyneighbors0,function(x) length(x)==0 | is.null(x))

    
    pixel.unrest=pixel.unrest[!unrest.noenemy]
    unrest.belongnation=as.numeric(pixel.belong$nation.fid[pixel.unrest])
    unrest.enemyneighbors=unrest.enemyneighbors0[which(unrest.noenemy==FALSE)]
    unrest.enemyneighbors.belongnation=lapply(unrest.enemyneighbors,function(x) pixel.belong$nation.fid[x])
    
    pixel.unrest.enemy.costY=mapply(costYs,unrest.enemyneighbors.belongnation,unrest.enemyneighbors)
    pixel.unrest.self.costY=mapply(costYs,unrest.belongnation,unrest.enemyneighbors)
    
    unrest.enemyneighbors.land=lapply(unrest.enemyneighbors,function(x) unique(unlist(pixel.borderM_land[x])))
    unrest.enemyneighbors.seavector=mapply(function(x,y) 1-x%in%y,unrest.enemyneighbors,unrest.enemyneighbors.land,SIMPLIFY=FALSE) 
    
    
    pixel.unrest.expgain=mapply(function(Yattack,Ydefense,target,seavector) (1-seavector*(1-alpha_sea))*unlist(Yattack)/(unlist(Yattack)+unlist(Ydefense))*pixel.Y[target]/(1+theta(target)),pixel.unrest.self.costY,pixel.unrest.enemy.costY,unrest.enemyneighbors,unrest.enemyneighbors.seavector)
    pixel.unrest.target=mapply(function(x,y) x[which.max(y)],unrest.enemyneighbors,pixel.unrest.expgain,SIMPLIFY = TRUE)
    
    seawar.id=which(mapply(function(x,y) 1-x%in%y,pixel.unrest.target,unrest.enemyneighbors.land,SIMPLIFY=TRUE)==1)
    #print(seawar.id)
    noseawar.id=seawar.id[which(rbinom(length(seawar.id), 1, (1-alpha_sea))==1)]
    
  }
  if(conflictmech=="neighborY"){
    pixel.unrest.neighbors.land=lapply(pixel.unrest,function(x) unique(unlist(pixel.borderM_land[x])))
    pixel.unrest.neighbors.sea=mapply(setdiff,pixel.unrest.neighbors,pixel.unrest.neighbors.land,SIMPLIFY=FALSE)
    pixel.unrest.landborder=lapply(pixel.unrest.neighbors.land,length)
    

    pixel.unrest.neighbors.land.Y=lapply(pixel.unrest.neighbors.land,function(x) pixel.Y[x])
    pixel.unrest.neighbors.land.Ysum=lapply(pixel.unrest.neighbors.land.Y,sum)
    pixel.unrest.neighbors.land.prob=lapply(pixel.unrest.neighbors.land.Y,function(x) x/sum(x))
    pixel.unrest.neighbors.land.prob=mapply(function(x,y) x*y/6,pixel.unrest.neighbors.land.prob,pixel.unrest.landborder,SIMPLIFY=FALSE)
  
    pixel.unrest.neighbors.sea.Y=lapply(pixel.unrest.neighbors.sea,function(x) pixel.Y[x])
    pixel.unrest.neighbors.sea.Ysum=lapply(pixel.unrest.neighbors.sea.Y,sum)
    pixel.unrest.neighbors.sea.prob=lapply(pixel.unrest.neighbors.sea.Y,function(x) x/sum(x))
    pixel.unrest.neighbors.sea.prob=mapply(function(x,y) (alpha_sea)*x*(6-y)/6,pixel.unrest.neighbors.sea.prob,pixel.unrest.landborder,SIMPLIFY=FALSE)
    
    
    pixel.unrest.neighbors=mapply(function(x,y) c(x,y),pixel.unrest.neighbors.land,pixel.unrest.neighbors.sea,SIMPLIFY=FALSE)
    pixel.unrest.neighbors.prob=mapply(function(x,y) c(x,y),pixel.unrest.neighbors.land.prob,pixel.unrest.neighbors.sea.prob,SIMPLIFY=FALSE)
    
    pixel.unrest.neighbors.prob=lapply(pixel.unrest.neighbors.prob,function(x) x*min(length(x),6)/6)
    pixel.unrest.target=mapply(function(x,y) ifelse(length(x)>1,sample(x,1,prob=y),x),pixel.unrest.neighbors,pixel.unrest.neighbors.prob,SIMPLIFY = TRUE)
    
    unrest.belongnation=as.numeric(pixel.belong$nation.fid[pixel.unrest])

    noseawar.prob=(6-unlist(pixel.unrest.landborder))/6*(1-alpha_sea)
    noseawar.id=which(apply(as.array(noseawar.prob),1,function(x) rbinom(1, 1, x))==1)
  }
  
  unrest.target.belongnation=as.numeric(pixel.belong$nation.fid[pixel.unrest.target])
  warid=which(unrest.belongnation!=unrest.target.belongnation)
  warid=setdiff(warid,noseawar.id)
  
  pixel.attack=pixel.unrest[warid]
  pixel.attack.target=pixel.unrest.target[warid]


  
  
  if(length(warid)>0){
    pixel.attack=pixel.unrest[warid]
    pixel.attack.target=pixel.unrest.target[warid]
    nation.attack=unrest.belongnation[warid]
    nation.defense=unrest.target.belongnation[warid]
    nation.warcnt=as.matrix(table(c(nation.defense,nation.attack)))


##############estimate resource###########
    defense.connectpixel=mapply(function(x,y) PixelLandConnected(x,y),pixel.attack,nation.territory[nation.defense], SIMPLIFY=FALSE)
    attack.connectpixel=mapply(function(x,y) PixelLandConnected(x,y),pixel.attack.target,nation.territory[nation.attack], SIMPLIFY=FALSE)
    seawar=mapply(function(x,y) length(x)==0 & length(y)==0,defense.connectpixel,attack.connectpixel)
    


    defense.connectY=sapply(defense.connectpixel,function(x) sum(pixel.Y[x]))
    defense.unconnectY=nation.resource[nation.defense]-defense.connectY
    defense.Y.est=defense.connectY+sigma*defense.unconnectY
    
    attack.connectY=sapply(attack.connectpixel,function(x) sum(pixel.Y[x]))
    attack.unconnectY=nation.resource[nation.attack]-attack.connectY
    attack.Y.est=attack.connectY+sigma*attack.unconnectY
    
    defender.oppos.Yest.sum=aggregate(attack.Y.est,list(nation.defense),sum,simplify = TRUE)
    attack.oppos.Yest.sum=aggregate(defense.Y.est,list(nation.attack),sum,simplify = TRUE)
    oppos.Yest=rbind(defender.oppos.Yest.sum,attack.oppos.Yest.sum)
    oppos.Yest.sum0=aggregate(oppos.Yest[,2],list(oppos.Yest[,1]),sum,simplify=TRUE)
    oppos.Yest.sum=rep(0,nationno)
    oppos.Yest.sum[oppos.Yest.sum0[,1]]=oppos.Yest.sum0[,2]
    

    
###############allocate resource################
    attack.Y.ratio=defense.Y.est/(oppos.Yest.sum[nation.attack])
    attack.Y.cost=attack.Y.est*attack.Y.ratio
    
    defense.Y.ratio=attack.Y.est/(oppos.Yest.sum[nation.defense])
    defense.Y.cost=defense.Y.est*defense.Y.ratio
    
    
##################################
    theta.attack=theta(pixel.attack)
    theta.defense=theta(pixel.attack.target)
    theta.max=pmax(theta.attack,theta.defense)

    
    if(river==1){
      anyriver=(pixel.riverdummy[pixel.attack]+pixel.riverdummy[pixel.attack.target])>0
      theta.max=theta.max-riverconnected(pixel.attack,pixel.attack.target)*theta.max+(1-riverconnected(pixel.attack,pixel.attack.target))*anyriver*2
    }
    
    if(Theta0==0){
      theta.max=theta.max*0
    }
    
    theta.max=theta.max+seawar*Theta_sea
    
    nation.attack.stepborder=(nation.orgpixel[nation.attack]%in%id.steppeeast)*(sapply(nation.territory[nation.attack],function(x) length(intersect0(x,id.steppeeast))>0))*(Psi_steppe-1)+1
    nation.defense.stepborder=(nation.orgpixel[nation.defense]%in%id.steppeeast)*(sapply(nation.territory[nation.defense],function(x) length(intersect0(x,id.steppeeast))>0))*(Psi_steppe-1)+1
    
    
    pixel.attack.winprob=nation.attack.stepborder*attack.Y.cost/(nation.attack.stepborder*attack.Y.cost+nation.defense.stepborder*defense.Y.cost)/(1+theta.max)
    pixel.attack.loseprob=nation.defense.stepborder*defense.Y.cost/(nation.attack.stepborder*attack.Y.cost+nation.defense.stepborder*defense.Y.cost)/(1+theta.max)
    pixel.attack.tieprob=1-pixel.attack.winprob-pixel.attack.loseprob
    outcome=mapply(function(x,y) rmultinom(1, size=1, prob=c(x,y,max(1-x-y,0))),pixel.attack.winprob,pixel.attack.loseprob)
    
    pixel.attack.win=pixel.attack[outcome[1,]==1]
    pixel.attack.lose=pixel.attack[outcome[2,]==1]
    pixel.attack.tie=pixel.attack[outcome[3,]==1]
    pixel.attack.win.fid=which(outcome[1,]==1)
    pixel.attack.lose.fid=which(outcome[2,]==1)
    pixel.attack.tie.fid=which(outcome[3,]==1)

    pixel.attack.win.ufid=pixel.attack.win.fid[duplicatedR(pixel.attack.target[pixel.attack.win.fid])==FALSE]
    if(length(pixel.attack.win.ufid)!=length(pixel.attack.win.fid)){
      outcome[1,setdiff(pixel.attack.win.fid,pixel.attack.win.ufid)]=0
      outcome[3,setdiff(pixel.attack.win.fid,pixel.attack.win.ufid)]=1
      pixel.attack.win=pixel.attack[outcome[1,]==1]
      pixel.attack.lose=pixel.attack[outcome[2,]==1]
      pixel.attack.tie=pixel.attack[outcome[3,]==1]
      pixel.attack.win.fid=which(outcome[1,]==1)
      pixel.attack.lose.fid=which(outcome[2,]==1)
      pixel.attack.tie.fid=which(outcome[3,]==1)
    }
    if(length(pixel.attack.lose.fid)>0 & length(pixel.attack.win.fid)>0){
      duplic=duplicatedR(c(pixel.attack[pixel.attack.lose.fid],pixel.attack.target[pixel.attack.win.fid]))
      pixel.attack.lose.ufid=pixel.attack.lose.fid[(duplic==FALSE)[1:length(pixel.attack.lose.fid)]]
      pixel.attack.win.ufid=pixel.attack.win.fid[(duplic==FALSE)[(length(pixel.attack.lose.fid)+1):(length(pixel.attack.lose.fid)+length(pixel.attack.win.fid))]]
      outcome[2,setdiff(pixel.attack.lose.fid,pixel.attack.lose.ufid)]=0
      outcome[3,setdiff(pixel.attack.lose.fid,pixel.attack.lose.ufid)]=1
      outcome[1,setdiff(pixel.attack.win.fid,pixel.attack.win.ufid)]=0
      outcome[3,setdiff(pixel.attack.win.fid,pixel.attack.win.ufid)]=1
      pixel.attack.win=pixel.attack[outcome[1,]==1]
      pixel.attack.lose=pixel.attack[outcome[2,]==1]
      pixel.attack.tie=pixel.attack[outcome[3,]==1]
      pixel.attack.win.fid=which(outcome[1,]==1)
      pixel.attack.lose.fid=which(outcome[2,]==1)
      pixel.attack.tie.fid=which(outcome[3,]==1)
    }

    if(lengthNA(pixel.attack.win)>0 | lengthNA(pixel.attack.lose)>0){
      pixellost.merge=split(c(pixel.attack.target[pixel.attack.win.fid],pixel.attack[pixel.attack.lose.fid]),factor(c(nation.defense[pixel.attack.win.fid],nation.attack[pixel.attack.lose.fid])))
      nation.lost.mergeid=as.numeric(names(pixellost.merge))
      pixelwin.merge=split(c(pixel.attack.target[pixel.attack.win.fid],pixel.attack[pixel.attack.lose.fid]),factor(c(nation.attack[pixel.attack.win.fid],nation.defense[pixel.attack.lose.fid])))
      nation.win.mergeid=as.numeric(names(pixelwin.merge))
      nation.territory[nation.lost.mergeid]=mapply(function(x,y) x[!(x%in%y)], nation.territory[nation.lost.mergeid], pixellost.merge, SIMPLIFY=FALSE)
      nation.territory[nation.win.mergeid]=mapply(function(x,y) c(x,y), nation.territory[nation.win.mergeid], pixelwin.merge, SIMPLIFY=FALSE)

      nationdied.fid=which(sapply(nation.territory,lengthNA)==0)
      if(length(nationdied.fid)>0){
        nation.territory=nation.territory[-nationdied.fid]
        nation.orgpixel=nation.orgpixel[-nationdied.fid]
      }
    }
    
  }


  ##########Secession#################

  if(Theta0!=0){
  pixel.belong=stack(setNames(nation.territory, seq_along(nation.territory)))
  colnames(pixel.belong)[1]="pixelid"
  colnames(pixel.belong)[2]="nation.fid"
  pixel.belong=pixel.belong[order(pixel.belong$pixelid),]
  nation.pixelborder=mapply(setdiff,lapply(split(pixel.borderM,pixel.belong$nation.fid),function(x) unique(unlist(x))),nation.territory,SIMPLIFY = FALSE)
  nation.pixelborderIn=mapply(intersect0,lapply(nation.pixelborder,function(x) unique(unlist(pixel.borderM[x]))),nation.territory,SIMPLIFY = FALSE)
  nation.borderlength=sapply(nation.pixelborderIn,length)
  nation.size2L=which(lapply(nation.territory,function(x) lengthNA(x)>1)==1)
  nation.territoryS2=nation.pixelborderIn[nation.size2L]

  

  nation.dellandborder=mapply(setdiff,nation.territory,nation.pixelborderIn,SIMPLIFY = FALSE)
  nation.coastL=sapply(nation.dellandborder,function(x) sum(pixel.coast[x]))
  nation.medcoastL=sapply(nation.dellandborder,function(x) sum(pixel.medcoast[x]))
  nation.allborderL=nation.borderlength+nation.coastL-nation.medcoastL

  pixel.seperate=lapply(nation.territoryS2,function(pixelvec) pixelvec[which(rbinom(length(pixelvec),1, sep.prob(pixelvec))==1)] )
  nation.territory[nation.size2L]=mapply(function(x,y) x[!(x%in%y)], nation.territory[nation.size2L], pixel.seperate, SIMPLIFY=FALSE)
  nation.territory=c(nation.territory,as.vector(unlist(pixel.seperate)))
  pixel.seperate.orgpixel=nation.orgpixel[as.numeric(names(unlist(pixel.seperate)))]
  nation.orgpixel=c(nation.orgpixel,pixel.seperate.orgpixel)
  nationdied.fid=which(sapply(nation.territory,lengthNA)==0)
  if(length(nationdied.fid)>0){
    nation.territory=nation.territory[-nationdied.fid]
    nation.orgpixel=nation.orgpixel[-nationdied.fid]
  }
  }
  
  ######## Shock ############
  
  if(shockprob_regime>0){
    nation.break=which(rbinom(length(nation.territory),1,shockprob_regime)==1)
    if(length(nation.break)>0){
      nation.territory=c(nation.territory[-nation.break],as.vector(unlist(nation.territory[nation.break])))
    }
  }

  if(shockprob_general>0){
    if(rbinom(1, 1, shockprob_general)==1){
      nation.territory=as.list(seq(1:pixelno))
    }
  }




  
  ############split disconnected parts##########
  nation.split=lapply(nation.territory,function(x) SplitDisconnect(x))

   
   splitid=which(lapply(nation.split,length)>1)
   if(length(splitid)>0){
     nation.split.territory=c()
     nation.split.orgpixel=c()
     for(i in splitid){
       nation.split.territory=c(nation.split.territory,nation.split[[i]])
       nation.split.orgpixel=c(nation.split.orgpixel,rep(nation.orgpixel[i],length(nation.split[[i]])))
     }
     nation.territory=nation.territory[-splitid]
     nation.territory=c(nation.territory,nation.split.territory)
     nation.orgpixel=nation.orgpixel[-splitid]
     nation.orgpixel=c(nation.orgpixel,nation.split.orgpixel)
   }
  
   
  ##########Calculate Herfindahl Index and Empire No.##############
  nationno=length(nation.territory)
  nation.territoryEU=lapply(nation.territory,function(x) x[x %in% pixel.EU])
  pixelno.EU=length(pixel.EU)
  Hihis.EU[timer]=sum((unlist(sapply(nation.territoryEU,length))/pixelno.EU)^2)
  regimeno.EU[timer]=sum(unlist(sapply(nation.territoryEU,length))>0)
  rangemin=1:min(10,regimeno.EU[timer])
  regimesize.EU[rangemin,timer]=sort(sapply(nation.territoryEU,length),TRUE)[rangemin]
  
  nation.territoryCN=lapply(nation.territory,function(x) x[x %in% pixel.CN])
  pixelno.CN=length(pixel.CN)
  Hihis.CN[timer]=sum((unlist(sapply(nation.territoryCN,length))/pixelno.CN)^2)
  regimeno.CN[timer]=sum(unlist(sapply(nation.territoryCN,length))>0)
  rangemin=1:min(10,regimeno.CN[timer])
  regimesize.CN[rangemin,timer]=sort(sapply(nation.territoryCN,length),TRUE)[rangemin]
  
  
  nation.territoryMed=lapply(nation.territory,function(x) x[x %in% medsea.region])
  pixelno.Med=length(medsea.region)
  Hihis.Med[timer]=sum((unlist(sapply(nation.territoryMed,length))/pixelno.Med)^2)
  regimeno.Med[timer]=sum(unlist(sapply(nation.territoryMed,length))>0)
  rangemin=1:min(10,regimeno.Med[timer])
  regimesize.Med[rangemin,timer]=sort(sapply(nation.territoryMed,length),TRUE)[rangemin]
  
  
  nation.territoryMideast=lapply(nation.territory,function(x) x[x %in% pixel.Mideast])
  pixelno.Mideast=length(pixel.Mideast)
  Hihis.Mideast[timer]=sum((unlist(sapply(nation.territoryMideast,length))/pixelno.Mideast)^2)
  regimeno.Mideast[timer]=sum(unlist(sapply(nation.territoryMideast,length))>0)
  rangemin=1:min(10,regimeno.Mideast[timer])
  regimesize.Mideast[rangemin,timer]=sort(sapply(nation.territoryMideast,length),TRUE)[rangemin]
  
  nation.territoryIndia=lapply(nation.territory,function(x) x[x %in% pixel.India])
  pixelno.India=length(pixel.India)
  Hihis.India[timer]=sum((unlist(sapply(nation.territoryIndia,length))/pixelno.India)^2)
  regimeno.India[timer]=sum(unlist(sapply(nation.territoryIndia,length))>0)
  rangemin=1:min(10,regimeno.India[timer])
  regimesize.India[rangemin,timer]=sort(sapply(nation.territoryIndia,length),TRUE)[rangemin]
  
  nation.territorySEAsia=lapply(nation.territory,function(x) x[x %in% pixel.SoutheastAsia])
  pixelno.SEAsia=length(pixel.SoutheastAsia)
  Hihis.SEAsia[timer]=sum((unlist(sapply(nation.territorySEAsia,length))/pixelno.SEAsia)^2)
  regimeno.SEAsia[timer]=sum(unlist(sapply(nation.territorySEAsia,length))>0)
  rangemin=1:min(10,regimeno.SEAsia[timer])
  regimesize.SEAsia[rangemin,timer]=sort(sapply(nation.territorySEAsia,length),TRUE)[rangemin]
  
  nation.territoryEuropaX=lapply(nation.territory,function(x) x[x %in% pixel.EuropaX])
  pixelno.EuropaX=length(pixel.EuropaX)
  Hihis.EuropaX[timer]=sum((unlist(sapply(nation.territoryEuropaX,length))/pixelno.EuropaX)^2)
  regimeno.EuropaX[timer]=sum(unlist(sapply(nation.territoryEuropaX,length))>0)
  rangemin=1:min(10,regimeno.EuropaX[timer])
  regimesize.EuropaX[rangemin,timer]=sort(sapply(nation.territoryEuropaX,length),TRUE)[rangemin]
  
  
  nation.territoryAfricaET=lapply(nation.territory,function(x) x[x %in% pixel.AfricaET])
  pixelno.AfricaET=length(pixel.AfricaET)
  regimeno.AfricaET[timer]=sum(unlist(sapply(nation.territoryAfricaET,length))>0)
  rangemin=1:min(10,regimeno.AfricaET[timer])
  regimesize.AfricaET[rangemin,timer]=sort(sapply(nation.territoryAfricaET,length),TRUE)[rangemin]
  
  nation.territoryAfricaWT=lapply(nation.territory,function(x) x[x %in% pixel.AfricaWT])
  pixelno.AfricaWT=length(pixel.AfricaWT)
  regimeno.AfricaWT[timer]=sum(unlist(sapply(nation.territoryAfricaWT,length))>0)
  rangemin=1:min(10,regimeno.AfricaWT[timer])
  regimesize.AfricaWT[rangemin,timer]=sort(sapply(nation.territoryAfricaWT,length),TRUE)[rangemin]
  
  nation.territoryAfricaME=lapply(nation.territory,function(x) x[x %in% pixel.AfricaME])
  pixelno.AfricaME=length(pixel.AfricaME)
  regimeno.AfricaME[timer]=sum(unlist(sapply(nation.territoryAfricaME,length))>0)
  rangemin=1:min(10,regimeno.AfricaME[timer])
  regimesize.AfricaME[rangemin,timer]=sort(sapply(nation.territoryAfricaME,length),TRUE)[rangemin]
  
  nation.territoryAfricaNH=lapply(nation.territory,function(x) x[x %in% pixel.AfricaNH])
  pixelno.AfricaNH=length(pixel.AfricaNH)
  regimeno.AfricaNH[timer]=sum(unlist(sapply(nation.territoryAfricaNH,length))>0)
  rangemin=1:min(10,regimeno.AfricaNH[timer])
  regimesize.AfricaNH[rangemin,timer]=sort(sapply(nation.territoryAfricaNH,length),TRUE)[rangemin]
  
  nation.territoryAfricaSH=lapply(nation.territory,function(x) x[x %in% pixel.AfricaSH])
  pixelno.AfricaSH=length(pixel.AfricaSH)
  regimeno.AfricaSH[timer]=sum(unlist(sapply(nation.territoryAfricaSH,length))>0)
  rangemin=1:min(10,regimeno.AfricaSH[timer])
  regimesize.AfricaSH[rangemin,timer]=sort(sapply(nation.territoryAfricaSH,length),TRUE)[rangemin]
  
  nation.territoryAmericaNH=lapply(nation.territory,function(x) x[x %in% pixel.AmericaNH])
  pixelno.AmericaNH=length(pixel.AmericaNH)
  regimeno.AmericaNH[timer]=sum(unlist(sapply(nation.territoryAmericaNH,length))>0)
  rangemin=1:min(10,regimeno.AmericaNH[timer])
  regimesize.AmericaNH[rangemin,timer]=sort(sapply(nation.territoryAmericaNH,length),TRUE)[rangemin]
  
  nation.territoryAmericaSH=lapply(nation.territory,function(x) x[x %in% pixel.AmericaSH])
  pixelno.AmericaSH=length(pixel.AmericaSH)
  regimeno.AmericaSH[timer]=sum(unlist(sapply(nation.territoryAmericaSH,length))>0)
  rangemin=1:min(10,regimeno.AmericaSH[timer])
  regimesize.AmericaSH[rangemin,timer]=sort(sapply(nation.territoryAmericaSH,length),TRUE)[rangemin]
  
  nation.territoryAmericaCL=lapply(nation.territory,function(x) x[x %in% pixel.AmericaCL])
  pixelno.AmericaCL=length(pixel.AmericaCL)
  regimeno.AmericaCL[timer]=sum(unlist(sapply(nation.territoryAmericaCL,length))>0)
  rangemin=1:min(10,regimeno.AmericaCL[timer])
  regimesize.AmericaCL[rangemin,timer]=sort(sapply(nation.territoryAmericaCL,length),TRUE)[rangemin]
  
  
  
  EUmax=max(sapply(nation.territoryEU,length))
  CNmax=max(sapply(nation.territoryCN,length))
  
  EUmax.id=which.max(sapply(nation.territoryEU,length))
  CNmax.id=which.max(sapply(nation.territoryCN,length))
  
  IntersectCN[timer]=length(unlist(nation.territoryCN[CNmax.id]));
  UnionCN[timer]=length(union(pixel.CN,unlist(nation.territory[CNmax.id])));
  
  IntersectEU[timer]=length(unlist(nation.territoryEU[EUmax.id]));
  UnionEU[timer]=length(union(pixel.EU,unlist(nation.territory[EUmax.id])));

  EUmax.orgpixel=nation.orgpixel[EUmax.id]
  CNmax.orgpixel=nation.orgpixel[CNmax.id]

}


return(list(Hihis.EU=Hihis.EU,Hihis.CN=Hihis.CN,Hihis.Med=Hihis.Med,Hihis.Mideast=Hihis.Mideast,Hihis.India=Hihis.India,Hihis.SEAsia=Hihis.SEAsia,Hihis.EuropaX=Hihis.EuropaX
            ,Empireno.CN=Empireno.CN,Empireno.EU=Empireno.EU,Empireno.India=Empireno.India,Empireno.SoutheastAsia=Empireno.SoutheastAsia,Empireno.EuropaX=Empireno.EuropaX, Empireno.Mideast=Empireno.Mideast
            ,Empireno.AfricaET=Empireno.AfricaET,Empireno.AfricaWT=Empireno.AfricaWT,Empireno.AfricaME=Empireno.AfricaME,Empireno.AfricaNH=Empireno.AfricaNH,Empireno.AfricaSH=Empireno.AfricaSH
            ,Empireno.AmericaNH=Empireno.AmericaNH,Empireno.AmericaSH=Empireno.AmericaSH,Empireno.AmericaCL=Empireno.AmericaCL
            ,EUmax.orgpixel=EUmax.orgpixel,CNmax.orgpixel=CNmax.orgpixel,EUmax=EUmax,CNmax=CNmax
            ,UnionEU=UnionEU,UnionCN=UnionCN,IntersectCN=IntersectCN,IntersectEU=IntersectEU
            ,regimeno.EU=regimeno.EU
            ,regimesize.EU=regimesize.EU
            ,regimeno.CN=regimeno.CN
            ,regimesize.CN=regimesize.CN
            ,regimeno.Med=regimeno.Med
            ,regimesize.Med=regimesize.Med
            ,regimeno.Mideast=regimeno.Mideast
            ,regimesize.Mideast=regimesize.Mideast
            ,regimeno.India=regimeno.India
            ,regimesize.India=regimesize.India
            ,regimeno.SEAsia=regimeno.SEAsia
            ,regimesize.SEAsia=regimesize.SEAsia
            ,regimeno.EuropaX=regimeno.EuropaX
            ,regimesize.EuropaX=regimesize.EuropaX
            ,regimesize.AfricaET=regimesize.AfricaET
            ,regimesize.AfricaWT=regimesize.AfricaWT
            ,regimesize.AfricaME=regimesize.AfricaME
            ,regimesize.AfricaNH=regimesize.AfricaNH
            ,regimesize.AfricaSH=regimesize.AfricaSH
            ,regimesize.AmericaNH=regimesize.AmericaNH
            ,regimesize.AmericaSH=regimesize.AmericaSH
            ,regimesize.AmericaCL=regimesize.AmericaCL
            ,regimeno.AfricaET=regimeno.AfricaET
            ,regimeno.AfricaWT=regimeno.AfricaWT
            ,regimeno.AfricaME=regimeno.AfricaME
            ,regimeno.AfricaNH=regimeno.AfricaNH
            ,regimeno.AfricaSH=regimeno.AfricaSH
            ,regimeno.AmericaNH=regimeno.AmericaNH
            ,regimeno.AmericaSH=regimeno.AmericaSH
            ,regimeno.AmericaCL=regimeno.AmericaCL
            ))




