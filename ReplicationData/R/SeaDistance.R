library(foreach)    #added
library(doParallel) #added
library(parallel)   #added
library(fanplot)
library(grid)
library(dplyr)
library(ggplot2)
library(foreign)
library(reshape2)
library(mapproj)
library(maptools)
library(rgeos)
library(plyr)
library(rgdal)
library(RColorBrewer)
library(spdep)

#start=coastid[6]
#territory=coastnb_list[[6]]

seadist <- function(start,territory){
  graph=lapply(seaneigbor[territory],function(x) intersect(x,territory))
  #territory=territory[sapply(graph,length)>0]
  #graph=lapply(seaneigbor[territory],function(x) intersect(x,territory))
  
  distances = rep(Inf, length(graph))
  visited = rep(FALSE, length(graph))
  distances[territory==start] = 0
  
  
  
  repeat{
    
    shortest_distance = Inf
    shortest_index = -1
    for(i in seq_along(distances)) {
      if(distances[i] < shortest_distance && !visited[i]){
        shortest_distance = distances[i]
        shortest_index = i
      }
    }
    #cat("Visiting node ", shortest_index, " with current distance ", shortest_distance, "\n")
    
    if(shortest_index == -1){
      
      othercoast=intersect(setdiff(coastid,start),territory)
      distances[territory %in% othercoast]=sapply(graph[territory %in% othercoast],function(x) min(distances[territory %in% x])+1)
      outdist=distances[territory %in% othercoast]
      outnbid=territory[territory %in% othercoast]
      
      return (list(outdist,outnbid))
    }
    if(length(graph[[shortest_index]])>0){
      for (i in graph[[shortest_index]]) {
        if(distances[territory==i] > distances[shortest_index] + 1){
          distances[territory==i] = distances[shortest_index] + 1
          #cat("Updating distance of node ", i, " to ", distances[i], "\n")
        }
        visited[shortest_index] = TRUE
        #cat("Visited nodes: ", visited, "\n")
        #cat("Currently lowest distances: ", distances, "\n")
      }
    }else{
      visited[shortest_index] = TRUE
    }



  }
}




#mapfile=paste0(mainfolder,"map/SeaCoast")



#coastid=which(map$coastal==1)

coastnb500km=read.dbf("E:/Development Economics/Geography and the Size of Nation/GIS/coast500km.dbf")
coastnb500km$sfid_1=coastnb500km$sfid_1+1
coastnb500km$sfid=coastnb500km$sfid+1
coastnb_list=split(coastnb500km$sfid_1,coastnb500km$sfid)
coastid=as.numeric(names(coastnb_list))

seaneigbor=poly2nb(readOGR(dsn="E:/Development Economics/Geography and the Size of Nation/GIS/seacoast.shp"), queen=TRUE,snap = 1)
#sapply(seaneigbor,length)
#seaneigbor[start]
seaneigbor=sapply(seaneigbor,function(x) if(length(x)==1 && x==0){return(integer(0))}else{return(x)})
seaneigbor=lapply(seaneigbor,function(x) setdiff(x,coastid))

#start=coastid[1000]
#territory=coastnb_list[[1000]]
#seadist(coastid[1000],coastnb_list[[1000]])

#seaneigbor[coastid[1]]
#coastid[1000] %in% coastnb_list[[1000]]
#start=coastid[1]
#territory=coastnb_list[[1]]



#seadist(coastid[1],coastnb_list[[1]])
#seaneigbor[28712]


DistMatrix=mapply(seadist,as.list(coastid),coastnb_list,SIMPLIFY=FALSE)
#NbMatrix=mapply(function(start,territory) territory[territory %in% intersect(setdiff(coastid,start),territory)],as.list(coastid),coastnb_list,SIMPLIFY=FALSE)
NbList=lapply(DistMatrix,function(x) unlist(x[2]))
DistList=lapply(DistMatrix,function(x) unlist(x[1]))


coastmatch=read.dbf("E:/Development Economics/Geography and the Size of Nation/GIS/coastal_mapfid_world.dbf")
coastmatch=coastmatch[coastmatch$mapfid!=0,]
coastmatch$sfid=coastmatch$sfid+1
coastmatch$mapfid=coastmatch$mapfid+1

coastid_map=coastmatch$mapfid[match(coastid,coastmatch$sfid)]
NbList_map=lapply(NbList,function(x) coastmatch$mapfid[match(x,coastmatch$sfid)])

#write.csv(data.frame(from=coastid_map[1000]-1,fid=NbMatrix_map[[1000]]-1,dist=DistList[1000]),"Distsea.csv")
outdist=data.frame(coastid=unlist(mapply(function(x,y) rep(x,length(y)),coastid_map,NbList_map)),NbList=unlist(NbList_map),Dist=unlist(DistList))

outdist=outdist[!(is.na(outdist$coastid) | is.na(outdist$NbList)),]
outdist=outdist[outdist$Dist<Inf,]

write.csv(outdist,"E:/Development Economics/Geography and the Size of Nation/GIS/DistMatrix_World.csv")


