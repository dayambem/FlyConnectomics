library(natverse)
library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
library(reticulate)
library(R.matlab)
library(dplyr)
# library(tidyverse)
library(pheatmap)
library(igraph)
library(pracma)
library(visNetwork)
library(zoo)

conn = neuprint_login(server= "https://neuprint.janelia.org",
                      token= "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImRheWFtYmUxQHN3YXJ0aG1vcmUuZWR1IiwibGV2ZWwiOiJub2F1dGgiLCJpbWFnZS11cmwiOiJodHRwczovL2xoMy5nb29nbGV1c2VyY29udGVudC5jb20vYS0vQU9oMTRHaEh4c1hpcVFPRlZMNVlUMF92QmRUMEt0YUR1Sm0yZnZPOVBMQ0M0UT1zOTYtYz9zej01MD9zej01MCIsImV4cCI6MTgzMzQ5MTQ3MX0.FUef2y8ICD62LqIM0paoxIpebQVHiL1ZK6EqsNi5IXU")

neuprint_connection(
  server = "https://neuprint.janelia.org",
  token = token,
  dataset = 'manc:v1.0',
  conn = conn,
  config = httr::config()
)

getDownstreamsNR <- function(ids=c(), nLayers=2, connStrength=5) {
  s_ConnT=c(list())
  s_Conn=data.frame()
  counter=c()
  for (i in 1:nLayers) { # Loop through Number of desired layers
    conn_N = neuprint_connection_table(
      bodyids=c(ids),
      partners = "outputs",
      roi = NULL,
      by.roi = FALSE,
      threshold = connStrength,
      summary = FALSE,
      details = FALSE,
      superLevel = FALSE,
      progress = FALSE,
      dataset = NULL,
      chunk = TRUE,
      all_segments = FALSE,
      conn = NULL,
    )
    nconn_N=conn_N['partner']

    if (i==1) {
      counterT=rep(i,length(unique(nconn_N[[1]])))
      counter=c(counter,counterT)
      s_Conn=data.frame(rbind(s_Conn,(unique(nconn_N))))
      ids=unique(nconn_N)# Repeat process if i ~= nLayers
    } else {

      P=setdiff(unique(nconn_N),s_Conn) # Extract connections unique to this layer with respect to previous layers
      # s_ConnT=c(s_ConnT,nconn_N) # Add unique connections to total list
      counterT=rep(i,length(P[[1]]))
      counter=c(counter,counterT)
      s_Conn=data.frame(rbind(s_Conn,P))
      ids=P# Repeat process if i ~= nLayers
    }
  }
  # s_ConnT=unique(s_ConnT)
  s_Conn=data.frame(s_Conn,counter)
  # col.names(s_Conn)=c('All Conns',1:nLayers)
  return(s_Conn)
}

pathsByEffic<-function(inputIds,outputIds,pL,ThrS=11,intermedConn=c()) {
  Thr=ThrS
  pathsn=data.frame()
  pathft=c()
  idsPc1=intermedConn
  if (length(intermedConn)>0){
    if (length(pathft)==0 & length(pathsn)!=0){
      break
    }else {
      while ((length(pathft)==0 & length(pathsn)==0)|(length(pathft)!=0 & length(pathsn)!=0)){
        paths=neuprint_get_shortest_paths(inputIds,outputIds,weightT = Thr)
        pathft=c()
        pathtt=c()
        pathwt=c()
        if (length(paths)==0) {
          break
        } else {
        for (j in 1:nrow(paths)) {
          if (paths['to'][j,1]%in%outputIds) {
            if (paths['depth'][j,1]==pL) {
              if (TRUE%in%(idsPc1%in%paths['to'][(j-paths['depth'][j,1]+1):(j-1),1])) {
                if (paths['from'][(j-paths['depth'][j,1]+1),1]%in%inputIds){

                pathft=paths['from'][(j-paths['depth'][j,1]+1):(j),1]
                pathtt=paths['to'][(j-paths['depth'][j,1]+1):(j),1]
                pathwt=paths['weight'][(j-paths['depth'][j,1]+1):(j),1]
                count=rep(Thr,length(pathft))
                pathsn=rbind(pathsn,data.frame(pathft,pathtt,pathwt,count))
                }
              }
            }
          }
        }
        }
        Thr=Thr+1

        }
    }
    pathsn=data.frame(pathsn,c(neuprint_get_neuron_names(pathsn['pathft'][[1]])),c(neuprint_get_neuron_names(pathsn['pathtt'][[1]])))
    colnames(pathsn)=c('from','to','weight','Threshold','name from','name to')
  } else {
    if (length(pathft)==0 & length(pathsn)!=0){
      break
    }else {
      while ((length(pathft)==0 & length(pathsn)==0)|(length(pathft)!=0 & length(pathsn)!=0)){
        paths=neuprint_get_shortest_paths(inputIds,outputIds,weightT = Thr)
        pathft=c()
        pathtt=c()
        pathwt=c()
        if (length(paths)==0) {
          break
        } else {
        for (j in 1:nrow(paths)) {
          if (paths['to'][j,1]%in%outputIds) {
            if (paths['depth'][j,1]==pL) {
              if (paths['from'][(j-paths['depth'][j,1]+1),1]%in%inputIds){
              pathft=paths['from'][(j-paths['depth'][j,1]+1):(j),1]
              pathtt=paths['to'][(j-paths['depth'][j,1]+1):(j),1]
              pathwt=paths['weight'][(j-paths['depth'][j,1]+1):(j),1]
              count=rep(Thr,length(pathft))
              count2=(1:length(count))
              pathsn=rbind(pathsn,data.frame(pathft,pathtt,pathwt,count,count2))
              }
            }
          }
        }
       }

        Thr=Thr+1
      }
    }
    pathsn=data.frame(pathsn,c(neuprint_get_neuron_names(pathsn['pathft'][[1]])),c(neuprint_get_neuron_names(pathsn['pathtt'][[1]])))
    colnames(pathsn)=c('from','to','weight','Threshold','depth','name from','name to')
  }

  uniquePath <- function(l){
    idxs <- 1:length(l)
    tmp <- lapply(l, table, useNA='always')
    l2 <- lapply(idxs, function(i){
      res <- l[[i]]
      for(j in idxs[-i]){
        if ( all(res %in% l[[j]]) & all(tmp[[i]] <= tmp[[j]])){
          res <- NULL; break
        }
      }
      res
    })
    Filter(Negate(is.null), l2)
  }

  nPaths=0
  ftL=list()
  tL=list()
  for (i in 1:nrow(pathsn)){
    if (pathsn['to'][i,1]%in%outputIds){
      if(pathsn['depth'][i,1]==pL) {
        nPaths=nPaths+1
        ftL[[nPaths]]=c((pathsn['from'][(i+1-pathsn['depth'][i,1]):i,1]),pathsn['to'][(i+1-pathsn['depth'][i,1]):i,1])
       # tL[[nPaths]]=(pathsn['to'][(i+1-pathsn['depth'][i,1]):i,1])
      }
    }
  }


  wghts=c()
  dpth=c()
  itr2=1
  itr3=1

  # ufL=list()
  # utL=list()
#  uLi=intersect(which(uniquePath(fL)%in%fL),which(uniquePath(tL)%in%tL))
  # uLi=(which(ftL%in%uniquePath(ftL)))
  #
  # for (i in uLi) {
  #   ufL[[itr3]]=ftL[[i]]
  #   #utL[[itr3]]=tL[[i]]
  #   itr3=itr3+1
  # }
  ufL=unique(ftL)

  connAdj=grouped_adjacency(c(unlist(ufL)),c(unlist(ufL)),ingroup = NULL,outgroup = NULL,)

  from=c()
  to=c()

  for (i in ufL) {
      itr1=1
      for (n in 1:(length(i)/2)) {
        wghts = c(wghts,connAdj[which(rownames(connAdj)==i[n]),which(colnames(connAdj)==i[n+pL])])
        #dpth =c(dpth,itr1)
        from=c(from,i[n])
        to=c(to,i[n+pL])
        itr1=itr1+1
      }
      itr2=itr2+1
      dpth=c(dpth,1:(length(i)/2))
  }

  pathsn=data.frame(from,to,wghts,dpth,neuprint_get_neuron_names((from)),neuprint_get_neuron_names((to)))
  colnames(pathsn)=c('from','to','weight','depth','name from','name to')


  nPaths=0
  count3=c()
  for (i in 1:nrow(pathsn)){
    if (pathsn['to'][i,1]%in%outputIds){
      if(pathsn['depth'][i,1]==pL) {
        nPaths=nPaths+1
        count3[nPaths]=sum(1/(pathsn['weight'][(i+1-pathsn['depth'][i,1]):i,1]))
      }
    }
  }


  ind1=1
  ind2=pL
  L=list()
  for (i in 1:nPaths) {
    L[[i]]=pathsn[ind1:ind2,1:ncol(pathsn)]
    ind1=ind1+pL
    ind2=ind2+pL
  }

  find_subseq = function(x,vec) {
    p=match(x[1],vec)
    if(is.na(p)||length(x)==1){ p }
    else { c(p,p+find_subseq(x[-1],vec[-seq_len(p)])) }
  }
  "%seq_somewhere_in%" = function(b,a) all(!is.na(find_subseq(b,a)))

  Ln=t(data.frame(c(1:ncol(pathsn))))
  colnames(Ln)<-names(pathsn)
  for (i in order(count3)){
    Lt=bind_rows(L[[i]])
   # if (!all(Lt['from'][[1]]%in%Ln['from'][[1]])&!all(Lt['to'][[1]]%in%Ln['to'][[1]])){
   # if (!(length(c(which(Lt['from'][[1]]==Ln['from'][[1]]),which(Lt['to'][[1]]==Ln['to'][[1]])))==length(c(Lt['from'][[1]],Lt['to'][[1]])))){
  # if (length(which(Ln['from'][[1]]%in%Lt['from'][[1]]))>1) {
   # if (!(round(mean(diff(which(Ln['from'][[1]]%in%Lt['from'][[1]])[1:length(Lt['from'][[1]])])),0)==1)){
    #  if (!(round(mean(diff(which(Ln['to'][[1]]%in%Lt['to'][[1]])[1:length(Lt['to'][[1]])])),0)==1)){
    # if (!(find_subseq(c(1:nrow(Ln)),which(Ln['from'][[1]]%in%Lt['from'][[1]]))[-length(find_subseq(c(1:nrow(Ln)),
         #                         (c(which(Ln['from'][[1]]%in%Lt['from'][[1]])))))]%seq_somewhere_in%1:nrow(Ln))){
      #if (!all(Lt['to'][[1]]%in%Ln['to'][[1]])){
        Ln=rbind(Ln,bind_rows(L[[i]]))
     # }
  #  }
 #  } else if (length(which(Lt['from'][[1]]%in%Ln['from'][[1]]))<length(Lt['from'][[1]])|
           #   length(which(Lt['to'][[1]]%in%Ln['to'][[1]]))<length(Lt['to'][[1]])) {
    # Ln=rbind(Ln,bind_rows(L[[i]]))

  # }
  }

  pathsn=Ln[-1,1:ncol(Ln)]


  # N=pathsn['weight'][[1]]
  # L=list()
  # A=data.frame()
  # An=list()
  # S=c()
  # j2=1
  # isN=FALSE
  # dpth=pL
  # while (!isN) {
  # j=1
  # L2=list()
  # for (i in 1:nrow(pathsn)) {
  #   if (pathsn['to'][i,1]%in%outputIds) {
  #     At=t(pathsn[(i+1-pathsn['depth'][i,1]):i,1:ncol(pathsn)])
  #     An[[j]]=1/(pathsn['weight'][(i+1-pathsn['depth'][i,1]):i,1])
  #     S[j]=sum(An[[i]])
  #     A=data.frame(A,At)
  #
  #     j=j+1
  #   }
  # }
  #

  # j=1
  # for (i in 1:nPaths){
  # L[[j]]=An[[which(S==min(S))]]
  # A=A[[-which(S==min(S))]]
  # S=S[-which(S==min(S))]
  # j=j+1
  # }
  #
  #   isN=TRUE
  # }

  return(pathsn)
}

getTypeList <- function(conns=s_ConnPV6) {
  cNames = neuprint_get_neuron_names(
    bodyids=c(conns),
    dataset = NULL,
    all_segments = TRUE,
    conn = NULL,
  )
  ttNames=c()
  tNames=c()

  for (i in 1:length(cNames)){
    if (is.na(sub("^([[:alpha:]]*).*", "\\1", cNames[i]))) {
      ttNames='N/A'
      tNames=c(tNames,ttNames)
    } else if (sub("^([[:alpha:]]*).*", "\\1", cNames[i])=='') {
      ttNames=cNames[i]
      tNames=c(tNames,ttNames)
    } else {
      ttNames = sub("^([[:alpha:]]*).*", "\\1", cNames[i])
      tNames=c(tNames,ttNames)
    }
  }
  return(tNames)
}

# findElbowPoint <- function(numbers) {
#   # Calculate total number of points
#   nPoints <- length(numbers)
#
#   # Coordinates of all the points
#   allCoord <- cbind(1:nPoints, numbers)
#
#   # Coordinates of the line from first to last point
#   lineCoord <- rbind(allCoord[1, ], allCoord[nPoints, ])
#
#   # Calculate distances of each point from the line
#   lineVec <- lineCoord[2,] - lineCoord[1,]
#   pointVec <- t(apply(allCoord, 1, function(point) {point - lineCoord[1,]}))
#   distanceVec <- abs(lineVec[2]*pointVec[,1] - lineVec[1]*pointVec[,2]) / sqrt(sum(lineVec^2))
#
#   # Return the index of the point with max distance
#   return(which.max(distanceVec))
# }

findElbowPoint <- function(values) {
  n <- length(values)
  if (n <= 2) {
    stop("The input list must contain at least 3 values.")
  }

  # Calculate cumulative sum of the values
  cum_sum <- cumsum(values)

  # Calculate the distances from each point to a line connecting the first and last point
  distances <- sapply(1:n, function(i) {
    x <- c(1, n)
    y <- c(cum_sum[1], cum_sum[n])
    abs((n*x[i] - sum(x))*(cum_sum[i]) - (sum(x^2) - sum(x)^2/n)*y[i]) /
      sqrt((n*sum(x^2) - sum(x)^2) * (sum((y-mean(y))^2)))
  })

  # Find the index at which the maximum distance occurs (elbow point)
  elbow_index <- which.max(distances)

  # Return the elbow point value
  elbow_point <- values[elbow_index]
  return(elbow_point)
}

find_elbow_point <- function(values) {
  n <- length(values)
  if (n <= 2) {
    stop("The input list must contain at least 3 values.")
  }

  # Create a data matrix with a single column for values
  data_matrix <- as.matrix(values)

  # Determine the maximum number of clusters to consider
  max_k <- min(10, n - 1)  # You can adjust the maximum number of clusters as needed

  # Calculate the Gap Statistic for each k from 1 to max_k
  gap_stat <- sapply(1:max_k, function(k) {
    kmeans_result <- kmeans(data_matrix, centers = k, nstart = 10)
    return(log(sum(kmeans_result$betweenss) / kmeans_result$tot.withinss))
  })

  # Find the optimal number of clusters (elbow point)
  elbow_point <- which(diff(gap_stat) <= 0)[1] + 1

  return(elbow_point)
}

findElbow <- function(values) {
  # First we calculate the differences between consecutive elements
  differences <- diff(values)

  # Then we calculate the differences between these differences (second order difference)
  second_order_diff <- diff(differences)

  # The elbow point is where the second order difference is maximum.
  # We add 2 because diff() reduces the length of the vector by 1 each time it's applied
  elbow_index <- which.max(second_order_diff) + 2

  return(list(elbow_point = values[elbow_index], index = elbow_index))
}


DNsT=dn.ids[which(getTypeList(dn.ids)%in%c('DNp','DNb','DNa','DNd'))]

##### edit this block

## List of input neurons' ids (ex: list(aIPg neuron ids, LHPV6 neuron ids))
inpN=list()
## Desired layers to traverse for search overlapping downstream targets
nLayers=1
## Desired threshold
Thr=3
## Desired path Length
pL=c(1,2)
# pL=1
# pL=2

## Whether or not return all or only overlapping downstream targets
intersN=FALSE

## Whether to find paths to downstream network neurons or only specified outpN
DS=TRUE
#
# ## Whether to threshold global conductances
# ThrC=TRUE

## Target overlapping  neuron types
BCs = list(rn.ids, ## [[1]]
           orn.ids, ## [[2]]
           hrn.ids, ## [[3]]
           pn.ids, ## [[4]]
           upn.ids, ## [[5]]
           mpn.ids, ## [[6]]
           vppn.ids, ## [[7]]
           dan.ids, ## [[8]]
           mbon.ids, ## [[9]]
           alln.ids, ## [[10]]
           ton.ids, ## [[11]]
           lhn.ids, ## [[12]]
           DNsT, ## [[13]] Descending neurons
           kc.ids, ## [[14]]
           apl.ids, ## [[15]]
           cent.ids, ## [[16]]
           lc.ids ## [[17]]
           )
## Choose target overlapping neuron type
# DNs=BCs[[13]]
DNs=c()
for (i in neuprint_ROIs()){
DNs=c(DNs,neuprint_bodies_in_ROI(i)$bodyid)
}
DNs=as.character(DNs)
#####

nNames=c('LHPV6a1','LHPV6a3_a','LHPV6a3_b','LHAV4a2','LHAV4a6')
nNames2 = c('pC1d','pC1e','aIPg1','aIPg2','aIPg3')
nNames3= c('CL062_a','CL062_b')
nNames4= c('DNa04','DNb01','DNut038','DNut037')
nNames5=c('b1 MN','b2 MN','b3 MN','i1 MN','i2 MN',
          # 'iii2 MN',
          'iii3 MN',
         # 'iii4 MN',
          'hg1 MN','hg2 MN','hg3 MN','hg4 MN')
nNames6=c('DNa04_R','DNa04_L','DNb01_R','DNb01_L','DNut038_R','DNut038_L','DNut037_R','DNut037_L')
nNames7=c('b1 MN_R','b1 MN_L','b2 MN_R','b2 MN_L','b3 MN_R','b3 MN_L',
          'i1 MN_R','i1 MN_L','i2 MN_R','i2 MN_L',
                  # 'iii2 MN',
                  'iii3 MN_R','iii3 MN_L',
                  # 'iii4 MN',
                  'hg1 MN_R','hg1 MN_L','hg2 MN_R','hg2 MN_L',
                  'hg3 MN_R','hg3 MN_L','hg4 MN_R','hg4 MN_L'
          )
nNames8=c('DNa04_R','DNb01_R','DNut038_R','DNut037_R')
nNames9=c('DNa04_L','DNb01_L','DNut038_L','DNut037_L')


ids=c()
for (i in nNames){
  idsT=neuprint_search(i,field='type')
  ids =c(ids,idsT$bodyid)
}

ids2=c()
for (i in nNames2){
  idsT=neuprint_search(i,field='type')
  ids2 =c(ids2,idsT$bodyid)
}

idsPc1=ids2[1:3]
intermedConn=idsPc1

ids3=c()
for (i in nNames3){
  idsT=neuprint_search(i,field='type')
  ids3 =c(ids3,idsT$bodyid)
}

ids4=c()
for (i in nNames4){
  idsT=neuprint_search(i,field='type')
  ids4 =c(ids4,idsT$bodyid)
}

ids4L=c()
for (i in nNames4){
  idsT=neuprint_search(i,field='type')
  ids4L =c(ids4L,idsT$bodyid[2])
}

ids4R=c()
for (i in nNames4){
  idsT=neuprint_search(i,field='type')
  ids4R =c(ids4R,idsT$bodyid[1])
}

# ids4R[2]=neuprint_search(nNames4[2],'type')$bodyid[2]
# ids4L[2]=neuprint_search(nNames4[2],'type')$bodyid[1]

ids5=c()
coln=c()
for (i in nNames5){
  idsT=neuprint_search(i,field='type')
  ids5 =c(ids5,idsT$bodyid)
  coln=c(coln,rep(i,length(idsT$bodyid)))
}

rowN=nNames8[2]
#rowN=nNames9[2]
coln=nNames7

#inpN=list(idsPc1[1:2],ids2[4:14],ids3)
#inpN=list(ids3)
#inpN=list(idsPc1[1:2])
#inpN=list(ids2[4:14])
inpN=as.character(ids4R)[2]
#inpN=as.character(ids4L)[2]

outpN=ids5
outpNn=nNames7
allN=data.frame(c(nNames7,nNames8[2]))
rownames(allN)=c(outpN,inpN)
###### main code

if (DS==TRUE){
  DN=list()
itr=1
for (i in inpN) {
  DN[[itr]]=as.character(getDownstreamsNR(i,nLayers,Thr)$partner)
  itr=itr+1
}

if (intersN==TRUE){
outN=DN[[1]]

for(i in DN){
  outN=intersect(outN[which(outN%in%DNs)],
                 i[which(i%in%DNs)])
  }
} else {

 # outN=unique(unlist(DN)[which(unlist(DN)%in%DNs)])
  outN=unique(c(unlist(DN)[which(unlist(DN)%in%DNs)],outpN))
  }
} else {
  outN=outpN
}
if (length(outN)>0){

if (length(outN)>295){
  outN=outN[1:295]
}

allPathsO=list()
itr=1
for (pLen in pL){
  allPathsO[[itr]]=pathsByEffic(unlist(inpN),outN,pLen,1)
  itr=itr+1
}
distComp=data.frame();

itr=1

for (j in allPathsO){
  distCompt=data.frame()
  ind1=1
  ind2=pL[itr]
  for (i in 1:(nrow(j)/pL[itr])) {
    distCompt[i,1]=sum(1/(j$weight[ind1:ind2]));
    # distCompt[i,2]=neuprint_get_neuron_names(j$to[ind2])
    distCompt[i,2]=(j$to[ind2])
    distCompt[i,3]=(j$from[ind1])
    ind1=ind1+pL[itr]
    ind2=ind2+pL[itr]

  }
distComp=rbind(distComp,distCompt)
itr=itr+1
}
# effDF=data.frame(matrix(nrow=length(inpN),ncol=(length(unique(neuprint_get_neuron_names(outN))))))
effDF=data.frame(matrix(nrow=length(inpN),ncol=(length(((outN))))))

# rownames(effDF)=c(unique(getTypeList(unlist(inpN))))
rownames(effDF)=as.character(inpN)
colnames(effDF)=outN

# distCompN=getTypeList(distComp$V3[1:nrow(distComp)])
distCompN=c(distComp$V3[1:nrow(distComp)])


for (i in rownames(effDF)){
  for (j in colnames(effDF)) {
    itr=1
    currDNs=c()
    while (itr <= nrow(distComp)){
    if (i==(distCompN[itr]) & (j==(distComp[itr,2]))) {
      currDNs=c(currDNs,1/distComp[itr,1])
      itr=itr+1
    } else {
        itr=itr+1
      }
    }
    effDF[which(rownames(effDF)==i),which(colnames(effDF)==j)]=1/sum(currDNs)
  }

}

rownames(effDF)=t(data.frame(rowN)[[1]])
colnames(effDF)[which(colnames(effDF)%in%outpN)]=coln
colnames(effDF)[which(!(outN%in%outpN))]=neuprint_get_neuron_names(outN[which(!(outN%in%outpN))])

# Function to calculate the cumulative variance explained up to a given percentile
cumulativeVarianceExplained <- function(values, percentile) {
  cutoff <- quantile(values, percentile)
  above_cutoff <- values[values >= cutoff]
  variance_explained <- sum((above_cutoff - mean(values))^2) / sum((values - mean(values))^2)
  return(variance_explained)
}

# Define your values
values <- sort(unlist(1/effDF),decreasing = TRUE)

# Define the range of percentiles to test
percentiles <- seq(0, 1, by = 0.01)

# Calculate the cumulative variance explained for each percentile
variance_explained <- sapply(percentiles, function(p) cumulativeVarianceExplained(values, p))

# Find the percentile at which the elbow point occurs
elbow_point <- findElbow(variance_explained)$elbow_point

findLargeCutoff <- function(values, percentile = 0.75) {
  # Compute the quantile
  cutoff <- quantile(values, percentile)

  # Return the cutoff
  return(cutoff)
}

#effThr <- findLargeCutoff(values, 0.85)


effThr=findElbow(sort(unlist(1/effDF),decreasing = TRUE))$elbow_point


from=c()
to=c()
label=c()
for (i in 1:nrow(effDF)){
  for (j in 1:ncol(effDF)) {
    if (effDF[i,j]=='Inf') {
      j=j+1
    } else if (colnames(effDF)[j]%in%outpNn){
      if (1/effDF[i,j]>=effThr) {
    from=c(from,rownames(effDF)[i])
    to=c(to,colnames(effDF)[j])
    label=c(label,effDF[i,j])
      }
    }
  }
}

nodes=data.frame(id=c(unique(from),unique(to)),
                 label=c(unique(from),unique(to)),
                 level=c(rep(1,length(unique(from))),rep(2,length(unique(to)))),
                 font.size=rep(20,length(c(unique(from),unique(to)))))

edges=data.frame(from=from,to=to,value=1/label)

#nodes=nodes[-which(unique(from)%in%c('DNut038_R','DNa04_R')),]

try(visNetwork(nodes, edges,width = '100%') %>%
  visEdges(shadow = TRUE,
           arrows ='to',
           selectionWidth = 50,
           color = list(color = "lightblue", highlight = "red"),
           length = 500)  %>%
  visOptions(selectedBy = "label",
             highlightNearest = TRUE,
             nodesIdSelection = TRUE) %>%
visHierarchicalLayout(levelSeparation = 250,nodeSpacing = 300) %>%
  #visNodes(font = list(size='200px',color='red',align='middle'))%>%
  visPhysics(stabilization = TRUE,repulsion = TRUE))

bThr=22
ex=rbind(allPathsO[[1]],allPathsO[[2]])
G=graph_from_edgelist(cbind(ex[,1],ex[,2]))
bC=data.frame(betweenness(G))
intN=rownames(bC)[which(bC[[1]]>bThr)]
count=c()
for (i in 1:nrow(bC)) {
  if (rownames(bC)[i]%in%c(inpN,outpN)){
    count=c(count,i)
  }
}
bCrown=rownames(bC)

bC=bC[-which(rownames(bC)%in%c(inpN,outpN)),1]
bCrown=bCrown[-count]
bC=data.frame(bC)
rownames(bC)=bCrown
#bThr=findElbowPoint(sort(bC[[1]],decreasing = TRUE))
#bThr=findElbow(sort(bC[[1]],decreasing = TRUE))$elbow_point
bThr=findLargeCutoff(sort(bC[[1]],decreasing = TRUE), 0.98)
#inMN=unique((rownames(bC)[which(bC[[1]]>bThr)])[1:length(inpN)])
# inMN=unique((rownames(bC)[which(bC[[1]]>bThr)]))[1:round(mean(c(length(ids5),length(ids4R))))]
inMN=unique((rownames(bC)[which(bC[[1]]>bThr)]))


condMat=1/effDF
cThr=as.integer(findElbowPoint(unlist(condMat)[which(unlist(condMat)!=0)]))

MN=c()

for (i in 1:nrow(effDF)){
  for (j in 1:ncol(effDF)){
    if ((condMat[i,j]>cThr)){
      MN=c(MN,colnames(condMat)[j])
    }
  }
}

MN=unique(MN[which(MN%in%outpNn)])

effDF2=grouped_adjacency(c(inMN),c(outpN),ingroup = 'name',outgroup = NULL)
# rownames(effDF)=c(unique(getTypeList(unlist(inpN))))
# rownames(effDF2)=unique(neuprint_get_neuron_names(inMN))
colnames(effDF2)=outpNn

#rowN=rowN[-which(rowN%in%c('DNut038_R','DNa04_R'))]
nodes=data.frame(id=c(rowN,unique(neuprint_get_neuron_names(inMN)),MN),
                 label=c(rowN,unique(neuprint_get_neuron_names(inMN)),MN),
                 level=c(rep(1,length(rowN)),rep(2,length(unique(neuprint_get_neuron_names(inMN)))),rep(3,length(MN)))
                 ,font.size=rep(15,length(c(rowN,unique(neuprint_get_neuron_names(inMN)),MN))))
from=c(rowN,unique(neuprint_get_neuron_names(inMN)))
to=c(unique(neuprint_get_neuron_names(inMN)),MN)
if ('TBD'%in%from){
from=from[-which(from=='TBD')]
}
if ('TBD'%in%to){
to=to[-which(to=='TBD')]
}

label=c()
fromt=c()
tot=c()

toiM=unique(neuprint_get_neuron_names(inMN))
#effThr2=findElbow(sort(unlist(effDF2),decreasing = TRUE))$elbow_point
#effThr2 <- findLargeCutoff(sort(unlist(effDF2),decreasing = TRUE), 0.75)
effThr2=Thr
for (i in (rownames(condMat))) {
  for (j in (to)){
    if (length(condMat[which(rownames(condMat)==i),which(colnames(condMat)==j)][[1]])>0){
     # if (condMat[which(rownames(effDF2)==i),which(colnames(effDF2)==j)][[1]]>=effThr2){
      label=c(label,condMat[which(rownames(condMat)==i),which(colnames(condMat)==j)][[1]])
      fromt=c(fromt,i)
      tot=c(tot,j)
     # }
    }
  }
}

for (i in (rownames(effDF2))) {
  for (j in (MN)){
     if (length(effDF2[which(rownames(effDF2)==i),which(colnames(effDF2)==j)][[1]])>0) {
       if (effDF2[which(rownames(effDF2)==i),which(colnames(effDF2)==j)][[1]]>effThr2){
      # } else if (from[i]%in%unique(neuprint_get_neuron_names(inMN))) {
      label=c(label,effDF2[which(rownames(effDF2)==i),which(colnames(effDF2)==j)][[1]])
      fromt=c(fromt,i)
      tot=c(tot,j)
       }
      }
  }
}

edges=data.frame(from=fromt,to=tot,value=label)

visNetwork(nodes, edges,width = '100%',height = ) %>%
  visEdges(shadow = TRUE,
           arrows ='to',
           selectionWidth = 20,
           color = list(color = "lightblue", highlight = "red"))  %>%
  visOptions(selectedBy = "label",
             highlightNearest = TRUE,
             nodesIdSelection = TRUE) %>%
  visHierarchicalLayout(levelSeparation = 200) %>%
  visPhysics(stabilization = FALSE,repulsion = FALSE)

edges2=edges
count=c()
for (i in 1:nrow(edges2)) {
  if (edges2$from[i]%in%rowN&edges2$to[i]%in%outpNn){
    count=c(count,i)
  }
}
edges2=edges2[-count,]

visNetwork(nodes, edges2,width = '100%',height = ) %>%
  visEdges(shadow = TRUE,
           arrows ='to',
           selectionWidth = 20,
           color = list(color = "lightblue", highlight = "red"),
           length=200)  %>%
  visOptions(selectedBy = "label",
             highlightNearest = TRUE,
             nodesIdSelection = TRUE) %>%
  visHierarchicalLayout(levelSeparation = 100) %>%
  visPhysics(stabilization = FALSE,repulsion = TRUE,)

# write.csv(effDF,"C:\\Users\\LabAdmin\\Desktop\\Connectomics\\Deven Connectomics\\condMatMN.csv", row.names = TRUE)

View(effDF)

RtoWrange<-colorRampPalette(c('white','red' ))
col=RtoWrange(100)
pheatmap(1/effDF,
         color = c(col),cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=15)


} else {
  paste('Specified nLayers does return intersecting neurons')
}
