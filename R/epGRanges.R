## Equipartition GRanges:
## Splits a GRanges object into a specified number of parts with equal (total) width.
## Each part is a GRanges and all parts are returned as GRangesList object. 

## Alexandre Kuhn (alexandre.m.kuhn@gmail.com)

epGRanges<-function(gr,partitionNb) {
        if (partitionNb<1) stop("The number of partition should be at least 1.")
        else if (partitionNb==1) return(GRangesList(gr))
        else
        {
        grLength<-length(gr)
        cs<-cumsum(width(gr))

        partitionWidth<-floor(cs[grLength]/partitionNb) #last partition is longer
        #determine partitions' upper bounds (ub)
        upperbound<-seq(partitionWidth,cs[grLength],partitionWidth)
        upperbound[partitionNb]<-cs[grLength]

        #for each ub, identify the gr element to which it locates
        ubIndex<-sapply(upperbound,function(x){which(cs>=x)[1]})

        #for each ub, determine its coord. relative to the gr element it locates to 
        ubCoord<-rep(NA,partitionNb)
        if (ubIndex[1]==1) ubCoord[1]<-start(gr[1])+partitionWidth-1
        else ubCoord[1]<-start(gr[ubIndex[1]])+partitionWidth-cs[ubIndex[1]-1]-1
        for (i in 2:partitionNb) {
                if (ubIndex[i]==ubIndex[i-1]) ubCoord[i]<-ubCoord[i-1]+partitionWidth
                else ubCoord[i]<-start(gr[ubIndex[i]])+i*partitionWidth-cs[ubIndex[i]-1]-1
        }
        ubCoord[partitionNb]<-end(gr[grLength]) #correct coord. of last upper bound

        #use ub index and relative coordinates to generate GRangesList object
        grl<-GRangesList()
        if (ubIndex[1]==1) grl[[1]]<-restrict(gr[1],end=ubCoord[1])
        else grl[[1]]<-c(gr[1:(ubIndex[1]-1)],restrict(gr[ubIndex[1]],end=ubCoord[1]))
        for (i in 2:partitionNb) {
                if (ubIndex[i]==ubIndex[i-1]) grl[[i]]<-restrict(gr[ubIndex[i]],start=ubCoord[i-1]+1,end=ubCoord[i])
                else if (ubIndex[i]-ubIndex[i-1]==1) grl[[i]]<-c(restrict(gr[ubIndex[i-1]],start=ubCoord[i-1]+1),restrict(gr[ubIndex[i]],end=ubCoord[i]))
                else if (ubIndex[i]-ubIndex[i-1]>1) grl[[i]]<-c(restrict(gr[ubIndex[i-1]],start=ubCoord[i-1]+1),gr[(ubIndex[i-1]+1):(ubIndex[i]-1)],restrict(gr[ubIndex[i]],end=ubCoord[i]))
                else stop("Decreasing ubIndex")
        }
        return(grl)
        }
}

