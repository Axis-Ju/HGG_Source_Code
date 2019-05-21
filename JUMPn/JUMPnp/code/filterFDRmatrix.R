# inputs:
# combineTable_for_heatmap
# FDRcutoff
# rankcutoff

#-----------------------------------------------------------------------------
# To parse R parameters
# Credit to Nurcan Tuncbag from MIT
#args=(commandArgs(TRUE))
for (e in commandArgs(T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  var = ta[[1]][1]
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp,',')[[1]]
    }
    assign(var,temp)
    cat("assigned ",var," the value of |",temp,"|\n")
  } else {
    var_fields = strsplit(var,'-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var,TRUE)
    cat("assigned ",var," the value of TRUE\n")
  }
}

#-----------------------------------------------------------------------------

FDRcutoff=as.numeric(FDRcutoff)
rankcutoff=as.numeric(rankcutoff)

# load input file (output from combineTables.pl)
tb=read.table(combineTable_for_heatmap,head=T,sep="\t")

# set which column to be used
keyCol="_BH"

# prepare p value matrix
pm=tb[,grep(keyCol,names(tb))]
rownames(pm)=tb[,1]

# assign  'NA' to 1
smallValue=1e-2
for (i in 1:nrow(pm)) {
	for (j in 1:ncol(pm)) {
		if (is.na(pm[i,j]))
			pm[i,j]=1
	}
}

# select significant rows
# criteria #1: top 10 of each list
pmrank=pm
for (i in 1:ncol(pm)) {
	pmrank[,i]=rank(pm[,i])
}
pmrank$bestrank=apply(pmrank[,grep(keyCol,names(pmrank))],1,min)

# criteria #2: FDR
tb$minFDR=apply(tb[,grep(keyCol,names(tb))],1,min, na.rm = T)
tb$bestrankFDR=rep(1,nrow(tb))
for (i in 1:nrow(tb)) {
	a=tb[i,grep(keyCol,names(tb))]
	b=pmrank[i,1:(ncol(pmrank)-1)]
	f=as.numeric(a[which(b==pmrank[i,ncol(pmrank)])])
	if (!is.na(f)) {
		tb$bestrankFDR[i]=f
	}
}

# assign 0 p value as e-10; 
smallValue=1e-2
for (i in 1:nrow(pm)) {
	for (j in 1:ncol(pm)) {
		if (!(is.na(pm[i,j])) & pm[i,j]<smallValue)
			pm[i,j]=smallValue
	}
}

# select significant ones
#FDRcutoff=0.05
#rankcutoff=10
if (method=='rank') {
	pms=pm[tb$bestrankFDR<=FDRcutoff & pmrank$bestrank<=rankcutoff,]
} else {
	pms=pm[tb$minFDR<=FDRcutoff,]
}
dim(pms)

write.table(pms,file="significant_FDR_matrix.txt",quote=F,sep="\t")

# -log10(P)
pms=-log10(pms)

# heatmap
library(gplots)
my_palette <- colorRampPalette(c( "white","blue"))(n = 299)

# 
m=as.matrix(pms)
colnames(m)=gsub('enrichAnalysis_BH','',colnames(m)) # shorten the colnames

#windows(8,8)
pdf(file='FDR_heatmap.pdf',width=10,height=10)
heatmap.2(m,scale='n' ,
#labRow = gn[mergedColors==group],
Colv = FALSE,
density.info='n',
trace="none",         # turns off trace lines inside the heat map
  margins =c(7,28),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 

  dendrogram ="r",
colsep=0:ncol(m),rowsep=0:nrow(m),sepwidth=c(0.001,0.001),sepcolor="grey",keysize =1
)
dev.off()
