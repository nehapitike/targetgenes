#Function to process chip sequence for a particular transcription factor and output target genes list
#Even gives the correlation values of the target genes with respect to TF for gastric, colon and breast cohorts.
chipseq <- function(file,bp,threshold)
{
peaksdata <- read.csv(file,header=TRUE)
use <- peaksdata[,7] > threshold
peaksdata <- peaksdata[use,]
refgene <- read.csv("refgenes.csv",header=TRUE)
targetgenes <- refgene[1,1:8]
peaksdata[,4] <- (peaksdata[,3] - peaksdata[,2])/2
peaksdata[,4] <- peaksdata[,2]+peaksdata[,4]
k<-1
targetgenes[,6] <- factor(targetgenes[,6],levels=levels(targetgenes[,1]))
#processing the chip seq data, using the genome data as reference
for ( i in 1:nrow(refgene))
{
	if (refgene[i,6] == "+")
	{
		startpos <- refgene[i,2]
		startpos <- startpos-bp
		endpos <- refgene[i,3]
		use <- peaksdata[,4] > startpos
		tempdata <- peaksdata[use,]
		use <- tempdata[,4] < endpos
		tempdata <- tempdata[use,]
		tempdata[,1] <- factor(tempdata[,1],levels=levels(refgene[,1]))
		use <- tempdata[,1] == refgene[i,1]
		tempdata <- tempdata[use,]
		if(nrow(tempdata) != 0)
		{
			for ( j in 1:nrow(tempdata))
			{
				targetgenes[k,1] <- refgene[i,1]
				targetgenes[k,2] <- refgene[i,2]
				targetgenes[k,3] <- refgene[i,3]
				targetgenes[k,4] <- refgene[i,4]
				targetgenes[k,5] <- tempdata[j,7]
				targetgenes[k,6] <- tempdata[j,1]
				targetgenes[k,7] <- tempdata[j,2]
				targetgenes[k,8] <- tempdata[j,3]
				k <- k+1
			}
		}
			
	}
	else
	{
		startpos <- refgene[i,2]
		endpos <- refgene[i,3]
		endpos <- endpos+bp
		use <- peaksdata[,4] > startpos
		tempdata <- peaksdata[use,]
		use <- tempdata[,4] < endpos
		tempdata <- tempdata[use,]
		tempdata[,1] <- factor(tempdata[,1],levels=levels(refgene[,1]))
		use <- tempdata[,1] == refgene[i,1]
		tempdata <- tempdata[use,]
		if(nrow(tempdata) != 0)
		{
			for ( j in 1:nrow(tempdata))
			{
				targetgenes[k,1] <- refgene[i,1]
				targetgenes[k,2] <- refgene[i,2]
				targetgenes[k,3] <- refgene[i,3]
				targetgenes[k,4] <- refgene[i,4]
				targetgenes[k,5] <- tempdata[j,7]
				targetgenes[k,6] <- tempdata[j,1]
				targetgenes[k,7] <- tempdata[j,2]
				targetgenes[k,8] <- tempdata[j,3]
				k <- k+1
			}
		}
			
	}

}
colnames(targetgenes) <- c("Chrom","genestart","geneend","name","signalvalue","Chrom1","tstart","tend")
write.csv(targetgenes, file = "TBX21targetgenes1.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

genenames <- unique(targetgenes[,4])
tg1 <- targetgenes[1,]
k <- 1
for (i in 1:length(genenames)){
	use <- targetgenes[,4] == genenames[i]
	temp <- targetgenes[use,]
	maxim <- max(temp[,5])
	use <- temp[,5] == maxim
	temp <- temp[use,]
	tg1[k,] <- temp[1,]
	k <- k+1	
}

write.csv(tg1, file = "TBX21targetgenes2.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#using ref2genesym for finding out the scientific gene name given gene symbol
genesym <- read.csv("ref2genesym.csv",header=TRUE)
tg1[,4] <- factor(tg1[,4],levels=levels(genesym[,1]))

for (i in 1:nrow(tg1))
{
	use <- tg1[i,4] == genesym[,1]
	temp <- genesym[use,]
	tg1[i,9] <- temp[1,2]

}
tg1 <- tg1[!duplicated(tg1[,9]),]
#writing the target gene list 
write.csv(tg1, file = "TBX21targetgenelist.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
finalans <- matrix(0,13,2)
finalans[1,1] <- c("ChipTG")
finalans[1,2] <- nrow(tg1)
#calculation correlation values for gastric cohort
gastric <- read.csv("gastricdatatot.csv",header=TRUE)
rown <- gastric[,1]
length(rown)
newgastric <- gastric[,2:ncol(gastric)]
new <- t(newgastric)
ncol(new)
gcorr <- matrix(0,ncol(new))
rownames(gcorr) <- rown
colnames(new) <- rown
for (i in 1:ncol(new))
{
        gcorr[i,1] <- cor(new[,i],new[,c("TBX21")])
}
head(gcorr)
#colnames(gcorr) <- c("GeneName","CorrVal")

write.csv(gcorr, file = "TBX21gastriccorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#calculation correlation values for breast cohort
breast <- read.csv("breastdatapro1.csv",header=TRUE)
rown <- breast[,1]
length(rown)
newbreast <- breast[,2:ncol(breast)]
new <- t(newbreast)
ncol(new)
gcorr <- matrix(0,ncol(new))
rownames(gcorr) <- rown
colnames(new) <- rown
for (i in 1:ncol(new))
{
        gcorr[i,1] <- cor(new[,i],new[,c("TBX21")])
}
head(gcorr)
#colnames(gcorr) <- c("GeneName","CorrVal")
write.csv(gcorr, file = "TBX21breastcorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#calculation correlation values for colon cohort
colon <- read.csv("colondatapro1.csv",header=TRUE)
rown <- colon[,1]
length(rown)
newcolon <- colon[,2:ncol(colon)]
new <- t(newcolon)
ncol(new)
gcorr <- matrix(0,ncol(new))
rownames(gcorr) <- rown
colnames(new) <- rown
for (i in 1:ncol(new))
{
        gcorr[i,1] <- cor(new[,i],new[,c("TBX21")])
}
head(gcorr)
#colnames(gcorr) <- c("GeneName","CorrVal")
write.csv(gcorr, file = "TBX21coloncorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#Finding out highly correlated genes using a 0.5 threshold and the common highly corr genes among all the cancer cohorts.
colon <- read.csv("TBX21coloncorr.csv",header=TRUE)
use <- is.na(colon[,2])
colon <- colon[!use,]
use <- colon[,2] > 0.5
temp <- colon[use,]
use <- colon[,2] < -0.3
temp1 <- colon[use,]
temp2 <- rbind(temp,temp1)
colnames(temp2) <- c("GeneName","CorrVal")
write.csv(temp2, file = "TBX21colonhighcorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
finalans[2,1] <- c("Colonhigh")
finalans[2,2] <- nrow(temp2)

comb <- temp2
breast <- read.csv("TBX21breastcorr.csv",header=TRUE)
use <- is.na(breast[,2])
breast <- breast[!use,]
use <- breast[,2] > 0.5
temp <- breast[use,]
use <- breast[,2] < -0.3
temp1 <- breast[use,]
temp2 <- rbind(temp,temp1)
colnames(temp2) <- c("GeneName","CorrVal")
write.csv(temp2, file = "TBX21breasthighcorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
finalans[3,1] <- c("Breasthigh")
finalans[3,2] <- nrow(temp2)

comb <- rbind(comb,temp2)
gastric <- read.csv("TBX21gastriccorr.csv",header=TRUE)
use <- is.na(gastric[,2])
gastric <- gastric[!use,]
use <- gastric[,2] > 0.5
temp <- gastric[use,]
use <- gastric[,2] < -0.3
temp1 <- gastric[use,]
temp2 <- rbind(temp,temp1)
colnames(temp2) <- c("GeneName","CorrVal")
write.csv(temp2, file = "TBX21gastrichighcorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
finalans[4,1] <- c("Gastrichigh")
finalans[4,2] <- nrow(temp2)

comb <- rbind(comb,temp2)
comb <- comb[!duplicated(comb[,1]),]

write.csv(comb, file = "TBX21combcancerhighcorr.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
finalans[5,1] <- c("Combhigh")
finalans[5,2] <- nrow(comb)

file <- matrix(c("TBX21colonhighcorr.csv","TBX21breasthighcorr.csv","TBX21gastrichighcorr.csv","TBX21combcancerhighcorr.csv"),1)
output <- matrix(c("TBX21colonchipintersect","TBX21breastchipintersect","TBX21gastricchipintersect","TBX21combcancerchipintersect"),1)
nchip <- read.csv("TBX21targetgenelist.csv",header=TRUE)
chip <- nchip[!duplicated(nchip[,9]),]
for ( j in 1:length(file))
{
corr <- read.csv(file[j],header=TRUE)
common <- intersect(chip[,9],corr[,1])
cm <- matrix(0,length(common),8)
for (i in 1:length(common))
{
        cm[i,1] <- common[i]
        use <- corr[,1] == common[i]
        cm[i,2] <- corr[use,2]
        use <- chip[,9] == common[i]
        cm[i,3] <- chip[use,5]
        cm[i,4] <- chip[use,1]
        cm[i,5] <- chip[use,2]
        cm[i,6] <- chip[use,3]
        cm[i,7] <- chip[use,7]
        cm[i,8] <- chip[use,8]
}
colnames(cm) <- c("GeneName","CorrVal","SignalVal","Chrm","genestart","geneend","transstart","transend")
write.csv(cm,paste(output[j],".csv",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

k <- j+5
finalans[k,1] <- c("intersect")
finalans[k,2] <- nrow(cm)


}
file <- matrix(c("TBX21colonchipintersect.csv","TBX21breastchipintersect.csv","TBX21gastricchipintersect.csv"),1)
colon <- read.csv(file[1],header=TRUE)
breast <- read.csv(file[2],header=TRUE)
gastric <- read.csv(file[3],header=TRUE)
cb <- intersect(colon[,1],breast[,1])
print("colon breast intersection:")
print(cb)

cg <-intersect(colon[,1],gastric[,1])
print("colon gastric intersection:")
print(cg)

bg <- intersect(gastric[,1],breast[,1])
print("breast gastric intersection:")
print(bg)

cbg <- intersect(cg,breast[,1])
print("colon breast gastric intersection:")
print(cbg)

inters <- cbind(cb,cg,bg,cbg)
colnames(inters) <- c("cb","cg","bg","cbg")


write.csv(inters,paste("TBX21intersectinggenes",".csv",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

finalans[10,1] <- c("cb")
finalans[10,2] <- length(cb)

finalans[11,1] <- c("cg")
finalans[11,2] <- length(cg)
finalans[12,1] <- c("bg")
finalans[12,2] <- length(bg)
finalans[13,1] <- c("cbg")
finalans[13,2] <- length(cbg)
print("Final Results:")
finalans
}
