require(readxl)
require(dplyr)
require(tidyr)
require(tibble)
require(ade4)
require(ggplot2)
set.seed(123)
setwd("D:/OneDrive/UCT/PhD Material/BTA")
#Read in the species matrix and get rid of zero rows, there shouldn't be any now but the code will remain in case JAdH: 2022/04/20
L = read_xlsx("spp_test.xlsx",sheet="Clean species matrix")
L = 
  L %>%
  remove_rownames() %>%
  column_to_rownames(var="Subsample")
#L=L[-grep("R1021_",row.names(L)),]
#L=L[-grep("R1124_",row.names(L)),]
L_bin=L
L_bin[L>0]=1
#These two checking fucntions are not needed now
#colSums(L_bin)
L=L[,which(colSums(L)>0)]
#colSums(L)
rm(L_bin)

#Read in the Enviromental matrix
R=read.csv("D:/OneDrive/UCT/PhD Material/Norway Data/Detailed Data/Env_all_V3_100m.csv")
#The values in the ReefID and Dist Cat are merged and converted to be the site names
row.names(R)=paste0(R$ReefID,"_",R$Distance_Cat)
#At site R1243 Madrepora was recorded as Clay.and.Spongespicules, this must be corrected
Madrepora=R[grepl("R1243_",rownames(R)),"Clay.and.Spongespicules"]
R[grepl("R1243_",rownames(R)),"Clay.and.Spongespicules"]=0
R$Madrepora=0
R[grepl("R1243_",rownames(R)),"Madrepora"]=Madrepora
#Since we are intrested in how cold-water coral reefs function all live coral is pooled
R$Live.Coral=R$Live.Lophelia+R$Madrepora
R$Coral.Rubble=R$Coral.rubble+R$lophelia.rubbel
#Drop the unsed columns from the Enviromental dataset
R=select(R,-all_of(c("X","ReefID","Distance_Cat","Area","Laser.point.distance","Observation_Duration","min_depth","Madrepora","Live.Lophelia","Coral.rubble","lophelia.rubbel","max_depth")))
rm(Madrepora)

#Read in the trait matrix and remove uneaded columns
Q = read_xlsx("traits_test.xlsx",sheet="Sheet2")
Q = 
  Q %>%
  dplyr::select(-contains(c("Valid_name","ID","Rank"))) %>%
  remove_rownames() %>%
  column_to_rownames(var="Taxon")

#Sort the species alphabetically and make sure the rows and columns of the tables match where required
L=L[,order(colnames(L))]
require(dplyr)
Q = 
  Q %>%
  arrange(row.names(Q))
R = 
  R %>%
  arrange(row.names(R))
obs_taxa = colnames(L)
obs_taxa[which(!(obs_taxa %in% row.names(Q)))]
L=L[,which((colnames(L) %in% row.names(Q)))]
Q=Q[which(row.names(Q) %in% colnames(L)),]
L=L[which(row.names(R) %in% row.names(L)),]
L_trans=L^(1/4)
L_trans=L_trans[,order(colnames(L_trans))]

#Standardise the R matrix, so that percentages and depth are centered at 0
Rstand=cbind(scale(R[,c(1,3,6:14)]),R[,c(4:5)])
#Perpare the fuzzy Q table
Q[is.na(Q)]=0
Q=Q[,colSums(Q,na.rm=T)>0]
#Removing social because so many taxa are colonial, it is causing tables L and Q to be unlinked
dropping=names(Q[,grepl("^SO[0-9]",colnames(Q))])
dropping=c(dropping,names(Q[,grepl("^LD[0-9]",colnames(Q))]))
Q=select(Q,-all_of(dropping))
#Automatically compute the block sizes for the traits
blocks=c(length(names(Q)[grep("^S[1-9]",names(Q))]),
         length(names(Q)[grep("^BF[1-9]",names(Q))]),
         length(names(Q)[grep("^SK[1-9]",names(Q))]),
#        length(names(Q)[grep("^SO[1-9]",names(Q))]),
         length(names(Q)[grep("^R[1-9]",names(Q))]),
#         length(names(Q)[grep("^LD[1-9]",names(Q))]),
         length(names(Q)[grep("^LH[1-9]",names(Q))]),
         length(names(Q)[grep("^MV[1-9]",names(Q))]),
         length(names(Q)[grep("^FH[1-9]",names(Q))]),
         length(names(Q)[grep("^SA[1-9]",names(Q))]))
#blocks=c(S=5,BF=7,SK=5,SO=3,R=4,LD=3,LH=6,MV=4,FH=6)
Q_Fuzzy=prep.fuzzy.var(Q,col.blocks=blocks)
Q_Fuzzy = 
  Q_Fuzzy %>%
  arrange(row.names(Q_Fuzzy))

#RLQ
Rstand=select(Rstand,all_of(c("mean_depth","grad2","Boulder","Cobble","Sand",
                              "Sandy.mud","Pebble","Dead.Lophelia","Live.Coral","Coral.Rubble")))
dudiL <- dudi.coa(L_trans, scannf = FALSE) ##densities
dudiR <- dudi.pca(Rstand, scannf = FALSE, nf = 2,row.w = dudiL$lw)##env
dudiQ <- dudi.pca(Q_Fuzzy, scannf = FALSE, nf = 2, row.w = dudiL$cw)##traits

biplot(dudiL,posieig="bottomright",clab.r  = 0.75,clab.c = 0.5) ##CA
biplot(dudiR,posieig="bottomright") ##PCA
biplot(dudiQ,posieig="bottomright")

vlt.rlq <- rlq(dudiR = dudiR, dudiL = dudiL, dudiQ = dudiQ, scannf = FALSE,nf=2)
plot(vlt.rlq)
randtest(vlt.rlq, nrepet = 4999) #If model 4 is not significant we can not proceed under the hypothesis tests as Q is not linked to L
#I have been able to make model 4 signficant after very carefully looking at the trait table
plot(randtest(vlt.rlq))
summary(vlt.rlq)

#Global significance
global_sig<-fourthcorner2(Rstand,L_trans, Q_Fuzzy,
                    modeltype = 6, p.adjust.method.G = "fdr", nrepet=4999) 
summary(global_sig)

four4c<-fourthcorner.rlq(vlt.rlq, modeltype=6,
                         nrepet = 4999,
                         p.adjust.method.G = "fdr",
                         p.adjust.method.D = "fdr",
                         typetest = "R.axes") 
plot(four4c, alpha=.05, stat="D2", type="biplot")
plot(four4c, alpha=.05, stat="D2")

rlq.env = as.data.frame(cbind((vlt.rlq$li),c("mean_depth","grad2","Boulder",
                                            "Cobble","Sand","Sandy.mud","Pebble",
                                            "Dead.Lophelia","Live.Coral","Coral.Rubble")))

names(rlq.env)[3]<-"Env"
head(rlq.env)

max_x=ceiling(max(abs(10*vlt.rlq$li["Axis1"])))/10+0.1
max_y=ceiling(max(abs(10*vlt.rlq$li["Axis2"])))/10+0.1

require(ggplot2)
require(ggrepel)
rlq.env.plot<-ggplot()+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2)+
  #The two geom_segment lines plot the arrows on the graph, the top command for the significant environmental variables and the second line for the grey non-significant ones.
  geom_segment(data=rlq.env[c(1,2,3,4,6,7,8,9,10),], aes(x=0, xend=Axis1, y=0, yend=Axis2),
               alpha=1, size=.4, color="black",
               arrow = arrow(length = unit(.2,"cm")))+
  geom_segment(data=rlq.env[c(5),], aes(x=0, xend=Axis1, y=0, yend=Axis2),
               alpha=1, size=.4, color="grey50",
               arrow = arrow(length = unit(.2,"cm")))+
  #Create labels for the variables (This code should be altered as I want the labels to not be on the arrow heads) This probably needs more work 20220425
  geom_label(data=rlq.env[c(1,2,3,4,6,7,8,9,10),], aes(x=Axis1, y=Axis2,
                                                  label=Env),color="black", size=5, vjust = "outward", hjust = "outward")+
  geom_label(data=rlq.env[c(5),], aes(x=Axis1, y=Axis2,
                                              label=Env),color="grey50",size=5, vjust = "outward", hjust = "outward")+
  theme_minimal()+
  theme(axis.title = element_text(size=18))+
  xlim(-max_x, max_x)+
  ylim(-max_y, max_y)+
  xlab(paste0("Projected inertia (",round((vlt.rlq$eig[1]/sum(vlt.rlq$eig))*100,2),"%): Axis 1"))+
  ylab(paste0("Projected inertia (",round((vlt.rlq$eig[2]/sum(vlt.rlq$eig))*100,2),"%): Axis 2"))
pdf("RLQ_envV4.pdf",width=8.27,height=5.83)
#windows(8.27,5.83)
rlq.env.plot
dev.off()

  
four2<-fourthcorner.rlq(vlt.rlq, modeltype=6,
                        nrepet = 4999,
                        p.adjust.method.G = "fdr",
                        p.adjust.method.D = "fdr",
                        typetest = "Q.axes") ##test modalities against environmental gradient
plot(four2, alpha=0.05, stat="D2") 
plot(four2, alpha=.05, stat="D2", type="biplot")

modsAbr=c("S1","S2", "S3","S4", "S5","S6",
           "BF1", "BF2", "BF3", "BF4", "BF5", "BF6","BF7",
           "SK1","SK2","SK3","SK4","SK5",
           "R1", "R2", "R3", "R4",
           "LH1", "LH2", "LH3","LH4", "LH5","LH6",
           "MV1", "MV2", "MV3", "MV4",
           "FH1", "FH3a","FH3p","FH4", "FH5", "FH6", "FH7",
           "SA1","SA2","SA3","SA4")
Traits=c(rep("Body size", 6), rep("Body form", 7),
         rep("Skeleton",5), rep("Reproduction", 4),
         rep("Living habit", 6), rep("Movement",4),
         rep("Feeding Habit", 7),rep("Substrate",4))
four2Coutput=as.data.frame(cbind(as.data.frame(four2$tabD2$obs), 
                                as.data.frame(four2$tabD2$names),
                                as.data.frame(four2$tabD2$adj.pvalue)))
names(four2Coutput)=c("Pearson","Relation", "adj.pvalue")
head(four2Coutput)
##make 2 data frames: for modality relationship with each axis
four2CSubsetaxis1 = as.data.frame(four2Coutput[grep("AxcR1", four2Coutput$Relation), ])
names(four2CSubsetaxis1)=c("PearsonAxis1", "RelationAxis1", "adj.pvalueAxis1")
four2CSubsetaxis2=as.data.frame(four2Coutput[grep("AxcR2", four2Coutput$Relation), ])
names(four2CSubsetaxis2)=c("PearsonAxis2", "RelationAxis2", "adj.pvalueAxis2")

rlq.mods = as.data.frame(vlt.rlq$c1)
four2CstatDoutput = as.data.frame(cbind(four2CSubsetaxis1[,c(1,3)],four2CSubsetaxis2[,c(1,3)], modsAbr,
                                        rlq.mods, Traits))

ModsSig1 = as.data.frame(subset(four2CstatDoutput, adj.pvalueAxis1<0.05))
ModsNotSig1 = as.data.frame(subset(four2CstatDoutput, adj.pvalueAxis1>0.05))
#Creating colour codes for the graph based on the modelity number, so all modality 1's are the same colour and so forth
ModsSig1$colour_code=gsub("[A-z]([0-9]*)","\\1",ModsSig1$modsAbr)
ModsNotSig1$colour_code=gsub("[A-z]([0-9]*)","\\1",ModsNotSig1$modsAbr)
#Manually changing the colour code for FH3a from 3 to 2, since FH3p is also 3
ModsNotSig1[which(ModsNotSig1$modsAbr=="FH3a"),"colour_code"]=2
ModsSig2 = as.data.frame(subset(four2CstatDoutput, adj.pvalueAxis2<0.05))
ModsNotSig2 = as.data.frame(subset(four2CstatDoutput, adj.pvalueAxis2>0.05))

CRLQ1<-ggplot() + 
  geom_boxplot(data=four2CstatDoutput, aes(y = Traits, x = CS1))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_label(data = ModsSig1,aes(y = Traits, x = CS1,
                                 label=modsAbr,
                                 size=abs(PearsonAxis2),fill=colour_code))+
  geom_point(data = ModsNotSig1,aes(y = Traits, x = CS1,colour=colour_code,),size=6)+
  scale_color_manual(values=c("#e69f00","#56b4e9", "#009e73","#F0E442","#d55e00","#cc79a7","#0072b2"))+
  scale_fill_manual(values=c("#F0E442"))+
  scale_size_continuous(range = c(4,6))+
  xlim(-0.5, 0.5)+
  theme_bw()+
  
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14, color="black"),
        legend.key = element_blank(),
        legend.position = "none")+
  xlab(paste0("RLQ1 (",round((vlt.rlq$eig[1]/sum(vlt.rlq$eig))*100,2),"%)"))+
  ylab("Trait")
pdf("RLQ_traitsV4.pdf",width=7.25,height=5.83)
CRLQ1
dev.off()
require(ggpubr)
C.rlq.outputs<-ggarrange(rlq.env.plot , CRLQ1,
                         ncol=1, nrow=2,align = "hv")
C.rlq.outputs
