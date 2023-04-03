require(readxl)
require(dplyr)
require(tidyr)
require(tibble)
require(ade4)
require(ggplot2)
require(vegan)
set.seed(123)
setwd("D:/OneDrive/UCT/PhD Material/BTA")
#Read in the species matrix and get rid of zero rows, there shouldn't be any now but the code will remain in case JAdH: 2022/04/20
L = read_xlsx("spp_test.xlsx",sheet="Clean species matrix")
L = 
  L %>%
  remove_rownames() %>%
  column_to_rownames(var="Subsample")
L_bin=L
L_bin[L>0]=1
#These two checking functions are not needed now
#colSums(L_bin)
L=L[,which(colSums(L)>0)]
#colSums(L)
rm(L_bin)

#Read in the Environmental matrix
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
R=dplyr::select(R,-all_of(c("X","ReefID","Distance_Cat","Area","Laser.point.distance","Observation_Duration","min_depth","Madrepora","Live.Lophelia","Coral.rubble","lophelia.rubbel","max_depth")))
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
Q=dplyr::select(Q,-all_of(dropping))
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

#The below code does the exact some as the function:
#LQ<-as.data.frame(FD::functcomp(as.matrix(Q_Fuzzy),as.matrix(L_trans)))
#summed_traits_by_station=matrix(0,dim(L_trans)[1],dim(Q_Fuzzy)[2])
#for(x in colnames(L_trans)){
#  print(x)
#  current_traits=as.matrix(Q_Fuzzy[x,])
#  current_species_obs=as.matrix(L_trans[,x])
#  colnames(current_species_obs)=x
#  rownames(current_species_obs)=rownames(L_trans)
#  temp=current_species_obs%*%current_traits
#  summed_traits_by_station = summed_traits_by_station+temp
#}
#summed_traits_by_station=summed_traits_by_station[,which(colSums(summed_traits_by_station)>0)]
#trans_summed_traits_by_station=prep.fuzzy.var(as.data.frame(summed_traits_by_station),col.blocks=blocks)

LQ<-as.data.frame(FD::functcomp(as.matrix(Q_Fuzzy),as.matrix(L_trans)))

rda_null=rda(LQ~1,data=Rstand[,c(1:3,5:11)])
plot(rda_null)
rda_full=rda(LQ~.+Condition(as.matrix(Rstand[,c(12,13)])),data=Rstand[,c(1:3,5:11)])
vif.cca(rda_full)
stepping=MASS::stepAIC(rda_null,
                       scope=list(upper=~mean_depth+grad2+Boulder+Cobble+Sand+Sandy.mud+Pebble+
                                    Dead.Lophelia+Live.Coral+Coral.Rubble+Condition(as.matrix(Rstand[,c(12,13)])),
                                  lower=~1),
                       direction = "both",
                       method="perm")
final_ord=rda(LQ ~ Live.Coral + Condition(as.matrix(Rstand[, c(12, 13)])) + mean_depth + Coral.Rubble + Pebble + Dead.Lophelia,data=Rstand[,c(1:11)])
vif.cca(final_ord)
temp=summary(final_ord)
temp$biplot
temp$cont[[1]][1:3,1:2]*100
ordiplot(rda_full,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6))
ordiplot(final_ord,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6))
text(final_ord,"species")

anova(final_ord) 
anova(final_ord, by="axis") ##First two axis are signficant
anova(final_ord, by="margin") ##All environmental variables are significant
RsquareAdj(final_ord)

RDA.site.scores = as.data.frame(scores(final_ord,display="sites",scaling = 2))
RDA.spp.scores = as.data.frame(scores(final_ord,display="species",scaling = 2))
RDA.env.scores = as.data.frame(scores(final_ord,display = "bp",scaling = 0)*ordiArrowMul(final_ord,display = "bp"))

require(ggrepel)
#It should be noted that if included, Sandy.mud which is mutuality exclusive with Sand in this dataset adds in the direction of LH4 (Burrowing, this makes logical sense)
pdf("CWM_RDAV3.pdf")
ggplot(RDA.site.scores,aes(RDA1,RDA2))+
  geom_point()+
  geom_point(data=RDA.spp.scores,aes(RDA1,RDA2),colour="Red",shape=8)+
  geom_text_repel(data=RDA.spp.scores,aes(RDA1,RDA2,label=rownames(RDA.spp.scores)),
                  colour="blue",cex=3,point.padding = NA,inherit.aes = FALSE,segment.colour = NA)+
  geom_segment(data=RDA.env.scores,aes(x=0,xend=RDA1,y=0,yend=RDA2),arrow=arrow(length=unit(0.2,"cm")),colour="Black")+
  geom_text_repel(data=RDA.env.scores,aes(x=1.1*RDA1,y=1.1*RDA2,label=rownames(RDA.env.scores)),
                  colour="Black",cex=4, segment.size = 0.25)+
  labs(x=paste0("RDA1\n(Proportion explained ",round(temp$cont[[1]][2,1]*100,2),"%)"),y=paste0("RDA2\n(Proportion explained ",round(temp$cont[[1]][2,2]*100,2),"%)"))+
  theme_minimal()+
  coord_fixed()
dev.off()
