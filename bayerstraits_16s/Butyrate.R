setwd("data/")
Blast=read.delim('Butyrate_filter.ws4.Chr.blastresult.txt',header=T)
#Blast=read.delim('all.butyrate.blast.txt',header=T)
library('ggplot2')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(Blast, aes(as.numeric(as.matrix(Blast$identity)))) + 
  geom_density(alpha=0.5,adjust=1,colour = cbPalette[8],fill=cbPalette[8],size=0.5) +###stat_density
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
#Blast=read.delim('all.butyrate.blast.80.txt',header=T)
#Meta=read.delim('Geba_meta.txt',header=T)
Meta=read.delim('assembly_summary_genbank_taxon.txt.notnormalized2.new.pathogen.habitat.8000.txt',header=T)
Gene=read.delim('Butyrate.pro.mapping.txt',header=T)
#Blast=merge(Blast,Meta,by='Genome',all.x=T)
#Blast=merge(Blast,Meta,by.x='Genome2',by.y='Genome',all.x=T,all.y=T)
Blast=merge(Blast,Meta,by.x='Genome',by.y='assembly_accession',all=T)
Blast=merge(Blast,Gene,by.x='Subject_id',by.y='Gene',all.x=T)

#Genome=unique(Blast$Genome2)
Genome=unique(Blast$Genome)
Genome_but=data.frame(Genome)
Genome_but$Butyrate=0
Genome_but=as.matrix(Genome_but)
for(i in 1:length(Genome))
{
  #temp=Blast[which(as.matrix(Blast$Genome2)==Genome[i]),] 
  temp=Blast[which(as.matrix(Blast$Genome)==Genome[i]),] 
  if('But' %in% as.matrix(temp$Function)){
    Genome_but[i,2]=1
  }else if(('Buk' %in% as.matrix(temp$Function))&&
           ('Ptb' %in% as.matrix(temp$Function))){
    Genome_but[i,2]=1
    # Ato not considered yet
  }
}
Genome_but=data.frame(Genome_but)
Genome_but=merge(Genome_but,Meta,by.x='Genome',
                 by.y='assembly_accession',all.x=T,all.y=T)

write.table(Blast,'all.butyrate.blast.meta.ws4.txt',quote=F,sep='\t',row.names = F)
write.table(Genome_but,'Genome_butyrate_noato.ws4.txt',quote=F,sep='\t',row.names = F)

Genome_but=data.frame(Genome)
Genome_but$Butyrate=0
Genome_but=as.matrix(Genome_but)
for(i in 1:length(Genome))
{
  temp=Blast[which(as.matrix(Blast$Genome2)==Genome[i]),] 
  if('But' %in% as.matrix(temp$Function)){
    Genome_but[i,2]=1
  }else if(('Buk' %in% as.matrix(temp$Function))&&
           ('Ptb' %in% as.matrix(temp$Function))){
    Genome_but[i,2]=1
  }else if(('AtoA' %in% as.matrix(temp$Function))&&
           ('AtoD' %in% as.matrix(temp$Function))){
    Genome_but[i,2]=1
  }
}
Genome_but=data.frame(Genome_but)
Genome_but=merge(Genome_but,Meta,by.x='Genome',by.y='Genome',all.x=T,all.y=T)
write.table(Genome_but,'Genome_butyrate_withato.txt',quote=F,sep='\t',row.names = F)

###butyrate hmm
Meta=read.delim('Geba_meta.txt',header=T)
Genome=read.delim('Butyrate.hmm',header=T)
Genome_but=data.frame(unique(Meta$Genome))
Genome_but$butyrate=0
Genome_but=as.matrix(Genome_but)
colnames(Genome_but)=c('Genome','Genome')

for(i in 1:nrow(Meta))
{
  temp=Genome[which(as.matrix(Genome$Genome)==Meta$Genome[i]),] 
  rownum=which(Genome_but[,1]==toString(as.matrix(Meta$Genome[i])))
  
  if('But' %in% as.matrix(temp$Gene)){
    Genome_but[rownum,2]=1
  }else if(('Buk' %in% as.matrix(temp$Gene))&&
           ('Ptb' %in% as.matrix(temp$Gene))){
    Genome_but[rownum,2]=1
  }else if(('AtoA' %in% as.matrix(temp$Gene))&&
           ('AtoD' %in% as.matrix(temp$Gene))){
    Genome_but[rownum,2]=1
  }
}
colnames(Genome_but)=c('Genome','Butyrate')
Genome_but=data.frame(Genome_but)
Genome_but=merge(Genome_but,Meta,by.x='Genome',by.y='Genome',all.x=T,all.y=T)
write.table(Genome_but,'Genome_butyrate_noato_hmm.txt',quote=F,sep='\t',row.names = F)

#merge blast and hmm
Genome_blast=read.delim('Genome_butyrate_noato.txt',header=T)
Genome_but=merge(Genome_but,Genome_blast,by='Genome',all=T)
write.table(Genome_but,'Genome_butyrate_noato_hmm_blast.txt',quote=F,sep='\t',row.names = F)

###NR

Meta=read.delim('Geba_meta.txt',header=T)
NR=read.delim('NR.hmm.result',header=T)
NR_for=data.frame(unique(Meta$Genome))
NR_for$napA=0
NR_for$narG=0
NR_for$nirS=0
NR_for$nirK=0
NR_for=as.matrix(NR_for)
for(i in 1:nrow(NR))
{
  colnum=which(colnames(NR_for)==toString(as.matrix(NR$Gene[i])))
  rownum=which(NR_for[,1]==toString(as.matrix(NR$Genome[i])))
  NR_for[rownum,colnum]=1
}
colnames(NR_for)[1]='Genome'
NR_for=data.frame(NR_for)

NR_for=merge(NR_for,Meta,by.x='Genome',by.y='Genome',all.x=T,all.y=T)
write.table(NR_for,'NR.txt',quote=F,sep='\t',row.names = F)

###NR
setwd("C:/Users/zhang/Desktop/News/project/Gut/models/butyrate/temp_file/ws4_chr/")
Meta=read.delim('assembly_summary_genbank_taxon.txt.notnormalized2.new.pathogen.habitat.8000.txt',header=T)
Genome_NR=read.delim('ws4.NR.HMM.filter',header=T)
Genome_NR=merge(Genome_NR,Meta,by='Genome',all=T)
CG_16S=read.delim('Genome_CG_16S.txt',header=T)
Genome_NR=merge(Genome_NR,CG_16S,by='Genome',all.y=T)

write.table(Genome_NR,'ws4.NR.HMM.meta.filter',quote=F,sep='\t',row.names = F)
#Genome_NR=read.delim('ws4.NR.HMM.meta.filter',header=T)
CG_16S=read.delim('Genome_CG_16S.txt',header=T)
CG_16S$NR=0
CG_16S$nap_nar=0
CG_16S$nir=0
CG_16S=as.matrix(CG_16S)

for(i in 1:nrow(CG_16S))
{
  temp=Genome_NR[which(as.matrix(Genome_NR$Genome)==CG_16S[i,1]),] 
  
  if(any(as.matrix(temp$NR) %in% c('napA','narG')))
    if(any(as.matrix(temp$NR) %in% c('nirS','nirK')))
      CG_16S[i,3]=1
    if(any(as.matrix(temp$NR) %in% c('napA','narG')))
      CG_16S[i,4]=1
    if(any(as.matrix(temp$NR) %in% c('nirS','nirK')))
      CG_16S[i,5]=1
}
CG_16S=data.frame(CG_16S)
CG_16S=merge(CG_16S,Meta,by.x='Genome',by.y='Genome',all.x=T)
write.table(CG_16S,'Genome_CG_16S_NR.txt',quote=F,sep='\t',row.names = F)

###SRB
setwd("C:/Users/zhang/Desktop/News/project/Gut/models/butyrate/temp_file/ws4_chr/")
Meta=read.delim('assembly_summary_genbank_taxon.txt.notnormalized2.new.pathogen.habitat.8000.txt',header=T)
Genome_SRB=read.delim('ws4.SRG.blastp.filter',header=T)
Genome_SRB=merge(Genome_SRB,Meta,by='Genome',all=T)
CG_16S=read.delim('Genome_CG_16S.txt',header=T)
Genome_SRB=merge(Genome_SRB,CG_16S,by='Genome',all.y=T)

write.table(Genome_SRB,'ws4.SRG.blastp.meta.filter',quote=F,sep='\t',row.names = F)
Genome_SRB=read.delim('ws4.SRG.blastp.meta.filter.flter',header=T)
CG_16S=read.delim('Genome_CG_16S.txt',header=T)
CG_16S$SRB=0
CG_16S$SAT=0
CG_16S$arp=0
CG_16S$dsr_dsv=0
CG_16S=as.matrix(CG_16S)

for(i in 1:nrow(CG_16S))
{
  temp=Genome_SRB[which(as.matrix(Genome_SRB$Genome)==CG_16S[i,1]),] 
  
  if('SAT' %in% as.matrix(temp$Subject))
    if(any(as.matrix(temp$Subject) %in% c('arpA','arpB')))
      if(any(as.matrix(temp$Subject) %in% c('dsrA','dsvC','dsvB','dsrB')))
        CG_16S[i,3]=1
      if('SAT' %in% as.matrix(temp$Subject))
        CG_16S[i,4]=1
      if(any(as.matrix(temp$Subject) %in% c('arpA','arpB')))
        CG_16S[i,5]=1
      if(any(as.matrix(temp$Subject) %in% c('dsrA','dsvC','dsvB','dsrB')))
        CG_16S[i,6]=1
}
CG_16S=data.frame(CG_16S)
CG_16S=merge(CG_16S,Meta,by.x='Genome',by.y='Genome',all.x=T)
write.table(CG_16S,'Genome_CG_16S_SRB.txt',quote=F,sep='\t',row.names = F)

###But
Genome_But=read.delim('ws4.But.blast.filter',header=T)
CG_16S=read.delim('Genome_CG_16S.txt',header=T)
CG_16S$Both=0
CG_16S$Either=0
CG_16S$But=0
CG_16S$Buk_Ptb=0
CG_16S$Ato=0
CG_16S=as.matrix(CG_16S)

for(i in 1:nrow(CG_16S))
{
  temp=Genome_But[which(as.matrix(Genome_But$Genome)==CG_16S[i,1]),] 
  
  if('But' %in% as.matrix(temp$Function))
    if(any(as.matrix(temp$Function) %in% c('Buk')))
      if(any(as.matrix(temp$Function) %in% c('Ptb')))
        CG_16S[i,3]=1
      if('But' %in% as.matrix(temp$Function) ||
         (any(as.matrix(temp$Function) %in% c('Buk')) 
          && any(as.matrix(temp$Function) %in% c('Ptb'))))
        CG_16S[i,4]=1
      if('But' %in% as.matrix(temp$Function))
        CG_16S[i,5]=1
      if(any(as.matrix(temp$Function) %in% c('Buk','Ptb')))
        CG_16S[i,6]=1
      if(any(as.matrix(temp$Function) %in% c('AtoA','AtoD')))
        CG_16S[i,7]=1
      
}
CG_16S=data.frame(CG_16S)
CG_16S=merge(CG_16S,Meta,by.x='Genome',by.y='Genome',all.x=T)
write.table(CG_16S,'Genome_CG_16S_But.txt',quote=F,sep='\t',row.names = F)

# ibd results merge
setwd("C:/Users/zhang/Desktop/News/project/Gut/models/butyrate/ibd/")
ID=read.delim('out',header=T)
Mapping=read.delim('Tree_mapping_file.txt',header=T)
Genome=read.delim('assembly_summary_genbank_taxon.txt.notnormalized2.new.pathogen.habitat.8000.txt',header=T)
ID=merge(ID,Mapping,by='ID',all.x=T)
ID=merge(ID,Genome,by.x='ID',by.y='Genome',all.x=T)
write.table(ID,'out.merge.txt',quote=F,sep='\t',row.names = F)

#ibd results
setwd("C:/Users/zhang/Desktop/News/project/Gut/models/butyrate/ibd/")
Traits=read.delim('ibd_huttenhower.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.abu',header=F)
Traits=data.frame(t(as.matrix(Traits[c(1,nrow(Traits)),])))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Metadata=read.delim('ibd_huttenhower.metadata.txt',header=T)
Traits=read.delim('ibd_huttenhower.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.t.abu',header=T)
Traits=merge(Traits,Metadata,by.x='OTU_ID',by.y='X.SampleID',all.x=T)
write.table(Traits,'ibd_huttenhower.metadata_abu.txt',quote=F,sep='\t',row.names = F)

Traits=read.delim('ibd_huttenhower.metadata_abu.txt',header=T)
Traits=Traits[order(Traits$DiseaseState),]

library(ggplot2)
dataplot=data.frame(
  x=seq(1,nrow(Traits),by=1),
  y=as.numeric(as.matrix(Traits$Traits)),
  colours=as.matrix(Traits$DiseaseState)
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y)),
                    color=as.matrix(dataplot$colours))) + ##or x2,y2, or x3,y3
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  # scale_y_continuous(limits = c(-5.2, 1),breaks = seq(-5, 1, by = 1))+
  scale_color_brewer(palette="Set1")

Traits=read.delim('ibd_engstrand_maxee.metadata_abu.txt',header=T)
unique(Traits$DiseaseState)
#ks test
healthy=as.matrix(Traits$Traits[which(Traits$DiseaseState=='nonIBD')])
unhealthy=as.matrix(Traits$Traits[which(Traits$DiseaseState %in% c('CD'))])
ks.test(healthy,
        unhealthy)
#wilcox test

#dataplot$colours[which(dataplot$colours %in% c('CD','UC','IBDundef'))]='CD'
dataplot=Traits[which(Traits$DiseaseState %in% c('CD','nonIBD')),]
dataplot=data.frame(
  Traits=as.numeric(as.matrix(dataplot$Traits)),
  DiseaseState=as.matrix(dataplot$DiseaseState)
)
wilcox.test(Traits ~ DiseaseState, dataplot) 

dataplot=data.frame(
  x=seq(1,nrow(dataplot),by=1),
  y=as.numeric(as.matrix(dataplot$Traits)),
  colours=as.matrix(dataplot$DiseaseState)
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y)),
                    color=as.matrix(dataplot$colours))) + ##or x2,y2, or x3,y3
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  # scale_y_continuous(limits = c(-5.2, 1),breaks = seq(-5, 1, by = 1))+
  scale_color_brewer(palette="Set1")

#test data
setwd("C:/Users/zhang/Desktop/News/project/Gut/models/butyrate/ibd/")
OTUs=read.delim('ibd_alm.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.txt',header=F)
Table=read.delim('ibd_alm.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.abu',header=F)
Table=merge(Table,OTUs,by='V1',all.x=T)
Table=rbind(as.matrix(Table),matrix(0,nrow=1,ncol=ncol(Table)))
Table[which(is.na(Table[,ncol(Table)])),ncol(Table)]=0
for(i in 1:(nrow(Table)-3))
{
  Table[nrow(Table),2:(ncol(Table)-1)]=
    as.numeric(Table[nrow(Table),2:(ncol(Table)-1)])+
    as.numeric(Table[i,2:(ncol(Table)-1)])*as.numeric(Table[i,ncol(Table)])
}

write.table(Table,'ibd_alm.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.R.abu',quote=F,sep='\t',row.names = F)
Traits=read.delim('ibd_alm.otu_seqs.100.fasta.filter.align.16S.nwk.infertraits.R.abu',header=F)
Traits=data.frame(t(as.matrix(Traits)))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Metadata=read.delim('ibd_alm.metadata.txt',header=T)
Traits=merge(Traits,Metadata,by.x='OTU_ID',by.y='X',all.x=T)
write.table(Traits,'ibd_alm.R.metadata_abu.txt',quote=F,sep='\t',row.names = F)

#Traits=read.delim('ibd_huttenhower.metadata_abu.txt',header=T)
Traits=Traits[order(Traits$DiseaseState),]

library(ggplot2)
dataplot=data.frame(
  x=seq(1,nrow(Traits),by=1),
  y=as.numeric(as.matrix(Traits$Traits_R)),
  colours=as.matrix(Traits$DiseaseState)
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y)),
                    color=as.matrix(dataplot$colours))) + ##or x2,y2, or x3,y3
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  # scale_y_continuous(limits = c(-5.2, 1),breaks = seq(-5, 1, by = 1))+
  scale_color_brewer(palette="Set1")


##SRB and NR
setwd("../16S/100wells/")
Metadata=read.delim('enigma.meta.core',header=T)
Traits=read.delim('eotus.fasta.filter.align.16S.nwk.infertraits.SRB.abu',header=F)
Traits=data.frame(t(as.matrix(Traits[c(1,nrow(Traits)),])))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Traits=merge(Traits,Metadata[,c(1,47,51,59)],by.x='OTU_ID',by.y='X',all.x=T)
#write.table(Traits,'eotus.fasta.filter.align.16S.nwk.infertraits.SRB.abu.meta',quote=F,sep='\t',row.names = F)
Traits=data.frame(Traits)

library(ggplot2)
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  y=as.numeric(as.matrix(Traits$SO4_mgL))
  #y=as.numeric(as.matrix(Traits$Sulfide_mgL))
  #x=as.numeric(as.matrix(Traits$S_all_mgL))
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y))))+
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_color_brewer(palette="Set1")+
  xlab('log10(Sulfate Reducting Bacteria Abundance)')+
  ylab('SO4_mgL(mgL)')
#+!scale_y_continuous(limits = c(-5.2, 1),breaks = seq(-5, 1, by = 1))+
# + scale_color_brewer(palette="Set1")
#ks test
Low=as.numeric(as.matrix(dataplot$x[which(dataplot$y>0.5)]))
High=as.numeric(as.matrix(dataplot$x[which(dataplot$y<=0.5)]))
ks.test(Low,
        High)
#wilcox test
dataplot=data.frame(
  y=log10(as.numeric(as.matrix(Traits$Traits))),
  x=as.numeric(as.matrix(Traits$Sulfide_mgL))
)
dataplot=as.matrix(dataplot)
dataplot[which(dataplot[,2]>0.5),2]='High'
dataplot[which(dataplot[,2]!='High'),2]='Low'
dataplot=dataplot[which(!(is.na(dataplot[,2]))),]
dataplot=data.frame(
  Traits=as.numeric(as.matrix(dataplot[,1])),
  Sulfide_mgL=as.matrix(dataplot[,2])
)
wilcox.test(Traits ~ Sulfide_mgL, dataplot) 

#---
#Metadata=read.delim('enigma.meta.core',header=T)
Traits=read.delim('eotus.fasta.filter.align.16S.nwk.infertraits.SRB.dsr.abu',header=F)
Traits=data.frame(t(as.matrix(Traits[c(1,nrow(Traits)),])))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Traits=merge(Traits,Metadata[,c(1,47,51,59)],by.x='OTU_ID',by.y='X',all.x=T)
#write.table(Traits,'eotus.fasta.filter.align.16S.nwk.infertraits.SRB.dsr.abu.meta',quote=F,sep='\t',row.names = F)
Traits=data.frame(Traits)

library(ggplot2)
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  y=as.numeric(as.matrix(Traits$SO4_mgL))
  #y=as.numeric(as.matrix(Traits$Sulfide_mgL))
  #x=as.numeric(as.matrix(Traits$S_all_mgL))
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y))))+
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_color_brewer(palette="Set1")+
  xlab('log10(Sulfate Reducting Bacteria Abundance (dsr))')+
  ylab('SO4(mgL)')
#+!scale_y_continuous(limits = c(-5.2, 1),breaks = seq(-5, 1, by = 1))+
# + scale_color_brewer(palette="Set1")


#ks test
Low=as.numeric(as.matrix(dataplot$x[which(dataplot$y>10)]))
High=as.numeric(as.matrix(dataplot$x[which(dataplot$y<=10)]))
ks.test(Low,
        High)
#wilcox test
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  y=as.numeric(as.matrix(Traits$SO4_mgL))
)
dataplot=as.matrix(dataplot)
dataplot[which(dataplot[,2]>1),2]='High'
dataplot[which(dataplot[,2]!='High'),2]='Low'
dataplot=dataplot[which(!(is.na(dataplot[,2]))),]
dataplot=data.frame(
  Traits=as.numeric(as.matrix(dataplot[,1])),
  SO4_mgL=as.matrix(dataplot[,2])
)
wilcox.test(Traits ~ SO4_mgL, dataplot) 

#NR
setwd("../16S/100wells/")
Metadata=read.delim('enigma.meta.core',header=T)
Traits=read.delim('eotus.fasta.filter.align.16S.nwk.infertraits.NR.abu',header=F)
Traits=data.frame(t(as.matrix(Traits[c(1,nrow(Traits)),])))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Traits=merge(Traits,Metadata[,c(1,39,40)],by.x='OTU_ID',by.y='X',all.x=T)
#write.table(Traits,'eotus.fasta.filter.align.16S.nwk.infertraits.NR.abu.meta',quote=F,sep='\t',row.names = F)
Traits=data.frame(Traits)

library(ggplot2)
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  #x=as.numeric(as.matrix(Traits$N2_mM))
  y=as.numeric(as.matrix(Traits$NO3_mgL))
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y))))+
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab('log10(Nitrate Reducting Bacteria Abundance)')+
  ylab('NO3(mgL)')


#ks test
Low=as.numeric(as.matrix(dataplot$x[which(dataplot$y>1000)]))
High=as.numeric(as.matrix(dataplot$x[which(dataplot$y<=1000)]))
ks.test(Low,
        High)
#wilcox test
dataplot=data.frame(
  y=log10(as.numeric(as.matrix(Traits$Traits))),
  #x=as.numeric(as.matrix(Traits$N2_mM))
  x=as.numeric(as.matrix(Traits$NO3_mgL))
)
dataplot=as.matrix(dataplot)
dataplot[which(dataplot[,2]>10),2]='High'
dataplot[which(dataplot[,2]!='High'),2]='Low'
dataplot=dataplot[which(!(is.na(dataplot[,2]))),]
dataplot=data.frame(
  Traits=as.numeric(as.matrix(dataplot[,1])),
  NO3_mgL=as.matrix(dataplot[,2])
)
wilcox.test(Traits ~ NO3_mgL, dataplot) 

---
  #Metadata=read.delim('enigma.meta.core',header=T)
  Traits=read.delim('eotus.fasta.filter.align.16S.nwk.infertraits.NR.nir.abu',header=F)
Traits=data.frame(t(as.matrix(Traits[c(1,nrow(Traits)),])))
colnames(Traits)=as.matrix(Traits[1,])
Traits=Traits[-1,]
Traits=merge(Traits,Metadata[,c(1,39,40)],by.x='OTU_ID',by.y='X',all.x=T)
#write.table(Traits,'eotus.fasta.filter.align.16S.nwk.infertraits.NR.NAP.NAR.abu.meta',quote=F,sep='\t',row.names = F)
Traits=data.frame(Traits)

library(ggplot2)
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  y=as.numeric(as.matrix(Traits$NO3_mgL))
)

ggplot(dataplot,aes(x=as.numeric(as.matrix(dataplot$x)),
                    y=as.numeric(as.matrix(dataplot$y))))+
  geom_point(shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab('log10(Nitrate Reducting Bacteria Abundance (nir))')+
  ylab('NO3(mgL)')


#ks test
Low=as.numeric(as.matrix(dataplot$x[which(dataplot$y>500)]))
High=as.numeric(as.matrix(dataplot$x[which(dataplot$y<=500)]))
ks.test(Low,
        High)
#wilcox test
dataplot=data.frame(
  x=log10(as.numeric(as.matrix(Traits$Traits))),
  y=as.numeric(as.matrix(Traits$NO3_mgL))
)
dataplot=as.matrix(dataplot)
dataplot[which(dataplot[,2]>1),2]='High'
dataplot[which(dataplot[,2]!='High'),2]='Low'
dataplot=dataplot[which(!(is.na(dataplot[,2]))),]
dataplot=data.frame(
  Traits=as.numeric(as.matrix(dataplot[,1])),
  NO3_mgL=as.matrix(dataplot[,2])
)
wilcox.test(Traits ~ NO3_mgL, dataplot) 

setwd("data/")
CG1 = read.delim('Genome_butyrate_both.txt', header = F)
CG2 = read.delim('Genome_butyrate_either.txt', header = F)
CG3 = read.delim('Genome_butyrate_Ato.txt', header = F)
CG4 = read.delim('Genome_butyrate_buk_ptb.txt', header = F)
CG5 = read.delim('Genome_butyrate_but.txt', header = F)
CG = merge(CG1,CG2,by='V1',all=T)
CG3 = merge(CG3,CG4,by='V1',all=T)
colnames(CG3)=c('V1','1','2')
CG3 = merge(CG,CG5,by='V1',all=T)
colnames(CG3)=c('Genome','Ato','buk_ptb','but')
#CG = merge(CG,CG3,by='V1',all=T)
colnames(CG)=c('Genome','both','either','Ato','buk_ptb','but')
CG$difference = CG$V2.x - CG$V2.y
table(CG$difference)
table(CG$both)
table(CG$either)