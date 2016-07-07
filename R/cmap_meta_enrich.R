
#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝
options(stringsAsFactors=F);library(colorout)
rm(list=ls());ls()

clicky.run<-function(module_genes_list,module_bkrnd,clicky_dir,dat_descr='',id_type='hsapiens__gene_symbol'){
  #hsapiens__ensembl_gene_stable_id
#cat('\tUSE : creates temp directory and runs clicky enrichments on a list of module gene names - list$ModuleName$ModuleGenes')
#  system(paste0("mkdir ",clicky_dir,"/tmp_in_files/"))

##  clunky, but saves having to paste the clickly path everywhere..
  if(class(module_genes_list)!='list'){
    stop('\n\tERROR :\trequire module_genes_list as list, each element -  vector or matrix of module genes')
  }



  current_dir=getwd()
  working_dir=paste0(clicky_dir,"/tmp_in_files/")
  setwd(clicky_dir)
#  getwd()
  
  write.table(module_bkrnd,file="tmp_in_files/tmp_bkgrnd.txt", quote=F, sep="\t", row.names=F, col.names=F)
  

  for(imod in 1:length(module_genes_list)){
    cat('\t',names(module_genes_list)[imod],dat_descr,'\t',length(module_genes_list[[imod]]),'\t',length(module_bkrnd),'\n')

    write.table(module_genes_list[[imod]],file="tmp_in_files/tmp_genelist.txt", quote=F, sep="\t", row.names=F, col.names=F)
  
  system(paste0("python ",clicky_dir,"/bin/runner_modified_commandline_osx.py "
#  system(paste0("python ",clicky_dir,"/bin/clicky_python_commandline_linux.py "
                ,working_dir,"tmp_genelist.txt "
                ,working_dir,"tmp_bkgrnd.txt "
                ,working_dir," "
                ,paste0(names(module_genes_list)[imod],'_',dat_descr)
                ," hsapiens ",id_type," ",id_type," 0.05"))
         
# rm(enrich_list)
# enrich_list=gestalt_read(enrich_type=c("GO","Commons","KEGG","WIKI"),
#                          mod_names= names(module_genes_list)[imod],
#                          in_path=paste0(outpath,"/"),
#                          out_path=paste0(outpath,"/graphics"),
#                          dat_descr="SpeDE_eachRegion")

# pdf(paste0(outpath,"/graphics/",names(module_genes_list)[imod],"_graph.pdf"),width=15, height=25)
# gestalt_plot(enrich_list,p_threshold=0.05)
# dev.off()

  }
  setwd(current_dir)
####  system(paste0("rm -rf ",clicky_dir,"/tmp_in_files/"))  ## if full plotting is done in a different directory, ie wrapping PDF around this function, with all 3 clicky functions enabled as above then can re-enable
}

idconvert.ensg<-function(dat_lis){
##  conveting list of vectors containing only official HUGO gene names (rownames) to ENSG

	for(ilis in names(dat_lis)){
		dat_lis[[ilis]]=idmap.ensg[idmap.ensg$gene%in%(dat_lis[[ilis]]),'ids']
	}
	return(dat_lis)
}


#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝
library(R.helper)


##  full data on differential expression using affy
Load('/Users/ks/Dropbox/PROJ/cmap/dtb/genes.difexpr.drug.cmap.Rdata')
	Head(degb) ## n.differential across all 1610 experiments


##  meta-analysis information
Load('/Users/ks/Dropbox/PROJ/cmap/out/table/Full_info.sig.randM.fisher.FDR.Rdata')
clusdat=make.numeric(clusdat[complete.cases(clusdat),])

head(metsum[[1]])

####  -----------------------------------------------------------------------------------------------------------
##     modules to use - ENSG format for now -----------------------------------------------------------------

#Load('/Users/ks/Dropbox/SHARED/liisi_modules/M5_drugs/genes.difexpr.drug.cmap.Rdata')
#Load('/Data/ks/dt/me2f1_m5_targets_up_down.Rdata')
Load('/Users/ks/Dropbox/SHARED/liisi_modules/M5_drugs/data.Rdata')
#[1,] "dif.up.hubs" 
#[2,] "E2F1.hubs.up"
#[3,] "E2F1.M5.down"
#[4,] "E2F1.M5.up"  
#[5,] "M5.down"     
#[6,] "M5.up"  



lmod=list()
	lmod$dif.hubs$up=idconvert.ensg(list(ensg=dif.up.hubs))$ensg
	lmod$E2F1.hubs$up=idconvert.ensg(list(ensg=E2F1.hubs.up))$ensg
	lmod$E2F1.M5$down=idconvert.ensg(list(ensg=E2F1.M5.down))$ensg
	lmod$E2F1.M5$up=idconvert.ensg(list(ensg=E2F1.M5.up))$ensg
	lmod$M5$down=idconvert.ensg(list(ensg=M5.down))$ensg
	lmod$M5$up=idconvert.ensg(list(ensg=M5.up))$ensg




bkg=read.delim('/Users/ks/Dropbox/SHARED/liisi_modules/M5_drugs/bkgrnd__power6.txt',header=F)
bkg=list(bkg=bkg$V1)
	bkg=idconvert.ensg(bkg)$bkg


####  -----------------------------------------------------------------------------------------------------------
##     enrichment stuffs  -----------------------------------------------------------------

dkey=gsub('DE_Drug_(.*)_Cell_.*_Array_.*_Conc_.*_Conc_.*_time_.*_res','\\1',names(degb))
	matst(names(metsum)%in%dkey)
ukey=unique(dkey)



de_thresh=0.01
n_genes=5
mstat=list()

k=1
for(idru in ukey){
	
	if(idru%in%names(metsum)){
		dummy=metsum[[idru]]
			rownames(dummy)=gsub('_at','',rownames(dummy))
		dummy=dummy[rownames(dummy)%in%bkg,]

		dummy$log.fc=dummy$logFC.mean
		dummy$adj.p=dummy$randmMod.fdr
		dummy=dummy[,c('log.fc','adj.p')]

	}
	if(!(idru%in%names(metsum))){
		dummy=degen[dkey==idru]
		if(length(dummy)>1){
			cat(idru,length(dummy))
		}

		if(length(dummy)==1){
			dummy=dummy[[1]]
				rownames(dummy)=gsub('_at','',rownames(dummy))
			dummy=dummy[rownames(dummy)%in%bkg,]

			dummy$log.fc=dummy$logFC
			dummy$adj.p=dummy$adj.P.Val
			dummy=dummy[,c('log.fc','adj.p')]
		}
	}

	holder=list()
#		holder$n.bg=nrow(dummy)
#		holder$n.bg.all=nrow(dummy)
		holder$n.bg.sig=rownames(dummy[dummy$adj.p<de_thresh,])
		holder$n.bg.sig_up=rownames(dummy[dummy$adj.p<de_thresh & dummy$log.fc>0,])
		holder$n.bg.sig_down=rownames(dummy[dummy$adj.p<de_thresh  & dummy$log.fc<0,])

	if(length(holder$n.bg.sig)>n_genes){

#	if(idru%in%names(metsum)){mstat$meta[[idru]]$bg=as.data.frame(humpty)}
#	if(!idru%in%names(metsum)){mstat$single[[idru]]$bg=as.data.frame(humpty)}

	for(imod in names(lmod)){
		if('up' %in% names(lmod[[imod]])){

			dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$up,]
			holder[[imod]]$up_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
			holder[[imod]]$up_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

			mstat[[idru]][[imod]]$up_downreg=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
			mstat[[idru]][[imod]]$up_upreg  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
		}

		if('down' %in% names(lmod[[imod]])){

			dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$down,]
			holder[[imod]]$down_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
			holder[[imod]]$down_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

			mstat[[idru]][[imod]]$down_downreg=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
			mstat[[idru]][[imod]]$down_upreg  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
		}
#		holder[[imod]]=as.data.frame(holder[[imod]])
	}
	rm(humpty)
	rm(dumpty)
	}
	k=lcount(k,length(ukey))
	rm(dummy)
	rm(holder)

}


#idru='trifluridine'
#mstat[['trifluridine']]








#msig=msig[unlist(lapply(msig,length))>n_genes]
#esig=esig[unlist(lapply(esig,length))>n_genes]




de_thresh=0.01
n_genes_thresh=5

ldrug=list()
lstat=list()
k=1
for(idru in names(mega)){
	for(imod in names(lmod)){
		dummy=mega[[idru]]
		dummy=dummy[rownames(dummy)%in%lmod[[imod]],]
		dummy=dummy[dummy[,grepl('adj.P.Val|randmMod.fdr',colnames(dummy)),drop=F]<de_thresh,]

		if(nrow(dummy)>n_genes_thresh){
			ldrug[[idru]][[imod]]=dummy


			lstat[[idru]][[imod]]$n.signif=nrow(dummy)
			lstat[[idru]]$up[[imod]]=nrow(dummy[dummy$logFC>0,])
			lstat[[idru]]$down[[imod]]=nrow(dummy[dummy$logFC<0,])

		}
	}
	k=lcount(k,length(mega))
}





names(signif)=paste0('Drug_',names(signif))
dover=list()
k=1
for(idru in names(signif)){
	holder=list()
	for(imod in names(lmod)){
		holder[[imod]]=unique(intersect(signif[[idru]],lmod[[imod]]))	##  getting some weird behaviour when saving straight to dover here
	}
	dover[[as.character(idru)]][[as.character(imod)]]=holder
	k=lcount(k,length(signif))
}


	holder=list()

	degen[[idru]]=degen[[idru]][grepl('ENSG',rownames(degen[[idru]])),]
	rownames(degen[[idru]])=gsub('_at','\\1',rownames(degen[[idru]]))

	dummy=degen[[idru]][rownames(degen[[idru]])%in%bkg$ids,]
	holder$n.bg=nrow(dummy)

	dummy=dummy[dummy$adj.P.Val<de_thresh,]


	if(nrow(dummy)>n_genes){

		holder$n.signif=nrow(dummy)
		holder$up.reg.all=nrow(dummy[dummy$logFC>0,])
		holder$down.reg.all=nrow(dummy[dummy$logFC<0,])		
		holder$dif.up.hubs=nrow(dummy[dummy$logFC<0 & rownames(dummy)%in%dhubu$ids,])
		holder$E2F1.hubs.up=nrow(dummy[dummy$logFC<0 & rownames(dummy)%in%ehubu$ids,])
		holder$E2F1.M5.up=nrow(dummy[dummy$logFC<0 & rownames(dummy)%in%emoup$ids,])
		holder$E2F1.M5.down=nrow(dummy[dummy$logFC>0 & rownames(dummy)%in%emodn$ids,])
		holder$M5.down=nrow(dummy[dummy$logFC>0 & rownames(dummy)%in%modup$ids,])
		holder$M5.up=nrow(dummy[dummy$logFC<0 & rownames(dummy)%in%moddn$ids,])

		change[[idru]]$stats=as.data.frame(holder)
		change[[idru]]$dif.up.hubs.genes=dummy[dummy$logFC<0 & rownames(dummy)%in%dhubu$ids,]
		change[[idru]]$E2F1.hubs.up.genes=dummy[dummy$logFC<0 & rownames(dummy)%in%ehubu$ids,]
		change[[idru]]$E2F1.M5.up.genes=dummy[dummy$logFC<0 & rownames(dummy)%in%emoup$ids,]
		change[[idru]]$E2F1.M5.down.genes=dummy[dummy$logFC>0 & rownames(dummy)%in%emodn$ids,]
		change[[idru]]$M5.down.genes=dummy[dummy$logFC>0 & rownames(dummy)%in%modup$ids,]
		change[[idru]]$M5.up.genes=dummy[dummy$logFC<0 & rownames(dummy)%in%moddn$ids,]

	}

	k=lcount(k,length(degen))
}


singleDrug_etoposide








#library(pvclust)

pdf('/Users/ks/Dropbox/PROJ/cmap/out/img/drug_clust_dist_ward.pdf',height=10,width=40)
	clust(clusdat,clust_method='ward.D2',k='dynamic',cor_method='dist',do_plots=T,plot_cex=0.8,dat_descr='',help=F)

##  sadly, pvclust is unlikely to ever finish.. try on 'blu'
#    fit <- pvclust(clusdat, method.hclust="ward.D2", method.dist="euclidean")
#    plot(fit) # dendogram with p values
##  add rectangles around groups highly supported by the data
#    pvrect(fit, alpha=.95) 

dev.off()



pdf('/Users/ks/Dropbox/PROJ/cmap/out/img/drug_clust_dist_average.pdf',height=10,width=40)
	clust(clusdat,clust_method='average',k='dynamic',cor_method='dist',do_plots=T,plot_cex=0.8,dat_descr='',help=F)

dev.off()




#cplot(-log10(dummy$randmMod.fdr), abs(dummy$scoreMixed))  ##  good concordance bewteen the two


##  counts of 'significantly' dif-expressed genes
##   5%
##   1%
	ngen=list()
#idru='wortmannin'
for(idru in ukey){
	humpty=degb[dkey==idru]

	holder1=list()
	holder5=list()
	for(idat in names(humpty)){
		holder1[[idat]]=sum(humpty[[idat]]$adj.P.Val<0.01,na.rm=T)
		holder5[[idat]]=sum(humpty[[idat]]$adj.P.Val<0.05,na.rm=T)
	}

	dummy=metsum[[idru]]
	ngen[[idru]]$mix.fdr1=sum(dummy$randmMod.fdr<0.01,na.rm=T)
	ngen[[idru]]$mix.fdr5=sum(dummy$randmMod.fdr<0.05,na.rm=T)
	ngen[[idru]]$fis.fdr1=sum(dummy$Fisher.fdr<0.01,na.rm=T)
	ngen[[idru]]$fis.fdr5=sum(dummy$Fisher.fdr<0.05,na.rm=T)

	ngen[[idru]]$mix.p1=sum(dummy$randomMod.P<0.01,na.rm=T)
	ngen[[idru]]$mix.p5=sum(dummy$randomMod.P<0.05,na.rm=T)
	ngen[[idru]]$fis.p1=sum(dummy$Fisher.P<0.01,na.rm=T)
	ngen[[idru]]$fis.p5=sum(dummy$Fisher.P<0.05,na.rm=T)

	ngen[[idru]]$difgn1=paste(sort(unlist(holder1)),collapse=';')
	ngen[[idru]]$difgn5=paste(sort(unlist(holder5)),collapse=';')


}

(ngen[[1]])

sumstats=(matrix(unlist(ngen),ncol=10,byrow=T))
	colnames(sumstats)=names(ngen[[1]])
	rownames(sumstats)=names(ngen)

#write.delim(sumstats,'/Users/ks/Dropbox/PROJ/cmap/out/table/sumstats_metaP_dif.genes.affy.txt')

#  astemizole
#  gossypol
#  monensin
#  5155877
#  F0447−0125
#  5194442
#  mefloquine
#  vorinostat



msig1=list()
msig5=list()
for(idru in names(metsum)){

	dummy=metsum[[idru]]
	msig1[[idru]]=gsub('_at','',rownames(dummy)[(dummy$randmMod.fdr<0.01)])
	msig5[[idru]]=gsub('_at','',rownames(dummy)[(dummy$randmMod.fdr<0.05)])
#	nmsig1gen[[idru]]$fis.fdr1=rownames(dummy)[sum(dummy$Fisher.fdr<0.01)]
#	msig5[[idru]]$fis.fdr5=rownames(dummy)[sum(dummy$Fisher.fdr<0.05)]

}

skey=ukey[!(ukey%in%names(metsum))]
	length(ukey)-length(skey)



esig1=list()
esig5=list()
k=1
for(idru in skey){
	humpty=degen[dkey==idru]

	if(length(humpty)!=1){
		print(idru)
	}
	if(length(humpty)==1){
		k=lcount(k,length(skey))
		humpty=humpty[[1]]
		esig1[[idru]]=gsub('_at','',rownames(humpty)[(humpty$adj.P.Val<0.01)])
		esig5[[idru]]=gsub('_at','',rownames(humpty)[(humpty$adj.P.Val<0.05)])
	}
}








query=c('pioglitazone','valproic acid','trifluridine','resveratrol','monobenzone','methotrexate','thioridazine')
	matst(query%in%names(msig1))

#
dfunc=(msig1[query[query%in%names(msig1)]])
#dfunc=(msig5[query[query%in%names(msig5)]])

 str(dfunc)
dfunc=dfunc[c('valproic acid','resveratrol','methotrexate','thioridazine')]
bfunc=gsub('_at','',rownames(degen[[1]]))




str(dfunc)


clicky.run(
#	module_genes_list=plm[!grepl('_bkgrnd',names(plm))]
	module_genes_list=dfunc
	,
#	module_bkrnd=plm[grepl('_bkgrnd',names(plm))][[1]]
	module_bkrnd=bfunc
	,
	clicky_dir='/Users/ks/Dropbox/PROJ/cmap/out/funct'
#	clicky_dir='/Users/ks/Dropbox/SHARED/tools/clicky'
	,dat_descr = '_fdr1_'
	,
#	id_type = "hsapiens__gene_symbol"
	id_type = "hsapiens__ensembl_gene_stable_id"
)


###====================================================================================================================================
###  functional enrichment plots ============================================================================================================
###====================================================================================================================================

library(clickyOSX)



enrich=gestalt_read(mod_names=paste0(names(dfunc),'__fdr5'),in_path='/Users/ks/Dropbox/PROJ/cmap/out/funct/tmp_in_files/')

pdf(paste0('/Users/ks/Dropbox/PROJ/cmap/out/funct/img/funct.',paste(names(dfunc),collapse='_'),'.mixedM.fdr5.pdf'),height=22,width=11)
	gestalt_plot(enrich,p_threshold=0.1)
dev.off()







####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     TODO     ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
## clean up databases
##  combine meta-analysed and single experiment - probably 2 separate lists
##  set up enrichment method script
##  functional enrichment for genes that are affected by drug within and outside module










































