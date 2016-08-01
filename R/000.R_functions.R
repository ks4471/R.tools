"
╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗
╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝
https://www.dropbox.com/s/4nhe1ukd7ee9b3h/000.R_functions.R?dl=0 

╔═╦╗╔╦╗╔═╦═╦╦╦╦╗╔═╗╔═╦╗╔╦╗╔═╦═╦╦╦╦╗╔═╗╔═╦╦╦╦╗╔═╗╔═╦╗╔═╦╗╔╦╗╔═╦═╦╦╦╦╗╔═╗╔═╗╔═╦╗╔═╦╗╔╔═╦╗╔═╦╗╔╗
╠╗║╚╝║║╠╗║╚╣║║║║║╚╣                   ╠╣║║║║║═╣║║╠╗║║╚╣╚╣╔╣╔╣╔╣╔╣║╚╣═╣║╚╣║║║╚╣╔╣╔╣║╚╣═╣║╗║╚╚╣
╚═╩══╩═╩═╩═╩╝╚╩═╩═╝╚═╩══╩═╩═╩═╩╝╚╩═╩═╩╝╚╩═╩═╝╚═╩══╩═╩═╩═╩╝╚╩═╩═╩╝╚╩═╩═╝═╩╝╚╩═╩═╩╝╩═╩═╩═╩╝╚╩═╩
"



#•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•#
##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•#- MOTHBALLS -#•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##
#•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•#
#•#
#•#		#git clone https://github.com/jalvesaq/colorout.git
#•#		#sudo R CMD INSTALL colorout
#•#
#•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•#
#•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•##•#

#install.packages('devtools')

#library(devtools)

#t1=Sys.time()
#library(colorout)
#devtools::install_github("ks471/R.helper")
#devtools::install_github("ks471/clickyOSX")
#devtools::install_github("ks471/clickyLinux")

#Sys.time()-t1

#library(R.helper)
#library(clickyOSX)
#library(clickyLinux)

#setwd('~/Dropbox/SHARED/tools/R_functions/R.helper/data/')
##Error: Could not find package root.		##  error possibly due to the fact that it checks for a R package structure around the folder where it is saving
#devtools::use_data(x)	# saves to working directory by default





####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
# nohup for R scrpts
# nohup Rscript foo.R > ~/R/logs/foo.out 2>&1 &
#http://www.r-bloggers.com/long-running-r-commands-unix-screen-nohup-and-r/
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■

## extracting files using cmd / mac os
# tar xvf file.tar
# gunzip file.gz
# gzip -d file.gz


####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■

####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#fnames=list.files("/Users/Shkura/Dropbox/Cognition/GWAS/magma/hrh.magma/out/netw.enrich", pattern = ".out", full.names = TRUE)
#basename
#dirname
#basename(fnames)

#load("~/Downloads/expr.maps.rda")
#attach("~/Downloads/expr.maps.rda")  # works in a way similar to load but checks workspaces for variables with the same names as the object and 'hides' them

#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╔╦╦╗
#source('~/Dropbox//000.R.functions.R')#╣║╚╣═╣║╚╣║║║╚╣╔╣╔╣║╚╣═╣╔╗║╚╚╣║╚╣
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩╩═╝

## useful concept : add a readme.object.in.saved.R.object when saving file create a "method" / "description" of the file for easier handover / re-analysis
##  readme.lires=
#            "\tlires - differentially expressed genes in the hippocampus compared to 3 other brain regions FC, OC, TC Hardy UKBEC dataset n=102\n
#            \tdesm = design matrix file for limma\n
#            limma used for this code.file : 003.DifferentiallyExpressedGenes.limma.R     |\n|     adj.P - calculated based using FDR, CX excluded because it is very differentially expressed compared to the other 3 - all genes were differential with CX"
##  save(lires,desm,readme,file="~/Dropbox/CapricaPrime/RU/dtb/secondary/HC.differentially.expressed.v.FC.OC.TC.R")

###  saves 'cat' output to file, file separator can be controlled using the "\t" in relevant places
#sink("~/Dropbox/bin/gsea/modules_bonn.gtf")
#  for(imod in 1:length(bgen)){
#   cat(as.vector(paste(c(names(bgen)[imod],bgen[[names(bgen)[imod]]]),collapse="\t")),"\n",sep="")
#  }
#sink()


## easier names for vars used by phyper and fisher.test: 
# phyper(success_in_sample, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail=TRUE)

#fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='less');
# Numerical parameters in order:
# (success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part).

#########################################################################################################
# ----------------------------------------------------------------------------------------------------- #
#
#library(colorout)

# alternative ?
#https://github.com/gaborcsardi/crayon
#devtools::install_github("gaborcsardi/crayon")
#library(crayon)


# ----------------------------------------------------------------------------------------------------- #
#########################################################################################################
# Alternate way to calculate PC1    ----------------------------------------------------
# removing the first Principal component from data and reconstructing full matrix
#SVD=svd(t(scale(t(Ccombi[2:(ncol(Ccombi)-3)]))))
#str(SVD)
#Expr=SVD$u%*%diag(SVD$d)%*%t(SVD$v)
#PC1=SVD$u[,1]%*%(SVD$d[1])%*%t(SVD$v[1,1])

#SVD=svd(t(scale(t(D1))))
#str(SVD)
#PC1=SVD$u[,1]%*%(SVD$d[1])%*%t(SVD$v[1,1])
#for(idat in 1:length(fnames)){
#cat("     ",dat_descr,"===========",which(fnames==fnames[idat]),"of",length(fnames),"\n")


#for(ireg in 1:length(names(list_expr))){
#   print(paste("--------------------",names(list_expr)[ireg],"---------------------",ireg,"of",length(names(list_expr))))

#   cat(round(j/ncol(Matrix),digits=2),"\r");flush.console()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
##  WARNING : MANUALLY select the file below based on covariates corrected (contains both rma and plier)
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ es in the folder OTHER than the ones matching the format below
##   i.e. "02.WGCNA.MODULES.hipp.plier-gcbgPEER.bb.sex.age.pmi.cod.merge.height=0.25.genes14646.samples102powr5.bicor.R"
###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#specify the name of the object in the R binary file you just loaded ||  can use this to load R binary files instead of 'load':

#file_path = "/User/dy/my_R_data/a_data_set.RData"
#attach(file_path, pos=2, name=choose_a_name, warn.conflict=T)
#'warn.conflicts=T' is the default option

#'pos=2' is also the default; "2" refers to the position in your search path. For instance, position 1 is ".GlobalEnv." To get the entire array of search paths, use search(). So you would access the search path for the new object by search()[2]

#use 'detach' to remove the object

"
╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗
╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝
"


####■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
###------------------------------------------------------------------------------------------------------------
###------------------------------            PEER data correction            ----------------------------------
###------------------------------------------------------------------------------------------------------------
##


#  potentially very useful plot  - eg factor plot
#par(mar=c(0,0,0,0))

# Set up the plotting area
#plot(NA, xlim=c(0,1), ylim=c(6.5, -0.5),
#    xaxt="n", yaxt="n",
#    xlab=NA, ylab=NA )

# Draw the lines
#for (i in 0:6) {
#    points(c(0.25,1), c(i,i), lty=i, lwd=2, type="l")
#}
# Add labels
#text(0, 0, "0. 'blank'"   ,  adj=c(0,.5))
#text(0, 1, "1. 'solid'"   ,  adj=c(0,.5))
#text(0, 2, "2. 'dashed'"  ,  adj=c(0,.5))
#text(0, 3, "3. 'dotted'"  ,  adj=c(0,.5))
#text(0, 4, "4. 'dotdash'" ,  adj=c(0,.5))
#text(0, 5, "5. 'longdash'",  adj=c(0,.5))
#text(0, 6, "6. 'twodash'" ,  adj=c(0,.5))


colmix=c(
"#0072B2"
,"#E69F00"
,"#009E73"
,"#56B4E9"
,"#D55E00"
,"#66A61E"
,"#7570B3"
,"#a50f15"
,"#A6761D"
,"#117733"
,"#332288"
,"#b15928"
,"#882255"
,"#999933"
,"#AA4499"
,"#1f78b4"
,"#F0E442"
,"#e31a1c"
,"#6a3d9a"
,"#b2df8a"
,"#08519c"
,"#ff7f00"
,"#fdbf6f"
,"#33a02c"
,"#b15928"
,"#f16913"
,"#238b45"
,"#807dba"
,"#d94801"
,"#41ab5d"
,"#fd8d3c"
,"#4292c6"
)

pastelcolmix=c(
"#cab2d6"
,"#ffff99"
,"#8dd3c7"
,"#ffffb3"
,"#bebada"
,"#fb8072"
,"#80b1d3"
,"#fdb462"
,"#b3de69"
,"#fccde5"
,"#d9d9d9"
,"#bc80bd"
,"#ccebc5"
,"#ffed6f"
,"#a6cee3"
,"#fb9a99"
)


colmixrb=c(colorRampPalette(c(
"#3f007d"
,"#313695"
,"#053061"
,"#08306b"
,"#045a8d"
,"#08519c"
,"#2171b5"
,"#6baed6"
))(8),rev(colorRampPalette(c(
"#67001f"
,"#67000d"
,"#a50f15"
,"#bd0026"
,"#cb181d"
,"#e31a1c"
,"#ef3b2c"
,"#f4a582"
))(8)))

colmixb=c(colorRampPalette(c(
"#3f007d"
,"#313695"
,"#053061"
,"#08306b"
,"#045a8d"
,"#08519c"
,"#2171b5"
,"#6baed6"
))(8))

colmixr=rev(colorRampPalette(c(
"#67000d"
,"#a50f15"
,"#bd0026"
,"#cb181d"
,"#e31a1c"
,"#ef3b2c"
,"#f4a582"
))(7))


colred=c('#67001f','#a50026','#f46d43','#fb9a99')
colblu=c('#2d004b','#053061','#2166ac','#4393c3')
colblus=c('#08306b','#08519c','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7','#f7fbff')
colreds=c('#67000d','#a50f15','#cb181d','#ef3b2c','#fb6a4a','#fc9272','#fcbba1','#fee0d2','#fff5f0')
colrbd  =c('#320000',"#800026","#bd0026","#e31a1c","#fc4e2a","#fd8d3c","#feb24c","#fed976"  # reds

         ,"#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58",'#49006a')  # blues # dark purple , -- need a corresponding smth for the darkest red
colrb  =c("#800026","#bd0026","#e31a1c","#fc4e2a","#fd8d3c","#feb24c","#fed976"  # reds
          ,'white'
         ,"#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58")  # blues # dark purple , -- need a corresponding smth for the darkest red


colbw=c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525","#000000")

#--#~~blu=c(
#blue- violet
#--#~~"#fff7fb"
#--#~~#--#~~,"#ece7f2"
#--#~~,"#d0d1e6"
#--#~~,"#a6bddb"
#--#~~,"#74a9cf"
#--#~~,"#3690c0"
#--#~~,"#0570b0"
#--#~~,"#045a8d"
#--#~~,"#023858"
# blue - violet
#--#~~,"#f7fcfd"
#--#~~,"#e0ecf4"
#--#~~,"#bfd3e6"
#--#~~,"#9ebcda"
#--#~~,"#8c96c6"
#--#~~,"#8c6bb1"
#--#~~,"#88419d"
#--#~~,"#810f7c"
#--#~~,"#4d004b"

# violet
#--#~~,"#fcfbfd"
#--#~~,"#efedf5"
#--#~~,"#dadaeb"
#--#~~,"#bcbddc"
#--#~~,"#9e9ac8"
#--#~~,"#807dba"
#--#~~,"#6a51a3"
#--#~~,"#54278f"
#--#~~,"#3f007d"
#--#~~,"#f7fbff"
# blue
#--#~~,"#deebf7"
#--#~~,"#c6dbef"
#--#~~,"#9ecae1"
#--#~~,"#6baed6"
#--#~~,"#4292c6"
#--#~~,"#2171b5"
#--#~~,"#08519c"
#--#~~,"#08306b"
# mix
#--#~~,"#d1e5f0"
#--#~~,"#e0f3f8"
#--#~~,"#92c5de"
#--#~~,"#abd9e9"

#--#~~,"#74add1"
#--#~~,"#4393c3"
#--#~~,"#4575b4"
#--#~~,"#2166ac"

#--#~~,"#053061"
#--#~~,"#313695"
#--#~~)






#--#~~red=c(
# mix
#--#~~"#67001f"
#--#~~,"#a50026"
#--#~~,"#b2182b"
#--#~~,"#d73027"
#--#~~,"#d6604d"
#--#~~,"#f46d43"

#--#~~,"#fdae61"
#--#~~,"#fee090"
#--#~~,"#f4a582"


#--#~~,"#ffffbf"
#--#~~,"#ffffcc"
#--#~~,"#ffeda0"
#--#~~,"#fed976"
#--#~~,"#feb24c"
#--#~~,"#fd8d3c"
#--#~~,"#fc4e2a"
#--#~~,"#e31a1c"
#--#~~,"#bd0026"
#--#~~,"#800026"
#--#~~,"#fff5eb"
#--#~~,"#fee6ce"
#--#~~,"#fdd0a2"
#--#~~,"#fdae6b"
#--#~~,"#fd8d3c"
#--#~~,"#f16913"
#--#~~,"#d94801"
#--#~~,"#a63603"
#--#~~,"#7f2704"
#--#~~,"#fff5f0"
#--#~~,"#fee0d2"
#--#~~,"#fcbba1"
#--#~~,"#fc9272"
#--#~~,"#fb6a4a"
#--#~~,"#ef3b2c"
#--#~~,"#cb181d"
#--#~~,"#a50f15"
#--#~~,"#67000d")





#--#~~tencol=c(
#--#~~"#0a0972"
##--#~~,"#313695"
#--#~~,"#08306b"
#--#~~,"#023858"



##--#~~,"#053061"

##--#~~,"#08519c"
##--#~~,"#045a8d"
#--#~~,"#2166ac"
#--#~~,"#2171b5"
#--#~~,"#4575b4"




#--#~~,"#4393c3"





#--#~~,"#3690c0"

#--#~~,"#4292c6"

#--#~~,"#6baed6"
#--#~~,"#74a9cf"
#--#~~)


#ramp=colorRampPalette(c("#0a0972"#--#~~,"#023858"#--#~~,"#2166ac"#--#~~,"#4292c6"))(10)

#--#~~rampb=colorRampPalette(c(
#--#~~"#3f007d"
#--#~~,"#313695"
#--#~~,"#053061"
#--#~~,"#08306b"
#--#~~,"#045a8d"
#--#~~,"#08519c"
#--#~~,"#2171b5"
#--#~~,"#6baed6"
#--#~~))(10)


#--#~~rampr=colorRampPalette(c(
#--#~~"#67001f"
#--#~~,"#67000d"
#--#~~,"#a50f15"
#--#~~,"#bd0026"
#--#~~,"#cb181d"
#--#~~,"#e31a1c"
#--#~~,"#ef3b2c"
#--#~~,"#f4a582"
#--#~~))(10)
#colmix=blu
#colmix=red
#colmix=tencol
#
#--#~~colmix=c(rampb,rev(rampr))

#--#~~  par(mar=c(0,0,0,0))
#--#~~  plot(NA, xlim=c(0,1), ylim=c(length(colmix), -0.5),
#--#~~      xaxt="n", yaxt="n",
#--#~~      xlab=NA, ylab=NA )

#--#~~  for (i in 1:length(colmix)){
#--#~~      points(c(0.25,1), c(i,i), lwd=15, type="l",col=colmix[i])
#--#~~    text(0, i, paste(i,colmix[i])   ,  adj=c(0,.5))
#--#~~}



#--#~~  colmix=c(
#--#~~  "white"
#--#~~  ,"#1f78b4"
#--#~~  ,"#E69F00"
#--#~~  ,"#882255"
#--#~~  ,"white"
#--#~~  ,"#56B4E9"
#--#~~  ,"#66A61E"
#--#~~  ,"#bd0026"
#--#~~  ,"white"
#--#~~  ,"#0072B2"
#--#~~  ,"#117733"
#--#~~  ,"#A6761D"
#--#~~  ,"white"
#--#~~  ,"#D55E00"
#--#~~  ,"#009E73"
#--#~~  ,"white"
#--#~~  ,"#AA4499"
#--#~~  ,"#6a3d9a"
#--#~~  ,"white"
#--#~~  ,"#7570B3"
#--#~~  ,"white"
#--#~~  ,"#332288"
#--#~~  ,"#999933"
#--#~~  ,"#b2df8a"
#--#~~  ,"#33a02c"
#--#~~  ,"#ff7f00"
#--#~~  ,"#e31a1c"
#--#~~  ,"#fdbf6f"
#--#~~  ,"#b15928"
#--#~~  ,"#fc4e2a"
#--#~~  ,"#e31a1c"
#--#~~  ,"#800026"
#--#~~  ,"#F0E442"
#--#~~  )


############################################################################################################
##############################=============     WGCNA Part 2      =============#############################
##############################-----         Clustering Analysis         -----###############################
############################################################################################################




#print("■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■")

####=================================================================================================
### Calculate consensus modules


## -------   1   -------
#peer="PEER."
#peer=""

## -------   2   -------
#pc1correct_data="pc1correct."
#
#pc1correct_data=""

####=================================================================================================
###-----------------------------------------------------------------------------------------------
# Plot the cut tree   ----------------------------------------------------------------------------




#=================================================================================================
#    3) Cluster again? Copy - paste this until there's no change. 
#------------------------------------------------------------------------------



####----------------------------------------------------------------------------------------------
##   calculate networks

#       May have to tinker with options to optimise tree cutting ---------------------------------
#       maxBlockSize - increased to cut as single tree (avoid merging)
#       deepSplit - reduced to 2 (from 3) 
#       an important parameter to be selected <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<***҉҉҉҉҉҉҉
#


#҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉҉


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



make.numeric<-function(Matrix,col_factor="",fac_legend=F,verbose=F,help=F,char_as_fac=F){ #char_as_fac=T,
if(help==T){
      cat("\tINPUTS:\tMatrix - rows : samples, columns - variables\n")
      cat("\tOPTIONS:\tfactors - column numbers to convert to factor instead of numeric  |  default='none'\n")
      cat("\tNOTE:\tcharacter columns can be converted to numeric 'categories' by as.numeric (default)\n\n")
      }
# for any correlation analysis like linear modeling non-numeric categories eg 'sex' should always be converted to factor, numeric 1 = male, 2=female is considered 2>1 vs 1!=2 if as.factor

# cat("converting",class(Matrix),"dimensions :",dim(Matrix),"\n")

#  print(head(Matrix))
#if(char_as_fac==T){
#  col_factor=names(which(lapply(Matrix,class)=="character"))
#}

legend=list()
  num.mat=matrix(NA,ncol=ncol(Matrix),nrow=nrow(Matrix))
    colnames(num.mat)=colnames(Matrix)
    rownames(num.mat)=rownames(Matrix)

#  if(char_as_fac){

#  }
col_numeric=1:ncol(Matrix)
col_numeric=col_numeric[!(col_numeric%in%col_factor)]

#print(col_numeric)
#print(col_factor)



    for(inum in col_numeric){
      if(verbose==T){
        cat("\t\t - column",inum,"coverted as.numeric\n")
      }
      num.mat[,inum]=as.numeric(Matrix[,inum])

      if(class(Matrix[,inum])=="factor" & fac_legend){
        descr=list(legend=Matrix[,inum])
        descr$values=num.mat[,inum,drop=F]
        descr=unique(as.data.frame(descr))
        legend[[colnames(Matrix)[inum]]]=descr
      }
#      print(head(num.mat))
    }
  if(sum(col_factor=="")!=1){
      if(verbose==T){
      cat("\n\tNOTE:\tmatrix can only hold 1 kind of object eg numeric, character, etc..\n\t\t\tdata.frame will be used instead\n")
      }
    num.mat=as.data.frame(num.mat)
    for(ifac in col_factor){
      if(verbose==T){
        cat("\t\t - column",ifac,"coverted as.factor\n")
      }
      num.mat[,ifac]=as.factor(Matrix[,ifac])
#      print(head(num.mat))
    }
  }

# cat("new",class(num.mat),"dimensions :", dim(num.mat))
  if(!fac_legend){
    return(num.mat)
  }
  if(fac_legend){
    return(list("numeric"=num.mat,"legend"=legend))
  }


      if(verbose==T){
        cat("\t---------------------   make.numeric - finished   ---------------------\n\n")
      }

}





### ----------------------------------------------------------------------
##   new version below is updated version to handle list objects more easily, and 
#Head<-function(Matrix){
#  if(is.vector(Matrix)==T & is.list(Matrix)==F){
#    print(as.matrix(Matrix[1:min(c(length(Matrix),10))]))
#    cat("\n vector length :",length(Matrix),"\n")
#    cat("\n       Numeric :",is.numeric(Matrix),"\n")
#  }
#  if(is.list(Matrix)==T & is.data.frame(Matrix)==F){
#    print(str(Matrix[1:min(c(length(Matrix),10))]))
#    print(as.matrix(Matrix[1:min(c(length(Matrix),10))]))
#    Matrix=as.data.frame(Matrix)
#    print(Matrix[1:min(c(nrow(Matrix),10)),1:min(c(ncol(Matrix),5))])
#    cat("\n list length :",length(Matrix),"\n")
#    cat("\n     Numeric :",is.numeric(Matrix),"\n")
#  }

#  if(is.matrix(Matrix)==T | is.data.frame(Matrix)==T){
#    print(Matrix[1:min(c(nrow(Matrix),10)),1:min(c(ncol(Matrix),5))])
#    cat("\n",class(Matrix),"dimensions :",dim(Matrix),"\n")
#    cat("\n           Numeric :",is.numeric(Matrix),"\n")
#  }
#}



### ----------------------------------------------------------------------
##   Updated version to handle list objects more easily
Head<-function(data_obj,nlist=1,ncol=1:5,nrow=1:10){
  cat("\n\tobject class : ",class(data_obj),"\n\n")
  prev=""

  if(class(data_obj)=="list"){

  	if(length(names(data_obj))<50){
	    cat("\t\tlist contains",length(names(data_obj)),"objects:\n")
	    cat("\t\t\t",as.matrix(names(data_obj)),"\n\n",sep="   ")
	    cat("\t\tlist[[",nlist,"]] contains ",class((data_obj[[1]]))," :",sep="")

	    prev=paste("list[[",nlist,"]] contains :",sep="")
	    data_obj=data_obj[[nlist]]
    }

  	if(length(names(data_obj))>50){
	    cat("\t\tlist contains",length(names(data_obj)),"objects, first 50:\n")
	    cat("\t\t\t",as.matrix(names(data_obj)[1:min(50,length(names(data_obj)))]),"...\n\n",sep="   ")
	    cat("\t\tlist[[",nlist,"]] contains ",class((data_obj[[1]]))," :",sep="")

	    prev=paste("list[[",nlist,"]] contains :",sep="")
	    data_obj=data_obj[[nlist]]
    }

  }

  if(class(data_obj)=="list"){
    cat("",length(names(data_obj)),"objects:\n")
    cat("\t\t\t",as.matrix(names(data_obj)),"\n\n",sep="   ")
  }

  if(class(data_obj)=="data.frame" | class(data_obj)=="matrix"){
    cat("\n\n")
    print(data_obj[min(1,min(nrow)):min(max(nrow),nrow(data_obj)),min(1,min(ncol)):min(max(ncol),ncol(data_obj)),drop=F])
    cat("\n\t",prev,class(data_obj),"dimensions : ",dim(data_obj),"\n")
    cat("\t\tis.numeric :",is.numeric(data_obj))
 
    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')
  }

  if(class(data_obj)=="vector"){
    cat("\n\n")
    print(data_obj[nrow])
    cat("\n\t",prev,class(data_obj),"length : ",length(data_obj),"\n")
    cat("\t\tis.numeric :",is.numeric(data_obj))

    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')

  }

  if(class(data_obj)!="list" & class(data_obj)!="data.frame" & class(data_obj)!="matrix" & class(data_obj)!="vector"){
    cat("\n\n")
    str(data_obj)
    cat("\n\tis.numeric :",is.numeric(data_obj))
    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')

  }
# cat("\t",class(data_obj),"dimensions : ",dim(data_obj),"\n")
}



#
sd.check<-function(dat_mat,check_rows=F,check_cols=T,verbose=T,help=F){
if(help==T){
  cat("\nUSE: check that values in rows/columns of matrix/dataframe vary: sd>0 or factors are informative: have >1 level / not ids n.levels=n.rows\n")
}
#start_time=Sys.time()

if(is.matrix(dat_mat)){dat_class=apply(dat_mat,2,class)}
if(is.data.frame(dat_mat)){dat_class=unlist(lapply(dat_mat,class))}
#if(verbose==T){
  cat("\n+++++++++++++++++ data qc check +++++++++++++++++\n")
  cat("dat_mat contains",ncol(dat_mat),"variables, of which :\n")
  print(table(dat_class))
#}
    dat_num=dat_mat[,dat_class==("numeric"),drop=F]
    dat_fac=dat_mat[,dat_class=="factor",drop=F]
    dat_otr=dat_mat[,!(dat_class%in%c("factor","numeric")),drop=F]

rowind=1:nrow(dat_mat)
colind=1:ncol(dat_mat)

if(check_cols==T){
  fac_col=apply(dat_fac,2,function(x) sum(table(x)!=0))
  num_col=apply(dat_num,2,sd)

  numvarcol=names(num_col)[num_col==0]
    if(length(numvarcol)){cat("\tcol - values do not vary (sd=0) :\t",paste(numvarcol,collapse=", "),"\n")}
  facvarcol=names(fac_col)[fac_col==1]
    if(length(facvarcol)){cat("\tcol - contains single value :\t\t",paste(facvarcol,collapse=", "),"\n")}
  facvaridc=names(fac_col)[fac_col==nrow(dat_fac)]
    if(length(facvaridc)){cat("\tcol - as many levels as rows :\t\t",paste(facvaridc,collapse=", "),"\n")}
    colind=!(colnames(dat_mat)%in% c(numvarcol,facvarcol,facvaridc))

#  cat("\ttotal n cols removed :\t",length(c(numvarrow,facvarrow,facvarids)),"\n")    # for this to work need to a track counter // or some count of n elements for each var
}

if(check_rows==T){
  fac_row=apply(dat_fac,1,function(x) sum(table(x)!=0))
  num_row=apply(dat_num,1,sd)

  numvarrow=names(num_row)[num_row==0]
    if(length(numvarrow)){cat("\trow - values do not vary (sd=0) :\t",paste(numvarrow,collapse=", "),"\n")}
  facvarrow=names(fac_row)[fac_row==1]
    if(length(facvarrow)){cat("\trow - contains single value :\t\t",paste(facvarrow,collapse=", "),"\n")}
  facvaridr=names(fac_row)[fac_row==nrow(dat_fac)]
    if(length(facvaridr)){cat("\trow - as many levels as rows :\t\t",paste(facvaridr,collapse=", "),"\n")}
  rowind=!(rownames(dat_mat)%in% c(numvarrow,facvarrow,facvaridr))

#  cat("\ttotal n rows removed :\t",length(c(numvarrow,facvarrow,facvarids)),"\n")    # for this to work need to a track counter // or some count of n elements for each var
}

if(ncol(dat_otr)>0){
    cat("\tcol - not factor nor numeric :\t",paste(names(dat_otr),collapse=", "),"\n")
}
  cat("\n")


  return(invisible(dat_mat[rowind,colind]))
#  print(Sys.time()-start_time)
}




psig<-function(dat_mat,p_col,sort=F){
  dat_mat$sig=""
  dat_mat[dat_mat[,p_col]<0.1,"sig"]="."
  dat_mat[dat_mat[,p_col]<0.05,"sig"]="+"
  dat_mat[dat_mat[,p_col]<0.01,"sig"]="*"
  dat_mat[dat_mat[,p_col]<0.001,"sig"]="**"
  dat_mat[dat_mat[,p_col]<0.0001,"sig"]="***"
  dat_mat[dat_mat[,p_col]<0.00001,"sig"]="****"

  if(sort){
    dat_mat=dat_mat[order(dat_mat[,p_col]),]
  }
  return(dat_mat)
}



# extract p-value form linear model
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f = summary(modelobject)$fstatistic
    p = pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) = NULL
    return(p)
}

mcols=  c(    # color pallete optimised for omitting 1, if plotting using all the alternatives for 4:5, 7:8 are more useful
  "#0072B2" 
  ,"#E69F00"
  ,"#009E73"
  ,"#56B4E9"  #,"#7570B3"
  ,"#D55E00"  #,"#882255"

  ,"#66A61E"
  ,"#7570B3"  #,"#D55E00" 
  ,"#882255"  #,"#D55E00"   ####,"#CC79A7"
  ,"#F0E442"
  ,"#A6761D"
  
  ,"#AA4499"# 
  ,"#117733"
  ,"#332288"
  ,"#999933")

lplot<-function(x,n,legend=T,xlab="",ylab="",main="",adjust=0,thresh=0.01){
  coords=max((adjust+1),1):min(11,(ncol(x)+(adjust)))
print(coords)
  par(mar=c(8.1,8.1,4.1,2.1),mgp=c(7,1,0))
  pchl=c(21:25,15:20)   # col - specifies color of 15:20 and border color in 21:25, bg = color of 21:25, nothing in 15:20
  coll=c(   # color pallete optimised for omitting 1, if plotting using all the alternatives for 4:5, 7:8 are more useful
  "#0072B2" 
  ,"#E69F00"
  ,"#009E73"
  ,"#56B4E9"  #,"#7570B3"
  ,"#D55E00"  #,"#882255"

  ,"#66A61E"
  ,"#7570B3"  #,"#D55E00" 
  ,"#882255"  #,"#D55E00"   ####,"#CC79A7"
  ,"#F0E442"
  ,"#A6761D"
  
  ,"#AA4499"# 
  ,"#117733"
  ,"#332288"
  ,"#999933")


# for testing colors
#x=1:20
#x=cbind(x,x)
#plot(x,main="color scale")
#for(i in 1:length(coll)){
# abline(h=i,col=coll[i],lwd=10)
#}
  if(ncol(x)>11){print("11 columns is the current maximum coded for plotting, columns 12 onwards ignored")}

plot(seq(floor(range(x[,(which(x==max(x,na.rm=T),arr=T))[1,2]])[1]),ceiling(range(x[,(which(x==max(x,na.rm=T),arr=T))[1,2]])[2]),length=nrow(x)),type="n",yaxt="n",xaxt="n",ylab=ylab,xlab=xlab,main=main,frame=F)    # (which(x==max(x,na.rm=T),arr=T))[2] determiines which column has the highest number, the lowest ==0 (assumed)
  abline(h=-log10(thresh),col="red",lty='dashed')
  axis(1, at=1:nrow(x), labels=rownames(x),las=2) # las - changes orientation of labels
  axis(2, at=0:ceiling(max(x,na.rm=T)),las=2) # las - changes orientation of labels
# axis(2, at=0:floor(max(x,na.rm=T)), labels=c("0","1","2","3","4","5","6","7")[0:(floor(max(x,na.rm=T))+1)],las=2) # las - changes orientation of labels
# axis(2, at=0:floor(max(x,na.rm=T)), labels=c("0","0.1","0.01","0.001","0.0001","0.00001","0.000001")[0:(floor(max(x,na.rm=T))+1)],las=2)  # las - changes orientation of labels
# axis(2, at=0:floor(max(x,na.rm=T)), labels=c("100 %","10 %","1 %","0.1 %","0.01 %","0.001 %","0.0001 %")[0:(floor(max(x,na.rm=T))+1)],las=2)  # las - changes orientation of labels
  

  for(i in coords){
    points(x[,i-(adjust)],pch=pchl[i],bg=coll[i],col=coll[i])
  }
  if(legend==T){
    legend(x="topright",pch=pchl[coords],box.lwd = 0,box.col = "white",col=coll[coords],pt.bg=coll[coords],legend=colnames(x)[1:min(11,ncol(x))],cex=1)
  }}




Heatmap<-function(cor.measures,min=-1,max=1,rowclust=F,colclust=F,ncols=101,dendrogram="none",main="",mode="cor",sig=T,cexrow=0.7,cexcol=0.7,margin=c(5,5),lheight=c(0.1,0.9),lwidth=c(0.7,0.3),lmatrix = rbind(c(4,3),c(1,2))){
  library(gplots)
# print("Heatmap( 'matrix to use' , 'min value for legend cols' , 'max ..' , 'clust by rows', 'clust by cols' , 'ncols to use for legend', 'dendrogram=c('none','row','column','both')' )")
# min=-1
# max=1

# default settings for testing plots
#min=-1;max=1;rowclust=F;colclust=F;ncols=101;dendrogram="none";main="";mode="cor"


if(mode=="cor"){
  print("plotting correlation based matrix")
#  heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=0.7,   cexRow=0.7,   lheight = lheight,symkey=T,main=main,lmat=lmatrix)#,lwid=c(0.5,0.5))
   heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,   symkey=T,main=main,lwid=lwidth,lmat=lmatrix)
#   heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#deebf7","#f7fbff","white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,   symkey=T,main=main,lwid=lwidth,lmat=lmatrix)

#"white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"
}
if(mode=="pval"){
  if(sig==T){
    print("plotting p-value based matrix with cellnotes")
    corm=round(cor.measures,digits=3)
    cor.match=matrix(as.numeric(NA),nrow=nrow(cor.measures),ncol=ncol(cor.measures))
      rownames(cor.match)=rownames(cor.measures)
      colnames(cor.match)=colnames(cor.measures)
    # Head(cor.match)
    for(i in 1:nrow(corm)){
      for(j in 1:ncol(corm)){
    #     cat(i,j,"\n")
        if(corm[i,j]>0.05){
          cor.match[i,j]=""
    #     print(">0.05")
        }
        if(corm[i,j]<=0.01){
          cor.match[i,j]=gsub(" ","",paste(rep("*",min(4,floor(-log10(corm[i,j]))-1)),collapse=" "))    # backup square == .
    #     print("<0.01")
        }
        if(corm[i,j]<0.05 & corm[i,j]>0.01){
          cor.match[i,j]="+"
    #     print("<0.05")
        }
      }
    }
#    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,
    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,   
#    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,   


  }
  if(sig==F){
    print("plotting p-value based matrix, no cellnotes")
#    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,    
    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)
#    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)

  }
}

#min(na.omit(as.vector(cor.measures))),max(na.omit(as.vector(cor.measures)))
}




Heat<-function(cor.measures,rowclust=T,colclust=T,ncols=101,cexrow=0.7,cexcol=0.7,margin=c(12,12),...){

  if(rowclust & colclust){dendrogram="both"}
  if(rowclust & !colclust){dendrogram="row"}
  if(!rowclust & colclust){dendrogram="column"}
# can feasibly include a parameter for dendrogram as well, for more plottting flexibility

    library(gplots)
    min=floor(min(cor.measures))
    max=ceiling(max(cor.measures))
    cat('\tmin=',min,', max=',max,'\n') 
    if(min<0 & max>=0){heat_colors =c(colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(ncols)); symmkey=T;ncols=ncols-1}
    if(min>=0 & max>=0){heat_colors=c('grey60',colorRampPalette(c("white","#F0E442","darkred"))(ncols));symmkey=F}
    if(min<=0 & max<=0){heat_colors=c('grey60',colorRampPalette(c("white","#56B4E9","#0072B2"))(ncols));symmkey=F}
  cor_heat=heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+2)),col = heat_colors,trace="none",dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,symkey=symmkey,hclustfun=function(x) hclust(x, method="ward.D2"),...)#,hclustfun=function(x) hclust(x, method="ward.D2"))
  return(invisible(cor_heat))
}



write.delim<-function(mat,file,row.names=T,col.names=T,missing.value.char="NA",sep="\t",...){
  if(col.names==T & row.names==T){
    col.names=NA
  }
  write.table(mat,file,row.names=row.names,col.names=col.names,sep=sep,quote=F,na=missing.value.char,...)
}


DAVID_enrich_list<-function(backg,clusters_list,idType= "ENSEMBL_GENE_ID", query_type="Gene"){
  library(RDAVIDWebService)
  david<-DAVIDWebService$new(email="aida.moreno-moral11@imperial.ac.uk")
 cat("\t",length(backg),"background genes -------------------\n")
  background <- addList(david, backg, idType= "ENSEMBL_GENE_ID", listName="backg", listType="Background")
  setCurrentBackgroundPosition(david, 1)

  DAVID_clusters=list()
  for (i in 1:length(clusters_list)){

    cat("\t",length(clusters_list[[i]]),"cluster genes\t\t",i," of ",length(clusters_list),"\n")
    result<-addList(david, inputIds=clusters_list[[i]], idType=idType, listName=names(clusters_list)[i],listType=query_type)
  
    setAnnotationCategories(david, "GOTERM_BP_ALL")
    GOBPchart <- getFunctionalAnnotationChart(david, threshold = 0.1)
    GOBPchart=GOBPchart[GOBPchart$FDR < 5,]
    setAnnotationCategories(david, "GOTERM_MF_ALL")
    GOMFchart <- getFunctionalAnnotationChart(david, threshold = 0.1)
    GOMFchart=GOMFchart[GOMFchart$FDR < 5,]
    setAnnotationCategories(david, "GOTERM_CC_ALL")
    GOCCchart <- getFunctionalAnnotationChart(david, threshold = 0.1)
    GOCCchart=GOCCchart[GOCCchart$FDR < 5,]
    setAnnotationCategories(david, "KEGG_PATHWAY")
    KEGGchart <- getFunctionalAnnotationChart(david, threshold = 0.1)
    KEGGchart=KEGGchart[KEGGchart$FDR < 5,]
    setAnnotationCategories(david, "OMIM_DISEASE")
    OMIMchart <- getFunctionalAnnotationChart(david, threshold = 0.1)
    OMIMchart=OMIMchart[OMIMchart$FDR < 5,]
#   DAVID_clusters[[i]]=list("GOBPchart"=GOBPchart,"GOMFchart"=GOMFchart,"GOCCchart"=GOCCchart,"KEGGchart"=KEGGchart,"OMIMchart"=OMIMchart)
    DAVID_clusters[[i]]=list(GOBPchart,GOMFchart,GOCCchart,KEGGchart,OMIMchart)
    names(DAVID_clusters[[i]])=c("GOBPchart","GOMFchart","GOCCchart","KEGGchart","OMIMchart")
  
  }
    names(DAVID_clusters)=names(clusters_list)
    return(DAVID_clusters)
}




DAVID_save_tables<-function(david_out,out_dir=getwd(),dataset_name=""){
  cat("\tsaving tables to :\t",paste(out_dir,dataset_name,names(david_out),".terms.txt",sep=""),"\n")

  # replace the david_out[[1]] with a higher level loop to support multiple module lists
  # - may need to add checks to determine if annotations are present (but works well enough for now)
  # - for proper pipeline may be worth editing out the un-necessary columns for easier formatting
    write.delim(names(david_out),paste(out_dir,dataset_name,names(david_out),".terms.txt",sep=""),rownm=F,colnm=F) 
  for(itrm in 1:length(david_out[[1]])){
    # spacers for easier overview
    write.delim("",paste(out_dir,dataset_name,names(david_out),".terms.txt",sep=""),rownm=F,colnm=F,append=T)
    write.delim("",paste(out_dir,dataset_name,names(david_out),".terms.txt",sep=""),rownm=F,colnm=F,append=T)
   

    annot_table=david_out[[1]][[itrm]][,c("Term","X.","FDR","Genes")]
    annot_table$Term=gsub(".*~","",annot_table$Term)          ##   remove trash from term names (should be fine to run both rather than loop KEGG/GO separately, neither have ':' or '~' in their names)
    annot_table$Term=gsub(".*:","",annot_table$Term)
    annot_table$FDR=annot_table$FDR/100                       ##   unfathonably the FDR is given as %
    annot_table$X.=annot_table$X./100                         ##   % are easier to handle in Excell as decimals
      colnames(annot_table)=c(names(david_out[[1]])[itrm],"Number of Genes (%)","Significance of Enrichment (FDR)","Ensembl Gene ID")
      colnames(annot_table)=gsub("GOBPchart","Gene ontology biologial process (BP)",colnames(annot_table))
      colnames(annot_table)=gsub("GOMFchart","Gene ontology molecular function (MF)",colnames(annot_table))
      colnames(annot_table)=gsub("GOCCchart","Gene ontology cell compartment (CC)",colnames(annot_table))
      colnames(annot_table)=gsub("KEGGchart","Gene ontology KEGG pathway",colnames(annot_table))
      colnames(annot_table)=gsub("OMIMchart","Gene ontology OMIM term",colnames(annot_table))

    
    write.delim(annot_table,paste(out_dir,dataset_name,names(david_out),".terms.txt",sep=""),append=T)

  }
  cat("\t•••••• NOTE : Multiple warnings about appending columns to file can be safely ignored ••••••\n")
}




plot_DAVID_enrich<-function(david_list_enrichments,out_path,name_module,module_info=""){
  library(gplots)
  library(ggplot2)  

if(module_info!=""){
  dir.create(file.path(out_path, module_info), showWarnings = F)   # suppresses warnings generated if folder already exists
  out_path=paste(out_path,module_info,sep="")

}

#    dir.create(file.path(out_path, "david.annot"), showWarnings = F)   # suppresses warnings generated if folder already exists
#    dir.create(file.path(paste(out_path,"/david.annot",sep=""), "plots_GO"), showWarnings = F)   # suppresses warnings generated if folder already exists
#    dir.create(file.path(paste(out_path,"/david.annot",sep=""), "plots_KEGG"), showWarnings = F)   # suppresses warnings generated if folder already exists
#    dir.create(file.path(paste(out_path,"/david.annot",sep=""), "plots_OMIM"), showWarnings = F)   # suppresses warnings generated if folder already exists

    dir.create(file.path(out_path, "david_GO"), showWarnings = F)   # suppresses warnings generated if folder already exists
    dir.create(file.path(out_path, "david_KEGG"), showWarnings = F)   # suppresses warnings generated if folder already exists
    dir.create(file.path(out_path, "david_OMIM"), showWarnings = F)   # suppresses warnings generated if folder already exists


  #combine GO terms results in single table
  GO_table=do.call(rbind,david_list_enrichments[1:3])
  if(nrow(GO_table)>0){
      GO_table=GO_table[with(GO_table, order(FDR,decreasing=T)), ]
      GO_table$minusLog10BH= -log((GO_table$FDR/100), base = 10)
      GO_table$Term=gsub("GO:[0-9]+~","",GO_table$Term)
      #add here if dim > 0, save pdf, add KEGG and OMIM
      if(nrow(GO_table) > 10){number_plot=10}else{number_plot=nrow(GO_table)}
        pdf(file=paste(out_path,"/david_GO/",module_info,"GO_functional_enrich_plot_",name_module,module_info,".pdf",sep=""),height=4,width=6)
          GOBPplot <- ggplot(GO_table[1:number_plot,c("Term","minusLog10BH")], aes(x=reorder(Term, -minusLog10BH), y = minusLog10BH)) + ylab("-Log10 (FDR)") +xlab("") +
          geom_bar(colour="black",stat= "identity", width=.4, fill = "darkred") + ggtitle(paste("GO ",name_module,module_info,sep="")) +
          theme_bw(base_size = 12, base_family = "") + coord_flip()
          print(GOBPplot)
        dev.off()

    write.delim(GO_table,paste(out_path,"/david_GO/",module_info,"GO_functional_enrich_plot_",name_module,module_info,".pdf",sep=""))
    }    


  
  if(nrow(david_list_enrichments$KEGGchart)>0){
      GO_table=david_list_enrichments$KEGGchart
      GO_table=GO_table[with(GO_table, order(FDR,decreasing=T)), ]
      GO_table$minusLog10BH= -log((GO_table$FDR/100), base = 10)
      GO_table$Term=unlist(lapply(strsplit(GO_table$Term,":"),function(x){x=x[2]}))
      if(nrow(GO_table) > 10){number_plot=10}else{number_plot=nrow(GO_table)}
        pdf(file=paste(out_path,"/david_KEGG/",module_info,"KEGG_functional_enrich_plot_",name_module,".pdf",sep=""),height=4,width=6)
          GOBPplot <- ggplot(GO_table[1:number_plot,c("Term","minusLog10BH")], aes(x=reorder(Term, -minusLog10BH), y = minusLog10BH)) + ylab("-Log10 (FDR)") +xlab("") +
          geom_bar(colour="black",stat= "identity", width=.4, fill = "darkgrey") + ggtitle(paste("KEGG ",name_module,module_info,sep="")) +
          theme_bw(base_size = 12, base_family = "") + coord_flip()
          print(GOBPplot)
        dev.off()
  
      write.delim(GO_table,paste(out_path,"/david_GO/",module_info,"KEGG_enrich_plot_",name_module,module_info,".pdf",sep=""))
    }    
  if(nrow(david_list_enrichments$OMIMchart)>0){
      GO_table=david_list_enrichments$OMIMchart
      GO_table=GO_table[with(GO_table, order(FDR,decreasing=T)), ]
      GO_table$minusLog10BH= -log((GO_table$FDR/100), base = 10)
      GO_table$Term=unlist(lapply(strsplit(GO_table$Term,":"),function(x){x=x[2]}))
      #add here if dim > 0, save pdf, add KEGG and OMIM
      if(nrow(GO_table) > 10){number_plot=10}else{number_plot=nrow(GO_table)}
        pdf(file=paste(out_path,"/david_OMIM/",module_info,"OMIM_functional_enrich_plot_",name_module,module_info,".pdf",sep=""),height=4,width=6)
          GOBPplot <- ggplot(GO_table[1:number_plot,c("Term","minusLog10BH")], aes(x=reorder(Term, -minusLog10BH), y = minusLog10BH)) + ylab("-Log10 (FDR)") +xlab("") +
          geom_bar(colour="black",stat= "identity", width=.4, fill = "darkgrey") + ggtitle(paste("OMIM ",name_module,module_info,sep="")) +
          theme_bw(base_size = 12, base_family = "") + coord_flip()
          print(GOBPplot)
        dev.off()

    write.delim(GO_table,paste(out_path,"/david_GO/",module_info,"OMIM_enrich_plot_",name_module,module_info,".pdf",sep=""))
    }    

}





heatmap_DAVID_enrich<-function(david_list_enrichments,out_path=getwd(),histH=20,histW=10,terms=c("KEGGchart","GOBPchart","GOMFchart","GOCCchart","OMIMchart"),points_sig=F,module_info=""){
  
  print("   NOTE : requires custom functions 'Heatmap' and 'make.numeric'")

####-----------------------------------------------------------------------------------------------------------
##   TERM
  for(iterm in 1:length(terms)){
    print(paste("============================",terms[iterm],"============================"))

  # create empty matrix to simplify the first merge
  TERM=as.data.frame(matrix("--",ncol=1,nrow=1))
    colnames(TERM)="Term"

  names(david_list_enrichments)=gsub("###","",names(david_list_enrichments))
  count=0
  for(imod in 1:length(david_list_enrichments)){
    if(nrow(david_list_enrichments[[names(david_list_enrichments)[imod]]][[terms[iterm]]])>=1){
      term=(david_list_enrichments[[names(david_list_enrichments)[imod]]][[terms[iterm]]][c("Term","FDR")])
      term$FDR=(term$FDR/100) # David FDR is in % => need to convert
        colnames(term)[2]=paste("M",names(david_list_enrichments)[imod],sep="")

      TERM=merge(TERM,term,by="Term",all=T)
        print(paste("module",names(david_list_enrichments)[imod],nrow(term),terms[iterm],"found"))
    }
    if(nrow(david_list_enrichments[[names(david_list_enrichments)[imod]]][[terms[iterm]]])<1){
      term=as.data.frame(matrix("--",ncol=2,nrow=1))
        colnames(term)=c("Term",paste("M",names(david_list_enrichments)[imod],sep=""))
      TERM=merge(TERM,term,by="Term",all=T)
        print(paste("module",names(david_list_enrichments)[imod],"no",terms[iterm],"found ***"))
    }
  }

  # clean up TERM matrix
    rownames(TERM)=TERM$Term
    TERM=(TERM[-which(TERM$Term=="--"),-which(colnames(TERM)=="Term")])  

# remove Pathway ID numbers
    if(terms[iterm] %in% c("GOBPchart","GOMFchart","GOCCchart")){
      rownames(TERM)=gsub("GO:.*~","",rownames(TERM))
    }
    if(terms[iterm]=="KEGGchart"){
      rownames(TERM)=gsub(".*:","",rownames(TERM))
    }

  TERM=TERM[do.call(order, as.data.frame(TERM)),]   #could be ..as.data.frame(mat[,index_vec])..
  write.table(TERM,file=paste(out_path,"combined.David.Annot.",module_info,terms[iterm],".txt",sep=""),sep="\t",col.names=NA)

    TERM[is.na(TERM)]=1 # the Heatmap function does not handle -log10(0) nor NA very well
    TERM=make.numeric(TERM)




  if(nrow(TERM)>=1){
    if(nrow(TERM<100)){   
      pdf(paste(out_path,module_info,terms[iterm],".pdf",sep=""),height=12,width=20)
      Heatmap(make.numeric(TERM),mode="pval",sig=points_sig,margin=c(5,1),cexrow=1.5,cexcol=1)  # Heatmap automatically does -log10() of the input matrix
      dev.off()
  }
  if(nrow(TERM)>100){
    pdf(paste(out_path,module_info,terms[iterm],".pdf",sep=""),height=histH,width=histW)
    itimes=floor(nrow(TERM)/100)
    for(iloop in 1:itimes){
      if(iloop==1){Heatmap(make.numeric(TERM[1:100,]),mode="pval",sig=points_sig,margin=c(5,1),cexcol=1)}
      if(iloop>1 & iloop<itimes){Heatmap(make.numeric(TERM[(iloop*100+1):(iloop*100+100),]),mode="pval",sig=points_sig,margin=c(5,1),cexcol=1)}
      if(iloop==itimes){
#        print("final loop")
        term1=TERM[(iloop*100):((iloop*100+100)-((iloop*100+100)-nrow(TERM))),]
        term2=matrix(1,nrow=100-nrow(term1),ncol=ncol(term1))
          colnames(term2)=colnames(term1)
        term=make.numeric(rbind(term1,term2))

        Heatmap(term,mode="pval",sig=points_sig,margin=c(5,1),cexcol=1)
      }
    }
    dev.off()
  }

  }
}

}






hist.norm<-function(x,normCurv=T,col='lightgray',points=F,density=F,prob=F,...){  #density=F, # need to add some sort of axis scaling is not very compatible
#  breaks=max(10,sqrt(length(x)+100))
  y=hist(x,prob=prob,...)


  if(density){  # requires prob=TRUE for probabilities not counts
     if(!prob){warning("\tWARNING: density plot is currently optimised for prob=T\n")}
    lines(density(x), col="darkred", lwd=2)                                 # add a density estimate with defaults
    lines(density(x, adjust=2), lty="dotted", col="darkgreen", lwd=2)    # add another "smoother" density


  }
  if(points){rug(x,)}   # add tick marks below histogram

  if(normCurv){
     if(prob){warning("\tWARNING: normal curve fit is currently optimised for prob=F\n")}
    a=min(x,na.rm=T)
    b=max(x,na.rm=T)
    xx=seq(a-(b-a)/10,b+(b-a)/10,length=100)
    lines(sort(xx),dnorm(sort(xx),median(x,na.rm=T),sd(x,na.rm=T))*sum(y$counts*diff(y$breaks)),col='dodgerblue')
    }
    return(invisible(y))

# to work need to construct the plot_cols and corresponding plot_lines & add a layout param dependent on legend=T to provide space for it on the left
#    legend(x="topright",pch=16,box.lwd = 0,box.col = "white",col=plot_cols,legend=plot_lines)
}



hist.dens<-function(x,normCurv=F,points=F,col='lightgray',density=T,prob=T,...){  #density=F, # need to add some sort of axis scaling is not very compatible
  y=hist(x,breaks=max(10,sqrt(length(x)+100)),prob=prob,...)


  if(density){  # requires prob=TRUE for probabilities not counts
     if(!prob){warning("\tWARNING: density plot is currently optimised for prob=T\n")}
    lines(density(x), col="darkred", lwd=2)                                 # add a density estimate with defaults
    lines(density(x, adjust=2), lty="dotted", col="darkgreen", lwd=2)    # add another "smoother" density


  }
  if(points){rug(x,)}   # add tick marks below histogram

  if(normCurv){
     if(prob){warning("\tWARNING: normal curve fit is currently optimised for prob=F\n")}
    a=min(x,na.rm=T)
    b=max(x,na.rm=T)
    xx=seq(a-(b-a)/10,b+(b-a)/10,length=100)
    lines(sort(xx),dnorm(sort(xx),median(x,na.rm=T),sd(x,na.rm=T))*sum(y$counts*diff(y$breaks)),col='dodgerblue')
    }
    return(invisible(y))

# to work need to construct the plot_cols and corresponding plot_lines & add a layout param dependent on legend=T to provide space for it on the left
#    legend(x="topright",pch=16,box.lwd = 0,box.col = "white",col=plot_cols,legend=plot_lines)
}


#hist(eruptions,,breaks=30, prob=TRUE, col="grey")# prob=TRUE for probabilities not counts
#lines(density(eruptions), col="blue", lwd=2) # add a density estimate with defaults
#lines(density(eruptions, adjust=2), lty="dotted", col="darkgreen", lwd=2) 




#### Andree Delahaye-Duriez #### 2014/01/12
      ###########################
      ### Fisher's exact test ###
      ###       pipeline      ###
      ###########################

############## to test whether each cluster of a list of clusters
############## has FET enrichment in EE DNMs compared to DNMs in NAFE controls

############## function EE_FET
EE_FET <- function(clusters_list){
  ### load data needed from the "dataToLoad_forFETinEE" folder
  load("~/Dropbox/LONDON-SINGAPORE/TOOLS/FETenrichment_in_EEDNMs_pipeline/input/dataToLoad_forFETinEE/NAFE_ctrlAllDNMs.Rdata") # path to be changed
  load("~/Dropbox/LONDON-SINGAPORE/TOOLS/FETenrichment_in_EEDNMs_pipeline/input/dataToLoad_forFETinEE/NAFE_ctrlNsDNMs.Rdata") # path to be changed
  load("~/Dropbox/LONDON-SINGAPORE/TOOLS/FETenrichment_in_EEDNMs_pipeline/input/dataToLoad_forFETinEE/EE_allDNMs.Rdata") # path to be changed
  load("~/Dropbox/LONDON-SINGAPORE/TOOLS/FETenrichment_in_EEDNMs_pipeline/input/dataToLoad_forFETinEE/EE_nsDNMs.Rdata") # path to be changed
  
  
  ### annotation of genes with ENS gene ID
  library(biomaRt)
  ensembl=useMart("ensembl")
  HUMensembl = useMart('ensembl',dataset = 'hsapiens_gene_ensembl')
  eeGene <- EE$Gene
  nseeGene <- nsEE$Gene
  nafeGene <-NAFE[,'CCDS_r14']
  nsnafeGene <-nsNAFE[,'CCDS_r14']
  
  eeENS <-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = eeGene, mart = HUMensembl)
  nseeENS <-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = nseeGene, mart = HUMensembl)
  nafeENS <-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = nafeGene, mart = HUMensembl)
  nsnafeENS <-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = nsnafeGene, mart = HUMensembl)
  
  ## keep in mind that there are sometimes several ENS gene id for only 1 external_gene_name
  #length(which(duplicated(eeENS$external_gene_name) == TRUE))
  #[1] 29
  
  ### create a matrix for results
  EE_FET_clusters <- matrix(nrow = length(clusters_list), ncol = 6)
  row.names(EE_FET_clusters) <- names (clusters_list)
  colnames(EE_FET_clusters) <- c("FET p.value all DNMs","OR all DNMs","[95% CI] all DNMs","FET Pvalue nsDNMs","OR nsDNMs","[95% CI] nsDNMs")  
  
  ### function to fill the matrix of results
  for (i in 1:length(clusters_list)){
    ## function to calculate the number Mc of DNMs in CTRL involving a gene of the cluster i
    y <- lapply(clusters_list[[i]],FUN = function(x) {nafeENS[which(nafeENS$ensembl_gene_id == x),'external_gene_name']})
    Mc <- sum(sapply(as.matrix(unique(y)),FUN = function(ym) {length(which(NAFE[,'CCDS_r14'] == ym ))}))
    
    ## number NMc of remaining DNMs in CTRL involving a gene not in the cluster i
    NMc <- nrow(NAFE)-Mc
    
    ##function to calculate the number Mee of DNMs in EE involving a gene of the cluster i
    z <- lapply(clusters_list[[i]],FUN = function(x) {eeENS[which(eeENS$ensembl_gene_id == x),'external_gene_name']})
    Mee <- sum(sapply(as.matrix(unique(z)),FUN = function(zm) {length(which(EE$Gene == zm ))}))
    
    ## number NMee of remaining DNMs in EE involving a gene not in the cluster i  
    NMee <- nrow(EE)-Mee
    
    ## function to calculate the number Mnsc of nsDNMs in CTRL NAFE involving a gene of the cluster i
    nsy <- lapply(clusters_list[[i]],FUN = function(x) {nsnafeENS[which(nsnafeENS$ensembl_gene_id == x),'external_gene_name']})
    Mnsc <- sum(sapply(as.matrix(unique(nsy)),FUN = function(ym) {length(which(nsNAFE[,'CCDS_r14'] == ym ))}))
    
    ## number NMnsc of remaining nsDNMs in CTRL involving a gene not in the cluster i
    NMnsc <- nrow(nsNAFE)-Mnsc
    
    ##function to calculate the number Mnsee of nsDNMs in EE involving a gene of the cluster i
    nsz <- lapply(clusters_list[[i]],FUN = function(x) {nseeENS[which(nseeENS$ensembl_gene_id == x),'external_gene_name']})
    Mnsee <- sum(sapply(as.matrix(unique(nsz)),FUN = function(zm) {length(which(nsEE$Gene == zm ))}))
    
    ## number NMnsee of remaining nsDNMs in EE involving a gene not in the cluster i  
    NMnsee <- nrow(nsEE)-Mnsee

    # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
    matrALL <- matrix(c(Mee,Mc,NMee,NMc), nrow = 2)
    matrNS <- matrix(c(Mnsee,Mnsc,NMnsee,NMnsc), nrow = 2)
    
    # FET
    FisherMEE <- fisher.test(matrALL)
    FisherMEEp <- FisherMEE$p.value
    FisherMEEor <- FisherMEE$estimate
    FisherMEEci <- paste(FisherMEE$conf.int[1],FisherMEE$conf.int[2], sep="-")
    FisherMnsEE <- fisher.test(matrNS)
    FisherMnsEEp <- FisherMnsEE$p.value
    FisherMnsEEor <- FisherMnsEE$estimate
    FisherMnsEEci <- paste(FisherMnsEE$conf.int[1],FisherMnsEE$conf.int[2], sep="-")
    
    EE_FET_clusters[i,]=c(FisherMEEp,FisherMEEor,FisherMEEci,FisherMnsEEp,FisherMnsEEor,FisherMnsEEci)
  }
  write.table(EE_FET_clusters, sep = '\t', file = 'EE_FET_clusters.txt', row.names = TRUE, quote = FALSE, col.names = NA)
  return(EE_FET_clusters)
}



gwas.enrich<-function(in_gwas,module,bkgrnd,nperm=100000,seed=0,type=1){
cat("\n\tFUNCTION : gwas.enrich function inputs : in_gwas, module, bkgrnd, nperm=100000, seed=0, type=1\n")
cat("\tINPUTS : in_gwas - link to a csv file, ENSG - pval  |  module - list of modules  |  bkgrnd - matrix first column contains the ENSG list\n\n")
  set.seed(seed)
  options(stringsAsFactors=FALSE)
  gwas<-read.csv(file=in_gwas,sep="\t",header=FALSE)
#  module<-read.table(file=in_module,sep=" ")
#  bkgrnd<-read.table(file=bkgrnd,sep=" ")
  
  enrich=as.data.frame(matrix(NA,nrow=length(names(module)),ncol=5))
    rownames(enrich)=names(module)
    colnames(enrich)=c("n.genes.module","pc.overlap.genes","GWAS.P.value","GWAS.FDR","GWAS.Bonferroni")


    for(imod in names(module)){
      idx<-match(bkgrnd[,1],gwas[,1])
      idx<-sort(idx)
      gwas_final<-gwas[idx,]

      idy<-match(as.matrix(module[[imod]])[,1],gwas_final[,1])
      idy<-sort(idy)
      module_final<-gwas_final[idy,]
      enrich[imod,"n.genes.module"]=length(module[[imod]])
      enrich[imod,"pc.overlap.genes"]=round(length(idy)/length(module[[imod]]),digits=3)

    ##   type 1 = z-test (USE THIS ONE)
      if(type==1){
        cat("\t'z-test' enrichment   ||   ",imod,"\t: ",which(names(module)==imod),"of",length(names(module)),"\n")
        module_p<-mean(-log10(module_final[,2]),na.rm=TRUE)  
      }
    ##  type 2 = fishers combined p-value
      else if(type==2){
        cat("  'fisher's combined p-value' enrichment   ||   ",imod," : ",which(names(module)==imod),"of",length(names(module)),"\n")
        df<-2*length(module_final[,2])
        temp1<-log(module_final[,2])
        temp2<- -2*sum(temp1)
        module_p<-pchisq(temp2,df,lower.tail=FALSE)
      }
    ##  type 3 = stouffer combined p-value
      else if(type==3){
        cat("  'stouffer combined p-value' enrichment   ||   ",imod," : ",which(names(module)==imod),"of",length(names(module)),"\n")
        temp1<-qnorm(1-module_final[,2])/sqrt(length(module_final[,2]))
        module_p<-sum(temp1[!is.infinite(temp1)])
      }
      output_stat<-vector()
      for(i in 1:nperm){
        rand1<-sample(gwas_final[,2],length(module_final[,2]),replace=FALSE)
        rand1<-rand1[!is.infinite(rand1)]

        if(type==1){
          output_stat[i]<-mean(-log10(rand1),na.rm=TRUE)  
        }else if(type==2){
          df<-2*length(rand1)
          temp1<-log(rand1)
          temp2<- -2*sum(temp1)
          output_stat[i]<-pchisq(temp2,df,lower.tail=FALSE)
        }else if(type==3){
          tempt<-qnorm(1-rand1)/sqrt(length(rand1))
          output_stat[i]<-sum(tempt[!is.infinite(tempt)])
        }
        
          cat(round(i/nperm,digits=2),"\r");flush.console()
      }

    output_stat_sd<-sd(output_stat)
    output_stat_mean<-mean(output_stat)

     Z<-(module_p - output_stat_mean)/(output_stat_sd)
    # print(Z)
    # enrich<-pnorm(abs(Z),low=FALSE)
    enrich[imod,"GWAS.P.value"]<-pnorm(module_p,output_stat_mean,output_stat_sd,lower.tail=FALSE)

    #enrich[imod,"n.genes.module"]=nrow(module_final)

    # testingWTF1<-t.test(output_stat,mu=module_p)
    # testingWTF2<-t.test(output_stat,mu=module_p,alternative="less")
    # testingWTF3<-t.test(output_stat,mu=module_p,alternative="greater")
    }

  enrich$GWAS.FDR=p.adjust(enrich$GWAS.P.value,method="fdr")
  enrich$GWAS.Bonferroni=p.adjust(enrich$GWAS.P.value,method="bonferroni")
  return(enrich)

}



p.adjust.mat<-function(p_mat,method='fdr',single_col=F,single_row=F,verbose=F){
##  p.adjust appears to be designed to handle NA values 

	if(single_row & single_col){
		 stop('\n\tERROR :\tto perform "p.adjust()" across whole matrix use "single_col=F" & "single_row=F", ie default\n\n')
	}
	if(!single_col & !single_row){
		if(verbose){cat('\n\t"p.adjust()" on whole matrix simultaneously using "method =',method,'"\n\n')}

		adj_mat=matrix(p.adjust(unlist(p_mat),method=method),nrow=nrow(p_mat))
		rownames(adj_mat)=rownames(p_mat)
		colnames(adj_mat)=colnames(p_mat)
		return(invisible(adj_mat))
	}

	if(single_col & ! single_row){
		 if(verbose){cat('\n\t"p.adjust()" individually on each column of supplied matrix "method =',method,'"\n\n')}
		adj_mat=matrix(NA,nrow=nrow(p_mat),ncol=ncol(p_mat))
		 rownames(adj_mat)=rownames(p_mat)
		 colnames(adj_mat)=colnames(p_mat)
		for(icol in colnames(p_mat)){
			adj_mat[,icol]=p.adjust(p_mat[,icol],method=method)
		}		
		return(invisible(adj_mat))	
	}

	if(single_row){
		 if(verbose){cat('\n\t"p.adjust()" individually on each row of supplied matrix "method =',method,'"\n\n')}
		adj_mat=matrix(NA,nrow=nrow(p_mat),ncol=ncol(p_mat))
		 rownames(adj_mat)=rownames(p_mat)
		 colnames(adj_mat)=colnames(p_mat)
		for(irow in colnames(p_mat)){
			adj_mat[irow,]=p.adjust(p_mat[irow,],method=method)
		}		
		return(invisible(adj_mat))	
	}

}


gwas_module_enrich<-function(in_gwas,module,bkgrnd,nperm=100000,seed=0,type=1){
cat("\tINPUTS : in_gwas - link to a csv file, ENSG - pval  |  module - list of modules  |  bkgrnd - matrix first column contains the ENSG list")
  set.seed(seed)
  options(stringsAsFactors=FALSE)
  gwas<-read.csv(file=in_gwas,sep=",",header=FALSE)
#  module<-read.table(file=in_module,sep=" ")
#  bkgrnd<-read.table(file=bkgrnd,sep=" ")
  
  enrich=as.data.frame(matrix(NA,nrow=length(names(module)),ncol=5))
    rownames(enrich)=names(module)
    colnames(enrich)=c("n.genes.module","pc.overlap.genes","GWAS.P.value","GWAS.FDR","GWAS.Bonferroni")


    for(imod in names(module)){
      idx<-match(bkgrnd[,1],gwas[,1])
      idx<-sort(idx)
      gwas_final<-gwas[idx,]

      idy<-match(as.matrix(module[[imod]])[,1],gwas_final[,1])
      idy<-sort(idy)
      module_final<-gwas_final[idy,]
      enrich[imod,"n.genes.module"]=length(module[[imod]])
      enrich[imod,"pc.overlap.genes"]=round(length(idy)/length(module[[imod]]),digits=3)

    ##   type 1 = z-test (USE THIS ONE)
      if(type==1){
        cat("  'z-test' enrichment     ||     ",imod," : ",which(names(module)==imod),"of",length(names(module)),"\n")
        module_p<-mean(-log10(module_final[,2]),na.rm=TRUE)  
      }
    ##  type 2 = fishers combined p-value
      else if(type==2){
        cat("  'fisher's combined p-value' enrichment   ||   ",imod," : ",which(names(module)==imod),"of",length(names(module)),"\n")
        df<-2*length(module_final[,2])
        temp1<-log(module_final[,2])
        temp2<- -2*sum(temp1)
        module_p<-pchisq(temp2,df,lower.tail=FALSE)
      }
    ##  type 3 = stouffer combined p-value
      else if(type==3){
        cat("  'stouffer combined p-value' enrichment   ||   ",imod," : ",which(names(module)==imod),"of",length(names(module)),"\n")
        temp1<-qnorm(1-module_final[,2])/sqrt(length(module_final[,2]))
        module_p<-sum(temp1[!is.infinite(temp1)])
      }
      output_stat<-vector()
      for(i in 1:nperm){
        rand1<-sample(gwas_final[,2],length(module_final[,2]),replace=FALSE)
        rand1<-rand1[!is.infinite(rand1)]

        if(type==1){
          output_stat[i]<-mean(-log10(rand1),na.rm=TRUE)  
        }else if(type==2){
          df<-2*length(rand1)
          temp1<-log(rand1)
          temp2<- -2*sum(temp1)
          output_stat[i]<-pchisq(temp2,df,lower.tail=FALSE)
        }else if(type==3){
          tempt<-qnorm(1-rand1)/sqrt(length(rand1))
          output_stat[i]<-sum(tempt[!is.infinite(tempt)])
        }
        
          cat(round(i/nperm,digits=2),"\r");flush.console()
      }

    output_stat_sd<-sd(output_stat)
    output_stat_mean<-mean(output_stat)

     Z<-(module_p - output_stat_mean)/(output_stat_sd)
    # print(Z)
    # enrich<-pnorm(abs(Z),low=FALSE)
    enrich[imod,"GWAS.P.value"]<-pnorm(module_p,output_stat_mean,output_stat_sd,lower.tail=FALSE)

    #enrich[imod,"n.genes.module"]=nrow(module_final)

    # testingWTF1<-t.test(output_stat,mu=module_p)
    # testingWTF2<-t.test(output_stat,mu=module_p,alternative="less")
    # testingWTF3<-t.test(output_stat,mu=module_p,alternative="greater")
    }

  enrich$GWAS.FDR=p.adjust(enrich$GWAS.P.value,method="fdr")
  enrich$GWAS.Bonferroni=p.adjust(enrich$GWAS.P.value,method="bonferroni")
  return(enrich)

}


sm<-function(x){
  return(as.matrix(sort(x,decreasing=T)))
}

smt<-function(x){
  return(as.matrix(sort(table(x),decreasing=T)))
}



wgcnaPickBeta<-function(list_expr,saveTable=F,savePlots=F,outDir=getwd(),datDescr=""){
print("  NOTE : Input data (list_expr) is expected as a list with each entry : rows = genes, columns = samples")
print("   - currently does not return the tables - information is printed // can also be saved as plots and tables")

  library('WGCNA')
    allowWGCNAThreads()
    set.seed(0)       # reproducibility 

  for(ireg in 1:length(names(list_expr))){
    run_time=Sys.time()
    print(paste("--------------------",names(list_expr)[ireg],"---------------------",ireg,"of",length(names(list_expr))))

  thChoice=pickSoftThreshold(as.matrix(t(list_expr[[names(list_expr)[ireg]]])),corFn='bicor')
  print(thChoice$power)

  if(saveTable==T){
    write.table(thChoice$fitIndices[order(thChoice$fitIndices$mean.k.,decreasing=T),],paste(outDir,"tables/WGCNA.PowerTABLE-(beta).",datDescr,".",names(list_expr)[ireg],".Power",thChoice$power,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }

  if(savePlots==T){
    pdf(file=paste(outDir,"plots/WGCNA.PowerPLOT-(beta).",datDescr,".",names(list_expr)[ireg],".Power",thChoice$power,".pdf",sep=""),width=15,height=5)

    #sizeGrWindow(9, 5)
    par(mfrow = c(1,3))
    cex1 = 0.9

    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(thChoice$fitIndices[2:10,"Power"], -sign(thChoice$fitIndices[2:10,"slope"])*thChoice$fitIndices[2:10,"SFT.R.sq"],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale independence"))
    text(thChoice$fitIndices[2:10,"Power"], -sign(thChoice$fitIndices[2:10,"slope"])*thChoice$fitIndices[2:10,"SFT.R.sq"],
      labels=thChoice$fitIndices[2:10,"Power"],cex=cex1,col="red")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.8,col="red")

    # Mean connectivity as a function of the soft-thresholding power
    plot(thChoice$fitIndices[2:10,"Power"], thChoice$fitIndices[2:10,"mean.k."],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean connectivity"))
    text(thChoice$fitIndices[2:10,"Power"], thChoice$fitIndices[2:10,"mean.k."],
      labels=thChoice$fitIndices[2:10,"Power"], cex=cex1,col="red")


    plot(thChoice$fitIndices[2:10,"mean.k."],-sign(thChoice$fitIndices[2:10,"slope"])*thChoice$fitIndices[2:10,"SFT.R.sq"],xlim=c(ceiling(max(thChoice$fitIndices[2:10,"mean.k."])),0)
      ,ylab="Mean Connectivity",xlab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))

    text(thChoice$fitIndices[2:10,"mean.k."],-sign(thChoice$fitIndices[2:10,"slope"])*thChoice$fitIndices[2:10,"SFT.R.sq"],labels=thChoice$fitIndices[c(2:10),"Power"],cex=1,col="red")
    abline(h=0.8,col="red")

    dev.off()

  }
  print(Sys.time()-run_time)
  }
#  return(pickSoftThreshold)  <= needs to be returned as a list of all the objects
}





wgcna.beta<-function(expr_mat,corFn='cor',networkType='signed',corOptions='spearman',...){
 cat('\tNOTE\t: use - apply WGCNA pickSoftThreshold() to matrix: samples = columns, genes = rows\n\n')
if(corFn=='cor'){cat('\t\t correlation function: ',corFn,' - ',corOptions,', network type: ',networkType,'\n\n',sep='')}
if(corFn=='bicor'){cat('\t\t correlation function',corFn,'network type:',networkType,'\n\n')}

  library('WGCNA')
    allowWGCNAThreads()
    set.seed(0)       # reproducibility 

  thChoice=pickSoftThreshold(t(expr_mat),...)
  cat('\n\tpower sugggested by WGCNA:',thChoice$power,'\n\n',verbose=1)

  return(invisible(thChoice))
}

wgcna.mods<-function(expr_mat,softPower=7,signType='signed',mergeCutHei=0.15,dat_descr='',...){
 cat('\tUSE:\tperform basic WGCNA clustering analysis for a given matrix, return full outputs as list\n')
  library(WGCNA)
  allowWGCNAThreads()
  t1=Sys.time()

  minModSize=min(40, nrow(expr_mat)/2 )
# mergeCutHei=0.15

     modules=blockwiseModules(t(expr_mat),   # Input data, expect: samples - rows, genes - columns
      # Data checking options ----------------------------------------------
         checkMissingData = TRUE,
       
      # Options for splitting data into blocks ----------------------------------------------
         blocks = NULL,
         maxBlockSize = 30000,
         blockSizePenaltyPower = 5,
         randomSeed = 12345,
       
      # if load TOM from previously saved file  ----------------------------------------------
      loadTOM = FALSE,
       
      # Network construction arguments: correlation options ----------------------------------------------
         corType = "bicor",            # pearson - default
         maxPOutliers = 1, 
         quickCor = 0,
         pearsonFallback = "individual",
         cosineCorrelation = FALSE,
       
      # Adjacency function options ----------------------------------------------
         power = softPower,
         networkType = signType,      # options 'signed', 'unsigned', default = 'signed'
       
      # Topological overlap options ----------------------------------------------
         TOMType = "signed",
         TOMDenom = "min",
       
      # Saving or returning TOM ----------------------------------------------
         getTOMs = NULL,
         saveTOMs = FALSE, 
         saveTOMFileBase = "blockwiseTOM",
   
    # Basic tree cut options ========================================================================
         deepSplit = 2,                # default=2 simplified version of tree-cutting  (may be worth looking at tree first and defining proper threshold)
         detectCutHeight = 0.995,    # dendrogram cut height for module detection. See cutreeDynamic for more details
         minModuleSize = minModSize,   # minimum module size for module detection. See cutreeDynamic for more details
        
      # Advanced tree cut options----------------------------------------------
         maxCoreScatter = NULL,         # maximum scatter of the core for a branch to be a cluster, given as the fraction of cutHeight relative to the 5th   percentile of joining heights. See cutreeDynamic for more details
         minGap = NULL,
         maxAbsCoreScatter = NULL, minAbsGap = NULL,
         minSplitHeight = NULL, minAbsSplitHeight = NULL,
       
         useBranchEigennodeDissim = FALSE,
         minBranchEigennodeDissim = mergeCutHei,
       
         stabilityLabels = NULL,
         minStabilityDissim = NULL,
       
         pamStage = TRUE, pamRespectsDendro = TRUE,
       
      # Gene reassignment, module trimming, and module "significance" criteria ----------------------------------------------
         reassignThreshold = 1e-6,
         minCoreKME = 0.5, 
#        minCoreKMESize = minModuleSize/3,    minModuleSize not found error..?
         minKMEtoStay = 0.3,
       
      # Module merging options ----------------------------------------------
          mergeCutHeight = mergeCutHei, 
         impute = TRUE, 
         trapErrors = FALSE, 
       
      # Output options ----------------------------------------------
         numericLabels = FALSE,
       
      # Options controlling behaviour ----------------------------------------------
         nThreads = 0,
         verbose = 1,      # options '0', '1' , '2' , '3'
         indent = 0,...)

  collectGarbage()

    print(Sys.time()-t1)


## generating mstat is a bit messy..
    mstat=as.data.frame(table(modules$colors))
    mstat=mstat[order(mstat[,2],decreasing=T),]
    msta0=mstat[(mstat[,1]=='grey'),]
    msta0$module='M0'
    mstat=mstat[!(mstat[,1]=='grey'),]
    mstat$module=paste0('M',1:nrow(mstat))
    mstat=rbind(mstat,msta0)
      colnames(mstat)=c('color','ngenes','module')
    mstat$color=as.character(mstat$color)
#    mstat$module=paste(mstat$module,mstat$color,sep="_")   ##  module name contains color
if(dat_descr!=''){mstat$module=paste(mstat$module,dat_descr,sep="_")}   ##  add info after module name

  mes=modules$MEs
    rownames(mes)=colnames(expr_mat)
    colnames(mes)=gsub('ME','',colnames(mes))
  mes=mes[,mstat$color]


  module_list=list()
  module_expr=list()
  for(imod in 1:nrow(mstat)){
    module_expr[[mstat$module[imod]]]=expr_mat[modules$colors==mstat$color[imod],]
    module_list[[mstat$module[imod]]]=rownames(module_expr[[mstat$module[imod]]])

  }

    mbg=as.data.frame('bkgrnd')
    mbg$length=nrow(expr_mat)
    mbg$name='bkgrnd'

      colnames(mbg)=colnames(mstat)
  mstat=rbind(mstat,mbg)


    module_expr[['bkgrnd']]=expr_mat
    module_list[['bkgrnd']]=rownames(expr_mat)


  cat('\n\tModules are named based on size M1 - biggest, M0 - unclustered, bkgrnd - all input genes, output contains :
    \t1.\tmodule_list - list containing names of genes in each module
    \t2.\tmodule_expr - expression matrix of all genes in module
    \t3.\tMEs         - module eigengenes (PC1) for each module
    \t4.\tmstat       - key used to name modules, includes module size
    \t5.\twgcna_out   - object containing full output of WGNCA blockwiseModules()
    \n')
  return(invisible(list(module_list=module_list,module_expr=module_expr,MEs=mes,mstat=mstat,wgcna_out=modules)))
}




wgcna.consensus<-function(list_expr,dat_descr='',corType='spearman',power = 6,signType = "signed",max_block_n=20000,...){
 cat("\tNOTE :\tInput data (list_expr) is expected as a list with each entry : rows = genes, columns = samples\n")
# cat("\tWARNING :\tthe checkMissingData option for WGCNA is disabled, set checkMissingData=T in the code to make robust to missing, alternatively use sd.check or is.missing to determine/remove non-varying or missing values\n")
  library(WGCNA)
  allowWGCNAThreads()

#  minModSize=min(40, nrow(expr_mat)/2 )

  for(ilis in 1:length(list_expr)){
    list_expr[[names(list_expr)[ilis]]]=list(data=t(list_expr[[names(list_expr)[ilis]]]))
  }

  t0=Sys.time()
  
  modules=blockwiseConsensusModules(
  list_expr          ##  lists with $data in each set containing expression: rows=samples, cols=genes

# Data checking options ------------------------------------------
  ,checkMissingData = T     ##  checks for missing and zero variance in expression.. 
     
# Blocking options ------------------------------------------
  ,blocks = NULL
  ,maxBlockSize = max_blocks        ##  ensure no separation is performed // if outdated machine ~4GB RAM, may need to change this to ~5,000, or better yet, consider upgrading
  ,blockSizePenaltyPower = 5
    ,randomSeed = 12345
     
# TOM precalculation arguments, if available ------------------------------------------
  ,individualTOMInfo = NULL         ##  pre-calculated TOM using blockwiseIndividualTOMs
     ,useIndivTOMSubset = NULL
     
# Network construction arguments: correlation options ------------------------------------------
# ,corType = "pearson"
     ,maxPOutliers = 1        ##  used only with bicor: Specifies the maximum percentile of data that can be considered outliers on either side of the median separately
     ,quickCor = 0            ##  handling of missing
     ,pearsonFallback = "individual"    ## using pearson when mean absolute deviation == zero -> cant perform bicor...
  ,cosineCorrelation = FALSE
     
# Adjacency function options ------------------------------------------
# ,power = 6                  ## moved to function inputs
  ,networkType = signType   ## moved to function inputs  ##  options: "unsigned"’, ‘"signed"’, ‘"signed hybrid" See ‘adjacency’
     ,checkPower = TRUE
# Topological overlap options ------------------------------------------
  ,TOMType = signType      ## moved to function inputs   ##  options: "none"’, ‘"unsigned"’, ‘"signed"’. If ‘"none"’, adjacency will be used for clustering
     ,TOMDenom = "min"          ## min - standard, mean - expreimental

# Save individual TOMs?  ------------------------------------------
  ,saveIndividualTOMs = F    ## TOM saved to disk
     ,individualTOMFileNames = "individualTOM-Set%s-Block%b.RData"
     
# Consensus calculation options: network calibration ------------------------------------------
  ,networkCalibration = c("single quantile","full quantile","none")     ## "single", "quantile", "full quantile", "none"
     
# Simple quantile calibration options ------------------------------------------
  ,calibrationQuantile = 0.95
     ,sampleForCalibration = TRUE
     ,sampleForCalibrationFactor = 1000
     ,getNetworkCalibrationSamples = FALSE
     
# Consensus definition ------------------------------------------
  ,consensusQuantile = 0       ## <<<<<<<<<<<<<<<<
     ,useMean = FALSE             ## <<<<<<<<<<<<<<<< use mean instead of consensusQuantile
     ,setWeights = NULL           ## <<<<<<<<<<<<<<<< for weighted mean
     
# Saving the consensus TOM ------------------------------------------
  ,saveConsensusTOMs = FALSE
     ,consensusTOMFileNames = "consensusTOM-block.%b.RData"
     
# Internal handling of TOMs ------------------------------------------
  ,useDiskCache = F         ## TRUE = slower but more RAM efficient
     ,chunkSize = NULL
     ,cacheBase = ".blockConsModsCache"
     ,cacheDir = "."
     
# Alternative consensus TOM input from a previous calculation  ------------------------------------------
  ,consensusTOMInfo = NULL
     
# Basic tree cut options  ------------------------------------------
  ,deepSplit = 2
  ,detectCutHeight = 0.995
     ,minModuleSize = 50           ##  default = 20
     ,checkMinModuleSize = TRUE
     
# Advanced tree cut options ------------------------------------------
  ,maxCoreScatter = NULL
     ,minGap = NULL
     ,maxAbsCoreScatter = NULL
     ,minAbsGap = NULL
     ,minSplitHeight = NULL
     ,minAbsSplitHeight = NULL
     ,useBranchEigennodeDissim = FALSE
     ,minBranchEigennodeDissim = mergeCutHeight
     ,stabilityLabels = NULL
     ,minStabilityDissim = NULL
  ,pamStage = TRUE
     ,pamRespectsDendro = TRUE

# Gene reassignment and trimming from a module, and module "significance" criteria ------------------------------------------
  ,reassignThresholdPS = 1e-4
#     ,trimmingConsensusQuantile = consensusQuantile    ##  breaks if specified here without outside variables
     ,minCoreKME = 0.5
#     ,minCoreKMESize = minModuleSize/3
     ,minKMEtoStay = 0.2
     
# Module eigengene calculation options ------------------------------------------
  ,impute = TRUE
     ,trapErrors = FALSE
     
  #Module merging options
  ,equalizeQuantilesForModuleMerging = FALSE
     ,quantileSummaryForModuleMerging = "mean"
     ,mergeCutHeight = 0.15
#  ,mergeConsensusQuantile = consensusQuantile    ##  breaks if specified here without outside variables
     
# Output options ------------------------------------------
  ,numericLabels = FALSE
# General options ------------------------------------------
  ,nThreads = 0
     ,verbose = 2
     ,indent = 1
     ,...
     )

     print(Sys.time()-t0)
     collectGarbage()


    mstat=as.data.frame(table(modules$colors))
    mstat=mstat[order(mstat[,2],decreasing=T),]
    msta0=mstat[(mstat[,1]=='grey'),]
    msta0$module='M0'
    mstat=mstat[!(mstat[,1]=='grey'),]
    mstat$module=paste0('M',1:nrow(mstat))
    mstat=rbind(mstat,msta0)
      colnames(mstat)=c('color','ngenes','module')
    mstat$color=as.character(mstat$color)

if(dat_descr!=''){mstat$module=paste(mstat$module,dat_descr,sep="_")}   ##  add info after module name

## module membership ------------------------------------------
  module_list=list()
    for(imod in 1:nrow(mstat)){
        module_list[[mstat$module[imod]]]=colnames(list_expr[[1]]$data)[modules$color==mstat$color[imod]]        ##  assuming all module gene names are in the same order, if not, modules are probably unreliable anyhow
    }
    module_list[['bkgrnd']]=colnames(list_expr[[1]]$data)


  mes=modules$multiMEs                              ##  comes as a list, one for each list input
      names(mes)=names(list_expr)

  module_expr=list()
  for(ilis in 1:length(list_expr)){

##  module eigene output ------------------------------------------
    rownames(mes[[names(list_expr)[ilis]]]$data)=rownames(list_expr[[names(list_expr)[ilis]]]$data)
    dummy=mes[[names(list_expr)[ilis]]]$data        ## § get rid of the $data sub-tag in every list..
    mes[[names(list_expr)[ilis]]]=dummy

##  re-name MEs to match module names
        colnames(mes[[names(list_expr)[ilis]]])=gsub('ME','',colnames(mes[[names(list_expr)[ilis]]]))
    mes[[names(list_expr)[ilis]]]=mes[[names(list_expr)[ilis]]][,mstat$color]
        colnames(mes[[names(list_expr)[ilis]]])=paste('me',mstat$module,sep='') # change the names to meM1_descr

## module expression ------------------------------------------
    for(imod in 1:length(module_list)){
        module_expr[[names(list_expr)[ilis]]][[names(module_list)[imod]]]=list_expr[[names(list_expr)[ilis]]]$data[,module_list[[names(module_list)[imod]]]]
    }
  }


## add bkgrnd info to mstat
    mbg=as.data.frame('bkgrnd')
    mbg$length=ncol(list_expr[[1]]$data)
    mbg$name='bkgrnd'

      colnames(mbg)=colnames(mstat)
  mstat=rbind(mstat,mbg)

  readme='\n\tModules are named based on size M1 - biggest, M0 - unclustered, bkgrnd - all input genes, output contains :
    \t1.\tmodule_list - list containing names of genes in each module
    \t2.\tmodule_expr - expression matrix of all genes in module / input dataset
    \t3.\tMEs         - module eigengenes (PC1) for each module / dataset
    \t4.\tmstat       - key used to name modules, includes module size
    \t5.\twgcna_out   - object containing full output of WGNCA blockwiseConsensusModules(), can be used for plots etc
    \n'
  cat(readme)
  return(invisible(list(module_list=module_list,module_expr=module_expr,MEs=mes,mstat=mstat,wgcna_out=modules,readme=readme)))
}



wgcnaBuildModules<-function(list_expr,pow_range,minModuleSize=40,mergeHeight=0.15,savePlots=F,outDir=getwd(),datDescr="",signType="signed"){
print("  NOTE : Input data (list_expr) is expected as a list with each entry : rows = genes, columns = samples")
if (savePlots==F){print("   - currently does not return the tables - information is printed // can also be saved as plots and tables")}

  library('WGCNA')
    allowWGCNAThreads()
    set.seed(0)       # reproducibility 

for(ireg in 1:length(names(list_expr))){
# pamStage=1            #  seems irrelevant for ult networks - ie cv, pval 0.5, log transformed etc

  print("■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■")
  cat("    range of power scan :",pow_range,"\n")
  print(paste(names(list_expr)[ireg],ireg,"of",length(names(list_expr))))
  
  system(paste("mkdir ",outDir,"/objects",sep=""))
  system(paste("mkdir ",outDir,"/modules",sep="")) 
  system(paste("mkdir ",outDir,"/plots",sep=""))
  #system(paste("mkdir ",outDir,"/plots/pw",softPower,sep="")) 

  for(ipow in pow_range){
      t0=Sys.time()
    softPower=ipow

  #merge.height.scan=c(0.15,0.2,0.25,0.3)
  #for(ihei in mergeHeight){
    ihei=mergeHeight

      print(softPower)
      print(names(list_expr)[ireg])


  wgcnaModules=blockwiseModules(t(as.matrix(list_expr[[names(list_expr)[ireg]]])), 
      maxBlockSize = 20000,
      power=softPower,
      #deepSplit=2,           # simplified version of tree-cutting  (may be worth looking at tree first and defining proper threshold)
      mergeCutHeight = ihei,
      #mergeCutHeight=0.25,
      #detectCutHeight = 0.995,     # dendrogram cut height for module detection. See cutreeDynamic for more details
      minModuleSize = minModuleSize,         # minimum module size for module detection. See cutreeDynamic for more details

      # maxCoreScatter          # maximum scatter of the core for a branch to be a cluster, given as the fraction of cutHeight relative to the 5th
                        # percentile of joining heights. See cutreeDynamic for more details
      corType="bicor",          # pearson - default
      #reassignThreshold = 0, 
      networkType=signType,
      #TOMtype="signed",  #is default
      numericLabels = FALSE,
      #pamStage = 1, pamRespectsDendro = TRUE,
      #TOMDenom="min",          # min=default/standard, mean - may produce better results but "experimental" atm
      #saveTOMs = TRUE,
      #saveTOMFileBase = "femaleMouseTOM",
      verbose = 1 # options '1' , '2' , '3'
    )

      mergedColor<-wgcnaModules$colors
      collectGarbage()


    print(Sys.time()-t0)

# save modules
  #module_ensg_expr=list()
  module_ensg_list=list()
  module_ensg_list[["M0_grey"]]=rownames(list_expr[[names(list_expr)[ireg]]])[mergedColor=="grey"]
  write.table(rownames(list_expr[[names(list_expr)[ireg]]])[mergedColor=="grey"],file=paste(outDir,"/modules/M0_",datDescr,"_power",softPower,".txt",sep=""),quote=F,row.names=F,col.names=F)
    n<-1
  for(imod in unique(mergedColor[-which(mergedColor =="grey")])) {
    module_ensg_list[[paste("M",n,"_",imod,sep="")]]= rownames(list_expr[[names(list_expr)[ireg]]])[mergedColor==imod]
    write.table(rownames(list_expr[[names(list_expr)[ireg]]])[mergedColor==imod],file=paste(outDir,"/modules/M",n,"_",imod,"_",datDescr,"_power",softPower,".txt",sep=""),quote=F,row.names=F,col.names=F)
    n<-n+1
  }
  module_ensg_list[["bkgrnd"]]=rownames(list_expr[[names(list_expr)[ireg]]])
  write.table(module_ensg_list[["bkgrnd"]],file=paste(outDir,"/modules/bkgrnd_",datDescr,"_power",softPower,".txt",sep=""),quote=F,row.names=F,col.names=F)    

  save(wgcnaModules,module_ensg_list,file=paste0(outDir,"/objects/WGCNA_MODULES_",
    names(list_expr)[ireg],"_",datDescr,
    ".merge-height=",ihei,
    ".genes=",nrow(list_expr[[names(list_expr)[ireg]]]),
    ".samples=",ncol(list_expr[[names(list_expr)[ireg]]]),
    "_pw",softPower,"_bicor.R"))

  if(savePlots==T){
    pdf(file=paste(outDir,"/plots/pw",softPower,"_WGCNA_dendo_",names(list_expr)[ireg],datDescr,
      #peer,
      #paste(covarnm$short,collapse="."),
      ".merge-height=",ihei,
      ".genes=",nrow(list_expr[[names(list_expr)[ireg]]]),
      ".samples=",ncol(list_expr[[names(list_expr)[ireg]]]),
      "_bicor.pdf",sep=""))
    ###--------------------------------------------------------------------------------------------------------------
    # Plot the cut tree   -----------------------------------------------------------------------------------------

    #sizeGrWindow(12, 9)    # open a graphics window

    plotDendroAndColors(wgcnaModules$dendrograms[[1]], wgcnaModules$colors[wgcnaModules$blockGenes[[1]]],
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main=paste("Cluster Dendrogram",names(list_expr)[ireg],"(power",softPower,")"))
    #mergedColors = labels2colors(wgcnaModules$colors)    # Convert labels to colors for plotting
    #plotDendroAndColors(wgcnaModules$dendrograms[[1]], colors[wgcnaModules$blockGenes[[1]]],
    # "Module colors",
    # dendroLabels = FALSE, hang = 0.03,
    # addGuide = TRUE, guideHang = 0.05) # Plot the dendrogram and the module colors underneath

    #quartz(height=7,width=7) # opens a new plot window of specified size

    #quartz()


    ##pdf(file=paste("./Dropbox/Cognition/WGCNA/finalised/plots/WGCNA.PDF.merge-height=",ihei,".pam=",pamStage,".",nrow(D1),"genes.",ncol(list_expr[[names(list_expr)[ireg]]]),"patients.powr",softPower,resid,quanNorm,logTransf,cv.filt,resid,".bicor.MEplot.only.pdf",sep=""))

    ## looking at the clustering of the MEs of the modules --------------------------------------------------------------------------

    #me.curr=(MEcurr[,-which(colnames(MEcurr)=="grey")])              # omit grey cluster
    #MEDiss = 1-cor(me.curr)                          # calculate the distances

    #MEcurr=wgcnaModules$MEs
    #rownames(MEcurr)=colnames(list_expr[[names(list_expr)[ireg]]])
    #colnames(MEcurr) = gsub("ME", "", colnames(MEcurr),perl=T)   # get rid of the ME_ prefix in colors
    # head(MEcurr)
    # dim(MEcurr)

    #write.table(MEcurr,paste("~/caprica/wgcna/tables/WGCNA.MEcurr.table.",bonn_background_restricted,names(list_expr)[ireg],peer,pc1correct_data,paste(covarnm$short,collapse="."),".merge-height=",ihei,".genes",nrow(list_expr[[names(list_expr)[ireg]]]),".samples",ncol(list_expr[[names(list_expr)[ireg]]]),"powr",softPower,".bicor.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

    #MEDiss = 1-cor(MEcurr)                         # calculate the distances
    #METree = hclust(as.dist(MEDiss), method = "average")           # Cluster module eigengenes



    #plot(METree, main = paste("Clustering of ",bonn_background_restricted,names(list_expr)[ireg],"module eigengenes (power",softPower,")"), xlab = "", sub = "",cex=0.7)   # Plot the result
    # abline(h=0.3,col="orange")
    # abline(h=0.25,col="red")
    # abline(h=0.20,col="blue")
    # abline(h=0.15,col="green")
    #

    dev.off()
  }
  }
}
}





mod.name.col <- function(module_ensg_list){
  mcol=module_ensg_list$bkgrnd
  for (i in 1:(length(module_ensg_list)-1)) {
    mcol[mcol%in%module_ensg_list[[i]]]=names(module_ensg_list[i])
    
  }
  return(mcol)
}



peer.correct<-function(expDat,covDat=matrix(),sanity=F,cov_match=T,verbose=T){
#  if(verbose){
#  cat("\tNOTE : exprDat expected as matrix/data.frame of expression columns=genes rows=samples\n")
#  cat("\tNOTE : covDat is optional covariate as a matrix columns=covariates, rows=samples, PEER requires numeric as.numeric(as.factor(cov))\n")
#  }

####-----------------------------------------------------------------------------------------------------------
###  PEER data correction for covariates and factors

# PEER package installation
# https://github.com/PMBio/peer/wiki/Installation-instructions
# R CMD INSTALL R_peer_source_1.3.tgz

####   PEER DATA CORRECTION   ---------------------------------------------------------------------------------
### Applying PEER in R amounts to creating the model and setting its parameters, followed by inference. The model object will then hold the posterior distributions of the parameters.

### PEER can automatically include an additional factor (covariate) to account for the mean expression. For most use cases, including the mean effect is likely to be a good choice. To active mean factors, use
### Run correlation analyses between the inferred variables and batch confounding effects. If several inferred factors correlated with batch effects/confounders, this can be indicative of
### a more complex, nonlinear effect of these known covariates on the mRNA levels. Scatter plots can help understand the nature of these dependencies.
  library(peer)

 # peercor=list()  # covariate adjusted PEER output
 # peerfac=list()  # list of factors generated by PEER and/or user-provided covariates used by PEER 

### sub-select covariate data based on exprDat -> potentially a problem if not in correct setup /&/ missing / non-matching names

    model = PEER()                          ##  create model object
    PEER_setPhenoMean(model, t(expDat))    ##  add observed data

    if(nrow(covDat)==1 & ncol(covDat)==1){
      cat("\tno covariates selected\n")
    }

    if(nrow(covDat)>1 & ncol(covDat)>=1){
        ### add in covariates, can be used to correct the data -> residuals retrieved at the end
        ### covariates - samples = rows, have to be numeric, :: as.numeric(as.factor())
    if(cov_match){
      cat('\tcovariate data matched to expression\n')
      covDat=covDat[rownames(covDat)%in%colnames(expDat),,drop=F]
    }


#        covDat
       cat("\t\tcovariates selected\t:",paste(colnames(covDat),collapse=", ")," : ", length(colnames(covDat)),"\n")
       cat("\t\tcovariates contain\t:", round(sum(rownames(covDat)==colnames(expDat))/nrow(covDat),digits=3)*100,"% matching samples \n")

      PEER_setCovariates(model, covDat)  ##  add covariates to model
    }

    #dim(PEER_getPhenoMean(model))

    ### Number of factors to use (ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/pdf/ukmss-48786.pdf)
    ###  If no prior information on the magnitude of confounding effects is available,
    ###  we recommend using 25% of the number of individuals contained in the study but no more than 100 factors.
#n_peer_fac=3
#if(ncol(expDat)<=30){n_peer_fac=3;cat('\tn samples </= 30, using',ncol(expDat)/n_peer_fac,'factors\n')}
#if(ncol(expDat)>30){n_peer_fac=4;cat('\tn samples > 30, using',ncol(expDat)/n_peer_fac,'factors\n')}
    n_peer_fac=4
    PEER_setNk(model,ceiling(ncol(expDat)/n_peer_fac))  ##  PEER recommends running n/4 factors => these should be checked to determine what they measure || ensure they are not correcting for biological variation
    PEER_getNk(model)


    PEER_update(model)    ###  perform the inference.

    ### The result is the model object with posterior distributions of the variables.
    ###  2. Observing output

    ### You can get the posterior mean of the inferred confounders (NxK matrix), their weights (GxK matrix), precision (inverse variance) of the weights (Kx1 matrix), and the residual dataset (NxG matrix):
    resid = t(PEER_getResiduals(model))    ##  the residual dataset (NxG matrix)  ## i.e. corrected data for all factors
    factr   = PEER_getX(model)            ##  the posterior mean of the inferred confounders (NxK matrix)
    weights   = PEER_getW(model)            ##  weights (GxK matrix) of the inferred confounders (NxK matrix)
    precision = PEER_getAlpha(model)        ##  precision (inverse variance) of the weights (Kx1 matrix)

#    print(cor(PEER_getX(model), PEER_getCovariates(model)))
    colnames(resid)=colnames(expDat)
    rownames(resid)=rownames(expDat)

    rownames(factr)=colnames(expDat)
    colnames(factr)=paste0('F',1:ncol(factr))

    colnames(weights)=colnames(factr)
    rownames(weights)=rownames(expDat)

    rownames(precision)=colnames(factr)
    colnames(precision)='inverse.var.weights'

## **************   need to re-name colnames and rownames in peer output since it does not save these automatically *************************

    if(nrow(covDat)>1 & ncol(covDat)>=1){
       readme="\n\tResults obtained via peerCovCorrect function using PEER PMC3398141 bayes PMC4158865
       \n\t\t1. resid - the residual dataset (NxG matrix)  ## i.e. corrected data for all factors
       \n\t\t2. factr - covaraites the data was corrected for - can be used as 'phenotypes'  ||  the posterior mean of the inferred confounders (NxK matrix)
       \n\t\t3. weights - weights (GxK matrix) of the inferred confounders (NxK matrix)
       \n\t\t4. precision - precision (inverse variance) of the weights (Kx1 matrix)  ||  plot(precision)
       \n\t\t5. model -  PEER model object can be used to perform diagnostic plot:  PEER_plotModel(model)
       \n\t\t6. expr - original expression matrix
       \n\t\t7. covar - original covariate matrix
       \n"

      return(list("expr"=expDat,"covar"=covDat,"resid"=resid,"factr"=factr,weights=weights,precision=precision,model=model,"readme"=readme))
    }

    if(nrow(covDat)==1 & ncol(covDat)==1){
       readme="\n\tResults obtained via peerCovCorrect function using PEER PMC3398141 bayes PMC4158865
       \n\t\t1. resid - the residual dataset (NxG matrix)  ## i.e. corrected data for all factors
       \n\t\t2. factr - covaraites the data was corrected for - can be used as 'phenotypes'  ||  the posterior mean of the inferred confounders (NxK matrix)
       \n\t\t3. weights - weights (GxK matrix) of the inferred confounders (NxK matrix)
       \n\t\t4. precision - precision (inverse variance) of the weights (Kx1 matrix)  ||  plot(precision)
       \n\t\t5. model -  PEER model object can be used to perform diagnostic plot:  PEER_plotModel(model)
       \n\t\t6. expr - original expression matrix
       \n"

      return(list("expr"=expDat,"resid"=resid,"factr"=factr,weights=weights,precision=precision,model=model,"readme"=readme))

    }



 cat("\tpeercov - peer corrected expression residuals\n")
 cat("\tpeerfac - covaraites the data was corrected for - can be used as 'phenotypes'\n")


}




peerCovCorrect<-function(listExpr,covDat=matrix(),sanity=F){
  cat("\tNOTE : listExpr expected as list of expression matrixes columns=genes rows=samples\n")
  cat("\tNOTE : covDat is optional covariate as a matrix columns=covariates, rows=samples\n")
####-----------------------------------------------------------------------------------------------------------
###  PEER data correction for covariates and factors

# PEER package installation
# https://github.com/PMBio/peer/wiki/Installation-instructions
# R CMD INSTALL R_peer_source_1.3.tgz

####   PEER DATA CORRECTION   ---------------------------------------------------------------------------------
### Applying PEER in R amounts to creating the model and setting its parameters, followed by inference. The model object will then hold the posterior distributions of the parameters.

### PEER can automatically include an additional factor (covariate) to account for the mean expression. For most use cases, including the mean effect is likely to be a good choice. To active mean factors, use
### Run correlation analyses between the inferred variables and batch confounding effects. If several inferred factors correlated with batch effects/confounders, this can be indicative of
### a more complex, nonlinear effect of these known covariates on the mRNA levels. Scatter plots can help understand the nature of these dependencies.
  library(peer)

  peercov=list()
  peerfac=list()

  

#pdf(paste("~/caprica/graphics/02.rma.plier-gcbg.",paste(covariates_selected,collapse="."),".peer.plots.pdf",sep=""))
# par(mfrow=c(1,2))

  for(i in 1:length(names(listExpr))){

    print(paste("--------------------",names(listExpr)[i],"---------------------",i,"of",length(names(listExpr))))
  # expr[[names(listExpr)[i]]]=read.delim(paste("~/dtb/primary/ensg.expression.",names(listExpr)[i],".txt",sep=""),row.names=1)

    # the whole of covarDat (below) can have no entries + comment out PEER_setCovariates

  if(nrow(covDat)==1 & ncol(covDat)==1){
      cat("no covariates selected\n")   
  }

    ### We run PEER three times: 
    ###  1. Using 8 factors and all covariates


    ### The data matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes.

    ### Now we can create the model object,

    model = PEER()
    PEER_setPhenoMean(model, t(listExpr[[names(listExpr)[i]]]))    ### set the observed data

    ### add in covariates, can be used to correct the data -> residuals retrieved at the end
    ### covariates - samples = rows, have to be numeric, :: as.numeric(as.factor())

    if(nrow(covDat)!=1 & ncol(covDat)!=1){

#     covarDat=as.matrix(covDat[colnames(listExpr[[names(listExpr)[i]]]),colnames(covDat) %in% covariates_selected,drop=F]) # ************* special provision required when ONLY ONE covariate is selected in covarDat -> to retain row and column names
      covarDat=as.matrix(covDat[intersect(rownames(covDat),colnames(listExpr[[names(listExpr)[i]]])),]) # ************* special provision required when ONLY ONE covariate is selected in covarDat -> to retain row and column names
    if(sanity){
      covarDat=sd.check(covarDat,T,F)
      
      }
       cat("\t\tcovariates selected\t:",paste(colnames(covarDat),collapse=", ")," : ", length(colnames(covarDat)),"\n")
       cat("\t\tcovariates contain\t:", round(sum(rownames(covarDat)==colnames(listExpr[[1]]))/nrow(covarDat),digits=3)*100,"% matching samples \n")
#      print(dim(covarDat))
        # making sure that this is selected only when there are covariates chosen
      PEER_setCovariates(model, covarDat)
    }
    # NULL
    dim(PEER_getPhenoMean(model))
    ## [1]    102 17410
    ## (NULL response means no error here), say we want to infer K=10 hidden confounders,

    ### Number of factors to use (ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/pdf/ukmss-48786.pdf)
    ###  If no prior information on the magnitude of confounding effects is available,
    ###  we recommend using 25% of the number of individuals contained in the study but no more than 100 factors.

    ### 102 samples / 4 => we use 25 factors
    PEER_setNk(model,ncol(listExpr[[i]])/4)
    #  NULL
    PEER_getNk(model)
    # [1] 25
    ### and perform the inference.
    PEER_update(model)
    ### The result is the model object with posterior distributions of the variables.
    ###  2. Observing output

    ### You can get the posterior mean of the inferred confounders (NxK matrix), their weights (GxK matrix), precision (inverse variance) of the weights (Kx1 matrix), and the residual dataset (NxG matrix):
    factors = PEER_getX(model)
      dim(factors)
    # [1] 102  25
    weights = PEER_getW(model)
      dim(weights)
    # [1] 17410    25
    precision = PEER_getAlpha(model)
      dim(precision)
    # [1] 25  1
    peercov[[names(listExpr)[i]]] = PEER_getResiduals(model)
    peerfac[[names(listExpr)[i]]] = factors
      dim(peercov$cov5peer_vst_filt)
    # [1]  102 17410


    cor (PEER_getX(model), PEER_getCovariates(model))

  ####   plotting metrics   ------------------------------ 
#      plot(precision,main=paste(names(listExpr)[i]))
#      PEER_plotModel(model)
#  }

#    dev.off()


  # peer does not save row/column names (also the matrixes are transposed so need to adjust for that)
  #peercor=list()


#  for(i in 1:length(names(listExpr))){
    print(paste(names(listExpr)[i],i,"of",length(names(listExpr))))
    peercorex=as.data.frame(t(peercov[[names(listExpr)[i]]]))

      rownames(peercorex)=rownames(listExpr[[names(listExpr)[i]]])
      colnames(peercorex)=colnames(listExpr[[names(listExpr)[i]]])
#    peercov[[names(listExpr)[i]]]=peercorex

      if(nrow(covDat)!=1 & ncol(covDat)!=1){
        peerfactr=as.data.frame(peerfac[[names(listExpr)[i]]])
          rownames(peerfactr)=rownames(covarDat)
          colnames(peerfactr)=c(colnames(covarDat),paste("F",1:(ncol(peerfactr)-ncol(covarDat)),sep=""))
 #       peerfac[[names(listExpr)[i]]]=peerfactr
      }
  # write.delim(peercov[[names(listExpr)[i]]],paste("~/dtb/secondary/peerCorrected.",names(listExpr)[i],".txt",sep=""),rownm=T,colnm=T)
  # write.delim(peercov[[names(listExpr)[i]]],paste("~/dtb/secondary/peerCorrected.cv.",names(listExpr)[i],".txt",sep=""),rownm=T,colnm=T)
  # write.delim(peercov[[names(listExpr)[i]]],paste("~/dtb/secondary/peerCorrected.bonn_restrict.",names(listExpr)[i],".txt",sep=""),rownm=T,colnm=T)
  }

 readme="\n\tResults obtained via peerCovCorrect function using PEER PMC3398141 bayes PMC4158865
 \n\t\t1. peercov - peer corrected expression residuals
 \n\t\t2. peerfac - covaraites the data was corrected for - can be used as 'phenotypes'\n"
# \n\t\t3. expr - original expression matrix
# \n\t\t4. covar - original covariate matrix\n")

return(list("peercov"=peercov,"peerfac"=peerfac,"readme"=readme))
#return(list("peercov"=peercov,"peerfac"=peerfac,"expr"=listExpr,"covar"=covDat,"readme"=readme))
 cat("\tpeercov - peer corrected expression residuals\n")
 cat("\tpeerfac - covaraites the data was corrected for - can be used as 'phenotypes'\n")
# cat("\texpr - original expression matrix\n")
# cat("\tcovar - original covariate matrix\n")
 cat("\texpr and covar - original inputs are no longer returned by peerCovCorrect()\n")
}






pcVarExpl<-function(cexpr){
  #
  corpc=list()
  pcstats=list()
  #pie.cols=c("#D55E00","#F0E442", "#009E73", "#56B4E9","#0072B2") 

  #pdf(paste("~/Dropbox/CapricaPrime/RU/out/plots/02.PC.ENSG.var.explained.",cv_correct,peer,paste(corrected,collapse="."),".pdf",sep=""),height=14,width=7)
  # par(mfrow=c(5,2),cex=0.6)

  #
  for(ireg in 1:length(names(cexpr))){
      print(paste("--------------------",names(cexpr)[ireg],"---------------------",ireg,"of",length(names(cexpr))))

    pcs=prcomp(t(cexpr[[names(cexpr)[ireg]]]),scale=T,center=T)
  #   plot(pcs$x[,1],pcs$x[,2],main=paste("PC1 v PC2",names(cexpr)[ireg],paste(corrected,collapse=".")))


    pcstat=as.matrix(summary(pcs)$importance["Proportion of Variance",1:5])
    # pie(pcstat, main=paste("Variance Explained (%) \n by PC1-5 of expr",names(cexpr)[ireg],paste(corrected,collapse=".")), radius =1,
    #   labels =paste("PC",1:5," ",round(pcstat*100),"%",sep=""), col=pie.cols)
    print(pcstat)

    pcstats[[names(cexpr)[ireg]]]=pcstat
    
  #
    corpc[[names(cexpr)[ireg]]]=as.data.frame(lm(as.matrix(t(cexpr[[names(cexpr)[ireg]]]))~pcs$x[,1])$residuals)


  # top and bottom 10% of genes in the loadings
  # octx pc1 - ubiquitin, proteolysis
  # octx pc1 - interesting things
  # correct for pc1 or use more stuff from *** http://www.nature.com/nprot/journal/v7/n3/full/nprot.2011.457.html?WT.ec_id=NPROT-201203

  # write.delim(names(pcs$rotation[(pcs$rotation[,1] > quantile(pcs$rotation[,1],0.9) | pcs$rotation[,1] < quantile(pcs$rotation[,1],0.1)),1]),paste("~/Dropbox/CapricaPrime/RU/out/tables/",bonn_background_restricted,cv_correct,peer,pc1correct_data,paste(corrected,collapse="."),".",names(cexpr)[ireg],".PC1.top.bottom.10pc.loadings.txt",sep=""),rownm=F,colnm=F)
  # write.delim(names(pcs$rotation[(pcs$rotation[,2] > quantile(pcs$rotation[,2],0.9) | pcs$rotation[,2] < quantile(pcs$rotation[,2],0.1)),2]),paste("~/Dropbox/CapricaPrime/RU/out/tables/",bonn_background_restricted,cv_correct,peer,pc1correct_data,paste(corrected,collapse="."),".",names(cexpr)[ireg],".PC2.top.bottom.10pc.loadings.txt",sep=""),rownm=F,colnm=F)
  # write.delim(names(pcs$rotation[(pcs$rotation[,2] > quantile(pcs$rotation[,2],0.9) | pcs$rotation[,2] < quantile(pcs$rotation[,2],0.1)),2]),paste("~/Dropbox/CapricaPrime/RU/out/tables/",bonn_background_restricted,cv_correct,peer,pc1correct_data,paste(corrected,collapse="."),".",names(cexpr)[ireg],".PC3.top.bottom.10pc.loadings.txt",sep=""),rownm=F,colnm=F)
  }

  #dev.off()

  return(invisible(corpc))
}







frac<-function(subset,full,num=T,sig_fig=2){
  if(!num){return(round(length(subset)/length(full),digits=sig_fig))}
  if(num){return(round(subset/full,digits=sig_fig))}

}


rmduplicates<-function(matrix,colName){
  clean=matrix[!duplicated(matrix[,which(colnames(matrix)==colName)]),]

      cat("     ",sum(duplicated(matrix[,which(colnames(matrix)==colName)]))," duplicates removed  || ",round(nrow(clean)/nrow(matrix),digits=3)*100,"% of data remaining || ",sum(duplicated(clean[,which(colnames(clean)==colName)])),"dupilicates remaining \n")
return(invisible(clean))
      }





list.as.df<-function(in_list){
  dflist=as.data.frame(matrix(NA,nrow=length(in_list),ncol=ncol(in_list[[1]])))
    colnames(dflist)=colnames(in_list[[1]])
    rownames(dflist)=names(in_list)

  for(idat in 1:length(in_list)){
    dflist[idat,]=in_list[[idat]]
  }
  return(unique(dflist))
}





gnameToENSG<-function(gene_names){
  hugo=read.delim("~/Dropbox/annotation/dtb/hgnc_complete_set.8.4.2015.txt")
#  sum(nur$Gene.Symbol %in% hugo$Approved.Symbol)
#    length(unique(nur$Gene.Symbol))
#  sum(esc$Gene.Symbol %in% hugo$Approved.Symbol)
#    length(unique(esc$Gene.Symbol))


  gene_names=gene_names[gene_names!=""]
  gene_names=as.data.frame(gene_names)

    colnames(gene_names)="Gene.Name"

    cat("\tinput : ",length(gene_names$Gene.Name)," genes   ||  ",length(unique(gene_names$Gene.Name))," unique\n")

  mapped=hugo[hugo$Approved.Symbol %in% gene_names$Gene.Name,c("Approved.Symbol","Ensembl.Gene.ID","Entrez.Gene.ID")]
  mapped$input.id=mapped$Approved.Symbol
  mapped$input.id.type="Approved.Symbol"


    cat("\t\t",nrow(mapped)," directly mapped   ||  ",length(unique(mapped$Approved.Symbol))," unique\n")

  unmapped=unique(gene_names$Gene.Name[!(gene_names$Gene.Name %in% hugo$Approved.Symbol)])
#  unmapped=unmapped[unmapped!=""]

  #


  remap=list()
  imatch=1

  for(igen in 1:nrow(hugo)){
  #   cat("\t\t\t",length(unlist(strsplit(hugo$Synonyms[igen],", "))),"     ",sum(unlist(strsplit(hugo$Synonyms[igen],", ")) %in% unmapped),"\n")
    if(sum(unlist(strsplit(hugo$Synonyms[igen],", ")) %in% unmapped)>0){
        cat("\t\t\t",igen,"  Synonyms               ",length(unlist(strsplit(hugo$Synonyms[igen],", "))),"     ",sum(unlist(strsplit(hugo$Synonyms[igen],", ")) %in% unmapped),"\n")
      remap[[imatch]]=(hugo[igen,c("Approved.Symbol","Ensembl.Gene.ID","Entrez.Gene.ID")])
      remap[[imatch]]$input.id=unlist(strsplit(hugo$Synonyms[igen],", "))[unlist(strsplit(hugo$Synonyms[igen],", ")) %in% unmapped]
      remap[[imatch]]$input.id.type="Synonym"
      imatch=imatch+1
    }

    if(sum(unlist(strsplit(hugo$Previous.Symbols[igen],", ")) %in% unmapped)>0){
        cat("\t\t\t",igen,"  Previous Symbols   ",length(unlist(strsplit(hugo$Previous.Symbols[igen],", "))),"     ",sum(unlist(strsplit(hugo$Previous.Symbols[igen],", ")) %in% unmapped),"\n")
      remap[[imatch]]=hugo[igen,c("Approved.Symbol","Ensembl.Gene.ID","Entrez.Gene.ID")]
      remap[[imatch]]$input.id=unlist(strsplit(hugo$Previous.Symbols[igen],", "))[unlist(strsplit(hugo$Previous.Symbols[igen],", ")) %in% unmapped]
      remap[[imatch]]$input.id.type="Prevous Symbol"
      imatch=imatch+1
    }

    if(sum(unlist(strsplit(hugo$Previous.Names[igen],", ")) %in% unmapped)>0){
        cat("\t\t\t",igen,"  Previoius Names    ",length(unlist(strsplit(hugo$Previous.Names[igen],", "))),"     ",sum(unlist(strsplit(hugo$Previous.Names[igen],", ")) %in% unmapped),"\n")
      remap[[imatch]]=hugo[igen,c("Approved.Symbol","Ensembl.Gene.ID","Entrez.Gene.ID")]
      remap[[imatch]]$input.id=unlist(strsplit(hugo$Previous.Names[igen],", "))[unlist(strsplit(hugo$Previous.Names[igen],", ")) %in% unmapped]
      remap[[imatch]]$input.id.type="Prevous Name"
      imatch=imatch+1
    }

    if(sum(unlist(strsplit(hugo$Name.Synonyms[igen],", ")) %in% unmapped)>0){
        cat("\t\t\t",igen,"  Name Synonyms      ",length(unlist(strsplit(hugo$Name.Synonyms[igen],", "))),"     ",sum(unlist(strsplit(hugo$Name.Synonyms[igen],", ")) %in% unmapped),"\n")
      remap[[imatch]]=hugo[igen,c("Approved.Symbol","Ensembl.Gene.ID","Entrez.Gene.ID")]
      remap[[imatch]]$input.id=unlist(strsplit(hugo$Name.Synonyms[igen],", "))[unlist(strsplit(hugo$Name.Synonyms[igen],", ")) %in% unmapped]
    remap[[imatch]]$input.id.type="Name Synonym"
      imatch=imatch+1
    }


    cat(round(igen/nrow(hugo),digits=2),"\r");flush.console()
  }


#  sum(nrow(mapped),length(remap))/nrow(gene_names)
#  nrow(gene_names)-sum(nrow(mapped),length(remap))

  length(remap)
  remap=list.as.df(remap)
    length(unique(remap$Approved.Symbol))

    cat("\n\t\t",length(unique(mapped$Approved.Symbol)),"remapped using synonyms\n\n")
    dim(remap)
  mapped=rbind(mapped,remap)
    dim(mapped)
#    dim(esc)


    # possible duplicates due to multiple 'old id type' for the same gene
    #erm[erm$old.id %in% erm$old.id[duplicated(erm$old.id)],]


qc=mapped
mapped=unique(mapped[,-which(colnames(mapped)=="input.id.type")])
mapped=mapped[!(mapped$input.id %in% mapped$input.id[(duplicated(mapped$input.id))]),]
unmapped=gene_names[!(gene_names$Gene.Name %in% mapped$input.id),]
ambigous=mapped[mapped$input.id %in% mapped$input.id[(duplicated(mapped$input.id))],]


print(mapped[mapped$input.id %in% mapped$input.id[(duplicated(mapped$input.id))],])
print(as.matrix(unmapped))

    return(list(mapped=mapped,unmapped=unmapped,ambigous=ambigous,qc=qc))

}




bgcommon<-function(list_dat,transform=F,dat_mat='',union=F,help=F,verbose=T){
if(help){
 cat("\n\tUSE\t: determine common background (union or intersect) across all entries in expression list\n")
 cat("\tNOTE\t: list_dat - list of expression matrices: row=genes, col=samples\n")
 cat("\tNOTE\t: transfrom=T also returns the list_dat processed for common bg\n")
 cat("\tNOTE\t: dat_mat required only if union=T & transform=T \n")
}

   bgcommon=rownames(list_dat[[1]])
  for(ilis in 2:length(list_dat)){
     if(!union){bgcommon=intersect(bgcommon,rownames(list_dat[[ilis]]))}
     if(union){bgcommon=union(bgcommon,rownames(list_dat[[ilis]]))}
  }

  if(!union&verbose){cat("\n\t\tintersect:",length(bgcommon),"genes common to all",length(list_dat),"datasets\n")}
  if(union&verbose){cat("\n\t\tunion",length(bgcommon),"genes common to all",length(list_dat),"datasets\n")}

  
    
  if(!transform){return(invisible(bgcommon))}

  if(transform & !union){
   listt=list()
    for(ilis in 1:length(list_dat)){
     listt[[names(list_dat)[ilis]]]=list_dat[[names(list_dat)[ilis]]][bgcommon,,drop=F]
    }

    return((listt))
  }

  if(transform & union){
   listt=list()
    for(ilis in 1:length(list_dat)){
     listt[[names(list_dat)[ilis]]]=dat_mat[bgcommon,colnames(list_dat[[names(list_dat)[ilis]]]),drop=F]
    }

    return((listt))
  }

}



bgfix<-function(mod_list,bg_vec){
  bgcommon=intersect(unique(unlist(mod_list)),bg_vec)
	cat('\n% genes remaining\n')
  mod_common=list()
  mod_common[[names(mod_list)[1]]]=intersect(mod_list[[names(mod_list)[1]]],bgcommon)
   cat('\n\t',names(mod_list)[1],'    \t',round(length(mod_common[[names(mod_list)[1]]])/length(mod_list[[names(mod_list)[1]]]),digits=2),'\n')
  for(ilis in 2:length(mod_list)){
    mod_common[[names(mod_list)[ilis]]]=intersect(mod_list[[names(mod_list)[ilis]]],bgcommon)
      cat('\t',names(mod_list)[ilis],'    \t',round(length(mod_common[[names(mod_list)[ilis]]])/length(mod_list[[names(mod_list)[ilis]]]),digits=2),'\n')
  }
  return(invisible(mod_common))
}




bgintersect<-function(list1,bg1){ #list2,,bg2
#bgintersect<-function(list1,bg.common){
    cat("\t\tUSE : background intersect to only retain genes present in bkgrnd provided\n")
#    cat("\t\tNOTE : list1 and list1 are expected as lists of IDs, bg1 and bg2 expected as vector, all are required\n\n")
#  bg.common=(intersect(bg1,bg2))
#    cat("\t",length(bg.common),"genes shared between datasets\t1:",length(bg1),pc(bg.common,bg1)*100,"%    ||     2:",length(bg2),pc(bg.common,bg2)*100,"%\n\n")

  alis=list()
  for(ilis in 1:length(list1)){
    alis[[names(list1)[ilis]]]=list1[[ilis]][list1[[ilis]] %in% bg1]
 #     cat("\t",names(list1)[ilis],"\t--------",length(alis[[ilis]]),"/",length(list1[[ilis]]),"\t:",pc((alis[[ilis]]),(list1[[ilis]])),"\t--------",ilis,"of",length(list1),"\n")
  }

#  blis=list()
#  for(ilis in 1:length(list2)){
#    blis[[names(list2)[ilis]]]=list2[[ilis]][list2[[ilis]] %in% bg.common]
#      cat("\t",names(list2)[ilis],"\t--------",length(blis[[ilis]]),"/",length(list2[[ilis]]),":\t",pc((blis[[ilis]]),(list2[[ilis]])),"\t--------",ilis,"of",length(list2),"\n")
#  }

return(invisible(alis))
#  return(list1=alis)
}


#bgintersect<-function(list1,list2,bg1,bg2){
#bgintersect<-function(list1,bg.common){
#    cat("\t\tUSE : background intersect to only retain genes present in both lists provided\n")
#    cat("\t\tNOTE : list1 and list1 are expected as lists of IDs, bg1 and bg2 expected as vector, all are required\n\n")
#  bg.common=(intersect(bg1,bg2))
#    cat("\t",length(bg.common),"genes shared between datasets\t1:",length(bg1),pc(bg.common,bg1)*100,"%    ||     2:",length(bg2),pc(bg.common,bg2)*100,"%\n\n")

#  alis=list()
#  for(ilis in 1:length(list1)){
#    alis[[names(list1)[ilis]]]=list1[[ilis]][list1[[ilis]] %in% bg.common]
 #     cat("\t",names(list1)[ilis],"\t--------",length(alis[[ilis]]),"/",length(list1[[ilis]]),"\t:",pc((alis[[ilis]]),(list1[[ilis]])),"\t--------",ilis,"of",length(list1),"\n")
#  }

#  blis=list()
#  for(ilis in 1:length(list2)){
#    blis[[names(list2)[ilis]]]=list2[[ilis]][list2[[ilis]] %in% bg.common]
#      cat("\t",names(list2)[ilis],"\t--------",length(blis[[ilis]]),"/",length(list2[[ilis]]),":\t",pc((blis[[ilis]]),(list2[[ilis]])),"\t--------",ilis,"of",length(list2),"\n")
#  }

#return(list(list1=alis,list2=blis,bg_common=bg.common))
#  return(list1=alis)
#}



empirical_mod_cons<-function(expr_cor_mat,mod_list,nperm=1000,do_plots=F,dat_descr=''){
  if(nrow(expr_cor_mat)!=ncol(expr_cor_mat)){
    cat('\n\tWARNING: expect a square matrix of gene correlations\n')
  }
  cat('\t')
  perm_stat=list()
  clus_stat=list()

  diag(expr_cor_mat)=NA
  expr_cor_mat=abs(expr_cor_mat)

  for(imod in 1:length(mod_list)){
    t2=Sys.time()
      cat("    --------- ",names(mod_list)[imod],length(mod_list[[imod]]),"genes","---------------------",imod,"of",length(mod_list)," \n")
      clus_stat[[names(mod_list)[imod]]]$mean.abscor_matrix=mean(expr_cor_mat[mod_list[[imod]],mod_list[[imod]]],na.rm=T)
#      clus_stat[[names(mod_list)[imod]]]$mean2abscor_matrix=mean(apply(expr_cor_mat[mod_list[[imod]],mod_list[[imod]]],1,mean,na.rm=T),na.rm=T)    # the results appear identical
      clus_stat[[names(mod_list)[imod]]]$clust.length=length(mod_list[[names(mod_list)[imod]]])

      for(iper in 1:nperm){
        perm.ids=sample(rownames(expr_cor_mat),size=clus_stat[[names(mod_list)[imod]]]$clust.length,replace=FALSE)
        perm.cor=as.matrix(expr_cor_mat[perm.ids,perm.ids])
        diag(perm.cor)=NA

      perm_stat[[names(mod_list)[imod]]][[as.character(iper)]]=mean(abs(perm.cor),na.rm=T)
        cat(round(iper/nperm,digits=2),"\r");flush.console()
      }

      clus_stat[[names(mod_list)[imod]]]$pval=(sum(perm_stat[[names(mod_list)[imod]]]>clus_stat[[names(mod_list)[imod]]]$mean.abscor_matrix)+1)/(nperm+1)
      cat('\temprirical P=',clus_stat[[names(mod_list)[imod]]]$pval,'\n')
      print(Sys.time()-t2)
  }

#### separating this out to optimise plotting axis across all plots

    if(do_plots){
#      x_lim=c((min(unlist(perm_stat),unlist(clus_stat))-0.2),(max(unlist(perm_stat),unlist(clus_stat))+0.2))
      x_lim=c(0,0.5)
      for(imod in 1:length(mod_list)){
        hist(perm_stat[[names(mod_list)[imod]]],breaks=30,xlim=x_lim,main=paste(names(mod_list)[imod],dat_descr,'\n',nperm,"permutations, p=",round(clus_stat[[names(mod_list)[imod]]]$pval,digits=2)),xlab="Average absolute correlation")
        abline(v=clus_stat[[names(mod_list)[imod]]]$mean.abscor_matrix,col="red")
      }
    }

  return(invisible(list(clust_stat=t(as.data.frame(lapply(clus_stat,unlist))),perm_stat=perm_stat)))
}





##### more on visualisation of PCs biplot() and alternative ggbiplot()
###http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
##### vital info if ever using ggplots
###http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/


#par(fig=c(0,0.8,0,1), new=TRUE)
#plot(ab$pcs,col=colmix[match(dat_vect,levels(dat_vect))],pch=16)
#par(fig=c(0.7,1,0,1),new=TRUE)
#legend("center",legend=unique(c1$SMTSD),col=1:length(c1$SMTSD),pch=1,cex=1)


pcstat<-function(expr_mat,plot_pcs=c(F,F,F),col="black",n_pcs=5,dat_descr="",pie_radius=1,help=F,verbose=T){
  if(help==T){cat("\tINPUTS\t:  expr_mat : rows - genes, columns - samples  |  n_pcs - number of PCs to use\n\n")}
# does not workpar???
#  if(sum(plot_pcs>1)){
#    par(mfrow=c(1,sum(plot_pcs)))
#     par(mfrow=c(2,2))
#  }
  pcs=prcomp(t(expr_mat),scale=T,center=T)
#   coll=c("#0072B2" ,"#E69F00","#009E73","#56B4E9","#D55E00","#66A61E","#7570B3","#882255","#F0E442","#A6761D","#AA4499","#117733","#332288","#999933")

  pcstat=round(as.matrix(summary(pcs)$importance["Proportion of Variance",1:n_pcs]),digits=2)*100
#  cat("\t\tPCs variance explained :\n\n")
    colnames(pcstat)="variance explained %"
    print(round(pcstat,digits=3))

if(plot_pcs[1]==T){
  pie(pcstat, main=paste("Variance Explained (%) \n by PC1-5 of expr",dat_descr), radius =pie_radius,
      labels =paste("PC",1:n_pcs," ",pcstat,"%",sep=""), col=colmix[1:n_pcs])
  }
if(plot_pcs[2]==T){
  plot(pcs$x[,1],pcs$x[,2],col=col,main=paste("PC1 v PC2",dat_descr),pch=16,frame.plot=F,xlab=paste("PC1 ",pcstat[1],"%",sep=""),ylab=paste("PC2 ",pcstat[2],"%",sep=""))


   }
#   plot(cbind())

  pcs=pcs$x[,1:n_pcs]

if(plot_pcs[3]==T){
  matplot(pcs,type="l",lty=1,lwd=2,col=colmix)
}
  return(invisible(pcs))
}



pcplot<-function(dat_list,colmix="",pch=16,dat_descr="",main="",legend_space=8,...){
   cat("\tUSE\t: use list of multiple datasets to plot PC1vPC2\n")
   cat("\tNOTE\t: list1 - expect a list of expression matrices, rows=genes, columns=samples |*| assumes same row order\n")
   cat("\tNOTE\t: if legend is clipping the plot, increase plot width or legend_space param\n")
   
   if(class(dat_list)!="list"){
    warning("\tWARNING\t: input is not a list\n")
   }

   if(colmix[1]==""){
    cat("\t\tusing default color mix\n")
    colmix=c("#0072B2","#E69F00","#009E73","#56B4E9","#D55E00","#66A61E","#7570B3","#a50f15","#A6761D","#117733","#332288","#b15928","#882255","#999933","#AA4499","#1f78b4","#F0E442","#e31a1c","#6a3d9a","#b2df8a","#08519c","#ff7f00","#fdbf6f","#33a02c","#b15928","#f16913","#238b45","#807dba","#d94801","#41ab5d","#fd8d3c","#4292c6")[1:length(dat_list)]
   }

   if(colmix[1]=="none"){
    cat("\t\tusing default color mix\n")
    colmix=rep('black',length(dat_list))
   }

   if(length(colmix)>length(dat_list)){
    warning("\tWARNING\t: colmix is longer than number of elements in list, only first ",length(dat_list)," colors will be used\n")
    colmix=colmix[1:length(dat_list)]
   }


   if(length(colmix)<length(dat_list)){
    warning("\tWARNING\t: colmix is shorter than number of elements in list, default mix will be used instead *(to make all same color, use colmix='none'\n")
      colmix=c("#0072B2","#E69F00","#009E73","#56B4E9","#D55E00","#66A61E","#7570B3","#a50f15","#A6761D","#117733","#332288","#b15928","#882255","#999933","#AA4499","#1f78b4","#F0E442","#e31a1c","#6a3d9a","#b2df8a","#08519c","#ff7f00","#fdbf6f","#33a02c","#b15928","#f16913","#238b45","#807dba","#d94801","#41ab5d","#fd8d3c","#4292c6")[1:length(dat_list)]
   }


   datleg=as.data.frame(names(dat_list))
   colnames(datleg)="name"
   datleg$color=colmix
   datleg$n=NA

     datleg$n[1]=ncol(dat_list[[1]])
    datcol=rep(datleg[1,"color"],datleg$n[1])

    for(idat in 2:length(dat_list)){
       datleg$n[idat]=ncol(dat_list[[idat]])
      datcol=c(datcol,rep(datleg[idat,"color"],datleg$n[idat]))
    }


    pcs=prcomp(t(as.data.frame(dat_list)))
    pcstat=round(as.matrix(summary(pcs)$importance["Proportion of Variance",1:2]),digits=3)*100
#  legend_space=8
  # creates a plot with a wide right margin
  cat('\t calculating PCs, will take time with big datasets\n')
      par(mar=c(5,4,4, legend_space))
    plot(pcs$x[,1:2],col=datcol,pch=pch,frame.plot=F,xlab=paste("PC1 ",pcstat[1],"%",sep=""),ylab=paste("PC2 ",pcstat[2],"%",sep=""),main=paste0(dat_descr,"\n",main),...)

  # create a new plot overlay (with no left margin) with legend on the topright
      par(fig = c(0,1,0,1), oma = c(0, 4, 0, 0), mar = c(0, 4, 0, 0), new = TRUE)
#      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  plot.new()
    legend(x="topright",pch=16,box.lwd = 0,box.col = "white",col=unique(datcol),legend=names(dat_list))

# would be nice to add biplot style arrows for key loadings (ie genes)
  return(invisible(list(legend=datleg,data=dat_list,pcs=pcs$x,pcobj=pcs,cols=datcol)))
}




dat.type<-function(dat_mat){
  if(is.matrix(dat_mat)){dat_class=apply(dat_mat,2,class)}
  if(is.data.frame(dat_mat)){dat_class=unlist(lapply(dat_mat,class))}
   print(matst(dat_class))

}



cov.impact.check<- function(expr_mat,cov_mat,plot_cov_distributions=F,single_gene_analysis=F,sanity=F,n_pcs=5,sd_check=T,dat_descr="",help=F){
if(help==T){ 
  cat("\tINPUTS:\texpr_mat : rows - genes, columns - samples   |   cov_mat rows - samples , columns - covariates\n")
# cat("\tNOTE\t:  cov_mat - non-numeric covariates will be used 'as.factor'\n")
  cat("\tNOTE:\tcov_mat - should contain numeric or factor values only\n")
  cat("\tNOTE:\tsd.check may be required if there are non-varying covariates or genes\n")
  cat("\tUSE:\tsanity=T can be used to correctly order colnames(expr_mat) and rownames(cov)\n")
  }
start_time=Sys.time()

if(!sanity){
  order_check=sum(colnames(expr_mat)!=rownames(cov_mat))
  if(order_check!=0){
    warning("\n\t\tWARNING\t: (colnames(expr_mat)==rownames(cov_mat)) shows :",order_check,"not in same order","\n")
  }
}
if(sanity){

#    matst(colnames(expr_mat)==rownames(cov_mat))
    expr_mat=expr_mat[,intersect(colnames(expr_mat),rownames(cov_mat))]
    cov_mat=cov_mat[intersect(colnames(expr_mat),rownames(cov_mat)),]
  cat("\t\tcolnames(expr_mat) n=",ncol(expr_mat)," rownames(cov_mat) n=",nrow(cov_mat)," same order",sum(colnames(expr_mat)==rownames(cov_mat)),"\n")

    expr_mat=sd.check(expr_mat,T,F)
    cov_mat=sd.check(cov_mat,F,T)
  }

#library(WGCNA)
# sampleTree2 = hclust(dist(expr_mat), method="average")

# traitColors = numbers2colors(cov_mat, signed = FALSE);
# plotDendroAndColors(sampleTree2, traitColors,groupLabels =names(datTraits),main ="Sample dendrogram and trait heatmap")
#print(dim(cov_mat))


#if(sd_check==T){
#  cov_mat=sd.check(cov_mat)
#}
#print(dim(cov_mat))
###------------------------------------------------------------------------------------------------------------
##  PC level analysis
#  cat("\t\tpca analyis of expression profiles\n")
  pcs=pcstat(expr_mat,plot_pcs=c(T,T,F),n_pcs=5,verbose=F)

  univpcP=matrix(NA,nrow=ncol(pcs),ncol=ncol(cov_mat))
    colnames(univpcP)=paste(colnames(cov_mat),"Pval",sep=".")
    rownames(univpcP)=colnames(pcs)

  univpcR=matrix(NA,nrow=ncol(pcs),ncol=ncol(cov_mat))
    colnames(univpcR)=paste(colnames(cov_mat),"Rsq",sep=".")
    rownames(univpcR)=colnames(pcs)



  for(icov in 1:ncol(cov_mat)){
    for(ipcs in 1:ncol(pcs)){
#     cat(icov,ipcs,"\n")
    
      univpcP[colnames(pcs)[ipcs],paste(colnames(cov_mat)[icov],"Pval",sep=".")]=lmp((lm(pcs[,ipcs]~cov_mat[,icov])))
      univpcR[colnames(pcs)[ipcs],paste(colnames(cov_mat)[icov],"Rsq",sep=".")]=summary(lm(pcs[,ipcs]~cov_mat[,icov]))$r.sq
    }
  }




# temp measure to avoid breaking due to 'missing' P-values
  univpcP[is.na(univpcP)]=1
# parmar=c(5.1,4.1,4.1,2.1) 
#  try(
#  Heatmap(univpcP,mode="pval",main=paste(dat_descr,"\n covariate effect on PCs ",nrow(cov_mat)," samples",sep=""))
# )
#  Heatmap(univpcR)
#  par(mar=parmar)
  library(corrplot)

  par(mfrow=c(1,1))
#  col2= colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7","darkgoldenrod1", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))#"#FFFFFF"
#  col2= colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","darkgoldenrod1", "#92C5DE", "#4393C3", "#2166AC", "#053061"))#"#FFFFFF"
#  corrplot(-log10(univpcP),method = "circle", is.corr=F,
#    p.mat = univpcP, insig = "blank",sig.level = 0.01,col=rev(col2(100)),cl.align="l",tl.col="black",
#    mar=c(0,0,4,0),tl.cex=0.7,
#    title="expression PCs correlation with covars")
    try(Heatmap(univpcP,mode="pval",margin=c(10,5),main="expression PCs correlation with covars"))
  par(mar=c(5,4,4,2)+0.1)
# distributions of covariates
  par(mfrow=c(1,1))

if(plot_cov_distributions){
  hist(make.numeric(expr_mat),breaks=100,col="grey",main="expression matrix")

#  hist(pcs[1,breaks=100,col="grey",main="expression matrix")
  par(mfrow=c(3,2))
  for(icov in 1:ncol(cov_mat)){

    hist(make.numeric(cov_mat[icov]),breaks=50,col="grey",main=colnames(cov_mat)[icov])
  }
  par(mfrow=c(1,1))
}
###------------------------------------------------------------------------------------------------------------
##  Single gene analysis
if(!single_gene_analysis){
   return(invisible(list("pcs"=pcs,"univpcP"=univpcP,"univpcR"=univpcR)))
}
if(single_gene_analysis){
cat("\t\tPerforming single gene analysis on",nrow(expr_mat),"genes\n")

  par(mfrow=c(5,2))
    univgenP=list()
    univgenFDR=list()
    univgenR=list()
    for(icov in 1:ncol(cov_mat)){
      cat("\t\t  ----------- ",colnames(cov_mat)[icov],"\t-----------",icov,"of",ncol(cov_mat),"\n")
        genp=list()
        genr=list()
#     for(igen in 1:100){
      for(igen in 1:nrow(expr_mat)){
  #     genp[[paste(rownames(expr_mat)[igen],colnames(cov_mat)[icov],"Pval",sep=".")]]=lmp((lm(pcs[,ipcs]~cov_mat[,icov])))

        test.pair=cbind(t(expr_mat[igen,,drop=F]),cov_mat[,icov,drop=F])
        test.pair=test.pair[complete.cases(test.pair),]
      lmobj=lm(test.pair[,1]~test.pair[,2])
#      lmobj=(lm(t(expr_mat[igen,,drop=F])~.,data=(cov_mat[,icov,drop=F])))       # does not handle NA values
        genp[[rownames(expr_mat)[igen]]]=lmp(lmobj)
  #     genr[[paste(rownames(expr_mat)[igen],colnames(cov_mat)[icov],"Rsq",sep=".")]]=summary(lm(pcs[,ipcs]~cov_mat[,icov]))$r.sq
        genr[[rownames(expr_mat)[igen]]]=summary(lmobj)$r.sq

        cat(round(igen/nrow(expr_mat),digits=2),"\r");flush.console()
      }
      univgenP[[paste(colnames(cov_mat)[icov],"Pval",sep=".")]]=unlist(genp)
      univgenR[[paste(colnames(cov_mat)[icov],"Rsq",sep=".")]]=unlist(genr)
      univgenFDR[[paste(colnames(cov_mat)[icov],"FDR",sep=".")]]=p.adjust(univgenP[[paste(colnames(cov_mat)[icov],"Pval",sep=".")]],method="fdr")

      
      hist((univgenP[[icov]]),main=paste(dat_descr,"n=",nrow(test.pair),"\nlm( gene ~",colnames(cov_mat)[icov],")",sep=" "),xlab="P-value")
  #   hist(-log10(univgenP[[icov]]),main=paste("gene~",colnames(cov_mat)[icov],"-log10(Pval)",sep=" "))
  #   boxplot((univgenP[[icov]]),main=paste("gene~",colnames(cov_mat)[icov],"Pval",sep=" "))
#      boxplot(-log10(univgenP[[icov]]),main=paste(dat_descr,"\nlm( gene ~",colnames(cov_mat)[icov],")",sep=" "),ylab="-log10(P-value)")
      plot.new()
      legend(x="center",box.lwd = 0,box.col = "white",xpd = TRUE,
        legend=c(paste(sum(univgenP[[icov]]<0.05,na.rm=T)," (",round(sum(univgenP[[icov]]<0.05,na.rm=T)/length(univgenP[[icov]]),digits=2)*100,"%)"," genes P<5%",sep="")
              ,paste(sum(univgenFDR[[icov]]<0.05,na.rm=T)," (",round(sum(univgenFDR[[icov]]<0.05,na.rm=T)/length(univgenFDR[[icov]]),digits=2)*100,"%)"," genes FDR<5%",sep="")
              ),cex=1)

  }
  par(mfrow=c(1,1))
#  par(mar = c(5,4,7,2) + 0.1)
  par(mar = c(15,4,4,2))
   boxplot(univgenR,main=paste(dat_descr,"\nlm( gene ~ single-covariate )",sep=" "),ylab="R-squared",ylim=c(0,1),las=2,frame=F,pch=16,cex=0.5)

  return(invisible(list("pcs"=pcs,"univpcP"=univpcP,"univpcR"=univpcR,"univP"=univgenP,"univR"=univgenR,"univFDR"=univgenFDR)))

}
  cat("\n\ttime taken : ",round(Sys.time()-start_time,digits=2),"mins\n\n")
}



Heatm<-function(cor.measures,min=-1,max=1,rowclust=F,colclust=F,ncols=101,dendrogram="none",main="",mode="cor",sig=T,cexrow=0.7,cexcol=0.7,margin=c(12,12)){
  library(gplots)
# print("Heatmap( 'matrix to use' , 'min value for legend cols' , 'max ..' , 'clust by rows', 'clust by cols' , 'ncols to use for legend', 'dendrogram=c('none','row','column','both')' )")
# min=-1
# max=1

# default settings for testing plots
#min=-1;max=1;rowclust=F;colclust=F;ncols=101;dendrogram="none";main="";mode="cor"


if(mode=="cor"){
  print("plotting correlation based matrix")
#  heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=0.7,   cexRow=0.7,   lheight = lheight,symkey=T,main=main,lmat=lmatrix)#,lwid=c(0.5,0.5))
   heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,symkey=T,main=main)
#   heatmap.2((cor.measures),breaks=seq(min,max,length=(ncols+1)),col = colorRampPalette(c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#deebf7","#f7fbff","white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,   symkey=T,main=main,lwid=lwidth,lmat=lmatrix)

#"white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"
}
if(mode=="pval"){
  if(sig==T){
    print("plotting p-value based matrix with cellnotes")
    corm=round(cor.measures,digits=3)
    cor.match=matrix(NA,nrow=nrow(cor.measures),ncol=ncol(cor.measures))
      rownames(cor.match)=rownames(cor.measures)
      colnames(cor.match)=colnames(cor.measures)
    # Head(cor.match)
    for(i in 1:nrow(corm)){
      for(j in 1:ncol(corm)){
    #     cat(i,j,"\n")
        if(corm[i,j]>0.05){
          cor.match[i,j]=""
    #     print(">0.05")
        }
        if(corm[i,j]<=0.01){
          cor.match[i,j]=gsub(" ","",paste(rep("*",min(4,floor(-log10(corm[i,j]))-1)),collapse=" "))    # backup square == .
    #     print("<0.01")
        }
        if(corm[i,j]<0.05 & corm[i,j]>0.01){
          cor.match[i,j]="+"
    #     print("<0.05")
        }
      }
    }
#    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,
    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,symkey=F,main=main)  #cexCol=1,cexRow=0.8,   
#    heatmap.2(-log10(cor.measures),cellnote=(cor.match),notecol="black",breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,   


  }
  if(sig==F){
    print("plotting p-value based matrix, no cellnotes")
#    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#F0E442","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)  #cexCol=1,cexRow=0.8,    
    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,symkey=F,main=main)
#    heatmap.2(-log10(cor.measures),breaks=seq(0,max(-log10(cor.measures)),length=(ncols+1)),col = colorRampPalette(c("white","#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))(ncols),tracecol=F,dendrogram=dendrogram,Rowv=rowclust,Colv=colclust,margins=margin,density.info="none",keysize=1,cexCol=cexcol,cexRow=cexrow,lhei = lheight,symkey=F,main=main,lwid=lwidth,lmat=lmatrix)

  }
}

#min(na.omit(as.vector(cor.measures))),max(na.omit(as.vector(cor.measures)))
}


































Load<-function(obj_path){
#  obj=ls()
#  obj=ls()
#  all_obj=load(obj_path)
 # print(as.matrix(sort(ls()[!(ls() %in% obj)]),decreasing=T))
#  print(parent.frame())
#  if(exists("readme")){cat(readme)}  # needs to be the one read-in
  return(as.matrix(sort(load(obj_path, parent.frame()),decreasing=F)))

#  return(load(obj_path, parent.frame()))
}







mat<-function(dat_mat,decreasing=T){
  return(as.matrix(dat_mat))
}

cmat<-function(dat_mat){
  return(as.matrix(colnames(dat_mat)))
}

rmat<-function(dat_mat,decreasing=T){
  return(as.matrix(rownames(dat_mat)))
}

mats<-function(dat_mat,decreasing=T){
  return(as.matrix(sort(dat_mat,decreasing=decreasing)))
}



matst<-function(dat_mat,sort=T,decreasing=T){
  if(sort){
    dummy=as.data.frame(sort(table(dat_mat),decreasing=decreasing))
  }
  if(!sort){
    dummy=as.data.frame((table(dat_mat)))
  }

  if(ncol(dummy==1)){
  	holder=as.data.frame(rownames(dummy))
  	holder$count=dummy[,1]
  	dummy=holder
  }
  colnames(dummy)=c('entry','count')
  dummy$percent=round(dummy$count/sum(as.numeric(dummy$count)),digits=3)
  return(dummy)

}





# writing lists to file / matrix with different lengths
#lapply(modg, write, "~/Dropbox/Cognition/list.writing.test.txt", append=TRUE, ncolumns=5850)


is.missing<-function(data_mat,make_plot=F,use_grid=F,dat_descr="",define_na=c(NA,NaN)){

# needs optimisation for large datasets (calc missn from the start and use that for cat) <<<<<<<§§§§§§§§§§§§§§§§§§§§§§§
  n_overlap=nrow(data_mat[complete.cases(data_mat),])

#    cat("\t----------------  Total Missing : ",sum(is.na(data_mat))," (",round(sum(is.na(data_mat))/length(is.na(data_mat)),digits=3)*100,"%)  ----------------\n",sep="")
    cat("\t--------------------  Total Missing : ",sum((data_mat%in%define_na))," (",round(sum((data_mat%in%define_na))/length((data_mat%in%define_na)),digits=3)*100,"%)  --------------------\n",sep="")
    cat("\t   ------ ",n_overlap,"samples remaining if 'complete.cases()'  ------\n\n")
#    missn=as.data.frame(mats(apply(data_mat,2,function(x)sum(is.na(x))),T))

    missn=as.data.frame(mats(apply(data_mat,2,function(x)sum(x%in%define_na)),T))
    colnames(missn)="n.missing"
    missn$pc.missing=round(missn$n.missing/nrow(data_mat),digits=3)

  if(sum(missn)<0 & make_plot){cat('\n\t no plot produced : missing detected\n')}
  if(sum(missn)>0){

#merge_dist[merge_dist%in%c(NA,'NA','NaN','Inf','-Inf','')]
#  data_mat[data_mat%in%c("NA","NaN")]=0
#  data_mat[(data_mat%in%define_na)]=1
#  data_mat[!(data_mat%in%define_na)]=0
    dat_mat=apply(data_mat,2,function(x)(x%in%define_na))
      rownames(dat_mat)=rownames(data_mat)
    data_mat=dat_mat*1


#    Head(data_mat)
 #   summary(as.numeric(as.matrix(data_mat)))
#  data_mat[data_mat%in%c("NA","NaN")]=1

#    table(as.numeric(as.matrix(data_mat)))

#Heat((matrix),min=0,max=1,rowclust=T,colclust=T,dendrogram='both')
# intended to remove rows with all missing, but these should be clustered in 1 corner and not cause issues anyhow
#  if(!full.row.missing){
#    apply(data_mat,1,function(x){sum(is.na(x))==length(x)})


#  }

if(dat_descr!=""){
  paste0(dat_descr,"\n")
}

#if(make_plot==T & use_grid==F & max(missn)!=0){
if(make_plot & !use_grid){
    library(gplots)
  heatmap.2(make.numeric(data_mat),breaks=seq(0,1,length=(3)),col = c("#9ebcda","#e6550d"),tracecol=F,dendrogram='both',Rowv=T,Colv=T,margins=c(12,12),density.info="none",keysize=1,cexCol=0.7,cexRow=0.7,symkey=F)
   
   mtext(dat_descr,adj=1,side=3,line=2)
   mtext(paste0("n.complete : ",n_overlap),adj=1,side=3)
  }

#if(make_plot==T & use_grid==T & max(missn)!=0){
if(make_plot & use_grid){
    library(gplots)
  heatmap.2(make.numeric(data_mat),breaks=seq(0,1,length=(3)),col = c("#9ebcda","#e6550d"),tracecol=F,dendrogram='both',Rowv=T,Colv=T,margins=c(12,12),density.info="none",keysize=1,cexCol=0.7,cexRow=0.7,symkey=F,colsep=0:ncol(data_mat),rowsep=0:nrow(data_mat),sepcolor="white",sepwidth=c(0.05,0.05))

   mtext(dat_descr,adj=1,side=3,line=2)
   mtext(paste0("n.complete : ",n_overlap),adj=1,side=3)

  }

    return(missn)
  
#    print(missn)
#  return(invisible(data_mat[complete.cases(data.mat),]))

}
}



cplot<-function(x,y,dat_descr="",legend.pos="topright",...){
#  cat(min(x),max(x),min(y),max(y),"\n")

#  dummy=cbind(seq(floor(min(x)),(ceiling(max(x))+log10(ceiling(max(x)))),length.out=10),seq(floor(min(y)),(ceiling(max(y))+log10(ceiling(max(y)))),length.out=10))
#  print(dummy)
  plot(x=x,y=y,pch=16,frame.plot=FALSE,...
    ,main=paste(dat_descr
                ,"spearman  P =",signif(cor.test(x,y,method="spearman")$p.val,digits=2)
                ,"  R-sq =",round(cor(x,y,method="spearman"),digits=3)
                ,"\nkendall      P =",signif(cor.test(x,y,method="kendall")$p.val,digits=2)
                ,"  R-sq =",round(cor(x,y,method="kendall"),digits=3)
                ,"\nlm           P =",signif(lmp(lm(x~y)),digits=2)
                ,"     R-sq =",signif(summary(lm(x~y))$r.sq,digits=3))
    )
#  points(x,y,pch=16)
  abline(lm(x~y),col="dodgerblue")
  
#  legend(legend.pos,legend=paste(
 #               "spearman P =",signif(cor.test(x,y,method="spearman")$p.val,digits=2)
  #              ,"  R-sq =",round(cor(x,y,method="spearman"),digits=3)
   #             ,"\nkendall      P =",signif(cor.test(x,y,method="kendall")$p.val,digits=2)
    #            ,"  R-sq =",round(cor(x,y,method="kendall"),digits=3)
     #           ,"\nlm             P =",signif(lmp(lm(x~y)),digits=2)
      #          ,"  R-sq =",signif(summary(lm(x~y))$r.sq,digits=3)
#
#                )
 #               ,box.lwd = 0,box.col = "white")
}

eplot<-function(x,xlab="",ylab="",main="",legend.pos="topright"){
#  cat(min(x),max(x),"\n")
  dummy=matrix(
    seq(floor(min(x)),(ceiling(max(x))+log10(ceiling(max(x)))+1),length.out=length(as.numeric(as.matrix(x))))
      ,nrow=nrow(x),ncol=ncol(x))

matplot(dummy,type="n",frame.plot=FALSE,xlab=xlab,ylab=ylab,main=main)
}




clust.analyse<-function(cov_mat,do_plots=c(T,T,T,T),sanity=F,box_sig=0.95,clust.method="ward.D2",cor.method="kendall",dat_descr=""){
    cat("\tINPUTS\t:   cov_mat rows - samples , columns - covariates (numeric or factor)\n")
    cat("\t\t\t - clustering method :",clust.method,"   |   correlation :",cor.method,"   |   n.samples :",nrow(cov_mat),"\n\n")
    # ward clustering : http://stats.stackexchange.com/questions/109949/what-algorithm-does-ward-d-in-hclust-implement-if-it-is-not-wards-criteria



    if(sanity){
      cov_mat=sd.check(cov_mat)
    }
    cor_mat=cor(make.numeric(cov_mat),method=cor.method)

###>> # attempt at introducing p-values of correlation (matrix solution rather than trying to loop thru manually)
###>>
#    library("psych")
###>>
#   cor_test=corr.test(cov_mat,method="kendall")

  if(do_plots[1]){
    plot(hclust(as.dist(1-cor_mat),method=clust.method),main=paste(dat_descr,"\n",clust.method,"clutering on",cor.method,"correlations, n=",nrow(cov_mat)))#\nline - 99th percentile dist"))
#    abline(h=quantile(clust1$height,0.99),col="red")
    plot(hclust(dist(t(cov_mat)),method=clust.method),main=paste(dat_descr,"\n",clust.method,"clutering on distance matrix (not correl), n=",nrow(cov_mat)))#\nline - 99th percentile dist"))
#    abline(h=quantile(clust1$height,0.99),col="red")
  }
  if(do_plots[2]){
   library(pvclust)
   # Ward Hierarchical Clustering with Bootstrapped p values
   cat("\n\tpvclust provides two types of p-values:\n\t\tBP (Bootstrap Probability)\t- computed by normal bootstrap resampling\n\t\tAU (Approximately Unbiased)\t- computed by multiscale bootstrap resampling (better approximation to unbiased p-value than BP)\n\n")
    fit=pvclust((cor_mat), method.hclust=clust.method,nboot=10000)
    plot(fit) # dendogram with p values
   # add rectangles around groups highly supported by the data
#    pvrect(fit, alpha=box_sig)
  }
  if(do_plots[3]){
   library(corrplot)
#   corrplot(colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols))
#  col2=colorRampPalette(c("#0072B2","#56B4E9","white","#F0E442","darkred"))(20)
###>>   corrplot.mixed((cor_test$r)^2,p.mat=cor_test$p,insig="blank",col=colmixrb,cl.align="l",tl.col="black",order = "hclust",main=paste(cor.method,"samples=",nrow(cov_mat),"ordered using hclust"))
  corrplot.mixed(cor_mat,col=colmixrb,cl.align="l",tl.col="black",order = "hclust",main=paste(dat_descr,cor.method,"samples=",nrow(cov_mat),"ordered using hclust"))
#  p.mat = res1[[1]], insig = "blank"
#    title="Maximum % of overlap with p hyper > 0.01"


#   library(PerformanceAnalytics)
#    chart.Correlation(cor_mat, pch=20,main=paste(cor.method,"samples=",nrow(cov_mat)))#bg=c("blue","red","yellow"), 

#  library(psych)
#    pairs.panels(cor_mat) # covarate style plot
  }
  if(do_plots[4]){
   library(gplots) 
#    cor_mat=cor(make.numeric(cov_mat),method=cor.method)
#
#    Heat(cor_mat,rowclust=T,colclust=T,dendrogram="both",main=paste(cor.method,"samples=",nrow(cov_mat)))
      Heat(cor_mat,main=paste(dat_descr,"\n",cor.method,"samples=",nrow(cov_mat)))
#    Heatmap(cor_mat,rowclust=T,colclust=T,dendrogram="both",main=paste(cor.method,"samples=",nrow(cov_mat)))

  }
#  if(do_plots[4]==T){
#library(psych)
#  biplot(pcs)
#  pairs.panels(expr_ma) # covarate style plot
#   }
}


#pclust<-function(cov_mat){
#    cat("\tINPUTS\t:   cov_mat rows - samples , columns - covariates (numeric or factor)\n")

#  library(pvclust)
  # Ward Hierarchical Clustering with Bootstrapped p values
#    fit <- pvclust((cov_mat), method.hclust="ward.D2", method.dist="euclidean")
#    plot(fit) # dendogram with p values
    # add rectangles around groups highly supported by the data
#    pvrect(fit, alpha=.95) 
#}




clustplot<-function(data,clust.method="ward.D2",cor.method="kendall",verbose=T){
if(verbose==T){  
  cat("\tHierarchical clustering on 1) correlation and 2) distance\n")
  cat("\t\tclust.method :",clust.method,"  |   cor.method :",cor.method,"  |   n :",nrow(data),"\n\n")
  }
#  library(WGCNA)
  par(mfrow=c(2,1))
    plot(hclust(as.dist(1-cor(data,method=cor.method)),method=clust.method),main=paste(clust.method,"clutering on",cor.method,"correlations, n=",nrow(data)))#\nline - 99th percentile dist"))
#    abline(h=quantile(clust1$height,0.99),col="red")
    
    plot(hclust(dist(t(data)),method=clust.method),main=paste(clust.method,"clutering on distance matrix, n=",nrow(data)))#\nline - 99th percentile dist"))
#    abline(h=quantile(clust1$height,0.99),col="red")

}






#rmerge<-function(data_mat1,data_mat2,all=T,...){
#  dat_out=merge(data_mat1,data_mat2,by="row.names",...)
#    rownames(dat_out)=dat_out$Row.names
#  dat_out=dat_out[,-which(colnames(dat_out)=="Row.names")]
#  return(dat_out)
#}


rmerge<-function(data_mat1,data_mat2,all=T){
  dat_out=merge(data_mat1,data_mat2,by="row.names",all=all)
    rownames(dat_out)=dat_out$Row.names
  dat_out=dat_out[,-which(colnames(dat_out)=="Row.names")]
  return(dat_out)
}

rmerge.list<-function(dat_lis,all=F){
  lis_merge=rmerge(dat_lis[[1]],dat_lis[[2]])
   cat('\tnrows inputs:',nrow(dat_lis[[1]]),nrow(dat_lis[[1]])==nrow(dat_lis[[2]]))   ##  keeping track of nrows in each list
  if(length(dat_lis)>2){
    for(ilis in 3:length(dat_lis)){
       cat(' ',nrow(dat_lis[[1]])==nrow(dat_lis[[ilis]]))               ##  keeping track of nrows in each list
      lis_merge=rmerge(lis_merge,dat_lis[[ilis]])
    }
  }
   cat('\tnrows output:',nrow(dat_lis[[1]])==nrow(lis_merge),'\n')                    ##  keeping track of nrows final output
  return(lis_merge)
}


univarlm<-function(data_mat,colname_y_var){
 cat("\tUSE:  univariate lm for",ncol(data_mat),"variables, y =",colname_y_var,"\n")
 cat("\t\tNOTE: missing values can be used, only relevant variable n will be affected\n")
 cat("\t\tNOTE: for factor variables, min P is used\n\n")
  vdr=make.numeric(data_mat[,which(colnames(data_mat)==colname_y_var),drop=F])
  data_mat=(data_mat[,-which(colnames(data_mat)==colname_y_var),drop=F])

  unistat=as.data.frame(matrix(NA,nrow=ncol(data_mat),ncol=4))
    colnames(unistat)=c("lmP","lmRsq","n","sig")
    rownames(unistat)=colnames(data_mat)
#     Head(unistat)


  for(icov in 1:ncol(data_mat)){
    tester=row.merge(vdr,data_mat[,colnames(data_mat)[icov],drop=F])
  # tester=merge(vdr,data_mat[,colnames(data_mat)[icov],drop=F],by="row.names")
  # tester=tester[,-which(colnames(tester)=="Row.names")]
    tester=tester[complete.cases(tester),]

    unistat[icov,"n"]=nrow(tester)
    #### using min for factors - 1 p-value per factor ==> only interested in the lowest one
    unistat[icov,"lmP"]=min(summary(lm(tester[,1]~tester[,2]))$coefficients[,"Pr(>|t|)"][-1])  # -1 removes the intercept
#    unistat[icov,"lmRsq"]=summary(lm(tester[,1]~tester[,2]))$r.sq*sign(min(summary(lm(tester[,1]~tester[,2]))$coefficients[,"Estimate"][-1]))
    unistat[icov,"lmRsq"]=summary(lm(tester[,1]~tester[,2]))$r.sq
  }

  unistat[unistat$lmP>0.1,"sig"]=""
  unistat[unistat$lmP<0.1,"sig"]="."
  unistat[unistat$lmP<0.05,"sig"]="+"
  unistat[unistat$lmP<0.01,"sig"]="*"


  return(unistat[order(unistat$lmP),])

}




multivarlm<-function(data_mat,colname_y_var,verbose=T){
if(verbose==T){
 cat("\tUSE:  multivariate lm  :  y =",colname_y_var,"~",ncol(data_mat),"variables, \n\n")
}
  lm_out=(lm(data_mat[,colname_y_var,drop=F]~.,data=as.data.frame(data_mat[,-which(colnames(data_mat)==colname_y_var)])))

if(verbose==T){
  lm_stat=as.data.frame(sort(summary(lm_out)$coefficients[,"Pr(>|t|)"][-1]),decreasing=T)
  colnames(lm_stat)="P"
  lm_stat$FDR=p.adjust(lm_stat$P,method="fdr")
  lm_stat$sig=""
  lm_stat[lm_stat$P<0.1,"sig"]="."
  lm_stat[lm_stat$P<0.05,"sig"]="+"
  lm_stat[lm_stat$P<0.01,"sig"]="*"
  print(lm_stat)
  cat("\n\tR-squared :\t\t",summary(lm_out)$r.sq,"\n\tadjusted R-squared :\t",summary(lm_out)$adj.r.sq,"\n\n")
}
  return((lm(data_mat[,colname_y_var,drop=F]~.,data=as.data.frame(data_mat[,-which(colnames(data_mat)==colname_y_var)]))))
}




meplot<-function(expr_mat,do_plots=T){
  
# expr_mat=bexpr$M
  pca=prcomp(t(expr_mat),scale=T,center=T)
  expr_mat=rmerge(pca$x[,1]*sign(pca$rotation[which.max(abs(pca$rotation[,1])),1]),t(expr_mat),all=F)
  Head(expr_mat)

  expr_mat=(as.data.frame(scale(expr_mat,center=T)))

  eplot(expr_mat,xlab="samples",ylab="gene experssion",main=paste("PC1 explains",round(summary(pca)$importance["Proportion of Variance",1],digits=2)*100,"% variation"))
  for(igen in 2:ncol(expr_mat)){
    lines((expr_mat[,igen]),col=rgb(1, 0, 0, 0.5),)
  }
  
  lines(expr_mat$x,col="midnightblue",lwd=3)
  legend("topright", legend = c("gene expression", "ME 'expression'"), bty = "n",lwd = c(2,4), cex = 1, col = c("darkred", "midnightblue"), lty = c(1, 1), pch = c(NA, NA))

}




net.me<-function(expr_list,nme=1){
  mes=list()
  for(imod in 1:length(expr_list)){
    cat('\t',names(expr_list)[imod])
    pcs=prcomp((expr_list[[names(expr_list)[imod]]]),scale=T,center=T)
    for(imes in 1:nme){
     cat(' ',paste0('ME',imes))
#     mes[[paste0('ME',imes)]][[names(expr_list)[imod]]]=pcs$x[,imes]*sign(pcs$rotation[which.max(abs(pcs$rotation[,imes])),imes])
      mes[[names(expr_list)[imod]]][[paste0('ME',imes)]]=pcs$x[,imes]*sign(pcs$rotation[which.max(abs(pcs$rotation[,imes])),imes])
    }
      mes[[names(expr_list)[imod]]]=as.data.frame(mes[[names(expr_list)[imod]]])
    cat('\n')
  }
  return(invisible(mes))
}


list.as.df<-function(in_list){
  cat("\tUSE: convert 'uneven' list to data frame (missing introduced at the end)\n")
  #'adapted' from :  http://stackoverflow.com/questions/27153979/converting-nested-list-unequal-length-to-data-frame          # N.B. 'adapted' is defined as 'shamelessly ripped off', which includes direct copy/paste of explanation below
  #  get the length of list element ('indx') by looping with sapply
  #  recent version of R - can use lengths to replace the sapply(.., length) step
  #  change the length of each element to the max length from the 'indx' (length<-) and thereby pad NA values at the end of the list elements with length less than the max length
  #  rbind the list elements, convert to data.frame and change the column names.

  indx=sapply(in_list, length)
  res=as.data.frame(do.call(rbind,lapply(in_list, `length<-`,max(indx))))

  return(res)
}


duplicates<-function(dat_vec,...){
  return(dat_vec[duplicated(dat_vec)])
}


get.duplicates<-function(dat_mat,col_dup,...){
  ndup=sort(table(dat_mat[,col_dup]),decreasing=T)
#   Head(ndup)
  ndup=ndup[ndup>1]

  if(length(ndup)<1){
    cat("\tno duplicates found in column :",col_dup,"\n")
  }

  if(length(ndup)>1){
  #   Head(ndup)
    cat("\ttop duplicates\n")
    print(as.matrix(ndup[1:min(15,length(ndup))]))
  
    dat_dup=dat_mat[(dat_mat[,col_dup]%in%names(ndup)),]
    dat_dup=unique(dat_dup[order(dat_dup[,col_dup]),])
  
    dat_unq=unique(dat_mat[!(dat_mat[,col_dup]%in%names(ndup)),])
  
    
  return(invisible(list(n.duplicated=ndup,duplicates=dat_dup,unique=dat_unq)))
  }
}



Boxplot<-function(dat_mat,pch=16,cex=0.5,las=2,frame=F,varwidth=T,...){
  boxplot(dat_mat,pch=pch,cex=cex,las=las,frame=frame,varwidth=varwidth,...)
}



Plot<-function(x,pch=16,cex=0.5,frame.plot=F,...){
  plot(x,pch=pch,cex=cex,frame.plot=frame.plot,...)
}



sva.fac<-function(expr_mat,non_adjust="",adjust="",nsv=""){
  ## installing sva
  #http://www.bioconductor.org/packages/release/bioc/html/sva.html
  #source("https://bioconductor.org/biocLite.R")      ##  try http:// if https:// URLs are not supported
  #biocLite("sva")

  # data for sva tutorial ## https://www.bioconductor.org/packages/3.3/bioc/vignettes/sva/inst/doc/sva.pdf
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("bladderbatch")

  library(sva)
  #library(bladderbatch)
  #data(bladderdata)

  #source("https://bioconductor.org/biocLite.R")
  #biocLite("limma")
  library(limma)

  # install.packages('pamr')
  library(pamr)

  dummy=as.data.frame(rep(1,ncol(expr_mat)))

    if(non_adjust!=""){
      cat('\tusing specified covars to keep :',ncol(adjust),'variables\n')
    mod = model.matrix(~.,data=non_adjust)
    }
    if(non_adjust==""){
      cat('\tno variables specified to keep\n')
    mod = model.matrix(~1,data=dummy)
    }

    if(adjust!=""){
      cat('\tusing specified covars to adjust for :',ncol(adjust),'variables\n')
    mod0 = model.matrix(~.,data=adjust)
    }
    if(adjust==""){
      cat('\tno variables specified to adjust for\n')
    mod0 = model.matrix(~1,data=dummy)
    }

    if(nsv!=""){
      nsv
      cat('\tusing pre-specified n.factors :',nsv,'\n')
    }
    if(nsv==""){
    nsv = num.sv(expr_mat,mod,method="leek")
      cat('\tn factors suggested by sva :',nsv,'\n')
    }
    

    svobj = sva(expr_mat,mod,mod0,n.sv=nsv)

}




## wrapper code to run multlm for multiple MEs 1 at a time.., same applicable for genes 
#mlm=list()
#for(imod in 1:nrow(mes)){
#  mlm[[rownames(mes)[imod]]]=multlm(c1[,'Verbal.delayed.recall',drop=F],cbind(c1[,!(colnames(c1)=='Verbal.delayed.recall')],(mes[imod,])),T)
#}
multlm<-function(y_mat,x_mat,leave1out=F,verbose=F,help=F){
  if(help){
    cat('\tUSE\t: perform multivariate linear model analysis : lm(y_mat[,i]~.,data=x_mat)\n')
    cat('\tNOTE\t: NA acceptable in x_mat and y_mat - for each model complete.cases is used\n')

    cat('\tINPUTS\t: y_mat - y variables, uses 1 column at a time (eg phenotype), rows=samples, cols=variables\n')
    cat('\tINPUTS\t: x_mat - x variables, all are used simultaneously (eg clinical covariates / genes etc), rows=samples, cols=variables\n')
  }
     sample_match=round(sum(rownames(y_mat)==rownames(x_mat))/nrow(x_mat),digits=3)
  if(verbose){
 
    cat('\t\t% sample names matching : ',sample_match,'\n')
  }

  if(sample_match<1){warning(paste0('\t\t% sample names matching : ',sample_match,'\n'))}

  lmstat=list()


  for(ivar in 1:ncol(y_mat)){
  if(verbose){
      cat('\t\t',colnames(y_mat)[ivar],'\t----------------------------\n')
  }
    covars=(cbind(y_mat[,ivar,drop=F],x_mat)) # cbind appears to remove factors => go bk to original input
    covars=covars[complete.cases(covars),]

    lmod=lm(as.matrix(y_mat[rownames(covars),ivar,drop=F])~.,data=x_mat[rownames(covars),])
#   lmod=lm(covars[,1]~.,data=covars[,-1])
# save linear model stuffs  --------------------------------------------------------------------------------
      lmstat[['n']][[colnames(y_mat)[ivar]]]=nrow(covars)
      lmstat[['mlm.var_P']][[colnames(y_mat)[ivar]]]=(summary(lmod)$coefficients[,"Pr(>|t|)"][-1])
      lmstat[['mlm.var_Est']][[colnames(y_mat)[ivar]]]=summary(lmod)$coefficients[,"Estimate"][-1]
      lmstat[['mlm.Rsq']][[colnames(y_mat)[ivar]]]=summary(lmod)$r.sq
      lmstat[['mlm.adjRsq']][[colnames(y_mat)[ivar]]]=summary(lmod)$adj.r.sq
      lmstat[['mlm.Pval']][[colnames(y_mat)[ivar]]]=lmp(lmod)

#      lmstat[['mlm.FDR']]=p.adjust(lmstat[['mlm.Pval']],method='fdr')

  if(leave1out){
# rough estimate of the magnitude of effect of covariate on whole model using leave1out principle  ---------
  if(verbose){
      cat('\t\tperform leave-one-out for each column in x_mat\n')
  }
    leave1=list()
    for(icov in 1:ncol(x_mat)){
      lmod=lm(as.matrix(y_mat[rownames(covars),ivar,drop=F])~.,data=x_mat[rownames(covars),-icov])
        leave1[['l1o.Rsq']][[colnames(x_mat)[icov]]]=summary(lmod)$r.sq
        leave1[['l1o.adjRsq']][[colnames(x_mat)[icov]]]=summary(lmod)$adj.r.sq
    }
  
  leave1[['l1o.adjRsq']][leave1[['l1o.adjRsq']]<0]=0        #  negative adj.R.sq ~ 0
    lmstat[['l1o.Rsq']][[colnames(y_mat)[ivar]]]=lmstat[['mlm.Rsq']]-leave1[['l1o.Rsq']]    #  calculate the difference in Rsquared between full model and leave-1-out
    lmstat[['l1o.adjRsq']][[colnames(y_mat)[ivar]]]=lmstat[['mlm.adjRsq']]-leave1[['l1o.adjRsq']]   #  calculate the difference in Rsquared between full model and leave-1-out
  }
  }
  if(verbose){
  cat('\n\toutput contains :
  \t\t1. n           - number of samples in model
  \t\t2. mlm.var_P   - significance of contribution of single variable to model
  \t\t3. mlm.var_Est - estimate of single variable (coeff)
  \t\t4. mlm.Rsq     - overall model fit R-squared
  \t\t5. mlm.adjRsq  - overall model fit adjusted R-squared : (1 - (1 - R_squared) * ((n - df.int)/error.in.df)) corresponds Wherry Formula-1
  \t\t\tmore info : http://stats.stackexchange.com/questions/48703/what-is-the-adjusted-r-squared-formula-in-lm-in-r-and-how-should-it-be-interpret
  \t\t6. mlm.Pval    - overall model significance P-value
  \t\t7. l1o.Rsq     - if leave1out=T, rough estimate of the magnitude of effect of covariate on whole model using leave1out principle on R-squared
  \t\t8. l1o.adjRsq  - if leave1out=T, rough estimate of the magnitude of effect of covariate on whole model using leave1out principle on adjusted R-squared
    ')
  }
  return(lmstat)

}




univlm<-function(ymat,xmat){
# univariate lm for each column in ymat against each column in xmat, i.e. ymat[,1]~xmat[,1]
 lisp=list()
 lisr=list()
 lisd=list()

 lmstat=list()
 lmpval=list()
  for(vary in 1:ncol(ymat)){
    for(varx in 1:ncol(xmat)){
      lmod=lm((ymat[,vary])~(as.matrix(xmat[,varx])))
      lisp[[colnames(xmat)[varx]]]=lmp(lmod)
      lisr[[colnames(xmat)[varx]]]=summary(lmod)$adj.r.squared
      ##  extract the direction of corrrelation for the variable with lowest p-val (only an issue for factors)
      coeff=summary(lmod)$coefficients[-1,,drop=F]    ##  disregard the intercept line
      lisd[[colnames(xmat)[varx]]]=sign(coeff[coeff[,'Pr(>|t|)']==min(coeff[,'Pr(>|t|)']),'Estimate'])

    }
    lmstat[[colnames(ymat)[vary]]]=as.data.frame(list(p=unlist(lisp),fdr=p.adjust(unlist(lisp),method='fdr'),r.sq=round(unlist(lisr),digits=3),dir=unlist(lisd)))
    lmstat[[colnames(ymat)[vary]]]=psig(lmstat[[colnames(ymat)[vary]]],'p')

    lmpval[[colnames(ymat)[vary]]]=unlist(lisp)     ##  matrix of just p-values for plotting etc
   }
  return(list(lmstat=lmstat,lmpval=as.data.frame(lmpval)))
}



## lm for all genes explained by cov_mat
#multlm<-function(dat_mat,cov_mat,leave1out=F,verbose=T){
#  if(verbose){
#    cat('\tINPUTS\t: dat_mat - y variables, used 1 at a time, rows=variables, cols=samples')
#    cat('\tINPUTS\t: cov_mat - x variables, all are used simultaneously, rows=samples, cols=variables')

#    sample_match=round(sum(colnames(dat_mat)==rownames(cov_mat))/ncol(dat_mat),digits=3)
#    cat('\t\t% sample names matching : ',sample_match,'\n')
#  }
#  if(sample_match<1){warning(paste0('\t\t% sample names matching : ',sample_match,'\n'))}

#  lmstat=list()

#  for(ivar in 1:nrow(dat_mat)){
#    covars=(cbind(t(dat_mat[ivar,,drop=F]),c1))
#    covars=covars[complete.cases(covars),]
#    lmod=lm(t(dat_mat[ivar,rownames(covars),drop=F])~.,data=as.data.frame(make.numeric(cov_mat[rownames(covars),])))
#   lmod=lm(covars[,1]~.,data=covars[,-1])
# parse to save linear model stuffs
#      lmstat[['n']][[rownames(dat_mat)[ivar]]]=nrow(covars)
#      lmstat[['mlm.var_P']][[rownames(dat_mat)[ivar]]]=(summary(lmod)$coefficients[,"Pr(>|t|)"][-1])
#      lmstat[['mlm.var_Est']][[rownames(dat_mat)[ivar]]]=summary(lmod)$coefficients[,"Estimate"][-1]
#      lmstat[['mlm.Rsq']][[rownames(dat_mat)[ivar]]]=summary(lmod)$r.sq
#      lmstat[['mlm.adjRsq']][[rownames(dat_mat)[ivar]]]=summary(lmod)$adj.r.sq
#      lmstat[['mlm.Pval']][[rownames(dat_mat)[ivar]]]=lmp(lmod)

#  if(leave1out){
#    for(icov in 1:ncol(cov_mat){

#    }
#  }

#  lmstat[['mlm.FDR']]=p.adjust(lmstat[['mlm.Pval']],method='fdr')

#  }
#}





## alternative way to run the two below for ppi only within the query
ppi<-function(query,dtb_loc='~/Dropbox/PROJ/ppi/dtb/processed/allDTB_collate_macro.GeneMania.Hippie.iRefWeb.Rdata'){
### NOTE : loaded dtb name is 'macro', first interactor column name='a', second interactor column name='b'
  cat('\n\tloading database file\n')
Load(dtb_loc)
  neta=macro[(macro$a%in%query & macro$b%in%query),]
#   Head(neta)
    cat('\t number of unique genes interacting :',length(unique(c(neta$a,neta$b))),pc(length(unique(c(neta$a,neta$b))),length(query)),'\n')

    cat('\tperforming filtering step, takes a long time based on number of interactors\n')
#  neta=ppi.dbfilt(neta)
  return(invisible(neta)) 
}



ppi.net<-function(dtb,query,genex,geney,degree=1){

  ppi=list()

  for(icon in 1:degree){
#   print(icon)
    
    if(icon==1){ppi[[paste('d',icon,sep='')]]=dtb[dtb[,which(colnames(dtb)%in%genex)]%in%query | dtb[,which(colnames(dtb)%in%geney)]%in%query,]}
    
    if(icon>1){
      query=unique(c(ppi[[paste('d',icon-1,sep='')]][,which(colnames(dtb)%in%c(genex))],ppi[[paste('d',icon-1,sep='')]][,which(colnames(dtb)%in%c(geney))]))
      ppi[[paste('d',icon,sep='')]]=dtb[dtb[,which(colnames(dtb)%in%genex)]%in%query | dtb[,which(colnames(dtb)%in%geney)]%in%query,]
    }

    cat('\tdegree =',icon,'  unique genes =',length(unique(c(ppi[[paste('d',icon,sep='')]][,which(colnames(dtb)%in%c(genex))],ppi[[paste('d',icon,sep='')]][,which(colnames(dtb)%in%c(geney))]))),'connections=',nrow(ppi[[paste('d',icon,sep='')]]),'\n')
  }
  return(invisible(ppi))
}



ppi.dbfilt<-function(dtb){
t1=Sys.time()
  dumpty=list() 
  n.entries=nrow(dtb)
  dtb$x=paste(dtb$a,dtb$b)
  dtb$y=paste(dtb$b,dtb$a)
while(nrow(dtb)>0){
  humpty=dtb[which(dtb$x==dtb$x[1] | dtb$x==dtb$y[1] | dtb$y==dtb$y[1] | dtb$y==dtb$y[1]), ]
  dtb=dtb[-which(dtb$x==dtb$x[1] | dtb$x==dtb$y[1] | dtb$y==dtb$y[1] | dtb$y==dtb$y[1]), ]

  dumpty[[paste(sort(unique(c(humpty$a,humpty$b))),collapse='_')]]=t(as.data.frame(list(
      a=sort(unique(c(humpty$a,humpty$b)))[1]
      ,b=sort(unique(c(humpty$a,humpty$b)))[2]
      ,weight=paste(sort(unique(unlist(strsplit(as.character(humpty$weight),'\\|')))),collapse=";")
      ,method=paste(sort(unique(unlist(strsplit(as.character(humpty$method),'\\|')))),collapse=";")
      ,type=paste(sort(unique(unlist(strsplit(as.character(humpty$type),'\\|')))),collapse=";")
      ,dtb=paste(sort(unique(unlist(strsplit(as.character(humpty$dtb),'\\|')))),collapse=";")
      ,ref=paste(sort(unique(unlist(strsplit(as.character(humpty$ref),'\\|')))),collapse=";")
      ,pmid=paste(sort(unique(unlist(strsplit(as.character(humpty$pmid),'\\|')))),collapse=";")
      ,mark=paste(sort(unique(unlist(strsplit(as.character(humpty$mark),'\\|')))),collapse=";")
  )))
  cat(1-round(nrow(dtb)/n.entries,digits=2),'\r');flush.console()
}

  print(Sys.time()-t1)
  return(invisible(as.data.frame(t(as.data.frame(dumpty)))))
}



deconv<-function(exp_mat,do_plots=F,verbose=F,...){
	set.seed(0)	##  required for reproducibility -> in a single + limited test the gene numbers in upper and lower dont change, but not fully tested
if(verbose){
  cat('\tUSE : applies deconvolution method to separate bimodal distirbution into two, returns the greater normal\n')
  cat('\tNOTE: expr_mat should be the log transformed matrix of gene expression\n')
}

  library(mixtools)

  M=apply(exp_mat,1,mean)
  S=apply(exp_mat,1,sd)
  yy=normalmixEM(M,mu=quantile(M,c(0.25,0.75),lambda=c(0.5,0.5)))
  Mth=qnorm(0.95,min(yy$mu),yy$sigma[which.min(yy$mu)])

  upper=exp_mat[M>Mth,,drop=F]
  lower=exp_mat[M<Mth,,drop=F]

  if(do_plots){
    hist(exp_mat,col='darkgrey',breaks=100,prob=T,...)
#    hist(lower,col=rgb(1,0,0,0.5),breaks=max(30,100*(nrow(lower)/nrow(exp_mat))),add=T)
#    hist(upper,col=rgb(0,0,1,0.5),breaks=max(30,100*(nrow(upper)/nrow(exp_mat))),add=T)
	lines(density(M, na.rm=T), lty=2, lwd=2)
	par(new=T)
	lines(density(upper, na.rm=T), lty=1, lwd=2,col='dodgerblue')
	par(new=T)
	lines(density(lower, na.rm=T), lty=1, lwd=2,col='darkred')
#    legend(x="topright",pch=15,box.lwd = 0,box.col = "white",col=c('darkgrey',rgb(1,0,0,0.5),rgb(0,0,1,0.5)),pt.bg=c('darkgrey',rgb(1,0,0,0.5),rgb(0,0,1,0.5)),legend=c(paste0('n=',nrow(exp_mat)),paste0('n=',nrow(lower)),paste0('n=',nrow(upper))),cex=1)
    legend(x="topright",pch=15,box.lwd = 0,box.col = "white",col=c('darkgrey','dodgerblue','darkred'),pt.bg=c('darkgrey',rgb(1,0,0,0.5),rgb(0,0,1,0.5)),legend=c(paste0('density n=',nrow(exp_mat)),paste0('lower n=',nrow(lower)),paste0('upper n=',nrow(upper))),cex=1)
    cat('\tdataset n=',nrow(exp_mat),'lower n=',nrow(lower),'upper n=',nrow(upper),'\n')
  }
    
  return(invisible(list(upper=upper,lower=lower)))
}


#cat('\tstuffs\n')


overlap <- function(A,B){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B
    return(table(inA,inB))
}

#set.seed(1)
#A <- sample(letters[1:20],10,replace=TRUE)
#B <- sample(letters[1:20],10,replace=TRUE)
#xtab_set(A,B)

rowstat<-function(data_mat,add=F,round_digits=2,diag_na=T){
  cat('\n\tcalculate mean and median for each row, na.rm=T\n\n')
# round(cbind(apply(pcmf,1,mean,na.rm=T),apply(pcmf,1,median,na.rm=T)),digits=2)
  data_mat_calc=data_mat
  if(diag_na){diag(data_mat_calc)=NA}
  if(add){
    data_mat$mean=round(apply(data_mat_calc,1,mean,na.rm=T),digits=round_digits)
    data_mat$median=round(apply(data_mat_calc,1,median,na.rm=T),digits=round_digits)
    return(data_mat)
  }
  if(!add){
    dat_stat=round(cbind(apply(data_mat_calc,1,mean,na.rm=T),apply(data_mat_calc,1,median,na.rm=T)),digits=round_digits)
      colnames(dat_stat)=c('mean','median')
    return(dat_stat)
  }
}



clust<-function(dat_mat,clust_method='ward.D2',k=1,cor_method='dist',do_plots=F,plot_cex=0.8,dat_descr='',help=F){
 if(help){
  cat('\n\tINPUT :\tdat_mat - matrix, rows=genes, cols=samples\n')
  cat('\tNOTE  :\tcor_method - correlation method to calculate distance see "cor" for options, option: "dist" - distance on raw data, no correlation calculated\n\n')
  cat('\tNOTE  :\tk- specify number of clusters for cutree, can use range i.e. 2:5, if k="dynamic", WGCNA function "cutreeHybrid" is used to cut the tree\n\n')
 }

    if(do_plots){
    	library(WGCNA)	##  moved it here to avoid waiting for clustering only to find the plots failed 
    }

  if(cor_method=='dist'){
   cat('\t- calculating eucledian distance\n')
        distmat=dist(t(dat_mat))
  }

  if(cor_method!='dist'){
   cat('\t- calculating distance based on',cor_method,'correlation\n')
        distmat=as.dist(1-abs(t(cor(dat_mat,method=cor_method))))
    }

    clustat=hclust(distmat,method=clust_method)

  if(k!='dynamic'){
    trestat=cutree(clustat,k=k)
    main_text=paste(dat_descr,'\ndistance=',cor_method,'cluster.method=',clust_method,'k =',paste(k,collapse=', '))
  }
  if(k=='dynamic'){
    trestat=cutreeHybrid(clustat,as.matrix(distmat),minClusterSize=2,deepSplit=0)$labels
    main_text=paste(dat_descr,'\ndistance=',cor_method,'cluster.method=',clust_method,'clusters determined by WGCNA: "cutreeHybrid"')
  }
#   cutree$labels

  clustnm=(unique(trestat))
  clust=list()
  for(iclust in 1:length(clustnm)){
    clust[[as.character(clustnm)[iclust]]]=names(trestat)[trestat==clustnm[iclust]]
  }
    if(do_plots){
		plotDendroAndColors(
		clustat
		,trestat
		#,cutHeight=300
		,hang=0.03
		#,addGuide=TRUE
		,guideHang=0.05,
		,main=main_text
		,cex.colorLabels=plot_cex
		,cex.dendroLabels=plot_cex
		,cex.rowText=plot_cex
	)
	}

    readme='\n\toutput contains :
            \t1. clust - members of modules based on k
            \t2. distmat - distance matrix used to perfom hierarchical clustering
            \t3. clustat - ouput of hcust(distmat)
            \t4. trestat - output of cutree(clustat)
            \n'
#  cat(readme)
  return(invisible(list(clust=clust,distmat=distmat,clustat=clustat,trestat=trestat)))
}





clust.list<-function(dat_lis,clust_method='ward.D2',k=1,cor_method='dist',do_plots=F,do_indiv_plots=F,plot_cex=0.8,dat_descr='',help=F){
 if(help){
  cat('\tINPUT :\tdat_lis - expect list() of matrixes, rows=genes, cols=samples')
  cat('\tNOTE  :\tcor_method - correlation method to calculate distance see "cor" for options, option: "dist" - distance on raw data, no correlation calculated\n\n')
  cat('\tNOTE  :\tk- specify number of clusters for cutree, can use range i.e. 2:5, if k="dynamic", WGCNA function "cutreeHybrid" is used to cut the tree\n\n')
 }

distmat=list()
clustat=list()
trestat=list()
  for(ilis in 1:length(dat_lis)){
   cat('\t',names(dat_lis)[ilis],'\t',ilis,'of',length(dat_lis))
    dummy=clust(dat_lis[[names(dat_lis)[ilis]]],clust_method=clust_method,k=k,cor_method=cor_method,do_plots=do_indiv_plots,dat_descr=paste(names(dat_lis)[ilis],dat_descr))

    distmat[[names(dat_lis)[ilis]]]=dummy$distmat
    clustat[[names(dat_lis)[ilis]]]=dummy$clustat
    trestat[[names(dat_lis)[ilis]]]=dummy$trestat

  }

options(warn=-1)

  trem=as.data.frame(trestat[[1]])
  for(ireg in 2:length(trestat)){
      trem=rmerge(trem,as.data.frame(trestat[[ireg]]),all=T)
  }
  colnames(trem)=names(trestat)


  trem$mean=apply(trem,1,mean,na.rm=T)
  trem$median=apply(trem,1,median,na.rm=T)

  trem[is.na(trem)]=0
#  if(do_plots){
    Heat(make.numeric(as.matrix(trem)),rowclust=F)
#  }

options(warn=0)
  return(invisible(list(distmat=distmat,clustat=clustat,trestat=trestat,trem=trem)))
}













clicky.run<-function(module_genes_list,module_bkrnd,clicky_dir,dat_descr='',id_type='hsapiens__gene_symbol',os_type='osx'){
##  os_type='osx' ; os_type='linux'
##  hsapiens__ensembl_gene_stable_id
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

	if(os_type=='osx'){
	  system(paste0("python ",clicky_dir,"/bin/runner_modified_commandline_osx.py "
	                ,working_dir,"tmp_genelist.txt "
	                ,working_dir,"tmp_bkgrnd.txt "
	                ,working_dir," "
	                ,paste0(names(module_genes_list)[imod],'_',dat_descr)
	                ," hsapiens ",id_type," ",id_type," 0.05"))

	}

	if(os_type=='linux'){
		  system(paste0("python ",clicky_dir,"/bin/clicky_python_commandline_linux.py "
		                ,working_dir,"tmp_genelist.txt "
		                ,working_dir,"tmp_bkgrnd.txt "
		                ,working_dir," "
		                ,paste0(names(module_genes_list)[imod],'_',dat_descr)
		                ," hsapiens ",id_type," ",id_type," 0.05"))

		}

         
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



Heatp<-function(pval_mat,sig=T){
  library(gplots)
#  if(breaks=="default"){breaks=seq(0,max(-log10(pval_mat)),length=(102))}
#  if(col=="default"){col=colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(length(breaks)-1)}

margin=c(max(13,(125-nrow(pval_mat))),12)
lheight=c(0.06,0.94)

  if(sig){
      dat_sig=as.data.frame(pval_mat)
      dummy=pval_mat<=1e-5
      dat_sig[dummy]="****"
      dummy=pval_mat<=1e-4 & pval_mat>1e-5
      dat_sig[dummy]="***"
      dummy=pval_mat<=1e-3 & pval_mat>1e-4
      dat_sig[dummy]="**"
      dummy=pval_mat<=1e-2 & pval_mat>1e-3
      dat_sig[dummy]="*"
      dummy=pval_mat<=0.05 & pval_mat>1e-2
      dat_sig[dummy]="+"
      dummy=pval_mat>0.05
      dat_sig[dummy]=""
    heatmap.2(-log10(pval_mat),margins=margin,lhei=lheight,cellnote=(dat_sig),col=colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(101),breaks=seq(0,max(-log10(pval_mat)),length=(102)),notecol="black",tracecol=F,dendrogram="none",Rowv=F,Colv=F,density.info="none",keysize=1,cexCol=0.7,cexRow=0.7,symkey=F,lwid=c(0.7,0.3),lmat=rbind(c(4,3),c(1,2)))
  }
 if(!sig){
    heatmap.2(-log10(pval_mat),margins=margin,lhei=lheight,col=colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(101),breaks=seq(0,max(-log10(pval_mat)),length=(102)),tracecol=F,dendrogram="none",Rowv=F,Colv=F,density.info="none",keysize=1,cexCol=0.7,cexRow=0.7,symkey=F,lwid=c(0.7,0.3),lmat=rbind(c(4,3),c(1,2)))
  }

}





gestaltheat<-function(dat_mat,descr,multi_page=F){
  if(multi_page){
    if(nrow(dat_mat>100)){
      k=0
      remainder=(nrow(dat_mat)-(k+100))
      cat("\tmulti-page plot :\n")
      while(remainder>0){
        cat("\t",k+1,k+100,"\n")
        dat_plot=dat_mat[(k+1):(k+100),]
            k=k+100
            remainder=(nrow(dat_mat)-(k+100))

            Heatp(dat_plot,sig=T)
            par(xpd = NA)
            mtext(descr, adj = 1, side = 3)
            par(xpd = F)

        }
      if(remainder<0){
        cat("\t",k+1,k+(100+remainder),"\n")
        dat_plot=dat_mat[(k+1):(k+(100+remainder)),]
        
        Heatp(dat_plot,sig=T)
          par(xpd = NA)
          mtext(descr, adj = 1, side = 3)
          par(xpd = F)
      }
    }
  }
}




clicky.read<-function(enrich_type=c('KEGG','WIKI','Commons','GO'),mod_names="",p_threshold=0.1,in_path=getwd(),out_path="",dat_descr="",verbose=T){
if(verbose==T){
 cat("\tUSE\t: read and process files saved from WebGestalt from specified folder, default: in_path=getwd()\n")
# cat("\tNOTE\t: if out_path is provided, the new tables will also be saved at the provided location\n")
 cat("\tINPUTS\t: enrich_type currently supported  : 'KEGG', 'WIKI', 'Commons', 'GO'\n")
 cat("\tOUTPUT\t: list of module enrichments for use with gestalt_plot() function\n\n")
}

  dtb=list()
  types="na"

mod_names=unique(mod_names)
###-------------------
# get file names matching input criteria
  inf=list.files(in_path)#,pattern=paste(enrich_type,"&",mod_names,collapse=""))
  inf=inf[grepl(paste("(?=.*",paste(mod_names,collapse="|"),")",sep=""),inf,perl=T)]
  inf=inf[grepl(paste("(?=.*",paste(enrich_type,collapse="|"),")",sep=""),inf,perl=T)]
  if(length(inf)==0){
    stop(paste("\tERROR\t: no files matching criteria",enrich_type,mod_names,"found at",in_path,"\n"))
  }

  cat("\t\t",length(inf),"files match criteria :\n")
  enrich=list()
    enrich_type=gsub("GO","GO:",enrich_type)

  for(ifle in 1:length(inf)){
    path=list()
#ifle=1
    cat("\t\t  ",(inf)[ifle])
    mod=readLines(paste(in_path,inf[ifle],sep=""))
    mod=mod[11:length(mod)]
    mod=mod[mod!=""]
  #
    pos=grep(paste("(?=.*",paste(enrich_type,collapse="|"),")",sep=""),mod,perl=T)



    if(length(pos)==0){
    cat("\tno significant terms\n")
    path=list(type="na",term="na")
    path$adjPval=1
#   if(out_path!=""){write.table(paste("No Significant results for",enrich_type,"at 0.1 FDR"),paste(out_path,"transform.",inf[imod],".txt",sep=""),sep="\t",row.names=F,quote=F)}
    }



    if(length(pos)>0){
    path=list(type=paste(gsub("\t.*","",mod[pos])))
    path$term=as.vector(t(as.data.frame(strsplit(mod[pos],"\t"))[2,]))
      cat("\tterms found:",length(path$term),"\n")
    path$nGenes=gsub(".*O=|;.*","",mod[pos+1])
    path$rawPval=gsub(".*rawP=|;.*","",mod[pos+1])
    path$adjPval=gsub(".*adjP=","",mod[pos+1])

    endc=c((pos[-1]-1),length(mod))
    ensg_all=list()
    name_all=list()

# add in gene names using overly elaborate method
    for(ipat in 1:length(pos)){
      gen=(mod[(pos[ipat]+2):(endc[ipat])])
      ensg=rep(NA,length(gen))
      name=ensg
      for(igen in 1:length(ensg)){
        ensg[igen]=strsplit(gen[igen],"\t")[[1]][1]
        name[igen]=strsplit(gen[igen],"\t")[[1]][3]
      }
      ensg_all[[ipat]]=paste(ensg,collapse=" ")
      name_all[[ipat]]=paste(name,collapse=" ")
    }
    path$gene.id=unlist(ensg_all)
    path$gene.name=unlist(name_all)
    }


#   for(ityp in 1:length(table(path$types))){
#     dtb[[table(path$types)[ityp]]][[inf]]
#   }

    for(imod in 1:length(mod_names)){
      if((grepl(mod_names[imod],inf[ifle]))){
      modnm=mod_names[imod]
#     print(modnm)
      }
    }


  types=c(types,path$type)
    path=as.data.frame(path)
    path$module=modnm
  enrich[[inf[ifle]]]=as.data.frame(path)

  }


###-------------------
# sort the data for plotting
  types=types[types!="na"]
  types=names(table(types))
  path=list()

  for(ityp in 1:length(types)){
      path[[types[ityp]]]=vector("list",length(mod_names))
      names(path[[types[ityp]]])=mod_names
  }

  for(ifle in 1:length(enrich)){

#   print(names(enrich)[ifle])
    enri=enrich[[ifle]]
    etyp=names(table(enri$type))
    for(ityp in 1:length(etyp)){
      if(etyp[ityp]%in%types){
#     cat(ifle,ityp,"<<\n")
#     print(etyp[ityp])
#     print(unique(enrich[[ifle]]$module))
      enrtyp=enri[enri$type==etyp[ityp],]
        rownames(enrtyp)=enrtyp$term
      enrtyp=enrtyp[,-(which(!(colnames(enrtyp)%in%"adjPval"))),drop=F]
        colnames(enrtyp)=unique(enrich[[ifle]]$module)

#     if(nrow(enrtyp)==0){
#       print("no enrichment terms")
#       enrtyp=as.data.frame(1)
#       rownames(enrtyp)="na"
#       colnames(enrtyp)=unique(enrich[[ifle]]$module)
#     }
      path[[etyp[ityp]]][[unique(enrich[[ifle]]$module)]]=enrtyp
    }
    }
  }

  names(path)=gsub("[ ]","_",names(path))

  for(ityp in 1:length((path))){
    for(imod in 1:length(path[[ityp]])){
      if(class(path[[ityp]][[imod]])=="NULL"){
        dummy=as.data.frame(NA)
        rownames(dummy)="no_significant_enrichment"
        colnames(dummy)=names(path[[ityp]][imod])
#       cat(colnames(dummy),"erk\n")
        path[[ityp]][[imod]]=dummy
      }
    }
  }

  return(path)
}







clicky.plot<-function(enrich_list,dat_descr="",p_threshold=0.1,height=22,width=11,margins=c(5, 20, 4, 2),points_sig=T){

for(ityp in 1:length(enrich_list)){
    cat("\t",names(enrich_list)[ityp],"\n")

  for(imod in 1:length(enrich_list[[ityp]])){
    print(imod)
      cat("\t\t",names(enrich_list[[ityp]])[imod],"\n")
      print(nrow(enrich_list[[ityp]][[imod]]))

    if(imod==1){enrich=as.data.frame(enrich_list[[ityp]][[imod]])}
    if(imod>1){enrich=rmerge(enrich,make.numeric(enrich_list[[ityp]][[imod]]))} 
  }
  enrich=make.numeric(enrich)

  if("no_significant_enrichment"%in%rownames(enrich)){
    enrich=(enrich[-which(rownames(enrich)=="no_significant_enrichment"),,drop=F])
  }
  enrich[is.na(enrich)]=1
  enrich[enrich>p_threshold]=1
  enrich=(make.numeric(enrich))
  dim(enrich)
# remove rows where all values == 1   ||    enrich[-0,] seems to break, hence the logic check first
  if(sum(apply(enrich,1,sum)==ncol(enrich))>0){
    enrich=(enrich[!(apply(enrich,1,sum)==ncol(enrich)),,drop=F])
  }
# print(str(enrich))
# cat("\t\tdim enrich :",str(enrich),"\n")
# Head(enrich)
  print(head(enrich))
  print(dim(enrich))


# pdf(paste("~/Dropbox/bin/clicky/graphics/dummy.",names(enrich_list)[ityp],".pdf",sep=""),width=width,height=height)
  if(nrow(enrich)==0){
    print("no plot")
    plot.new()
    legend("top",legend=paste(dat_descr,names(enrich_list)[ityp],colnames(enrich),"\n","NO SIGNIFICANT TERMS at P<",p_threshold,"\n"),cex=2)

  }

  if(nrow(enrich)>0){

  if(ncol(enrich)==1){
    print("single plot")
    print(dim(enrich))
    enrich=-log10(make.numeric(enrich))
    enrich=enrich[order(enrich[,1]),,drop=F]
    margins=c(100-(nrow(enrich)*1),(max(nchar(rownames(enrich)))/1.8),2.5,1)
    par(mar=margins)
    barplot(t(enrich),horiz=T,las=1,main=paste(names(enrich_list)[ityp],"\n",dat_descr,colnames(enrich)))
    abline(v=-log10(0.01),col="dodgerblue")
    abline(v=-log10(0.05),col="red")
    margins=c(5.1,4.1,4.1,2.1)  # reset default margins
  }

  if(ncol(enrich)>1){
    print("multi plot")
    library(corrplot)
    enrich=enrich[do.call(order, (lapply(1:NCOL(enrich), function(i) enrich[, i]))), ]
    gestaltheat(enrich,names(enrich_list)[ityp],T)
    enrich=-log10(make.numeric(enrich))
  }
# dev.off()
  }
# return(invisible(enrich))
  rm(enrich)
}

}





FUNDNMmap<-function(clusters_list,runname="",inpath="~/Dropbox/SHARED/tools/FUNDNMmap/",outpath=getwd(),selection=F){
  cat('\tNOTE:\tclusters_list only ENSG ids currently supported')
library(MetaDE)
library('parallel')

  Load(paste(inpath,"/PathoGeneENS.Rdata",sep="")) # PathoGene & Patho_ENSgeneID
  Load(paste(inpath,"/ctr_GeneENS.Rdata",sep="")) # ctrGene & ctr_ENSgeneID
  map <- read.table(file = paste(inpath,'/functional_mut_rate.bias_corrected.local.canonical_tx_only.bed.txt',sep=""),
                   header = TRUE, sep = '\t', blank.lines.skip = TRUE)

  ### create a matrix for results
  FBET <- matrix(nrow = length(clusters_list), ncol = 23)
  row.names(FBET) <- names (clusters_list)
  colnames(FBET) <- c("module size","patho","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                        "module DNMs in patients","module DNMs in controls",
                        "non module DNMs in patients","non module DNMs in controls",
                        "gene names of modules DNMs in patients",
                        "gene names of modules DNMs in controls",
                      "BET p.value","BET FDR","Theo Ps (mutation rate of map.M)","Estimated Ps","Ratio Obs/Exp",
                        "[95% EsPs CI] inf","[95% EsPs CI] sup",
                        "x = nb of DNM in map.M (success)","n = nb of DNM in map/mutation rate all genes of map (trials)",
                        #"ENSgeneID of M not in map",
                        "nb ENSgeneID of M not in map"
                      )
  
  npDNM<- list()
  for (pat in 1:length(PathoGene)){
    pathoENS<-Patho_ENSgeneID[[pat]]
    DNM_Gene<-PathoGene[[pat]]
    ctrlENS<-ctr_ENSgeneID[['lgd']]
    ctrl_Gene <- ctrGene[['lgd']]
    
    ### function to fill the matrix of results
    #for (i in 1:length(clusters_list)){
    FUNC = function(i){
      Ms<-length(clusters_list[[i]])
      #### FET
      cat('\t',names(clusters_list)[i],'   \tpatho: ',names(PathoGene)[pat],'\t',pat,'\tof ',length(PathoGene),'\n',sep='')
      ## function to calculate the number Mc of DNMs in CTRL involving a gene of the cluster i
      y <- lapply(clusters_list[[i]],FUN = function(x) {ctrlENS[which(ctrlENS$ensembl_gene_id == x),'external_gene_name']})
      Mc <- sum(sapply(as.matrix(unique(y)), FUN = function(ym) {length(which(ctrl_Gene == ym ))}))
      McID <- paste(unlist(y),collapse=", ") 
      
      ## number NMc of remaining DNMs in CTRL involving a gene not in the cluster i
      NMc <- length(ctrl_Gene)-Mc
      
      ##function to calculate the number Mee of DNMs in patho involving a gene of the cluster i
      z <- lapply(clusters_list[[i]],FUN = function(x) {pathoENS[which(pathoENS$ensembl_gene_id == x),'external_gene_name']})
      Mp <- sum(sapply(as.matrix(unique(z)), FUN = function(zm) {length(which(DNM_Gene == zm ))}))
      MpID <- paste(unlist(z),collapse=", ") 
      
      ## number NMee of remaining DNMs in EE involving a gene not in the cluster i  
      NMp <- length(DNM_Gene)-Mp
      
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr <- matrix(c(Mp,Mc,NMp,NMc), nrow = 2)
      
      # FET
      #    FisherM <- fisher.test(matr,alternative = "greater")
      FisherM <- fisher.test(matr)
      Fisher.p <- FisherM$p.value
      Fisher.or <- FisherM$estimate
      Fisher.cinf <- FisherM$conf.int[1]
      Fisher.cis <- FisherM$conf.int[2]
      
      #### BET
            
      # theorical Ps = Theorical probablity of sucess based on mutation rate map
      CL.map.ens <- intersect(clusters_list[[i]],map$Gene)
      #NnID<-paste(setdiff(clusters_list[[i]],map$Gene), collapse=", ")
      nID<-length(setdiff(clusters_list[[i]],map$Gene))
      # if lgd (nonsens + missense)
      ThPs = sum(sapply(CL.map.ens,FUN=function(x){sum(map[which(map$Gene == x),c("Missense_rate","Nonsense_rate")])}))

      # nb of trials = nb of DNM falling in all map genes divided by the mutation rate on all map genes
      patho.map.ens <- intersect(pathoENS$ensembl_gene_id,map$Gene)
      y <- lapply(patho.map.ens,FUN = function(x) {pathoENS[which(pathoENS$ensembl_gene_id == x),'external_gene_name']})
      n <- round(sum(sapply(unique(y), FUN = function(x) {length(which(DNM_Gene == x ))}))/sum(map[,c("Missense_rate","Nonsense_rate")]))

      # nb of sucess = nb of DNM falling in map of the module
      z <- lapply(CL.map.ens,FUN = function(x) {pathoENS[which(pathoENS$ensembl_gene_id == x),'external_gene_name']})
      xz <- sum(sapply(z, FUN = function(x) {length(which(DNM_Gene == x ))}))

      BET<- binom.test(xz,n,ThPs)
      Binomial.p<-BET$p.value
      EsPs<-BET$estimate
      RobsE = xz /(n*ThPs)
      CI.inf<-BET$conf.int [1]
      CI.sup<-BET$conf.int [2]

      #FBET[i,]=c(Ms,names(PathoGene[pat]),Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mp,Mc,NMp,NMc,MpID,McID,Binomial.p,NA,ThPs,EsPs,RobsE,CI.inf,CI.sup,xz,n,NnID,nID)
      FBET[i,]=c(Ms,names(PathoGene[pat]),Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mp,Mc,NMp,NMc,MpID,McID,Binomial.p,NA,ThPs,EsPs,RobsE,CI.inf,CI.sup,xz,n,nID)
    }

    fbet <- mclapply(1:length(clusters_list),FUNC,mc.cores=detectCores())
    #The fbet output object of the mclapply function is a list of 44 vectors FBET[i,] in the good order

    for (i in 1:length(clusters_list)){
      FBET[i,]<-fbet[[i]]
    }

    FBET[,"FET FDR"]<- p.adjust(FBET[,"FET p.value"],method="fdr")
    FBET[,"BET FDR"]<- p.adjust(FBET[,"BET p.value"],method="fdr")
    write.table(FBET, sep = '\t', file = paste(outpath,"/",names(PathoGene)[pat],"_FBET_",runname,".txt",sep=""), row.names = TRUE, quote = FALSE, col.names = NA)
    npDNM[[pat]]<-FBET
  }  
  names(npDNM) <- names(PathoGene)

  if(selection==T){
    select<- intersect(which(as.numeric(npDNM[[1]][,"FET FDR"]) < 0.2),which(as.numeric(npDNM[[1]][,"BET FDR"]) < 0.2))
     cat("\tnumber of selected modules for\t\t", names(npDNM)[1]," :",length(select),'\n')
    if (length(select)==1){
        SignifT<- npDNM[[1]][c(select,NA),]
        SignifT<- SignifT[-which(is.na(rownames(SignifT)) ==T),]
      }else{
        SignifT<- npDNM[[1]][select,]
    }
    for (pat in 2:length(npDNM)){
      select<- intersect(which(as.numeric(npDNM[[pat]][,"FET FDR"]) < 0.2),which(as.numeric(npDNM[[pat]][,"BET FDR"]) < 0.2))
       cat("\tnumber of selected modules for\t\t", names(npDNM)[pat]," :",length(select),'\n')
      if (length(select)==1){
          SignifT<- rbind(SignifT,npDNM[[pat]][c(select,NA),])
          SignifT<- SignifT[-which(is.na(rownames(SignifT)) ==T),]
        }else{
          SignifT<- rbind(SignifT,npDNM[[pat]][select,])
      }          
    }
    write.table(SignifT, sep = '\t', file =paste(outpath,'/significantFBET_',runname,'.txt',sep=''), row.names = TRUE, quote = FALSE, col.names = NA) 
  }else{
    allT<- npDNM[[1]]
    for (pat in 2:length(npDNM)){
      allT<- rbind(allT,npDNM[[pat]])
    }
    write.table(allT, sep = '\t', file =paste(outpath,'/ALL_FBET_',runname,'.txt',sep=''), row.names = TRUE, quote = FALSE, col.names = NA) 
  }
  return(invisible(npDNM))
}





fet<-function(sampl,bkgrnd,success,counts=F,samp.success,bkgrnd.success,samp.fail,bkgrnd.fail,...){
# alternative ='greater'
# phyper(success_in_sample, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail=TRUE)

#fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='less');
# Numerical parameters in order:
# (success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part).

	if(!counts){
	    bkgrnd=bkgrnd[!(bkgrnd%in%sampl)]
	    fet=list(samp.success=sum(sampl%in%success),bkgrnd.success=sum(bkgrnd%in%success),samp.fail=sum(!(sampl%in%success)),bkgrnd.fail=sum(!(bkgrnd%in%success)))
	    test_mat=matrix(unlist(fet),nrow=2,dimnames=list(c('samp','bkgrnd'),c('success','fail')))
	    test_out=fisher.test(test_mat,...)
	
#		print(test_mat)

	    fet$n.genes=length(sampl)
	    fet$FETp=(test_out$p.value)
	    fet$fetOR=round(test_out$estimate) #,digits=3 
	    fet$lowerCI=round(test_out$conf.int[1])#,digits=3
	    fet$upperCI=round(test_out$conf.int[2])#,digits=3
#		fet$samp.success=paste(sampl[sampl%in%success],collapse=' ')
	    return(invisible((fet)))
	}

	if(counts){
		fet=list(samp.success=samp.success,bkgrnd.success=bkgrnd.success,samp.fail=samp.fail,bkgrnd.fail=bkgrnd.fail)
		test_mat=matrix(unlist(fet),nrow=2,dimnames=list(c('samp','bkgrnd'),c('success','fail')))
		test_out=fisher.test(test_mat,...)
    
#    	 print(test_mat)

		fet$n.genes=sum(samp.success,samp.fail)
		fet$FETp=(test_out$p.value)
		fet$fetOR=round(test_out$estimate,digits=3)
		fet$lowerCI=round(test_out$conf.int[1],digits=3)
		fet$upperCI=round(test_out$conf.int[2],digits=3)
#		fet$samp.success=paste(sampl[sampl%in%success],collapse=' ')
		return(invisible((fet)))

	}
}



msCellFET<-function(clusters_list,runname="",inpath = "~/Dropbox/tools/Data_to_load_CellFET",outpath=getwd(),selection=F){
cat('\tNOTE:\tclusters_list only ENSG ids currently supported')
  library(MetaDE)
  library('parallel')

  ###load data
  Load(paste(inpath,"/HUMQb.Rdata",sep="")) #"HUMQb" human ENSid orthologous genes of mice background genes
  Load(paste(inpath,"/hmscDF.Rdata",sep="")) #"hmscDF" human ENSid orthologous of mice single cell enriched by class dataframe
  
  ### create a matrix for results
  cFET <- matrix(nrow = length(clusters_list), ncol = 14)
  row.names(cFET) <- names (clusters_list)
  colnames(cFET) <- c("cell class","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                      "module cell enriched genes","module out cell enriched genes",
                      "non module cell enriched genes","non module out cell enriched genes",
                      "gene names of in modules cell enriched genes","module size","cell enriched genes size",
                      "% of genes in module in cell background"
  )
  print("Cell background is the genes with one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015 (Science)")
  
  resMsc<- list()
  for (ccl in 1:length(hmscDF)){ # ccl: cell class
    cclENS<-hmscDF[[ccl]]
    ### function to fill the matrix of results
    #for (i in 1:length(clusters_list)){
    FUNC = function(i){
      Ms<-length(clusters_list[[i]]) #Ms: module size
      CB<- HUMQb[,'hsapiens_homolog_ensembl_gene'] #CB: cell background
      Cs<-length(cclENS) #Cs: cell enriched genes size
      MCBp<-length(intersect(CB,clusters_list[[i]]))/Ms #MCBp: % of genes in module in cell background

      #cFET
      print(paste(names(clusters_list[i]),", cell class:",names(hmscDF)[ccl]))
      #calculate the number Mc of module i cell enriched genes (Mc: in module AND in cell class)
      Mc <- length(intersect(cclENS,clusters_list[[i]]))
      McID<- paste(unlist(HUMQb[which(CB %in% intersect(cclENS,clusters_list[[i]])),'external_gene_name']),collapse=", ")
      #calculate the number NMc of remaining genes not in module but in cell class
      NMc<- length(cclENS)-Mc
      #calculate the number Mnc of genes in module but not in cell class
      Mnc<- length(intersect(CB,clusters_list[[i]]))-Mc
      #calculate the number NMnc of genes out of module AND not in cell class
      NMnc<- length(CB)-(Mc+NMc+Mnc)
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr <- matrix(c(Mc,NMc,Mnc,NMnc), nrow = 2)
      #FET
      #FisherM <- fisher.test(matr,alternative = "greater")
      FisherM <- fisher.test(matr)
      Fisher.p <- FisherM$p.value
      Fisher.or <- FisherM$estimate
      Fisher.cinf <- FisherM$conf.int[1]
      Fisher.cis <- FisherM$conf.int[2]
      cFET[i,]=c(names(hmscDF)[ccl],Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mc,Mnc,NMc,NMnc,McID,Ms,Cs,MCBp)
    }
    cfet <- mclapply(1:length(clusters_list),FUNC,mc.cores=detectCores())
    #The cfet output object of the mclapply function is a list of n vectors cFET[i,] in the good order
    for (i in 1:length(clusters_list)){
      cFET[i,]<-cfet[[i]]
    }
    cFET[,"FET FDR"]<- p.adjust(cFET[,"FET p.value"],method="fdr")
    write.table(cFET, sep = '\t', file = paste(outpath,"/",names(hmscDF)[ccl],"_cFET_",runname,".txt",sep=""), row.names = TRUE, quote = FALSE, col.names = NA)
    resMsc[[ccl]]<-cFET
  }  
  names(resMsc) <- names(hmscDF)
  if(selection==T){
    select<- which(as.numeric(resMsc[[1]][,"FET FDR"]) < 0.2)
    print(paste("number of selected modules for ", names(resMsc)[1]," :",length(select)))
    if (length(select)==1){
      SignifT<- resMsc[[1]][c(select,NA),]
      SignifT<- SignifT[-which(is.na(rownames(SignifT)) ==T),]
    }else{
      SignifT<- resMsc[[1]][select,]
    }
    for (ccl in 2:length(resMsc)){
      select<- which(as.numeric(resMsc[[ccl]][,"FET FDR"]) < 0.2)
      print(paste("number of selected modules for ", names(resMsc)[ccl]," :",length(select)))
      if (length(select)==1){
        SignifT<- rbind(SignifT,resMsc[[ccl]][c(select,NA),])
        SignifT<- SignifT[-which(is.na(rownames(SignifT)) ==T),]
      }else{
        SignifT<- rbind(SignifT,resMsc[[ccl]][select,])
      }          
    }
    write.table(SignifT, sep = '\t', file =paste(outpath,'/significant_cFET_',runname,'.txt',sep=''), row.names = TRUE, quote = FALSE, col.names = NA) 
  }else{
    allT<- resMsc[[1]]
    for (ccl in 2:length(resMsc)){
      allT<- rbind(allT,resMsc[[ccl]])
    }
    write.table(allT, sep = '\t', file =paste(outpath,'/ALL_cFET_',runname,'.txt',sep=''), row.names = TRUE, quote = FALSE, col.names = NA) 
  }
  return(resMsc)      
}











forest.plot<-function(dat_lis,mode='single',x_lim=c(0,10),use_log=NA,phen=NA,module=NA,p_thresh=NA,point_width_scale=4,line_width=0.1,ci_bar_height=0.35){
  library(ggplot2)
  cat('\tWARNING : currently does not work for a single module\n\n')

#  colvec=c('red','darkblue','darkorange','magenta','orange')
   colvec=colmix
   cat('\tphenotypes detected :',paste(names(dat_lis),collapse=', '),'\n')
dat_lis=dat_lis[rev(names(dat_lis))]
options(warn=-1)
  if(!is.na(phen)){
     cat('\tdataset selected:',paste(phen,collapse=', '),'\n')
    dat_lis=dat_lis[names(dat_lis)%in%phen]
  }
    if(!is.na(module)){
       cat('\tmodules selected:',paste(module,collapse=', '),'\n')
  }
  plot_dat=list()
  for(ilis in 1:length(dat_lis)){
    if(!is.na(module)){
      dat_lis[[names(dat_lis)[ilis]]]=dat_lis[[names(dat_lis)[ilis]]][rownames(dat_lis[[names(dat_lis)[ilis]]])%in%module,]
   }

    plot_dat[[names(dat_lis)[ilis]]]=as.data.frame(make.numeric(dat_lis[[names(dat_lis)[ilis]]][,c('module size','FET p.value','OR','[95% OR CI] inf','OR [95% OR CI] sup')]))
      colnames(plot_dat[[names(dat_lis)[ilis]]])=c('n.genes','FETp','fetOR','lowerCI','upperCI')

    plot_dat[[names(dat_lis)[ilis]]]$module=as.factor(paste(rownames(plot_dat[[names(dat_lis)[ilis]]]),names(dat_lis)[ilis],sep='_'))
    plot_dat[[names(dat_lis)[ilis]]]$color=colvec[ilis]#c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
#    plot_dat[[names(dat_lis)[ilis]]]$color=c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
  plot_dat[[names(dat_lis)[ilis]]]$point_width=(plot_dat[[names(dat_lis)[ilis]]]$fetOR*3 )+point_width_scale
  }


   if(!is.na(p_thresh)){
    dummy=list()
    for(ilis in 1:length(plot_dat)){
      tempr=plot_dat[[names(dat_lis)[ilis]]][plot_dat[[names(dat_lis)[ilis]]]$FETp<p_thresh,]
      if(nrow(tempr)>0){
        dummy[[names(plot_dat)[ilis]]]=plot_dat[[names(plot_dat)[ilis]]][plot_dat[[names(plot_dat)[ilis]]]$FETp<p_thresh,]
      }
    }
    plot_dat=dummy
   }
options(warn=0)

#  if(!is.na(use_log)){
#    plot_dat$fetOR=log(as.numeric(plot_dat$fetOR),base=use_log)
#    plot_dat$lowerCI=log(as.numeric(plot_dat$lowerCI),base=use_log)
#    plot_dat$upperCI=log(as.numeric(plot_dat$upperCI),base=use_log)

#  }

####### since gplots insists on reversing all the phoenotypes, pre-empt it to keep original rownames
  for(ilis in 1:length(plot_dat)){
      plot_dat[[names(plot_dat)[ilis]]]=plot_dat[[names(plot_dat)[ilis]]][rev(rownames(plot_dat[[names(plot_dat)[ilis]]])),]
   }


gplots_dat=list()
  if(mode=='single'){
    for(idat in 1:length(plot_dat)){
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_dat[[idat]], aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data = plot_dat[[idat]], aes(y = module, x = fetOR), shape=15, cex=plot_dat[[idat]]$point_width, colour=plot_dat[[idat]]$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:9,seq(10,2000,by=10))) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))

    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat[[names(plot_dat)[idat]]]=(myplot + theme_bw())
    }
  }

  if(mode=='combined'){
    plot_merge=plot_dat[[1]]
    if(length(plot_dat)>1){
        for(ilis in 2:length(plot_dat)){
          plot_merge=rbind(plot_merge,plot_dat[[ilis]])
        }}
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_merge, aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data=plot_merge, aes(y = module, x = fetOR), shape=15, cex=plot_merge$point_width, colour=plot_merge$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:9,seq(10,2000,by=10))) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))



    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat=(myplot + theme_bw())
  }
  return(gplots_dat)
}






spatiotemp<-function(dat_lis,mod_gene,scale_data=F,row_clust=T,use_cols=NA,ncols=10,...){
#cat('\n\tNOTE: this function requires an object "name_here", available from:\n https://\n\n')	##  not yet but needs to be implemented


#### spatiotemproral transcriptome of human brain data:
#Load('~/Dropbox/PROJ/spatem/dtb/expr/GSE25219_GPL5175.exon_array.HUGO_GENE_Mapped.ngene13830.log26.nsamp1192.rin6.5_pmi24.covar_pmi_rin.Rdata')

## built to handle data_list with incomplete sub-list categories (ie missing)
# HIP contains all dev stages -> spatemp() works with all stages in correct order
#ballg=ballg[c('HIP',names(ballg)[names(ballg)!='HIP'])]
#ballc=ballc[names(ballg)]
#save(ballc,ballg,bfolc,bfold,bstats,gexpc,gexpr,ids,readme,sampl,file='~/Dropbox/PROJ/spatem/dtb/expr/GSE25219_GPL5175.exon_array.HUGO_GENE_Mapped.ngene13830.log26.nsamp1192.rin6.5_pmi24.covar_pmi_rin.Rdata')

  library(gplots)

  if(sum(mod_gene%in%rownames(dat_lis[[1]][[1]]))==0){
    stop('\tno query genes not found in dataset')
  }

  cat('\tnumber of query genes in dataset:\n')
  print(matst(mod_gene%in%rownames(dat_lis[[1]][[1]])))

  dev_stage=names(dat_lis[[1]])
  for(ilis in 2:length(dat_lis)){
    dev_stage=unique(c(dev_stage,names(dat_lis[[ilis]])))
  }


plot_dat=matrix(NA,nrow=length(dat_lis),ncol=length(dev_stage))
  colnames(plot_dat)=dev_stage
  rownames(plot_dat)=names(dat_lis)
 
#   Head(plot_dat)
 cat('\tcalculate median expression for query genes\n')
  for(ireg in 1:length(dat_lis)){
    for(idev in 1:length(dat_lis[[ireg]])){
      dummy=dat_lis[[names(dat_lis)[ireg]]][[names(dat_lis[[ireg]])[idev]]]
#      Head(dummy)
      plot_dat[names(dat_lis)[ireg],names(dat_lis[[ireg]])[idev]]=mean(apply(dummy[rownames(dummy)%in%mod_gene,,drop=F],1,mean))

    }
    cat(round(ireg/length(dat_lis),digits=2),'\r');flush.console()
  }

# cat(min(plot_dat,na.rm=T),max(plot_dat,na.rm=T),'\n')
  if(!scale_data){
    na_val_loc=is.na(plot_dat)
    plot_dat[na_val_loc]=(min(plot_dat,na.rm=T))-0.2

  }
  if(scale_data){
    plot_dat=t(make.numeric(scale(t(plot_dat),center=T,scale=F)))
      na_val_loc=is.na(plot_dat)
    plot_dat[na_val_loc]=(min(plot_dat,na.rm=T))-0.1
  }

#  Head(plot_dat)
#### check for user defined colors
 options(warn=-1)
#  if(is.na(use_cols)){heat_cols=c('grey60',colorRampPalette(c("#0072B2","#56B4E9","#abd9e9"))(ncols/2),colorRampPalette(c("#ffffbf","#F0E442","darkred"))(ncols/2))} #'#e31a1c'   ,
  if(is.na(use_cols)){heat_cols=c('grey60',colorRampPalette(c("#053061","#56B4E9","#abd9e9"))(ncols/2),colorRampPalette(c("#ffffbf","#F0E442","darkred"))(ncols/2))} #'#e31a1c'   ,

  if(!is.na(use_cols)){
    heat_cols=use_cols
  }

      heatmap.2(
        plot_dat
#        ,main= dat_descr
        ,Colv=F
        ,Rowv=row_clust
#        ,col=c('grey60',colorRampPalette(c("white","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026","darkred"))(ncols))
#        ,col=c('grey60',bluered(ncols))
        ,col=heat_cols
        ,key.title = NA
        ,key.xlab = NA
        ,key.ylab = NA
        ,trace='none'
        ,srtCol=0
        ,keysize = 1.1
        ,density.info='none'
        ,colsep=1:(ncol(plot_dat)-1),rowsep=1:(nrow(plot_dat)-1),sepcolor="white",sepwidth=c(0.001,0.001)
#        ,margins=c(5,5)
        ,...
        )
#  plot.new()
        par(new = TRUE)
      mtext('grey - no samples available',adj=1,side=3,line=1)
#      mtext(dat_descr,adj=1,side=3,line=2)
 #     mtext(paste0("n.complete : ",n_overlap),adj=1,side=3)

 options(warn=0)

 plot_dat[na_val_loc]=NA
 return(invisible(plot_dat))
}




Venn<-function(dat_lis,main='',...){
 library(gplots)
  venn(dat_lis)
  mtext(main,side=3,...)
}




gsea.run<-function(
  path_GSEA="/Users/adelahay/Documents/Andree/projets/iGENEE_ICL/data_and_scripts/R_scripts/GSEA/",
  file_gene_sets_with_path,rank_file_with_path,chip_file_with_path,nperm=10000,name_run,max_clust_size=5000,min_clust_size=10,path_out_folder){
  system(paste("mkdir ",path_out_folder,sep=""))
  options(scipen=999)
  ##-Xmx5000mm flags the amount of memory available to java. The default is -Xmx512m
  system(paste("java -cp ",path_GSEA, "gsea2-2.2.1.jar -Xmx11000mm xtools.gsea.GseaPreranked -gmx ",file_gene_sets_with_path," -collapse false -mode Max_probe -norm meandiv -nperm ",nperm,
               " -rnk ",rank_file_with_path, " -scoring_scheme classic -rpt_label ",name_run, " -chip ",chip_file_with_path, " -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max ",
               max_clust_size, " -set_min ",min_clust_size," -zip_report false -out ",path_out_folder, " -gui false",sep=""))
}










dnm.enrich<-function(clusters_list,runname="",selection=F){ 				#,inpath="~/Dropbox/SHARED/tools/FUNDNMmap/"
  cat('\tNOTE:\tclusters_list only ENSG ids currently supported\n\n')
library(MetaDE)
library('parallel')

#  Load(paste(inpath,"/PathoGeneENS.Rdata",sep="")) # PathoGene & Patho_ENSgeneID
#  Load(paste(inpath,"/ctr_GeneENS.Rdata",sep="")) # ctrGene & ctr_ENSgeneID
#  gene_map_dat=read.table(file=paste(inpath,'/functional_mut_rate.bias_corrected.local.canonical_tx_only.bed.txt',sep=""),header=TRUE, sep='\t', blank.lines.skip=TRUE)

  ### create a matrix for results
  FBET=matrix(nrow=length(clusters_list), ncol=23)
  row.names(FBET)=names(clusters_list)
  colnames(FBET)=c("module size","patho","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                        "module DNMs in patients","module DNMs in controls",
                        "non module DNMs in patients","non module DNMs in controls",
                        "gene names of modules DNMs in patients",
                        "gene names of modules DNMs in controls",
                      "BET p.value","BET FDR","Theo Ps(mutation rate of map.M)","Estimated Ps","Ratio Obs/Exp",
                        "[95% EsPs CI] inf","[95% EsPs CI] sup",
                        "x=nb of DNM in map.M(success)","n=nb of DNM in map/mutation rate all genes of map(trials)",
                        #"ENSgeneID of M not in map",
                        "nb ENSgeneID of M not in map"
                      )
  
  npDNM=list()
  for(pat in 1:length(PathoGene)){
    cat('\t=====================',names(PathoGene)[pat],'=====================',pat,' of ',length(PathoGene),'\n')
    pathoENS=Patho_ENSgeneID[[pat]]
    DNM_Gene=PathoGene[[pat]]
    ctrlENS=ctr_ENSgeneID[['lgd']]
    ctrl_Gene=ctrGene[['lgd']]
    
    ### function to fill the matrix of results
    #for(i in 1:length(clusters_list)){
    FUNC=function(i){
      Ms=length(clusters_list[[i]])
      #### FET
      cat('\t\t',names(clusters_list)[i],'\n')
      ## function to calculate the number Mc of DNMs in CTRL involving a gene of the cluster i
      y=lapply(clusters_list[[i]],FUN=function(x) {ctrlENS[which(ctrlENS$ensembl_gene_id==x),'external_gene_name']})
      Mc=sum(sapply(as.matrix(unique(y)), FUN=function(ym) {length(which(ctrl_Gene==ym ))}))
      McID=paste(unlist(y),collapse=", ") 
      
      ## number NMc of remaining DNMs in CTRL involving a gene not in the cluster i
      NMc=length(ctrl_Gene)-Mc
      
      ##function to calculate the number Mee of DNMs in patho involving a gene of the cluster i
      z=lapply(clusters_list[[i]],FUN=function(x) {pathoENS[which(pathoENS$ensembl_gene_id==x),'external_gene_name']})
      Mp=sum(sapply(as.matrix(unique(z)), FUN=function(zm) {length(which(DNM_Gene==zm ))}))
      MpID=paste(unlist(z),collapse=", ") 
      
      ## number NMee of remaining DNMs in EE involving a gene not in the cluster i  
      NMp=length(DNM_Gene)-Mp
      
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr=matrix(c(Mp,Mc,NMp,NMc), nrow=2)
      
      # FET
      #    FisherM=fisher.test(matr,alternative="greater")
      FisherM=fisher.test(matr)
      Fisher.p=FisherM$p.value
      Fisher.or=FisherM$estimate
      Fisher.cinf=FisherM$conf.int[1]
      Fisher.cis=FisherM$conf.int[2]
      
      #### BET
            
      # theorical Ps=Theorical probablity of sucess based on mutation rate map
      CL.map.ens=intersect(clusters_list[[i]],gene_map_dat$Gene)
      #NnID=paste(setdiff(clusters_list[[i]],map$Gene), collapse=", ")
      nID=length(setdiff(clusters_list[[i]],gene_map_dat$Gene))
      # if lgd(nonsens + missense)
      ThPs=sum(sapply(CL.map.ens,FUN=function(x){sum(gene_map_dat[which(gene_map_dat$Gene==x),c("Missense_rate","Nonsense_rate")])}))

      # nb of trials=nb of DNM falling in all map genes divided by the mutation rate on all map genes
      patho.map.ens=intersect(pathoENS$ensembl_gene_id,gene_map_dat$Gene)
      y=lapply(patho.map.ens,FUN=function(x) {pathoENS[which(pathoENS$ensembl_gene_id==x),'external_gene_name']})
      n=round(sum(sapply(unique(y), FUN=function(x) {length(which(DNM_Gene==x ))}))/sum(gene_map_dat[,c("Missense_rate","Nonsense_rate")]))

      # nb of sucess=nb of DNM falling in map of the module
      z=lapply(CL.map.ens,FUN=function(x) {pathoENS[which(pathoENS$ensembl_gene_id==x),'external_gene_name']})
      xz=sum(sapply(z, FUN=function(x) {length(which(DNM_Gene==x ))}))

      BET=binom.test(xz,n,ThPs)
      Binomial.p=BET$p.value
      EsPs=BET$estimate
      RobsE=xz /(n*ThPs)
      CI.inf=BET$conf.int [1]
      CI.sup=BET$conf.int [2]

      #FBET[i,]=c(Ms,names(PathoGene[pat]),Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mp,Mc,NMp,NMc,MpID,McID,Binomial.p,NA,ThPs,EsPs,RobsE,CI.inf,CI.sup,xz,n,NnID,nID)
      FBET[i,]=c(Ms,names(PathoGene[pat]),Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mp,Mc,NMp,NMc,MpID,McID,Binomial.p,NA,ThPs,EsPs,RobsE,CI.inf,CI.sup,xz,n,nID)
    }

    fbet=mclapply(1:length(clusters_list),FUNC,mc.cores=detectCores())
    #The fbet output object of the mclapply function is a list of 44 vectors FBET[i,] in the good order

    for(i in 1:length(clusters_list)){
      FBET[i,]=fbet[[i]]
    }

    FBET[,"FET FDR"]=p.adjust(FBET[,"FET p.value"],method="fdr")
    FBET[,"BET FDR"]=p.adjust(FBET[,"BET p.value"],method="fdr")
#    write.table(FBET, sep='\t', file=paste(outpath,"/",names(PathoGene)[pat],"_FBET_",runname,".txt",sep=""), row.names=TRUE, quote=FALSE, col.names=NA)
    npDNM[[pat]]=FBET
  }  
  names(npDNM)=names(PathoGene)

  if(selection==T){
    select=intersect(which(as.numeric(npDNM[[1]][,"FET FDR"]) < 0.2),which(as.numeric(npDNM[[1]][,"BET FDR"]) < 0.2))
     cat("\tnumber of selected modules for\t\t", names(npDNM)[1]," :",length(select),'\n')
    if(length(select)==1){
        SignifT=npDNM[[1]][c(select,NA),]
        SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
      }else{
        SignifT=npDNM[[1]][select,]
    }
    for(pat in 2:length(npDNM)){
      select=intersect(which(as.numeric(npDNM[[pat]][,"FET FDR"]) < 0.2),which(as.numeric(npDNM[[pat]][,"BET FDR"]) < 0.2))
       cat("\tnumber of selected modules for\t\t", names(npDNM)[pat]," :",length(select),'\n')
      if(length(select)==1){
          SignifT=rbind(SignifT,npDNM[[pat]][c(select,NA),])
          SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
        }else{
          SignifT=rbind(SignifT,npDNM[[pat]][select,])
      }          
    }
#    write.table(SignifT, sep='\t', file=paste(outpath,'/significantFBET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }else{
    allT=npDNM[[1]]
    for(pat in 2:length(npDNM)){
      allT=rbind(allT,npDNM[[pat]])
    }
#    write.table(allT, sep='\t', file=paste(outpath,'/ALL_FBET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }
  return(invisible(npDNM))
}


cellt.enrich<-function(clusters_list,runname="",outpath=getwd(),selection=F){	#inpath="~/Dropbox/SHARED/tools/Data_to_load_CellFET",
cat('\tNOTE:\tclusters_list only ENSG ids currently supported')
cat("Cell background is the genes with one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)")
  library(MetaDE)
  library('parallel')

  ###load data
#  Load(paste(inpath,"/HUMQb.Rdata",sep="")) #"HUMQb" human ENSid orthologous genes of mice background genes
#  Load(paste(inpath,"/hmscDF.Rdata",sep="")) #"hmscDF" human ENSid orthologous of mice single cell enriched by class dataframe
  
  ### create a matrix for results
  cFET=matrix(nrow=length(clusters_list), ncol=14)
  row.names(cFET)=names(clusters_list)
  colnames(cFET)=c("cell class","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                      "module cell enriched genes","module out cell enriched genes",
                      "non module cell enriched genes","non module out cell enriched genes",
                      "gene names of in modules cell enriched genes","module size","cell enriched genes size",
                      "% of genes in module in cell background"
  )

  
  resMsc=list()
  for(ccl in 1:length(hmscDF)){ # ccl: cell class
    cat('\t=====================',names(hmscDF)[ccl],'=====================',ccl,' of ',length(hmscDF),'\n')
    cclENS=hmscDF[[ccl]]
    ### function to fill the matrix of results
    #for(i in 1:length(clusters_list)){
    FUNC=function(i){
      Ms=length(clusters_list[[i]]) #Ms: module size
      CB=HUMQb[,'hsapiens_homolog_ensembl_gene'] #CB: cell background
      Cs=length(cclENS) #Cs: cell enriched genes size
      MCBp=length(intersect(CB,clusters_list[[i]]))/Ms #MCBp: % of genes in module in cell background

      #cFET
      cat('\t\t',names(clusters_list)[i],'\n')
      #calculate the number Mc of module i cell enriched genes(Mc: in module AND in cell class)
      Mc=length(intersect(cclENS,clusters_list[[i]]))
      McID=paste(unlist(HUMQb[which(CB %in% intersect(cclENS,clusters_list[[i]])),'external_gene_name']),collapse=", ")
      #calculate the number NMc of remaining genes not in module but in cell class
      NMc=length(cclENS)-Mc
      #calculate the number Mnc of genes in module but not in cell class
      Mnc=length(intersect(CB,clusters_list[[i]]))-Mc
      #calculate the number NMnc of genes out of module AND not in cell class
      NMnc=length(CB)-(Mc+NMc+Mnc)
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr=matrix(c(Mc,NMc,Mnc,NMnc), nrow=2)
      #FET
      #FisherM=fisher.test(matr,alternative="greater")
      FisherM=fisher.test(matr)
      Fisher.p=FisherM$p.value
      Fisher.or=FisherM$estimate
      Fisher.cinf=FisherM$conf.int[1]
      Fisher.cis=FisherM$conf.int[2]
      cFET[i,]=c(names(hmscDF)[ccl],Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mc,Mnc,NMc,NMnc,McID,Ms,Cs,MCBp)
    }
#    cfet=mclapply(1:length(clusters_list),FUNC,mc.cores=detectCores())
    cfet=lapply(1:length(clusters_list),FUNC)
    #The cfet output object of the mclapply function is a list of n vectors cFET[i,] in the good order
    for(i in 1:length(clusters_list)){
      cFET[i,]=cfet[[i]]
    }
    cFET[,"FET FDR"]=p.adjust(cFET[,"FET p.value"],method="fdr")
#    write.table(cFET, sep='\t', file=paste(outpath,"/",names(hmscDF)[ccl],"_cFET_",runname,".txt",sep=""), row.names=TRUE, quote=FALSE, col.names=NA)
    resMsc[[ccl]]=cFET
  }  
  names(resMsc)=names(hmscDF)
  if(selection==T){
    select=which(as.numeric(resMsc[[1]][,"FET FDR"]) < 0.2)
    cat("number of selected modules for ", names(resMsc)[1]," :",length(select),'\n')
    if(length(select)==1){
      SignifT=resMsc[[1]][c(select,NA),]
      SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
    }else{
      SignifT=resMsc[[1]][select,]
    }
    for(ccl in 2:length(resMsc)){
      select=which(as.numeric(resMsc[[ccl]][,"FET FDR"]) < 0.2)
      cat("number of selected modules for ", names(resMsc)[ccl]," :",length(select),'\n')
      if(length(select)==1){
        SignifT=rbind(SignifT,resMsc[[ccl]][c(select,NA),])
        SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
      }else{
        SignifT=rbind(SignifT,resMsc[[ccl]][select,])
      }          
    }
#    write.table(SignifT, sep='\t', file=paste(outpath,'/significant_cFET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }else{
    allT=resMsc[[1]]
    for(ccl in 2:length(resMsc)){
      allT=rbind(allT,resMsc[[ccl]])
    }
#    write.table(allT, sep='\t', file=paste(outpath,'/ALL_cFET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }
  return(resMsc)      
}



dnm.enrich.bg<-function(module_list,bg_vect=NA,use_counts=F){
#cat('\n\tNOTE:\tid_type - options: "name" - gene name / HUGO name; "ensg" - ensembl gene id;"mouse.ensg" - mouse ensembl gene id - one2one orthologs to human only\n')
cat('\tNOTE:\tbg_vect - options: "NA" - one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n\n')
##  based on bg input - if data.frame or matrix -> create a list of identical bg length=length(module_list)
##  + alternatively a list of modules with backgrounds in the same order -> flexibility option for unique backgrounds for each module
##  + added benefit of allowing different color coding for each bg (eg darker or different style of line / OR dot)

##  option for supporting multiple ID types -> ensg /+/ gene name
##  + auto-detect ID option - a possibility but seems unnecessary given only 2 types


##===================================================================================================================================
## generate background as per options, if list provided (same length as modules)=================================
options(warn=-1)

#  if(!is.na(bg_vect)){
  if(class(bg_vect)=='list'){
    if(length(bg_vect)==length(module_list)){
      cat('\tusing multiple bg list as provided\n')
    }
    if(length(bg_vect)!=length(module_list)){
      stop('\tlength of bg_vect provided not the same length as module_list. if same bkgrnd for all modules, do not provide as.list()\n')
    }
  }

  if(class(bg_vect)!='list' & !is.na(bg_vect)){
    dnmid=lapply(dnmid,function(x){unique(x[x[,'ids']%in%bg_vect|x[,'gene']%in%bg_vect,c('gene','count','ids')])})  ## limit ids to just those in bg_vect
    cat('\tusing the provided bkgrnd, same for all modules\n')
  }

  if(is.na(bg_vect)){
    cat('\tno background list provided, using default\n')
#    bg_vect=unique(unlist(lapply(dnmid[c('EE','ASD','ID','SCZ','DDD','lgd')],function(x){unique(x[,c('gene')])})))  ## basically all the genes assayed
  }
options(warn=0)

##===================================================================================================================================
## perform the fisher.test() using fet() function=================================

  plot_dat=list()
  for(idat in 1:(length(dnmid)-4)){
     cat('\t=====================',names(dnmid)[idat],'=====================',idat,' of ',(length(dnmid)-4),'\n')
     dumpty=list()
    for(imod in 1:length(module_list)){

       cat('\t\t',names(module_list)[imod],'\n')
       humpty=list()

## case DNM in module
       ms=dnmid[[idat]]
       ms=unique(ms[ms$ids%in%module_list[[names(module_list)[imod]]] | ms$gene%in%module_list[[names(module_list)[imod]]],c('gene','count')])

## control DNM in module
       mf=dnmid[['lgd']]
       mf=unique(mf[mf$ids%in%module_list[[names(module_list)[imod]]] | mf$gene%in%module_list[[names(module_list)[imod]]],c('gene','count')])

## case DNM not in module
       bs=dnmid[[idat]]
       bs=unique(bs[!(bs$ids%in%module_list[[names(module_list)[imod]]] | bs$gene%in%module_list[[names(module_list)[imod]]]),c('gene','count')])

## control DNM not in module
       bf=dnmid[['lgd']]
       bf=unique(bf[!(bf$ids%in%module_list[[names(module_list)[imod]]] | bf$gene%in%module_list[[names(module_list)[imod]]]),c('gene','count')])


       if(!use_counts){
        ms=length(unique(ms$gene))
        mf=length(unique(mf$gene))
        bs=length(unique(bs$gene))
        bf=length(unique(bf$gene))
        
        humpty=fetc(samp.sucess=ms,bkgrnd.success=bs,samp.fail=mf,bkgrnd.fail=bf)
       }
       if(use_counts){
        ms=sum(ms$count)
        mf=sum(mf$count)
        bs=sum(bs$count)
        bf=sum(bf$count)
        
        humpty=fetc(samp.sucess=ms,bkgrnd.success=bs,samp.fail=mf,bkgrnd.fail=bf)

        }
        humpty$n.genes=length(module_list[[names(module_list)[imod]]])    
#      if(is.infinite(humpty$upperCI)){humpty$upperCI=1}
#      if(is.infinite(humpty$lowerCI)){humpty$upperCI=0}

#      if(do_plots){
#        colvec=colmix
#        humpty$module=as.factor(paste(names(module_list)[imod],names(dnmid)[imod],sep='_'))
#        humpty$color=colvec[imod]#c(rep(colvec[imod],nrow(plot_dat[[names(module_list)[imod]]])-1),'black')
#        humpty$point_width=(humpty[[names(module_list)[imod]]]$fetOR*3)+point_width_scale
 #       }
      dumpty[[names(module_list)[[imod]]]]=(unlist(humpty))
    }

    plot_dat[[names(dnmid)[idat]]]=(t(as.data.frame(dumpty)))

  }

#  if(!do_plots){
    return(plot_dat)
#  }

#  if(do_plots){
  ##===================================================================================================================================
  ## generate plots=================================

}




#### this enrichment test is only for DNM in module compared to genome/background
##<<>>dnm.enrich.bg<-function(module_list,bg_vect=NA,use_counts=F){
##<<>>#cat('\n\tNOTE:\tid_type - options: "name" - gene name / HUGO name; "ensg" - ensembl gene id;"mouse.ensg" - mouse ensembl gene id - one2one orthologs to human only\n')
##<<>>cat('\tNOTE:\tbg_vect - options: "NA" - one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n\n')
##<<>>##  based on bg input - if data.frame or matrix -> create a list of identical bg length=length(module_list)
##<<>>##  + alternatively a list of modules with backgrounds in the same order -> flexibility option for unique backgrounds for each module
##<<>>##  + added benefit of allowing different color coding for each bg (eg darker or different style of line / OR dot)
##<<>>
##<<>>##  option for supporting multiple ID types -> ensg /+/ gene name
##<<>>##  + auto-detect ID option - a possibility but seems unnecessary given only 2 types
##<<>>
##<<>>
##<<>>##===================================================================================================================================
##<<>>## generate background as per options, if list provided (same length as modules)=================================
##<<>>options(warn=-1)
##<<>>
##<<>>#  if(!is.na(bg_vect)){
##<<>>  if(class(bg_vect)=='list'){
##<<>>    if(length(bg_vect)==length(module_list)){
##<<>>      cat('\tusing multiple bg list as provided\n')
##<<>>    }
##<<>>    if(length(bg_vect)!=length(module_list)){
##<<>>      stop('\tlength of bg_vect provided not the same length as module_list. if same bkgrnd for all modules, do not provide as.list()\n')
##<<>>    }
##<<>>  }
##<<>>
##<<>>  if(class(bg_vect)!='list' & !is.na(bg_vect)){
##<<>>    dnmid=lapply(dnmid,function(x){unique(x[x[,'ids']%in%bg_vect|x[,'gene']%in%bg_vect,c('gene','count','ids')])})  ## limit ids to just those in bg_vect
##<<>>    cat('\tusing the provided bkgrnd, same for all modules\n')
##<<>>  }
##<<>>
##<<>>  if(is.na(bg_vect)){
##<<>>    cat('\tno background list provided, using default\n')
##<<>>    bg_vect=unique(unlist(lapply(dnmid[c('EE','ASD','ID','SCZ','DDD','lgd')],function(x){unique(x[,c('gene')])})))  ## basically all the genes assayed
##<<>>  }
##<<>>options(warn=0)
##<<>>
##<<>>##===================================================================================================================================
##<<>>## perform the fisher.test() using fet() function=================================
##<<>>
##<<>>  plot_dat=list()
##<<>>  for(idat in 1:(length(dnmid)-4)){
##<<>>     cat('\t=====================',names(dnmid)[idat],'=====================',idat,' of ',(length(dnmid)-4),'\n')
##<<>>     dumpty=list()
##<<>>    for(imod in 1:length(module_list)){
##<<>>
##<<>>       cat('\t\t',names(module_list)[imod],'\n')
##<<>>       humpty=list()
##<<>>
##<<>>	   not_mod=bg_vect[!(bg_vect%in%module_list[[names(module_list)[imod]]])]
##<<>>## case dnm in module
##<<>>       ms=dnmid[[idat]]
##<<>>       ms=unique(ms[ms$ids%in%module_list[[names(module_list)[imod]]] | ms$gene%in%module_list[[names(module_list)[imod]]],c('gene','count')])
##<<>>
##<<>>## control dnm in module
##<<>>       mf=dnmid[['lgd']]
##<<>>       mf=unique(mf[mf$ids%in%module_list[[names(module_list)[imod]]] | mf$gene%in%module_list[[names(module_list)[imod]]],c('gene','count')])
##<<>>
##<<>>## case dnm outside module
##<<>>       bs=dnmid[[idat]]
##<<>>       bs=unique(bs[bs$ids%in%not_mod | bs$gene%in%not_mod,c('gene','count')])
##<<>>
##<<>>## control dnm outside module
##<<>>       bf=dnmid[['lgd']]
##<<>>       bf=unique(bf[bf$ids%in%not_mod | bf$gene%in%not_mod,c('gene','count')])
##<<>>
##<<>>
##<<>>       if(!use_counts){
##<<>>       	ms=length(unique(ms$gene))
##<<>>       	mf=length(unique(mf$gene))
##<<>>       	bs=length(unique(bs$gene))
##<<>>       	bf=length(unique(bf$gene))
##<<>>       	
##<<>>       	humpty=fetc(samp.sucess=ms,bkgrnd.success=bs,samp.fail=mf,bkgrnd.fail=bf)
##<<>>       }
##<<>>       if(use_counts){
##<<>>       	ms=sum(ms$count)
##<<>>       	mf=sum(mf$count)
##<<>>       	bs=sum(bs$count)
##<<>>       	bf=sum(bf$count)
##<<>>       	
##<<>>       	humpty=fetc(samp.sucess=ms,bkgrnd.success=bs,samp.fail=mf,bkgrnd.fail=bf)
##<<>>
##<<>>       	}
##<<>>       	humpty$n.genes=length(module_list[[names(module_list)[imod]]])    
##<<>>#      if(is.infinite(humpty$upperCI)){humpty$upperCI=1}
##<<>>#      if(is.infinite(humpty$lowerCI)){humpty$upperCI=0}
##<<>>
##<<>>#      if(do_plots){
##<<>>#        colvec=colmix
##<<>>#        humpty$module=as.factor(paste(names(module_list)[imod],names(dnmid)[imod],sep='_'))
##<<>>#        humpty$color=colvec[imod]#c(rep(colvec[imod],nrow(plot_dat[[names(module_list)[imod]]])-1),'black')
##<<>>#        humpty$point_width=(humpty[[names(module_list)[imod]]]$fetOR*3)+point_width_scale
##<<>> #       }
##<<>>      dumpty[[names(module_list)[[imod]]]]=(unlist(humpty))
##<<>>    }
##<<>>
##<<>>    plot_dat[[names(dnmid)[idat]]]=(t(as.data.frame(dumpty)))
##<<>>
##<<>>  }
##<<>>
##<<>>#  if(!do_plots){
##<<>>    return(plot_dat)
##<<>>#  }
##<<>>
##<<>>#  if(do_plots){
##<<>>  ##===================================================================================================================================
##<<>>  ## generate plots=================================
##<<>>
##<<>>}
##<<>>
##<<>>

### appears to work..
cellt.enrich.bg<-function(module_list,bg_list=NA,id_type='name'){
cat('\n\tNOTE:\tid_type - options: "name" - gene name / HUGO name; "ensg" - ensembl gene id;"mouse.ensg" - mouse ensembl gene id - one2one orthologs to human only\n')
cat('\tNOTE:\tbg_list - options: "NA" - one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n\n')
##  based on bg input - if data.frame or matrix -> create a list of identical bg length=length(module_list)
##  + alternatively a list of modules with backgrounds in the same order -> flexibility option for unique backgrounds for each module
##  + added benefit of allowing different color coding for each bg (eg darker or different style of line / OR dot)

##  option for supporting multiple ID types -> ensg /+/ gene name
##  + auto-detect ID option - a possibility but seems unnecessary given only 2 types
#dtb_path='/Users/ks/Dropbox/PROJ/annot/dtb/processed'
#load(paste0(dtb_path,'/cellt.enrich.zeizel.dtb.Rdata'))

##===================================================================================================================================
## generate background as per options, if list provided (same length as modules)=================================
options(warn=-1)


#  if(!is.na(bg_list)){
  if(class(bg_list)=='list'){
    if(length(bg_list)==length(module_list)){
      cat('\tusing multiple bg list as provided\n')
    }
    if(length(bg_list)!=length(module_list)){
      stop('\tlength of bg_list provided not the same length as module_list. if same bkgrnd for all modules, do not provide as.list()\n')
    }
  }

  if(class(bg_list)!='list' & !is.na(bg_list)){
  	cellid=lapply(cellid,function(x){x[,id_type][x[,id_type]%in%bg_list]})	## limit ids to just those in bg_list
    cat('\tusing the provided bkgrnd, same for all modules\n')
    bg_list=list()
    for(imod in 1:length(module_list)){
      bg_list[[names(module_list)[imod]]]=cellid$bkgrnd
    }

  }

  if(is.na(bg_list)){
    cat('\tno background list provided, using default\n')
    bg_list=list()
    for(imod in 1:length(module_list)){
      bg_list[[names(module_list)[imod]]]=cellid$bkgrnd
    }
  }
options(warn=0)

##===================================================================================================================================
## perform the fisher.test() using fet() function=================================

  plot_dat=list()
  for(idat in 1:(length(cellid)-1)){
     cat('\t=====================',names(cellid)[idat],'=====================',idat,' of ',(length(cellid)-1),'\n')
     dumpty=list()
    for(imod in 1:length(module_list)){
      humpty=list()
       cat('\t\t',names(module_list)[imod],'\n')
      humpty=fet(sampl=module_list[[names(module_list)[imod]]],bg_list[[imod]],success=cellid[[idat]])
#      if(is.infinite(humpty$upperCI)){humpty$upperCI=1}
#      if(is.infinite(humpty$lowerCI)){humpty$upperCI=0}

#      if(do_plots){
#        colvec=colmix
#        humpty$module=as.factor(paste(names(module_list)[imod],names(cellid)[imod],sep='_'))
#        humpty$color=colvec[imod]#c(rep(colvec[imod],nrow(plot_dat[[names(module_list)[imod]]])-1),'black')
#        humpty$point_width=(humpty[[names(module_list)[imod]]]$fetOR*3)+point_width_scale
 #       }
      dumpty[[names(module_list)[[imod]]]]=(unlist(humpty))
    }

    plot_dat[[names(cellid)[idat]]]=(t(as.data.frame(dumpty)))

  }

#  if(!do_plots){
    return(plot_dat)
#  }

#  if(do_plots){
  ##===================================================================================================================================
  ## generate plots=================================

}





forestp<-function(dat_lis,mode='single',x_lim=c(0,10),use_log=NA,phen=NA,module=NA,p_thresh=NA,point_width_scale=4,line_width=0.1,ci_bar_height=0.35){
  library(ggplot2)
  cat('\tWARNING : currently does not work for a single module\n\n')

#  colvec=c('red','darkblue','darkorange','magenta','orange')
   colvec=colmix
   cat('\tphenotypes detected :',paste(names(dat_lis),collapse=', '),'\n')

options(warn=-1)
  if(!is.na(phen)){
     cat('\tdataset selected:',paste(phen,collapse=', '),'\n')
    dat_lis=dat_lis[names(dat_lis)%in%phen]
  }
    if(!is.na(module)){
       cat('\tmodules selected:',paste(module,collapse=', '),'\n')
  }
  plot_dat=list()
  for(ilis in 1:length(dat_lis)){
    if(!is.na(module)){
      dat_lis[[names(dat_lis)[ilis]]]=dat_lis[[names(dat_lis)[ilis]]][rownames(dat_lis[[names(dat_lis)[ilis]]])%in%module,]
   }

    plot_dat[[names(dat_lis)[ilis]]]=as.data.frame(make.numeric(dat_lis[[names(dat_lis)[ilis]]][,c('n.genes','FETp','fetOR.odds ratio','lowerCI','upperCI')]))
      colnames(plot_dat[[names(dat_lis)[ilis]]])=c('n.genes','FETp','fetOR','lowerCI','upperCI')

    plot_dat[[names(dat_lis)[ilis]]]$module=as.factor(paste(rownames(plot_dat[[names(dat_lis)[ilis]]]),names(dat_lis)[ilis],sep='_'))
    plot_dat[[names(dat_lis)[ilis]]]$color=colvec[ilis]#c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
#    plot_dat[[names(dat_lis)[ilis]]]$color=c(rep(colvec[ilis],nrow(plot_dat[[names(dat_lis)[ilis]]])-1),'black')
  plot_dat[[names(dat_lis)[ilis]]]$point_width=(plot_dat[[names(dat_lis)[ilis]]]$fetOR*3 )+point_width_scale
  }


   if(!is.na(p_thresh)){
    dummy=list()
    for(ilis in 1:length(plot_dat)){
      tempr=plot_dat[[names(dat_lis)[ilis]]][plot_dat[[names(dat_lis)[ilis]]]$FETp<p_thresh,]
      if(nrow(tempr)>0){
        dummy[[names(plot_dat)[ilis]]]=plot_dat[[names(plot_dat)[ilis]]][plot_dat[[names(plot_dat)[ilis]]]$FETp<p_thresh,]
      }
    }
    plot_dat=dummy
   }
options(warn=0)

#  if(!is.na(use_log)){
#    plot_dat$fetOR=log(as.numeric(plot_dat$fetOR),base=use_log)
#    plot_dat$lowerCI=log(as.numeric(plot_dat$lowerCI),base=use_log)
#    plot_dat$upperCI=log(as.numeric(plot_dat$upperCI),base=use_log)

#  }

gplots_dat=list()
  if(mode=='single'){
    for(idat in 1:length(plot_dat)){
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_dat[[idat]], aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data = plot_dat[[idat]], aes(y = module, x = fetOR), shape=15, cex=plot_dat[[idat]]$point_width, colour=plot_dat[[idat]]$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:x_lim[2])) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))

    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat[[names(plot_dat)[idat]]]=(myplot + theme_bw())
    }
  }

  if(mode=='combined'){
    plot_merge=plot_dat[[1]]
    if(length(plot_dat)>1){
        for(ilis in 2:length(plot_dat)){
          plot_merge=rbind(plot_merge,plot_dat[[ilis]])
        }}
    myplot=ggplot() +                                                     # initiate the plot space, saved to myplot variable for later display
      #  - subsequent additions / modifications of the plot can be done by manipulating this space
      #  geom_point(data = df, aes(y = Module, x = OR),colour = 'red', size = 3) +
      #  geom_errorbar(data=pdat,aes(x=Module,y=OR,ymin=lower.95..CI,ymax=upper.95..CI),colour="grey60",width=0.2)
      geom_vline(xintercept=1, linetype="dashed")+                         # add vertical line
      #    geom_vline(xintercept=(1), linetype="dashed")+                         # add vertical line
      geom_errorbarh(data=plot_merge, aes(y=module,x=fetOR, xmin=lowerCI, xmax=upperCI), colour="black", height=ci_bar_height,size=2)+  # add HORISONTAL error bars, geom_errorbar used for vertical..
      geom_point(data=plot_merge, aes(y = module, x = fetOR), shape=15, cex=plot_merge$point_width, colour=plot_merge$color)+    # add plotting points
      #    xlim(0, 5)                                                         # removes datapoints outside the range == scale_x_continuous(limits = c(-5000, 5000))
      #
      #    coord_cartesian(xlim = c(0, 60)) +                                  # moves the window only
      coord_cartesian(xlim = x_lim) +    #        coord_cartesian(xlim = c(-0.5, 1))                        # moves the window only
      scale_x_continuous(breaks=c(0,1:x_lim[2])) +
      #    scale_x_continuous(breaks=c(0,1,5,seq(10,2000,by=10))) +
      #  scale_size_area() +
      xlab("FET odds ratio and CI") +
      ylab("module") 
    ##  plot line styles      ||  http://www.cookbook-r.com/Graphs/Shapes_and_line_types/  
    #  ggtitle(names(table(dat$serie.name)[idat]))



    ## essential if using the ugly mess that is ggplot2, for more options : http://felixfan.github.io/rstudy/2013/11/27/ggplot2-remove-grid-background-margin/ 
    #  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #
      gplots_dat=(myplot + theme_bw())
  }
  return(gplots_dat)
}






bg.list<-function(module_list,bg_list){
# cat('\tUSE:\tcheck and, if necessary, generate a list of backgrounds the same length as mod_list based on bg_list class\n')

  if(class(bg_list)=='list'){
    if(length(bg_list)==length(module_list)){
      cat('\tusing multiple bg list as provided\n')
    }
    if(length(bg_list)!=length(module_list)){
      stop('\tlength of bg_list provided not the same length as module_list. if same bkgrnd for all modules, do not provide as.list()\n')
    }
  }

  if(class(bg_list)!='list'){
    cat('\tusing same bkgrnd for all modules, as provided (assumes vector)\n')
    
    bg_nw=list()
    for(imod in 1:length(module_list)){
      bg_nw[[names(module_list)[imod]]]=bg_list
    }
  }
  return(bg_nw)
}





net.cons<-function(alis,blis,abg=NA,bbg=NA,do_plots=T,p_thresh=0.01,main="Maximum % of overlap with p hyper > 0.01"){
# cat('\tUSE:\tcalculte overlaps between two lists')
# cat('\tUSE:\toptional - add backgrounds for bg.common() - not implemented')
#cat('\tNOTE:\tbg_list - options: 'NA' - one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n\n')

##===================================================================================================================================
## generate background as per options, if list provided (same length as modules)=================================

 options(warn=-1)

   if(!is.na(abg)){
    abg=bg.list(alis,abg)

   }
   if(!is.na(bbg)){
    bbg=bg.list(blis,bbg)
   }

 options(warn=0)

    results_mat=as.data.frame(matrix(NA,nrow=length(alis),ncol=length(blis)))
      rownames(results_mat)=names(alis)
      colnames(results_mat)=names(blis)


  if(class(abg)=='list'){
    cons_stat=list(
      pc_test=results_mat
      ,phyper_test=results_mat
      ,pc_repl=results_mat
      ,phyper_repl=results_mat
      )
  }

  if(class(abg)!='list'){
    cons_stat=list(
      pc_test=results_mat
      ,pc_repl=results_mat
      )
  }


  for(amod in 1:length(alis)){
    for(bmod in 1:length(blis)){
#       cat(amod,'  |  ',bmod,'\n')

      if(class(abg)=='list'){
        bg_comn=intersect(abg[[amod]],bbg[[bmod]])
#          cat('\t\t background overlap: ',round(length(bg_comn)/length(abg[[amod]]),digits=2),round(length(bg_comn)/length(bbg[[bmod]]),digits=2),'\n')
        alis[[amod]]=alis[[amod]][alis[[amod]]%in%bg_comn]
        blis[[bmod]]=blis[[bmod]][blis[[bmod]]%in%bg_comn]
        comn=intersect(alis[[amod]],blis[[bmod]])
          cat('\t\t gene overlap: ',names(alis[amod]),round(length(comn)/length(alis[[amod]]),digits=2),'  |  ',names(blis[bmod]),round(length(comn)/length(blis[[bmod]]),digits=2),'\n')
      }

      if(class(abg)!='list'){
        comn=intersect(alis[[amod]],blis[[bmod]])
          cat('\t\t gene overlap: ',names(alis[amod]),round(length(comn)/length(alis[[amod]]),digits=2),'  |  ',names(blis[bmod]),round(length(comn)/length(blis[[bmod]]),digits=2),'\n')
      }


      cons_stat$pc_test[names(alis[amod]),names(blis[bmod])]=round(length(comn)/length(alis[[amod]]),digits=2)
      cons_stat$pc_repl[names(alis[amod]),names(blis[bmod])]=round(length(comn)/length(blis[[bmod]]),digits=2)

      if(class(abg)=='list'){
        cons_stat$phyper_test[names(alis[amod]),names(blis[bmod])]=phyper(
                                                q=length(comn)                      # n.success.sample
                                                ,m=sum(blis[[bmod]] %in% bg_comn)   # n.success.popn
                                                ,n=length(bg_comn)                  # n.population
                                                ,k=length(alis[[amod]])             # n.sample
                                                ,lower.tail=F)

        cons_stat$phyper_repl[names(alis[amod]),names(blis[bmod])]=phyper(
                                                q=length(comn)                      # n.success.sample
                                                ,m=sum(alis[[amod]] %in% bg_comn)   # n.success.popn
                                                ,n=length(bg_comn)                  # n.population
                                                ,k=length(blis[[bmod]])             # n.sample
                                                ,lower.tail=F)
      }
    }
  }

  cons_stat$pc_max=pmax(cons_stat$pc_test,cons_stat$pc_repl)
  if(class(abg)=='list'){
  cons_stat$phyper_min=pmin(cons_stat$phyper_test,cons_stat$phyper_repl)
  }

  if(do_plots){
    library(corrplot)

    if(class(abg)!='list'){
      corrplot(make.numeric(cons_stat$pc_max),method = "circle",is.corr=FALSE
        ,cl.lim=c(0, 1),col=rev(c(colrb)),cl.align="l",tl.col="black"
        ,mar=c(0,0,4,0)
        ,title=main
        )    
    }

    if(class(abg)=='list'){
      corrplot(make.numeric(cons_stat$pc_max),method = "circle", is.corr=FALSE
        ,p.mat = make.numeric(cons_stat$phyper_min), insig = "blank",sig.level = p_thresh
        ,cl.lim=c(0, 1),col=rev(c(colrb)),cl.align="l",tl.col="black"
        ,mar=c(0,0,4,0)
        ,title=main
        )
    }
  }
  return(invisible(cons_stat))
}




net.overlap<-function(alis,abg=NA,do_plots=T,...){
# cat('\t USE:\t calculte overlaps between two lists')
# cat('\t USE:\t optional - add backgrounds for bg.common() - not implemented')
#cat('\t NOTE:\t bg_list - options: 'NA' - one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n\n')

##===================================================================================================================================
## generate background as per options, if list provided (same length as modules)=================================

 options(warn=-1)
   if(!is.na(abg)){
    abg=bg.list(alis,abg)
   }
 options(warn=0)

    results_mat=as.data.frame(matrix(NA,nrow=length(alis),ncol=length(alis)))
      rownames(results_mat)=names(alis)
      colnames(results_mat)=names(alis)

    stats_mat=as.data.frame(matrix(NA,nrow=length(alis),ncol=3))
      rownames(stats_mat)=names(alis)
      colnames(stats_mat)=c('n.genes','n.overlap')

  if(class(abg)=='list'){
    cons_stat=list(
      pc=results_mat
      ,n.overlap=results_mat
      ,phyper=results_mat
      )
  }

  if(class(abg)!='list'){
    cons_stat=list(
      pc=results_mat
      ,n.overlap=results_mat
      )
  }

  for(amod in 1:length(alis)){
    for(bmod in 1:length(alis)){

      if(class(abg)=='list'){
        bg_comn=intersect(abg[[amod]],bbg[[bmod]])
          cat('\t\t background overlap: ',round(length(bg_comn)/length(abg[[amod]]),digits=2),round(length(bg_comn)/length(bbg[[bmod]]),digits=2),'\n')
        humpty=alis[[amod]][alis[[amod]]%in%bg_comn]
        dumpty=alis[[bmod]][alis[[bmod]]%in%bg_comn]

        comn=intersect(alis[[amod]][alis[[amod]]%in%bg_comn],alis[[amod]][alis[[amod]]%in%bg_comn])
      }

      if(class(abg)!='list'){
        comn=intersect(alis[[amod]],alis[[bmod]])
          cat('\t\t gene overlap: ',names(alis[amod]),round(length(comn)/length(alis[[amod]]),digits=2),'   \t|  ',names(alis[bmod]),round(length(comn)/length(alis[[bmod]]),digits=2),'\n')
        humpty=alis[[amod]]
        dumpty=alis[[bmod]]
      }


      cons_stat$pc[names(alis[amod]),names(alis[bmod])]=round(length(comn)/length(alis[[amod]]),digits=2)
      cons_stat$pc[names(alis[bmod]),names(alis[amod])]=round(length(comn)/length(alis[[bmod]]),digits=2)

      cons_stat$pc[names(alis[amod]),names(alis[bmod])]=round(length(comn)/length(alis[[amod]]),digits=2)
      cons_stat$pc[names(alis[bmod]),names(alis[amod])]=round(length(comn)/length(alis[[bmod]]),digits=2)


      if(class(abg)=='list'){
        cons_stat$phyper[names(alis[amod]),names(alis[bmod])]=phyper(
                                                q=length(comn)                              # n.success.sample
                                                ,m=sum(alis[[bmod]] %in% bg_comn)   # n.success.popn
                                                ,n=length(bg_comn)                  # n.population
                                                ,k=length(alis[[amod]])             # n.sample
                                                ,lower.tail=F)

        cons_stat$phyper[names(alis[bmod]),names(alis[amod])]=phyper(
                                                q=length(comn)                              # n.success.sample
                                                ,m=sum(alis[[amod]] %in% bg_comn)   # n.success.popn
                                                ,n=length(bg_comn)                  # n.population
                                                ,k=length(alis[[bmod]])             # n.sample
                                                ,lower.tail=F)
      }
    }
  }

  if(do_plots){
    library(corrplot)
    diag(cons_stat$pc)=0.000001
    corrplot(make.numeric(cons_stat$pc),p.mat=make.numeric(cons_stat$pc)*100,sig.level=0.001,col=rev(c(colrb)),cl.align="l",tl.col="black",method='circle',is.corr=F,insig='p-value',...)#,cl.lim=c(0,70)
    diag(cons_stat$pc)=1
  }
  return(invisible(cons_stat))
}



cyt.connect<-function(cor_mat,thresh=0.001,use_pcor=T){
#	USE: \t- convert correlation matrix to cytoscape connections
#	NOTE:\t- pcor=T - use Aracne can to calculate partial correlations
#	NOTE:\t- thresh=0.001, pcor=F - use correlation matrix to get connections that pass the R2 threshold
	if(nrow(cor_mat)!=ncol(cor_mat)){
		stop('\n\tERROR: expect a square correlation matrix \n\n')
	}

	if(use_pcor){
		 cat('\tcalculating partial correlations using Aracne PMID: 16723010\n')
		library(minet)
		cor_mat=aracne(abs(cor_mat))
	}

	for(irow in 1:nrow(cor_mat)){
		humpty=as.data.frame(cor_mat[irow,][cor_mat[irow,]>=thresh])
			colnames(humpty)='pcor'
			humpty$a=rownames(cor_mat)[irow]
			humpty$b=rownames(humpty)
#	cat(nrow(humpty),'\n')
#	print(irow)
		if(irow==1){dumpty=humpty}
		if(irow>1){dumpty=rbind(dumpty,humpty)}

#print(nrow(dumpty))
	}
	dumpty=dumpty[dumpty$a!=dumpty$b,]	##  remove genes interacting with themselves
	rownames(dumpty)=1:nrow(dumpty)
	cat('\t',length(unique(c(dumpty$a,dumpty$b))),'of',nrow(cor_mat),' : ',(length(unique(c(dumpty$a,dumpty$b)))/nrow(cor_mat))*100,'% of genes have connections','\n')

## filtering to remove duplicate connections ie a~b, b~a

	dtb=dumpty
	rm(humpty)
	rm(dumpty)
	## more elegant way to do this is first sort the genes then paste them..........
	  dtb$x=paste(dtb$a,dtb$b)
	  dtb$y=paste(dtb$b,dtb$a)
	  dumpty=list()
	 k=nrow(dtb)
	cat('\tremoving duplicate connectoins\n')
	while(nrow(dtb)>0){
	  humpty=dtb[which(dtb$x==dtb$x[1] | dtb$x==dtb$y[1] | dtb$y==dtb$y[1] | dtb$y==dtb$y[1]), ]
	  dtb=dtb[-which(dtb$x==dtb$x[1] | dtb$x==dtb$y[1] | dtb$y==dtb$y[1] | dtb$y==dtb$y[1]), ]

		dumpty[[paste(sort(unique(c(humpty$a,humpty$b))),collapse='_')]]=t(as.data.frame(list(
		      a=sort(unique(c(humpty$a,humpty$b)))[1]
		      ,b=sort(unique(c(humpty$a,humpty$b)))[2]
		      ,pcor=max(humpty$pcor)
		  )))
		  	cat(round(1-nrow(dtb)/k,digits=2),'\r');flush.console()
	}
	cat('\n')
	return(invisible(as.data.frame(t(as.data.frame(dumpty)))))

	## if partial correl does not give all genes connections, it is feasible to add a single connection to missing genes based on max of correl
}


install.bioc<-function(package_name){
	source("https://bioconductor.org/biocLite.R")
	biocLite(package_name)
}


install.dependencies<-function(){
##  USE : install all packages used by one of the functions (yes these can indeed be specified as dependencies, not implemented yet)
##  NOTE: example error below means the mirror used is not functional/unavailable, try another..

# Warning: unable to access index for repository https://mirrors.ebi.ac.uk/CRAN/src/contrib:
#   cannot download all files
# Warning: unable to access index for repository https://mirrors.ebi.ac.uk/CRAN/bin/macosx/mavericks/contrib/3.3:
#   cannot download all files
# Warning message:
# package ‘gplots’ is not available (for R version 3.3.0) 

	system('git clone https://github.com/jalvesaq/colorout.git')
	system('sudo R CMD INSTALL colorout')
	library(colorout)

	install.packages('devtools')
	install.packages('gplots')

##  installing WGCNA - requires dependencies first - simplest is to use instruction as per website :
##  http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall
	source("http://bioconductor.org/biocLite.R") 
	biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
	install.packages("WGCNA")

	install.packages('MetaDE')
	install.packages('mixtools')
	install.packages('pamr')
	install.packages("psych")
	install.packages('corrplot')

	install.packages('gplots')
	install.packages('ggplot2')
	install.packages('pvclust')

	install.bioc('minet')
	install.bioc('limma')
	install.bioc('biomaRt')
}



annot_combine<-function(expr_mat,annot_mat,annot_from,annot_to){
##  process expression matrix with corresponding annotation matrix
##  - re-mapping id types and taking median of non-uniquly mapping ids
##  - annot_from & annot_to - expect the name of colname containing current and new ids respectively


##  identify multiple id to gene mappings and take median expression
	dupli=get.duplicates(annot_mat,'name')
	duplg=names(dupli$n.duplicated)

	 cat('\n\tobserved correlations between "probes" mapping to the same "gene"\n')
	dup=list()
	for(igen in duplg){
		humpty=expr_mat[annot_mat[annot_mat[,annot_to]==igen,annot_from],]
		dummy=cor(t(humpty))
		cat('\t',igen,'        \t',round(dummy[1,2:ncol(dummy)],digits=3),'\n')

		dumpty=(as.data.frame(apply(humpty,2,median)))
#			colnames(dumpty)=igen
		dup[[igen]]=dumpty

	}

	dup=t(as.data.frame(dup))
		rownames(dup)=duplg

##  re-map ids for 'unique' ids to gene mapping
	ids=intersect(annot_mat[!(annot_mat[,annot_to]%in%duplg),annot_from],rownames(expr_mat))
	unc=expr_mat[ids,]
		rownames(annot_mat)=annot_mat[,annot_from]
	annot_mat=annot_mat[ids,]
		rownames(unc)=annot_mat[,annot_to]


		cat('\n\tuniquely mapping ids:',nrow(unc),', ',round(nrow(unc)/nrow(expr_mat),digits=2)*100,'%    multiple mapping ids :',nrow(dup),', ',round(nrow(dup)/nrow(expr_mat),digits=2)*100,'%\n')

	expr_out=rbind(unc,dup)
		cat('\tfinal output contains',nrow(expr_out),', ',round(nrow(expr_out)/nrow(expr_mat),digits=2)*100,'%\n')
	return(invisible(expr_out))

}





lcount<-function(x,length){
	cat(round(x/length,digits=2),"\r");flush.console()
	return(x+1)
}




overlap<-function (A, B){
    both = union(A, B)
    inA = both %in% A
    inB = both %in% B
    return(table(inA, inB))
}



Intersect <- function(a,b,...){
##  full credit to Abhishek K Baikady at http://stackoverflow.com/questions/3695677/how-to-find-common-elements-from-multiple-vectors
  Reduce(intersect, list(a,b,...))
}

Union <- function(a,b,...){
##  full credit to Abhishek K Baikady at http://stackoverflow.com/questions/3695677/how-to-find-common-elements-from-multiple-vectors
  Reduce(union, list(a,b,...))
}



idconvert<-function(ids,verbose=T){
	ids=toupper(ids)

	if(length(unique(ids))<length(ids)){
		warning('\tWARNING : some ids are not unique, ',length(unique(ids)),' of ',length(ids),' are unique after "toupper" case conversion\n')
	}
	if(verbose){
	cat('\tnum genes recognised : ',sum(ids%in%idmap$ids),', ',frac(sum(ids%in%idmap$ids),length(ids),num=T)*100,'%\n',sep='')
	}
	idmap$ids=toupper(idmap$ids)
#	idmap[idmap$ids%in%ids,c('gene','ids','name')]

	return(idmap[idmap$ids%in%ids,c('gene','ids','name')])
#	for(igen )

}


idconvert.ensg<-function(dat_lis){
##  conveting list of vectors containing only official HUGO gene names (rownames) to ENSG

	for(ilis in names(dat_lis)){
		dat_lis[[ilis]]=idmap.ensg[idmap.ensg$gene%in%(dat_lis[[ilis]]),'ids']
	}
	return(dat_lis)
}





check.data<-function(dat_lis){
  dummy=list()
  for(ilis in names(dat_lis)){
    dummy[[ilis]]=rownames(dat_lis[[ilis]])
    cat('\t',ilis,'\tnrow=',length(dummy[[ilis]]),'\tncol=',ncol(dat_lis[[ilis]]),'\n')
  }

  stata=list()
  for(ilis in 1:length(dummy)){
    if(ilis==1){
      holder=dummy[[names(dummy)[ilis]]]
    }
    if(ilis>1){
      stata[[paste(names(dummy)[ilis-1],names(dummy)[ilis],sep='_')]]=sum(holder==dummy[[names(dummy)[ilis]]])/length(holder)
      cat('\t\t% genes in same order between',names(dummy)[ilis-1],'&',names(dummy)[ilis],'wrt',names(dummy)[ilis-1],stata[[paste(names(dummy)[ilis-1],names(dummy)[ilis],sep='_')]],'\n')
    }
  }
}


applydiffcoex <- function(beta2,corr_mat_list,signtype=signType){
  correl=vector(mode="list", length=length(corr_mat_list));
  compDij=vector(mode="list", length=length(corr_mat_list));
  #compute the cij0 (for all the seven conditions)
  for (k in 1:length(corr_mat_list)){
    #we convert each element in the list in a matrix
    datmat=as.matrix(corr_mat_list[[k]]) 
    #compute sign(corr)*(corr)^2 and store it in a correlation vector
    correl[[k]]= sign(datmat)*(datmat)^2
    diag(correl[[k]])=0
  }
  
  #we compute cij0 (reduce applies the binary function to all the elements in the list).
  cij0 = Reduce("+",correl)/length(correl);
  #compute Dij
  for (element in 1:length(correl)){
    compDij[[element]]=abs(correl[[element]]- cij0)/2; 
  }
  Dij=(1/(length(corr_mat_list)-1))*Reduce("+", compDij)^(beta2/2)
  dissTOM=TOMdist(Dij, TOMType = signtype);
  collectGarbage()
  return(dissTOM)
}




wgcna.diffcoex<-function(list_expr,pow=5,minModuleSize=40,mergeHeight=0.15,datDescr="",signType="unsigned"){
  print("  NOTE : Input data (list_expr) is expected as a list with each entry : rows = genes, columns = samples")
  library('WGCNA')
  enableWGCNAThreads()
  print("Can deal with only one softpower (pow) at a time ")

  set.seed(0)       # reproducibility 
  bicorL <- list()
  
  
  for(ireg in 1:length(list_expr)){
    print("■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■")
    print(paste(names(list_expr)[ireg],ireg,"of",length(names(list_expr))))
      t0=Sys.time()
      COND<- list_expr[[ireg]]
      #add a line to substract mean of gene expression row by row in each condition
#      CONDav <- scale(COND,scale=F)		##  the data is likely scaled already (ie if removed a covariate etc)
      bicorL[[ireg]]<- bicor(t(as.matrix(COND))) # iteration to adapt to the seq chosen
  }
  
  names(bicorL)<-names(list_expr)
  collectGarbage()

    t0=Sys.time()
    softPower=pow
    print(softPower)

    dissTOM<-applydiffcoex(softPower,bicorL,signtype=signType)
    print(Sys.time()-t0) #41.17 min
    collectGarbage()
  geneTree = hclust(as.dist(dissTOM), method = "average")
  #print(Sys.time()-t0)
  collectGarbage()
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                        method="hybrid",cutHeight=.996,
                        deepSplit = T, pamRespectsDendro = FALSE,
                        minClusterSize = minModuleSize)
  dynamicColors = labels2colors(dynamicMods)
#  print(paste("dynamicColors =",length(unique(dynamicColors)),unique(dynamicColors)))
  #print(Sys.time()-t0)
  collectGarbage()

  Datall <- t(as.data.frame(list_expr))

  collectGarbage()
  mergedColor<-mergeCloseModules(Datall,dynamicColors,cutHeight=mergeHeight)$color
  print(paste("mergedColor =",length(unique(mergedColor)),unique(mergedColor)))
  print(Sys.time()-t0)
  collectGarbage()

  #' Plot the dendrogram and colors underneath
#  pdf(file=paste(outDir,"/plots/",datDescr,"_power",softPower,"_dendogram.pdf",sep=""),width=10,paper = 'a4r')
#  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColor), "Hybrid Tree Cut",
#    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, 
#    guideHang = 0.05,main = paste("Cluster Dendrogram ",datDescr," (power ",softPower,")",sep=""))
#  dev.off()

#  print(Sys.time()-t0)
    mstat=as.data.frame(table(mergedColor))
    mstat=mstat[order(mstat[,2],decreasing=T),]
    msta0=mstat[(mstat[,1]=='grey'),]
    msta0$module='M0'
    mstat=mstat[!(mstat[,1]=='grey'),]
    mstat$module=paste0('M',1:nrow(mstat))
    mstat=rbind(mstat,msta0)
      colnames(mstat)=c('color','ngenes','module')
    mstat$color=as.character(mstat$color)

if(datDescr!=''){mstat$module=paste(mstat$module,dat_descr,sep="_")}   ##  add info after module name


## module membership ------------------------------------------
  module_list=list()
    for(imod in 1:nrow(mstat)){
        module_list[[mstat$module[imod]]]=rownames(list_expr[[1]])[mergedColor==mstat$color[imod]]        ##  assuming all module gene names are in the same order, if not, modules are probably unreliable anyhow
    }
    module_list[['bkgrnd']]=rownames(list_expr[[1]])


  module_expr=list()
  for(ilis in 1:length(list_expr)){
## module expression ------------------------------------------
    for(imod in 1:length(module_list)){
        module_expr[[names(list_expr)[ilis]]][[names(module_list)[imod]]]=list_expr[[names(list_expr)[ilis]]][module_list[[names(module_list)[imod]]],]
    }
  }


## add bkgrnd info to mstat
    mbg=as.data.frame('bkgrnd')
    mbg$length=nrow(list_expr[[1]])
    mbg$name='bkgrnd'

      colnames(mbg)=colnames(mstat)
  mstat=rbind(mstat,mbg)

  readme='\n\tModules are named based on size M1 - biggest, M0 - unclustered, bkgrnd - all input genes, output contains :
    \t1. module_list - list containing names of genes in each module
    \t2. module_expr - expression matrix of all genes in module / input dataset
    \t3. mstat       - key used to name modules, includes module size
    \t4. GeneTree     - object to plot the WGCNA style dendrogram
    \n'
  cat(readme)
  return(invisible(list(module_list=module_list,module_expr=module_expr,mstat=mstat,plotobj=geneTree,readme=readme)))

}






cmap.meta<-function(lmod,bkg,de_thresh=0.01,n_genes=5){  ## combine the stuffs below to use with function rather than combined stuff as is atm
cat('\n\tNOTE: this function requires two objects:	"metsum" & "degen", available from:\nhttps://www.dropbox.com/s/xjg3xpyxwjgodyw/DTB.full_info.sig.randM.fisher.DE_genes_single.Rdata?dl=0\n\n')
####   input format :
##> str(bkg)
# chr [1:13210] "ENSG00000121410" "ENSG00000175899" "ENSG00000166535" ...
#> str(lmod)
#List of 4
# $ dif.hubs :List of 1
#  ..$ up: chr [1:130] "ENSG00000198826" "ENSG00000105011" "ENSG00000066279" "ENSG00000087586" ...
# $ E2F1.hubs:List of 1
#  ..$ up: chr [1:66] "ENSG00000198826" "ENSG00000066279" "ENSG00000087586" "ENSG00000178999" ...
# $ E2F1.M5  :List of 2
#  ..$ down: chr [1:12] "ENSG00000127837" "ENSG00000159842" "ENSG00000149925" "ENSG00000100307" ...
#  ..$ up  : chr [1:101] "ENSG00000198826" "ENSG00000066279" "ENSG00000156802" "ENSG00000176208" ...
# $ M5       :List of 2
#  ..$ down: chr [1:216] "ENSG00000127837" "ENSG00000165029" "ENSG00000127220" "ENSG00000159842" ...
#  ..$ up  : chr [1:262] "ENSG00000159251" "ENSG00000198826" "ENSG00000105011" "ENSG00000066279" ...

##  - lmod - list of modules - for each module 2 lists of gene ids "ENSG" - "up" & "down" - genes up/down-regulated between treatment and control
#de_thresh=0.01
#n_genes=5
mstat=list()
sigen=list()
k=1

dkey=gsub('DE_Drug_(.*)_Cell_.*_Array_.*_Conc_.*_Conc_.*_time_.*_res','\\1',names(degen))
#	matst(names(metsum)%in%dkey)
ukey=unique(dkey)

#idru='wortmannin'
cat('\n\tperform enrichment test for module genes in cmap\n')
	for(idru in ukey){

		if(idru%in%names(metsum)){
			dummy=metsum[[idru]]
				rownames(dummy)=gsub('_at','',rownames(dummy))
			dummy=dummy[rownames(dummy)%in%bkg,]

			dummy$log.fc=dummy$logFC.mean
			dummy$adj.p=dummy$randmMod.fdr
			dummy=dummy[,c('log.fc','adj.p')]



		holder=list()
	#		holder$n.bg=nrow(dummy)
	#		holder$n.bg.all=nrow(dummy)
			holder$n.bg.sig=rownames(dummy[dummy$adj.p<de_thresh,])
			holder$n.bg.sig_up=rownames(dummy[dummy$adj.p<de_thresh & dummy$log.fc>0,])
			holder$n.bg.sig_down=rownames(dummy[dummy$adj.p<de_thresh  & dummy$log.fc<0,])

			sigen[[paste0('mixedM_',idru)]]$bkg=(holder)

		if(length(holder$n.bg.sig)>n_genes){

	#	if(idru%in%names(metsum)){mstat$meta[[idru]]$bg=as.data.frame(humpty)}
	#	if(!idru%in%names(metsum)){mstat$single[[idru]]$bg=as.data.frame(humpty)}

		for(imod in names(lmod)){
			if('up' %in% names(lmod[[imod]])){

				dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$up,]
				holder[[imod]]$up_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
				holder[[imod]]$up_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

				mstat[[paste0('mixedM_',idru)]][[paste(imod,'up_downreg',sep='.')]]=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
				mstat[[paste0('mixedM_',idru)]][[paste(imod,'up_upreg',sep='.')]]  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
			}

			if('down' %in% names(lmod[[imod]])){

				dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$down,]
				holder[[imod]]$down_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
				holder[[imod]]$down_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

				mstat[[paste0('mixedM_',idru)]][[paste(imod,'down_downreg',sep='.')]]=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
				mstat[[paste0('mixedM_',idru)]][[paste(imod,'down_upreg',sep='.')]]  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
			}
			sigen[[paste0('mixedM_',idru)]][[paste(imod,'down_upreg',sep='.')]]=(holder[[imod]])
		}
		rm(dumpty)
		}
		k=lcount(k,length(ukey))
		rm(dummy)
		rm(holder)

		}
		if(!(idru%in%names(metsum))){
			dummy=degen[dkey==idru]
			if(length(dummy)>1){		##  if removing the above constraint will need to introduce another loop for each experiment within same drug
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

		holder=list()
	#		holder$n.bg=nrow(dummy)
	#		holder$n.bg.all=nrow(dummy)
			holder$n.bg.sig=rownames(dummy[dummy$adj.p<de_thresh,])
			holder$n.bg.sig_up=rownames(dummy[dummy$adj.p<de_thresh & dummy$log.fc>0,])
			holder$n.bg.sig_down=rownames(dummy[dummy$adj.p<de_thresh  & dummy$log.fc<0,])

			sigen[[paste0('single_',idru)]]$bkg=(holder)


		if(length(holder$n.bg.sig)>n_genes){

	#	if(idru%in%names(metsum)){mstat$meta[[idru]]$bg=as.data.frame(humpty)}
	#	if(!idru%in%names(metsum)){mstat$single[[idru]]$bg=as.data.frame(humpty)}

		for(imod in names(lmod)){
			if('up' %in% names(lmod[[imod]])){

				dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$up,]
				holder[[imod]]$up_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
				holder[[imod]]$up_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

				mstat[[paste0('single_',idru)]][[paste(imod,'up_downreg',sep='.')]]=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
				mstat[[paste0('single_',idru)]][[paste(imod,'up_upreg',sep='.')]]  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
			}

#sampl=rownames(dumpty)
#bkgrnd=rownames(dummy)
#success=holder$n.bg.sig_down

#success=holder$n.bg.sig_up

#table(sampl%in%bkgrnd)
#table(success%in%sampl)
#table(success%in%bkgrnd)


			if('down' %in% names(lmod[[imod]])){

				dumpty=dummy[rownames(dummy)%in%lmod[[imod]]$down,]
				holder[[imod]]$down_downreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc<0,]
				holder[[imod]]$down_upreg=dumpty[dumpty$adj.p<de_thresh & dumpty$log.fc>0,]

				mstat[[paste0('single_',idru)]][[paste(imod,'down_downreg',sep='.')]]=as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_down,counts=F))
				mstat[[paste0('single_',idru)]][[paste(imod,'down_upreg',sep='.')]]  =as.data.frame(fet(sampl=rownames(dumpty),bkgrnd=rownames(dummy),success=holder$n.bg.sig_up,counts=F))
			}
			sigen[[paste0('single_',idru)]][[paste(imod,'down_upreg',sep='.')]]=(holder[[imod]])
		}
#		rm(dumpty)
		}
		k=lcount(k,length(ukey))
#		rm(dummy)
#		rm(holder)

		}
	}


#	str(mstat[['mixedM_wortmannin']])
#	str(sigen[['mixedM_wortmannin']])
#	(sigen[['mixedM_wortmannin']][['E2F1.M5.down_upreg']]$down_upreg)


	sigdru=list()
	sigsum=list()
	sigpc1=list()
	sigpc2=list()
cat('\n\tcompiling results..\n\n')
	for(idru in names(mstat)){
		humpty=mstat[[idru]]
		dumpty=as.data.frame(matrix(unlist(humpty),nrow=length(humpty),byrow=T))
			rownames(dumpty)=names(humpty)
			colnames(dumpty)=colnames(humpty[[1]])

		dumpty$n_hits_pc=dumpty$samp.success/(dumpty$samp.success+dumpty$bkgrnd.success)
		dumpty$n_mod_pc=dumpty$samp.success/(dumpty$n.genes)

	#	if(sum(dumpty$FETp<0.05)>=1){
			sigdru[[idru]]=dumpty
			sigsum[[idru]]=(dumpty[,c('FETp'),drop=F])
			sigpc1[[idru]]=(dumpty[,c('n_hits_pc'),drop=F])
				colnames(sigsum[[idru]])=paste(as.character(idru))#,colnames(sigsum[[idru]]))
				colnames(sigpc1[[idru]])=paste(as.character(idru))#,colnames(sigpcs[[idru]]))

			sigpc2[[idru]]=(dumpty[,c('n_mod_pc'),drop=F])
				colnames(sigpc2[[idru]])=paste(as.character(idru))#,colnames(sigpcs[[idru]]))



##				colnames(sigsum[[idru]])=paste(as.character(idru),colnames(sigsum[[idru]]),sep='.')	#paste0('Drug_',as.character(idru))
	#	}

	}

	sigsum=t(as.data.frame(sigsum))
	sigpc1=t(as.data.frame(sigpc1))
	sigpc2=t(as.data.frame(sigpc2))

#sigsum['single_etoposide',]
#sigsum['single_monobenzone',]
#sigsum['single_trifluridine',]

#sigpcs['single_etoposide',]
#sigpcs['single_monobenzone',]
#sigpcs['single_trifluridine',]


#mstat['single_dequalinium chloride']

####  M2.up_downreg
## inf## 'single_dequalinium.chloride'
## NaN## 'single_dicycloverine'
#		Head(sigsum)

	readme='\tcmap differentially expressed genes from drug treatment, raw and meta-analysis
		\t\tNOTE :  currently individual experiments used to create the meta-analysis results are ignored
		\t\t$sigsum   - summary results - p-values of enrichment only
		\t\t$sigpc1   - summary results - % of differentially expressed genes in module / background (n.success.module/n.success.bkg)
		\t\t$sigpc2   - summary results - % of differentially expressed genes in module / all module genes (n.success.module/n.genes.module)
		\t\t$sigdru   - summary results - summary fisher exact test statistics
		\t\t$sigen    - lists of genes used to calculate fisher exact test in mstat
'
cat(readme)
	return(invisible(list(sigsum=(sigsum),sigpc1=(sigpc1),sigpc2=(sigpc2),sigdru=sigdru,sigen=sigen,readme=readme,dumpty=dumpty)))
}


