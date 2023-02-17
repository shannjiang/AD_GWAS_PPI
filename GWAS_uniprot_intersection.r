library(GenomicRanges)
library(stringr)
library(tidyr)
library(rtracklayer)
uniprot_dir = '/home/shann/Documents/database/uniprot/'
gwas_dir = '/home/shann/Documents/GWAS_sumstats/AD/'
out_dir = '/home/shann/Documents/AD_GWAS_PPI_project/'
chain_file = '/home/shann/Documents/R_package/rtracklayer/liftover_chain/hg19ToHg38.over.chain'
chain = import.chain(chain_file)

uniprot_var_df = read.table(file = paste0(uniprot_dir,'homo_sapiens_variation.txt.gz'), skip = 161, header = T, sep = '\t')
uniprot_var_df = uniprot_var_df[-1,]
#df[!grepl("REVERSE", df$Name),]
#uniprot_var_df = uniprot_var_df[!grepl('dup',uniprot_var_df$Chromosome.Coordinate),]
uniprot_var_coord_df = str_split_fixed(uniprot_var_df$Chromosome.Coordinate,'\\.',3)

uniprot_var_coord_df = as.data.frame(uniprot_var_coord_df)
uniprot_var_df$chr = uniprot_var_coord_df$V1
chrNC = c(paste0('NC_00000',c(1:9)),paste0('NC_0000',c(10:24)),'NC_012920')
uniprot_var_df = uniprot_var_df[uniprot_var_df$chr %in% chrNC,]

uniprot_var_coord_df = str_split_fixed(uniprot_var_df$Chromosome.Coordinate,'\\.',3)
uniprot_var_coord_df = as.data.frame(uniprot_var_coord_df)
uniprot_var_coord_df$V1 = gsub('NC_0+','',uniprot_var_coord_df$V1)
colnames(uniprot_var_coord_df)[1] = 'chr'
uniprot_var_coord_df$chr = gsub('23','X',uniprot_var_coord_df$chr)
uniprot_var_coord_df$chr = gsub('24','Y',uniprot_var_coord_df$chr)
uniprot_var_coord_df$chr = gsub('12920','M',uniprot_var_coord_df$chr)
uniprot_var_coord_df$V3 = gsub('_.*','',uniprot_var_coord_df$V3)
uniprot_var_coord_df$V3 = gsub('\\[.*','',uniprot_var_coord_df$V3)
uniprot_var_coord_df$pos = extract_numeric(uniprot_var_coord_df$V3)
uniprot_var_coord_df$A12 = gsub('^[0-9]+','',uniprot_var_coord_df$V3)
A1_A2_df = as.data.frame(str_split_fixed(uniprot_var_coord_df$A12,'>',2))
uniprot_var_df$chr = uniprot_var_coord_df$chr
uniprot_var_df$pos = uniprot_var_coord_df$pos
uniprot_var_df$refA = A1_A2_df$V1
uniprot_var_df$altA = A1_A2_df$V2
remove(uniprot_var_coord_df,A1_A2_df)
uniprot_var_df$chr = paste0('chr',uniprot_var_df$chr)
#length(which(uniprot_var_df$pos=='NA'))
#uniprot_var_df2 = uniprot_var_df[!uniprot_var_df$pos %in% 'NA',]
uniprot_var_df$pos = as.integer(uniprot_var_df$pos)
#z = uniprot_var_df[is.na(uniprot_var_df$pos),]
#uniprot_var_df = uniprot_var_df[!is.na(uniprot_var_df$pos),]
uniprot_var_gr = with(uniprot_var_df,GRanges(chr,IRanges(pos,pos,Source.DB.ID=Source.DB.ID,refA=refA,altA=altA,Gene.Name=Gene.Name,AC=AC,Variant.AA.Change=Variant.AA.Change,Consequence.Type=Consequence.Type,Cytogenetic.Band=Cytogenetic.Band,Chromosome.Coordinate=Chromosome.Coordinate,Ensembl.gene.ID=Ensembl.gene.ID,Ensembl.transcript.ID=Ensembl.transcript.ID,Ensembl.translation.ID=Ensembl.translation.ID,Evidence=Evidence)))

#jansen gwas
jansen_df = read.table(file = paste0(gwas_dir,'AD_Jansen_GWAS.assoc'), header = T, sep = ' ')
jansen_df$CHR = paste0('chr',jansen_df$CHR)
jansen_df$pos = with(jansen_df,ifelse(A1=='I' | A1=='D',BP+1, BP))
#convert gwas df to GRanges object
jansen_gr = with(jansen_df,GRanges(CHR,IRanges(pos,pos,SNP=SNP,A1=A1,A2=A2,P=P,OR=OR)))
