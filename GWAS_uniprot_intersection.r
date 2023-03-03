library(GenomicRanges)
library(stringr)
library(tidyr)
library(rtracklayer)
uniprot_dir = '/home/shann/Documents/database/uniprot/'
gwas_dir = '/home/shann/Documents/GWAS_sumstats/AD/'
out_dir = '/home/shann/Documents/AD_GWAS_PPI_project/'
chain_file = '/home/shann/Documents/references/liftover_chain/hg19ToHg38.over.chain'
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
jansen_uniprot_overlap_hit = findOverlaps(jansen_gr,uniprot_var_gr,type = 'equal')
jansen_uniprot_overlap_hit_df = as.data.frame(jansen_uniprot_overlap_hit)
jansen_overlap_part = jansen_df[jansen_uniprot_overlap_hit_df$queryHits,]
colnames(jansen_overlap_part) = paste0('gwas.',colnames(jansen_overlap_part))
uniprot_overlap_part = uniprot_var_df[jansen_uniprot_overlap_hit_df$subjectHits,]
colnames(uniprot_overlap_part) = paste0('uniprot.',colnames(uniprot_overlap_part))
jansen_uniprot_overlap_df = cbind(jansen_overlap_part,uniprot_overlap_part)
jansen_uniprot_overlap_df = jansen_uniprot_overlap_df[order(jansen_uniprot_overlap_df$gwas.P),]
write.csv(jansen_uniprot_overlap_df,file = paste0(out_dir,'jansen_uniprot_overlap_df.csv'))
jansen_uniprot_overlap_sugThre_df = jansen_uniprot_overlap_df[jansen_uniprot_overlap_df$gwas.P < 1e-5,]
write.csv(jansen_uniprot_overlap_sugThre_df,file = paste0(out_dir,'jansen_uniprot_overlap_sugThre_df.csv'))
jansen_uniprot_overlap_nomThre_df = jansen_uniprot_overlap_df[jansen_uniprot_overlap_df$gwas.P < 5e-2,]
write.csv(jansen_uniprot_overlap_nomThre_df,file = paste0(out_dir,'jansen_uniprot_overlap_nomThre_df.csv'))
#remove jansen related files to save memory
rm(list = ls(pattern = 'jansen'))

#kunkle gwas
kunkle_df = read.table(file = paste0(gwas_dir,'AD_Kunkle_GWAS.assoc'), header = T, sep = ' ')
kunkle_df$CHR = paste0('chr',kunkle_df$CHR)
#kunkle_df$pos = with(kunkle_df,ifelse(A1=='I' | A1=='D',BP+1, BP))
#convert gwas df to GRanges object
kunkle_gr = with(kunkle_df,GRanges(CHR,IRanges(BP,BP,SNP=SNP,A1=A1,A2=A2,P=P,OR=OR)))
kunkle_uniprot_overlap_hit = findOverlaps(kunkle_gr,uniprot_var_gr,type = 'equal')
kunkle_uniprot_overlap_hit_df = as.data.frame(kunkle_uniprot_overlap_hit)
kunkle_overlap_part = kunkle_df[kunkle_uniprot_overlap_hit_df$queryHits,]
colnames(kunkle_overlap_part) = paste0('gwas.',colnames(kunkle_overlap_part))
uniprot_overlap_part = uniprot_var_df[kunkle_uniprot_overlap_hit_df$subjectHits,]
colnames(uniprot_overlap_part) = paste0('uniprot.',colnames(uniprot_overlap_part))
kunkle_uniprot_overlap_df = cbind(kunkle_overlap_part,uniprot_overlap_part)
kunkle_uniprot_overlap_df = kunkle_uniprot_overlap_df[order(kunkle_uniprot_overlap_df$gwas.P),]
write.csv(kunkle_uniprot_overlap_df,file = paste0(out_dir,'kunkle_uniprot_overlap_df.csv'))
kunkle_uniprot_overlap_sugThre_df = kunkle_uniprot_overlap_df[kunkle_uniprot_overlap_df$gwas.P < 1e-5,]
write.csv(kunkle_uniprot_overlap_sugThre_df,file = paste0(out_dir,'kunkle_uniprot_overlap_sugThre_df.csv'))
kunkle_uniprot_overlap_nomThre_df = kunkle_uniprot_overlap_df[kunkle_uniprot_overlap_df$gwas.P < 5e-2,]
write.csv(kunkle_uniprot_overlap_nomThre_df,file = paste0(out_dir,'kunkle_uniprot_overlap_nomThre_df.csv'))

#Bellenguez
#bell_df = read.table(file = paste0(gwas_dir,'Bellenguez_GCST90027158_buildGRCh38_top10000.tsv.gz'), header = T, sep = '\t')
bell_df = read.table(file = paste0(gwas_dir,'Bellenguez_GCST90027158_buildGRCh38.tsv.gz'), header = T, sep = '\t')
bell_df$chromosome = paste0('chr',bell_df$chromosome)
bell_gr = with(bell_df,GRanges(chromosome,IRanges(base_pair_location,base_pair_location,SNP=variant_id,A1=other_allele,A2=effect_allele,P=p_value,OR=odds_ratio)))
bell_uniprot_overlap_hit = findOverlaps(bell_gr,uniprot_var_gr,type = 'equal')
bell_uniprot_overlap_hit_df = as.data.frame(bell_uniprot_overlap_hit)
bell_overlap_part = bell_df[bell_uniprot_overlap_hit_df$queryHits,]
colnames(bell_overlap_part) = paste0('gwas.',colnames(bell_overlap_part))
uniprot_overlap_part = uniprot_var_df[bell_uniprot_overlap_hit_df$subjectHits,]
colnames(uniprot_overlap_part) = paste0('uniprot.',colnames(uniprot_overlap_part))
bell_uniprot_overlap_df = cbind(bell_overlap_part,uniprot_overlap_part)
bell_uniprot_overlap_df = bell_uniprot_overlap_df[order(bell_uniprot_overlap_df$gwas.p_value),]
write.csv(bell_uniprot_overlap_df,file = paste0(out_dir,'bell_uniprot_overlap_df.csv'))
bell_uniprot_overlap_sugThre_df = bell_uniprot_overlap_df[bell_uniprot_overlap_df$gwas.P_value < 1e-5,]
write.csv(bell_uniprot_overlap_sugThre_df,file = paste0(out_dir,'bell_uniprot_overlap_sugThre_df.csv'))
bell_uniprot_overlap_nomThre_df = bell_uniprot_overlap_df[bell_uniprot_overlap_df$gwas.p_value < 5e-2,]
write.csv(bell_uniprot_overlap_nomThre_df,file = paste0(out_dir,'bell_uniprot_overlap_nomThre_df.csv'))

#wightman
wightman_df = read.table(file = paste0(gwas_dir,'Wightman_PGCALZ2sumstatsExcluding23andMe.txt.gz'), header = T, sep = '\t')
wightman_df$chr = paste0('chr',wightman_df$chr)
wightman_df$SNP = with(wightman_df,paste0(chr,'_',PosGRCh37,'_',otherAllele,'_',testedAllele))
wightman_gr = with(wightman_df,GRanges(chr,IRanges(PosGRCh37,PosGRCh37,SNP=SNP,A1=otherAllele,A2=testedAllele,P=p)))
#bim_hg38 = as.data.frame(liftOver(bim_gr,chain))
wightman_hg38_gr = liftOver(wightman_gr,chain)
wightman_hg38_gr = unlist(wightman_hg38_gr)
wightman_uniprot_overlap_hit = findOverlaps(wightman_gr,uniprot_var_gr,type = 'equal')
wightman_uniprot_overlap_hit_df = as.data.frame(wightman_uniprot_overlap_hit)
wightman_overlap_part = wightman_df[wightman_uniprot_overlap_hit_df$queryHits,]
colnames(wightman_overlap_part) = paste0('gwas.',colnames(wightman_overlap_part))
uniprot_overlap_part = uniprot_var_df[wightman_uniprot_overlap_hit_df$subjectHits,]
colnames(uniprot_overlap_part) = paste0('uniprot.',colnames(uniprot_overlap_part))
wightman_uniprot_overlap_df = cbind(wightman_overlap_part,uniprot_overlap_part)
wightman_uniprot_overlap_df = wightman_uniprot_overlap_df[order(wightman_uniprot_overlap_df$gwas.p),]
write.csv(wightman_uniprot_overlap_df,file = paste0(out_dir,'wightman_uniprot_overlap_df.csv'))
wightman_uniprot_overlap_sugThre_df = wightman_uniprot_overlap_df[wightman_uniprot_overlap_df$gwas.p < 1e-5,]
write.csv(wightman_uniprot_overlap_sugThre_df,file = paste0(out_dir,'wightman_uniprot_overlap_sugThre_df.csv'))
wightman_uniprot_overlap_nomThre_df = wightman_uniprot_overlap_df[wightman_uniprot_overlap_df$gwas.p < 5e-2,]
write.csv(wightman_uniprot_overlap_nomThre_df,file = paste0(out_dir,'wightman_uniprot_overlap_nomThre_df.csv'))
