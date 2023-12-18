#1. 整理出饮食数据，sas部分
#2. 做pca, 代码 diet_pca.R,改文件名，不能有缺失
#3. 做phewas, prepare icd10.csv, cov.csv and diet_pca.csv R & python
#4. 做pca和食物组的correlation
#5. 挑出关注的疾病，做pc和疾病的cox模型，sas部分
#6. 做pc和疾病prs的交互作用 目前没有发现有交互作用
#7. 做pc和全基因snps的交互作用

dos2unix xxx

/cluster/apps/R/4.2.2/bin/./Rscript --vanilla phewas.R pc11 pc1_all_results.csv pc1_p_results.csv
/cluster/apps/R/4.2.2/bin/./R
服务器的R路径: 
容器路径: singularity exec --bind /cluster:/cluster ~/gen_soft/MetaGenome.sif R 

#3. phewas
library(data.table)
library(dplyr)
library(PheWAS) 

# import pca file and icd code (ukb format)
nutrition <- read.csv("diet_pca.csv")
nutrition <- nutrition %>% rename(id = eid)

raw_data <- fread("raw.csv")
raw_data <- raw_data %>% rename(id = n_eid)

#merge and export，统一样本量
merge <- left_join(nutrition, raw_data, by = "id")
trans_icd <- merge %>% select(-pc1, -pc2, -pc3)
write.csv(trans_icd, file = "trans_icd.csv", row.names = FALSE)

#run.python, export icd1.csv transform 

#transfer icd
icd <- fread("icd1.csv")

icd2 <- icd %>% mutate(icd2 = case_when(icd == "icd10" ~ "ICD10CM", icd == "icd9" ~ "ICD9CM")) %>% rename( vocabulary_id = icd2, code = disease)%>% 
                filter(vocabulary_id == "ICD10CM")%>% select(id, vocabulary_id, code)  
   
icd2$code <- short_to_decimal(icd2$code)#change short to decimal
icd2$index <- 5

write.csv(icd2, file = "icd10.csv", row.names = FALSE)

#服务器上完成
library(data.table)
library(dplyr)
library(PheWAS)


icd4 <- fread("icd10.csv")
icd4$code <- as.character(icd4$code)
phenotypes <- createPhenotypes(icd4,min.code.count=2)


cov <- fread("cov.csv")
cov <- cov %>% rename(id = n_eid)

cov$sex<-as.factor(cov$sex)
cov$edu_level<-as.factor(cov$edu_level)
cov$race_group<-as.factor(cov$race_group)
cov$smoking<-as.factor(cov$smoking)
cov$drinking<-as.factor(cov$drinking)


diet_pca <- read.csv("diet_pca.csv")
diet_pca <- diet_pca %>% rename(id = eid)

diet_pca2 <-diet_pca %>% select(id,pc1)

data <- left_join(left_join(phenotypes,diet_pca2),cov)

results <- phewas(phenotypes=names(phenotypes)[-1], genotypes=names(diet_pca2)[-1], data = data, covariates=c("sex","age_attend","edu_level","bmi","race_group","townsend_score","drinking","smoking","sum_physical_met"), cores=20)
results_d <- addPhecodeInfo(results, alpha = 1)

write.csv(results_d, file = "test_pc1_results.csv", row.names = FALSE)

#filter
filter_results <-  results_d %>% filter(p < 0.0005 )


#4. correlation
nutrition <- read.csv("new_results/food_class.csv")
nutrition_pc <- read.csv("diet_pca.csv")
nutrition_pc <- nutrition_pc %>% rename(n_eid = eid)

merge <- left_join(nutrition, nutrition_pc, by = "n_eid")

test <- merge[, c(2:44, 45)]

result <- data.frame()
result_p <- data.frame()
for (i in 1:43) {
cor <- cor.test(x= test[, i], y = test[, 44]) 
cor_value_table <- as.data.frame(cor$estimate)
cor_p <- as.data.frame(cor$p.value)
result <- rbind(result, cor_value_table)
result_p <- rbind(result_p, cor_p)

result_final <- cbind(result, result_p)
}

result_final <- arrange(result_final, `cor$p.value`)

write.csv(result_final, file = "pc1_covalue.csv")


#5. 挑出关注的疾病，做pc和疾病的cox模型
#6. 做pc和疾病prs的交互作用 目前没有发现有交互作用，prs计算运行calculate_prs.sh

library("epiR")
library(interactionR)
library(dplyr)
library(data.table)

raw_data <- fread("inter_data.csv")

raw_data <- fread("inter_ffq.csv")


raw_data$sex <- factor(raw_data$sex) 
#raw_data$hyper_icd10 <- factor(raw_data$hyper_icd10) 
raw_data$drinking <- factor(raw_data$drinking) 
raw_data$edu_level <- factor(raw_data$edu_level) 
raw_data$smoking <- factor(raw_data$smoking) 
raw_data$race_group <- factor(raw_data$race_group)
raw_data$history_hyper <- factor(raw_data$history_hyper)

raw_data <- raw_data %>% mutate(pc1_adjust = ifelse(pc1 < -0.10, NA, pc1))
raw_data <- raw_data %>% mutate(sbp_274_adjust = ifelse(sbp_274 < -0.02 | sbp_274 > 0.02, NA, sbp_274))

raw_data <- raw_data %>% mutate(sbp_274_group = case_when(sbp_274 < -0.001044365 ~ 1, -0.001044365 <= sbp_274 & sbp_274 < 0.006343085 ~ 2, sbp_274 >= 0.006343085 ~ 3))
quantile(raw_data$sbp_274)

#pc
fit_pc <- glm(hyper_icd10 ~ pc11  + sex + drinking+edu_level+smoking+race_group+base_age+sum_physical_met,family= binomial,
           data = raw_data)
table <- as.data.frame(summary(fit_pc)$coefficients[c(2,3),]) 

#prs
fit_prs <- glm(hyper_icd10 ~  sbp_107 + sex+drinking+edu_level+smoking+race_group+base_age+sum_physical_met,family= binomial,
              data = raw_data)
table <- as.data.frame(summary(fit_prs)$coefficients[c(2,3),]) 



#pc*prs
fit <- glm(prevalent_diabetes ~ ffq_pc1_group + diabetes_prs + ffq_pc1_group:diabetes_prs + sex+drinking+edu_level+smoking+race_group+base_age+sum_physical_met+history_hyper,family= binomial,
            data = raw_data)


fit <- glm(prevalent_hyper ~ ffq_pc2_group + sbp_107 + ffq_pc2_group:sbp_107 + sex+drinking+edu_level+smoking+race_group+base_age+sum_physical_met+history_hyper,family= binomial,
           data = raw_data)


table <-as.data.frame(summary(fit)$coefficients[c(2,3,20),]) 
table(raw_data$prevalent_hyper)
epi.interaction(model = fit, param = "product", coef = c(2,3,20),
                conf.level = 0.95)

#7. 做pc和全基因snps的交互作用 plink
#记录：maf 0.005, 0.05,0.001都不行，没有显著的snp，插补分数大于0.6
for i in {1..4}; do
    plink2 --bgen ~/convert_new/ukb22828_c${i}_b0_v3.bgen \
           --sample ~/convert_new/ukb22828_c10_b0_v3_s487159.sample \
           --memory 15000 \
           --extract-col-cond ~/convert_new/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 \
           --maf 0.005 \
           --keep-fam new_id.txt \
           --export bgen-1.2 \
           --out chr${i}_qc
done

#test
    plink2 --bgen ~/convert/ukb_imp_chr21_v3.bgen \
           --sample ~/convert/ukb22828_c1_b0_v3_s487197.sample \
           --memory 15000 \
           --extract-col-cond ~/convert/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 \
           --maf 0.005 \
           --keep-fam sample.txt \
           --export bgen-1.2 \
           --out chr21_qc

#创建bgi文件
for i in {1..22}; do
  bgenix -g chr${i}_qc.bgen -index -clobber
done   

#设置哑变量，注意数parameter
for i in {1..22};do
plink2 --bgen chr${i}_qc.bgen ref-first --sample chr${i}_qc.sample --pheno pheno.txt --covar cov.txt --covar-variance-standardize --glm interaction --parameters 1-19 --out chr${i}_restuls
done 


#提取ADD, pc11, ADDxpc11
for i in {1..22}; do
grep 'ADD' chr${i}_results.mean_chole2.glm.linear > chr${i}_filter.txt

done

#计算每一个snp的maf
plink2 --bgen chr1_qc.bgen ref-first --sample chr1_qc.sample --freq --out chr1_qc_freq 

#做qc的时候keep-fam不能和SNP质控放到一起，因为是有顺序的

for i in {11..22};do
plink2 --bgen ~/convert_new/ukb22828_c${i}_b0_v3.bgen ref-first --sample ~/convert_new/ukb22828_c10_b0_v3_s487159.sample --memory 15000 --extract-col-cond ~/convert_new/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 --geno 0.02 --mind 0.02 --hwe 1e-6 --maf 0.05 --export bgen-1.2 --out chr${i}_qc 
done 


for i in {11..22};do
plink2 --bgen ../chr${i}_qc.bgen ref-first --sample ../chr${i}_qc.sample --memory 15000 --keep-fam ../subset_id.txt --export bgen-1.2 --out chr${i}_sample_qc 
done 
######################

#DP 做GWAS
#--geno 0.2 --mind 0.2 --hwe 1e-6 --maf 0.005 signficantly:23
#--hwe 1e-6 --maf 0.005 signficantly:25
#--geno 0.02 --mind 0.02 --hwe 1e-6 --maf 0.005 signficantly:12
#for i in {1..100};do echo "nohup plink --file b  --allow-no-sex --pheno mphe.txt --linear --out y_${i}_result --mpheno $i "|bash;done 多性状参数

for i in {9..22};do
plink2 --bgen ~/convert_new/ukb22828_c${i}_b0_v3.bgen ref-first --sample ~/convert_new/ukb22828_c10_b0_v3_s487159.sample --memory 15000 --extract-col-cond ~/convert_new/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 --geno 0.02 --mind 0.02 --hwe 1e-6 --maf 0.05 --keep-fam ../new_id.txt --pheno ../pc11.txt --covar ../cov.txt --covar-variance-standardize --glm  --out chr${i}_results 
done 



plink2 --bgen ~/convert_new/ukb22828_c14_b0_v3.bgen ref-first --sample ~/convert_new/ukb22828_c10_b0_v3_s487159.sample --memory 15000 --extract-col-cond ~/convert_new/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 --geno 0.02 --mind 0.02 --hwe 1e-6 --maf 0.005 --keep-fam new_id.txt --pheno pc31.txt --covar cov.txt --covar-variance-standardize --glm  --out chr14_results 


for i in {1..8};do
plink2 --bgen ~/convert_new/ukb22828_c${i}_b0_v3.bgen ref-first --sample ~/convert_new/ukb22828_c10_b0_v3_s487159.sample --memory 15000 --extract-col-cond ~/convert_new/ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 --geno 0.02 --mind 0.02 --hwe 1e-6 --maf 0.05 --keep-fam id_control.txt --pheno pc31.txt --covar cov.txt --covar-variance-standardize --glm  --out chr${i}_results 
done 
#提取ADD
for i in {1..22}; do
grep 'ADD' chr${i}_results.pc11.glm.linear > chr${i}_pc11.txt

done

mv chr*_normal.txt


#整理
/cluster/apps/R/4.2.2/bin/./Rscript --vanilla merge_file.R pc3_gene_race ./pc3_gene_race/ pc3_gene_race/pc3_p_gwas.csv

library(dplyr)
 


a = list.files("pc1")
dir = paste("./pc1/", a, sep="")
n = length(dir)

merge_data <- read.table(file = dir[1], head = FALSE)
for (i in 2:n){
  new_data <- read.table(file = dir[i], head = FALSE)
  merge_data <- rbind(merge_data, new_data)
 }

merge_data2 <- merge_data %>% rename(rsid = V3, chr = V1, BP = V2, P = V12, Beta=V9) %>% select(rsid,chr,BP,Beta, P)
write.csv(merge_data2, file = "pc1_vigious/pc1_all_gwas.csv", quote = FALSE, row.names = FALSE)
merge_adj2 <- merge_data2 %>% filter(P < 5e-08)

write.csv(merge_adj2, file = "pc1/pc1_gwas.csv", quote = FALSE, row.names = FALSE)

##转换rsid

library(MungeSumstats)
df <- read.csv("gwas_DP/pc21_gwas.csv")
df <- df %>% select(chr, BP, A1, A2, Beta, P)
reformatted <- format_sumstats(df,
                               ref_genome = "GRCh37",
                               nThread = 5,
                               return_data =T)
##############################################

#LD score
awk -F, '{ if (NR>1) { print sprintf("%02d", $2)":"$3"-"$3 }}' pc3_ld.csv > chrposlist.txt
bgenix -g ~/convert_new/ukb22828_c17_b0_v3.bgen -incl-range chrposlist.txt > pc3_ld.bgen

cmd=""
for i in {1..22}
do
  bgenix -g ukb22828_c${i}_b0_v3.bgen -incl-range chrposlist.txt > chr_${i}.bgen
  cmd=$cmd"chr_${i}.bgen "
done

# Combine the .bgen files for each chromosome into one
cat-bgen -g $cmd -og pc1_ld.bgen -clobber
# Write index file .bgen.bgi
bgenix -g pc1_ld.bgen -index -clobber


# Remove the individual chromosome files
for i in {1..22}
do
  rm chr_${i}.bgen
done


plink2 --bgen pc1_ld.bgen --sample ukb22828_c10_b0_v3_s487159.sample --indep-pairwise 500 50 0.2 --out input_pruned

#绘图部分
library(dplyr)
 


a = list.files("pc1_vigious")
dir = paste("./pc1_vigious/", a, sep="")
n = length(dir)

merge_data <- read.table(file = dir[1], head = FALSE)
for (i in 2:n){
  new_data <- read.table(file = dir[i], head = FALSE)
  merge_data <- rbind(merge_data, new_data)
 }

merge_data2 <- merge_data %>% rename(rsid = V3, chr = V1, BP = V2, P = V12) %>% select(rsid,chr,BP, P)

merge_adj2 <- merge_data2 %>% filter(P < 5e-08)

write.csv(merge_adj2, file = "pc1_vigious/final_p_snp.csv", quote = FALSE, row.names = FALSE)




#曼哈顿图
library(CMplot)

CMplot(merge_data2, plot.type = "m", threshold = c(5e-08, 1e-10), threshold.col = c('grey','black'),
       threshold.lty = c(1,2), threshold.lwd = c(1,1),amplify = T, 
       signal.cex = c(1,1), signal.pch = c(20,20), file = "tiff", file.output = TRUE, file.name = "test")


















