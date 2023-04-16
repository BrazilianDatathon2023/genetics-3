library(dplyr)

# subindo os arquivos
df_ld <- read.table(file = 'C:/Users/DRT50771/Downloads/ld-snps.tsv', sep = '\t', header = TRUE)
df_pheno <- read.table(file = 'C:/Users/DRT50771/Downloads/datathon-pheno.tsv', sep = '\t', header = TRUE)
df_snps <- read.table(file = 'C:/Users/DRT50771/Downloads/snps-functional-classification.tsv', sep = '\t', header = TRUE)
df_effect <- read.table(file = 'C:/Users/DRT50771/Downloads/PGS002296_hmPOS_GRCh38.tsv', sep = '\t', header = TRUE)


#####################
# criando as chaves para referenciar 
df_snps$column_key <- paste0(df_snps$chromosome,":",df_snps$pos_hg19,":",df_snps$effect_allele,":",df_snps$non_effect_allele)
df_snps$column_key_invert <- paste0(df_snps$chromosome,":",df_snps$pos_hg19,":",df_snps$non_effect_allele,":",df_snps$effect_allele)

df_effect$semchr <- gsub(pattern = "chr",replacement = "",x = df_effect$chr)
df_effect$key_19 <- paste0(df_effect$semchr,":",df_effect$pos_hg19,":",df_effect$effect_allele,":",df_effect$other_allele)
df_effect$key_19_invert <- paste0(df_effect$semchr,":",df_effect$pos_hg19,":",df_effect$other_allele,":",df_effect$effect_allele)

df_effect$key_38 <- paste0(df_effect$semchr,":",df_effect$pos_hg38,":",df_effect$effect_allele,":",df_effect$other_allele)
df_effect$key_38_invert <- paste0(df_effect$semchr,":",df_effect$pos_hg38,":",df_effect$other_allele,":",df_effect$effect_allele)

##############################

snp_to_effect <- lista_chr %>% 
  left_join(., df_effect,by=c("newchr2"="key_38"))

############################

effects_to_featsnp <- df_effect %>% 
  left_join(.,df_snps,by=c("key_19"="column_key"))


##### produção dos filtros

#################################################
# Criando os filtros que achamos necessários

df_ld_to_snp <- df_ld %>% 
  left_join(.,df_snps, by=c("rsID"="rs_id"))

df_ld_group_func <- df_ld %>% 
  group_by(IndSigSNP,func) %>% summarise(count_func = n())%>% 
  mutate(func = paste0("func_",func)) %>% 
  reshape2::dcast(.,IndSigSNP ~ func,fun.aggregate = sum,value.var = "count_func")

df_ld_group_rdb <- df_ld %>% 
  group_by(IndSigSNP,RDB) %>% summarise(count_rdb = n()) %>% 
  mutate(RDB = paste0("rdb_",RDB)) %>% 
  reshape2::dcast(.,IndSigSNP ~ RDB,fun.aggregate = sum,value.var = "count_rdb")

df_ld_group_eqtl_cimp <- df_ld %>% 
  group_by(IndSigSNP) %>% summarise(eqtlMapFilt = sum(eqtlMapFilt,na.rm = TRUE)
                                    ,ciMapFilt = sum(ciMapFilt,na.rm = TRUE))

df_ld_aggr <- df_ld_group_rdb %>% 
  left_join(.,df_ld_group_func,by="IndSigSNP") %>% 
  left_join(.,df_ld_group_eqtl_cimp,by="IndSigSNP") 


df_snp_to_ldagg <- df_snps %>% 
  left_join(.,df_ld_aggr,by=c("rs_id"="IndSigSNP")) 

# Filtro da função da variante
funcs_descount <- c("synonymous_variant","non_coding_transcript_exon_variant","3_prime_UTR_variant","mature_miRNA_variant","downstream_gene_variant")
funcs_descount_withintro <- c("intron_variant","synonymous_variant","non_coding_transcript_exon_variant","3_prime_UTR_variant","mature_miRNA_variant","downstream_gene_variant")

df_snp_to_ldagg$filter_funcs <- ifelse(!(df_snp_to_ldagg$consequence %in% funcs_descount),1,0)
df_snp_to_ldagg$filter_funcs_wointro <- ifelse(!(df_snp_to_ldagg$consequence %in% funcs_descount_withintro),1,0)

# filtro do RDB
df_snp_to_ldagg$filter_rdb_aux <- df_snp_to_ldagg$rdb_2a+df_snp_to_ldagg$rdb_2c+df_snp_to_ldagg$rdb_3b+df_snp_to_ldagg$rdb_4
df_snp_to_ldagg$filter_rdb <- ifelse(is.na(df_snp_to_ldagg$filter_rdb_aux),0,ifelse(df_snp_to_ldagg$filter_rdb_aux>0,1,0))

# Filtro de EqtlMapping
df_snp_to_ldagg$filter_eqtlmap <- ifelse(is.na(df_snp_to_ldagg$eqtlMapFilt),0,ifelse(df_snp_to_ldagg$eqtlMapFilt>0,1,0))

# Filtro de gene
genes_count <- c("ENSG00000144891","ENSG00000159640","ENSG00000143839","ENSG00000135744","ENSG00000135744","ENSG00000118762","ENSG00000081052","ENSG00000157388","ENSG00000123454","ENSG00000151623","ENSG00000087274","ENSG00000164867","ENSG00000169432","ENSG00000033867","ENSG00000070961","ENSG00000175206","ENSG00000120937","ENSG00000165995","ENSG00000169344")
df_snp_to_ldagg$filter_gene <- ifelse(df_snp_to_ldagg$gene %in% genes_count,1,0)

# Filtro de tecido
df_snp_to_ldagg$filter_aorta <- ifelse(is.na(df_snp_to_ldagg$Artery_Aorta),0,ifelse(df_snp_to_ldagg$Artery_Aorta=="high-expression",1,0))

# Criação do filtro unico
df_snp_to_ldagg$filters_wintro <- df_snp_to_ldagg$filter_funcs + df_snp_to_ldagg$filter_rdb + df_snp_to_ldagg$filter_eqtlmap + df_snp_to_ldagg$filter_gene + df_snp_to_ldagg$filter_aorta
df_snp_to_ldagg$filters_wintro <- ifelse(is.na(df_snp_to_ldagg$filters_wintro),0,df_snp_to_ldagg$filters_wintro)+1


df_snp_to_ldagg$filters_wointro <- df_snp_to_ldagg$filter_funcs_wointro + df_snp_to_ldagg$filter_rdb + df_snp_to_ldagg$filter_eqtlmap + df_snp_to_ldagg$filter_gene + df_snp_to_ldagg$filter_aorta
df_snp_to_ldagg$filters_wointro <- ifelse(is.na(df_snp_to_ldagg$filters_wointro),0,df_snp_to_ldagg$filters_wointro)+1

# Cnferencia
table(df_snp_to_ldagg$filters_wintro)
table(df_snp_to_ldagg$filters_wointro)

base_filter_wintro <- df_snp_to_ldagg %>% select(rs_id,filters_wintro,filter_funcs,filter_rdb,filter_eqtlmap,filter_gene,filter_aorta) %>% 
  filter(filters_wintro>0)


base_filter_wointro <- df_snp_to_ldagg %>% select(rs_id,filters_wointro,filter_funcs_wointro,filter_rdb,filter_eqtlmap,filter_gene,filter_aorta) %>% 
  filter(filters_wointro>0)

# salvando os arquivos de filtro
write.table(x = base_filter_wintro,file = "C:/Users/DRT50771/Desktop/filter_withintron.txt")
write.table(x = base_filter_wointro,file = "C:/Users/DRT50771/Desktop/filter_withoutintron.txt")



################################
# Subindo os resultados de prs
df_prs <- read.csv("C:/Users/DRT50771/Downloads/prs_scores.csv")
df_prs_withintro <- read.csv("C:/Users/DRT50771/Downloads/prs_withintron.csv")
df_prs_withoutintro <- read.csv("C:/Users/DRT50771/Downloads/prs_withoutintron.csv")
df_prs_withintro_sum1 <- read.csv("C:/Users/DRT50771/Downloads/prs_withintron_sum1.csv")
df_prs_withoutintro_sum1 <- read.csv("C:/Users/DRT50771/Downloads/prs_withoutintron_sum1.csv")


# juntando para fazer o modelo

df_pheno_prs <- df_pheno %>% 
  left_join(.,df_prs,by=c("s"="X")) %>% 
  mutate(prs = scale(X0)
         ,decil = ntile(prs,n = 10))

df_pheno_prswithintron <- df_pheno %>% 
  left_join(.,df_prs_withintro,by=c("s"="X")) %>% 
  mutate(prs = scale(X0)
         ,decil = ntile(prs,n = 10))

df_pheno_prswithoutintron <- df_pheno %>% 
  left_join(.,df_prs_withoutintro,by=c("s"="X")) %>% 
  mutate(prs = scale(X0)
         ,decil = ntile(prs,n = 10))

df_pheno_prswithintron_sum1 <- df_pheno %>% 
  left_join(.,df_prs_withintro_sum1,by=c("s"="X")) %>% 
  mutate(prs = scale(X0)
         ,decil = ntile(prs,n = 10))

df_pheno_prswithoutintron_sum1 <- df_pheno %>% 
  left_join(.,df_prs_withoutintro_sum1,by=c("s"="X")) %>% 
  mutate(prs = scale(X0)
         ,decil = ntile(prs,n = 10))

agg_full <- cbind(df_pheno_prs %>% group_by(decil) %>% summarise(avg_full = mean(target_phenotype))
                  ,df_pheno_prswithintron %>% group_by(decil) %>% summarise(avg_withintron = mean(target_phenotype)) %>% select(-decil)
                  ,df_pheno_prswithoutintron %>% group_by(decil) %>% summarise(avg_withintron = mean(target_phenotype)) %>% select(-decil)
                  ,df_pheno_prswithintron_sum1 %>% group_by(decil) %>% summarise(avg_withintron = mean(target_phenotype)) %>% select(-decil)
                  ,df_pheno_prswithoutintron_sum1 %>% group_by(decil) %>% summarise(avg_withintron = mean(target_phenotype)) %>% select(-decil)
)

library(ggplot2)

ggplot(df_pheno_prs, aes(x=as.factor(decil), y=prs, fill=as.factor(target_phenotype))) +
  geom_boxplot()

ggplot(df_pheno_prswithintron, aes(x=as.factor(decil), y=prs, fill=as.factor(target_phenotype))) +
  geom_boxplot()

ggplot(df_pheno_prswithoutintron, aes(x=as.factor(decil), y=prs, fill=as.factor(target_phenotype))) +
  geom_boxplot()

ggplot(df_pheno_prswithintron_sum1, aes(x=as.factor(decil), y=prs, fill=as.factor(target_phenotype))) +
  geom_boxplot()

ggplot(df_pheno_prswithoutintron_sum1, aes(x=as.factor(decil), y=prs, fill=as.factor(target_phenotype))) +
  geom_boxplot()



ggplot(df_pheno_prs,aes(x = prs, fill = as.factor(target_phenotype))) +
  geom_density(alpha = 0.7)

ggplot(df_pheno_prswithintron,aes(x = prs, fill = as.factor(target_phenotype))) +
  geom_density(alpha = 0.7)

ggplot(df_pheno_prswithoutintron,aes(x = prs, fill = as.factor(target_phenotype))) +
  geom_density(alpha = 0.7)

ggplot(df_pheno_prswithintron_sum1,aes(x = prs, fill = as.factor(target_phenotype))) +
  geom_density(alpha = 0.7)

ggplot(df_pheno_prswithoutintron_sum1,aes(x = prs, fill = as.factor(target_phenotype))) +
  geom_density(alpha = 0.7)



#### modelagem

glm_full <- glm(target_phenotype ~ prs + age + bmi + sex, data = df_pheno_prs, family = binomial("logit"))
glm_f_withintron <- glm(target_phenotype ~ prs + age + bmi + sex, data = df_pheno_prswithintron, family = binomial("logit"))
glm_f_withouintron <- glm(target_phenotype ~ prs + age + bmi + sex, data = df_pheno_prswithoutintron, family = binomial("logit"))
glm_s_withintron <- glm(target_phenotype ~ prs + age + bmi + sex, data = df_pheno_prswithintron_sum1, family = binomial("logit"))
glm_s_withoutintron <- glm(target_phenotype ~ prs + age + bmi + sex, data = df_pheno_prswithoutintron_sum1, family = binomial("logit"))

summary(glm_full)$coef %>% as.data.frame(.) %>% mutate(OR = exp(Estimate))
summary(glm_f_withintron)$coef %>% as.data.frame(.) %>% mutate(OR = exp(Estimate))
summary(glm_f_withouintron)$coef %>% as.data.frame(.) %>% mutate(OR = exp(Estimate))
summary(glm_s_withintron)$coef %>% as.data.frame(.) %>% mutate(OR = exp(Estimate))
summary(glm_s_withoutintron)$coef %>% as.data.frame(.) %>% mutate(OR = exp(Estimate))

