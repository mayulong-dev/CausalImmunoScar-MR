
library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(data.table)
library(MendelianRandomization)
library(MRPRESSO)
library(plinkbinr)
library(ieugwasr)
library(seqminer)
library(friendly2MR)
library(xlsx)

rm(list=ls())
setwd("D:/yanwork/IIT-I-231222663-MR")
load("D:/yanwork/IIT-I-231222663-MR/result/Hypertrophic scar.Rdata")
load("D:/yanwork/IIT-I-231222663-MR/Exposure inflammatory cytokines.RData")

###############1、读取暴露和结局数据###############

exp_CTACK <- extract_instruments(outcomes = "ebi-a-GCST004420",p1=5e-06)
exp_CTACK$exposure<-"CTACK/CCL27"
exp_Eotaxin <- extract_instruments(outcomes = "ebi-a-GCST004460",p1=5e-06)
exp_Eotaxin$exposure<-"Eotaxin"
exp_GROa <- extract_instruments(outcomes = "ebi-a-GCST004457",p1=5e-06)
exp_GROa$exposure<-"GROa/CXCL1"
exp_IP10 <- extract_instruments(outcomes = "ebi-a-GCST004440",p1=5e-06)
exp_IP10$exposure<-"IP10"
exp_MCP1 <- extract_instruments(outcomes = "ebi-a-GCST004438",p1=5e-06)
exp_MCP1$exposure<-"MCP1"
exp_MCP3 <- extract_instruments(outcomes = "ebi-a-GCST004437",p1=5e-06)
exp_MCP3$exposure<-"MCP3"
exp_MIG <- extract_instruments(outcomes = "ebi-a-GCST004435",p1=5e-06)
exp_MIG$exposure<-"MIG"
exp_MIP1a <- extract_instruments(outcomes = "ebi-a-GCST004434",p1=5e-06)
exp_MIP1a$exposure<-"MIP1a"
exp_MIP1b <- extract_instruments(outcomes = "ebi-a-GCST004433",p1=5e-08)
exp_MIP1b$exposure<-"MIP1b"
exp_RANTES <- extract_instruments(outcomes = "ebi-a-GCST004431",p1=5e-06)
exp_RANTES$exposure<-"RANTES/CCL5"
exp_Bngf <- extract_instruments(outcomes = "ebi-a-GCST004421",p1=5e-06)
exp_Bngf$exposure<-"Bngf"
exp_FGFBasic <- extract_instruments(outcomes = "ebi-a-GCST004459",p1=5e-06)
exp_FGFBasic$exposure<-"FGFBasic"
exp_GCSF <- extract_instruments(outcomes = "ebi-a-GCST004458",p1=5e-06)
exp_GCSF$exposure<-"GCSF"
exp_HGF <- extract_instruments(outcomes = "ebi-a-GCST004449",p1=5e-06)
exp_HGF$exposure<-"HGF"
exp_MCSF <- extract_instruments(outcomes = "ebi-a-GCST004436",p1=5e-06)
exp_MCSF$exposure<-"MCSF"
exp_PDGFbb <- extract_instruments(outcomes = "ebi-a-GCST004432",p1=5e-06)
exp_PDGFbb$exposure<-"PDGFbb"
exp_SCF <- extract_instruments(outcomes = "ebi-a-GCST004429",p1=5e-06)
exp_SCF$exposure<-"SCF"
exp_SCGFb <- extract_instruments(outcomes = "ebi-a-GCST004428",p1=5e-08)
exp_SCGFb$exposure<-"SCGFb"
exp_VEGF <- extract_instruments(outcomes = "ebi-a-GCST004422",p1=5e-08)
exp_VEGF$exposure<-"VEGF"
exp_IL10 <- extract_instruments(outcomes = "ebi-a-GCST004444",p1=5e-06)
exp_IL10$exposure<-"IL-10"
exp_IL12p70 <- extract_instruments(outcomes = "ebi-a-GCST004439",p1=5e-06)
exp_IL12p70$exposure<-"IL-12p70"
exp_IL13 <- extract_instruments(outcomes = "ebi-a-GCST004443",p1=5e-06)
exp_IL13$exposure<-"IL-13"
exp_IL16 <- extract_instruments(outcomes = "ebi-a-GCST004430",p1=5e-06)
exp_IL16$exposure<-"IL-16"
exp_IL17 <- extract_instruments(outcomes = "ebi-a-GCST004442",p1=5e-06)
exp_IL17$exposure<-"IL-17"
exp_IL18 <- extract_instruments(outcomes = "ebi-a-GCST004441",p1=5e-06)
exp_IL18$exposure<-"IL-18"
exp_IL1b <- extract_instruments(outcomes = "ebi-a-GCST004448",p1=5e-06)
exp_IL1b$exposure<-"IL-1b"
exp_IL1ra <- extract_instruments(outcomes = "ebi-a-GCST004447",p1=5e-06)
exp_IL1ra$exposure<-"IL-1ra"
exp_IL2 <- extract_instruments(outcomes = "ebi-a-GCST004455",p1=5e-06)
exp_IL2$exposure<-"IL-2"
exp_IL2ra <- extract_instruments(outcomes = "ebi-a-GCST004454",p1=5e-06)
exp_IL2ra$exposure<-"IL-2ra"
exp_IL4 <- extract_instruments(outcomes = "ebi-a-GCST004453",p1=5e-06)
exp_IL4$exposure<-"IL-4"
exp_IL5 <- extract_instruments(outcomes = "ebi-a-GCST004452",p1=5e-06)
exp_IL5$exposure<-"IL-5"
exp_IL6 <- extract_instruments(outcomes = "ebi-a-GCST004446",p1=5e-06)
exp_IL6$exposure<-"IL-6"
exp_IL7 <- extract_instruments(outcomes = "ebi-a-GCST004451",p1=5e-06)
exp_IL7$exposure<-"IL-7"
exp_IL8 <- extract_instruments(outcomes = "ebi-a-GCST004445",p1=5e-06)
exp_IL8$exposure<-"IL-8"
exp_IL9 <- extract_instruments(outcomes = "ebi-a-GCST004450",p1=5e-06)
exp_IL9$exposure<-"IL-9"
exp_IFNg<- extract_instruments(outcomes = "ebi-a-GCST004456",p1=5e-06)
exp_IFNg$exposure<-"IFNg"
exp_MIF<- extract_instruments(outcomes = "ebi-a-GCST004423",p1=5e-06)
exp_MIF$exposure<-"MIF"
exp_TNFa<- extract_instruments(outcomes = "ebi-a-GCST004426",p1=5e-06)
exp_TNFa$exposure<-"TNFa"
exp_TNFb<- extract_instruments(outcomes = "ebi-a-GCST004425",p1=5e-06)
exp_TNFb$exposure<-"TNFb"
exp_TRAIL<- extract_instruments(outcomes = "ebi-a-GCST004424",p1=5e-08)
exp_TRAIL$exposure<-"TRAIL"


exp_data<-rbind(exp_CTACK,exp_Eotaxin,exp_GROa,exp_IP10,exp_MCP1,exp_MCP3,exp_MIG,exp_MIP1a,exp_MIP1b,exp_RANTES,
                exp_Bngf,exp_FGFBasic,exp_GCSF,exp_HGF,exp_MCSF,exp_PDGFbb,exp_SCF,exp_SCGFb,exp_VEGF,exp_IL10,
                exp_IL12p70,exp_IL13,exp_IL16,exp_IL17,exp_IL18,exp_IL1b,exp_IL1ra,exp_IL2,exp_IL2ra,exp_IL4,exp_IL5,
                exp_IL6,exp_IL7,exp_IL8,exp_IL9,exp_IFNg,exp_MIF,exp_TNFa,exp_TNFb,exp_TRAIL)

colnames(exp_data)
#[1] "chr.exposure"           "samplesize.exposure"    "pos.exposure"           "pval.exposure"          "beta.exposure"         
#[6] "se.exposure"            "id.exposure"            "SNP"                    "effect_allele.exposure" "other_allele.exposure" 
#[11] "eaf.exposure"           "exposure"               "mr_keep.exposure"       "pval_origin.exposure"   "data_source.exposure"

snp_add_eaf <- function(dat, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(dat))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  snp_reverse_base <- function(x)
  {
    x <- stringr::str_to_upper(x)
    stopifnot(x %in% c("A","T","C","G"))
    switch(x,"A"="T","T"="A","C"="G","G"="C")
  }
  
  res_tab <- lapply(1:nrow(dat), function(i)
  {
    print(paste0("seaching for No.", i, " SNP"))
    dat_i <- dat[i,]
    
    ext <- paste0("/variation/Homo_sapiens/",dat_i$SNP, "?content-type=application/json;pops=1")
    url <- paste(server, ext, sep = "")
    res <- httr::GET(url)
    
    # Converts http errors to R errors or warnings
    httr::stop_for_status(res)
    
    # Convert R objects from JSON
    res <- httr::content(res)
    res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
    
    # Filter query results based on population set
    res_pop <- try(res_pop[res_pop$population == pop,])
    if("try-error" %in% class(res_pop))
    {
      print(paste0("There is not information for population ",pop))
      queried_effect_allele <- "NR"
      queried_other_allele <- "NR"
      queried_eaf <- -1
    }
    else
    {
      if(nrow(res_pop)==0)
      {
        print(paste0("There is not information for population ",pop))
        queried_effect_allele <- "NR"
        queried_other_allele <- "NR"
        queried_eaf <- -1
      }
      else
      {
        queried_effect_allele <- res_pop[1,"allele"][[1]]
        queried_other_allele <- res_pop[2,"allele"][[1]]
        queried_eaf <- res_pop[1,"frequency"][[1]]    
      }
    }
    
    effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                            dat_i$effect_allele.exposure,
                            dat_i$effect_allele)
    
    other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                           dat_i$other_allele.exposure,
                           dat_i$other_allele)
    
    if("effect_allele.exposure" %in% names(dat))
    {
      name_output <- unique(c(names(dat), "eaf.exposure","reliability.exposure"))
    }
    else
    {
      name_output <- unique(c(names(dat), "eaf","reliability.exposure"))
    }
    
    len_effect_allele <- nchar(effect_allele)
    len_other_allele <- nchar(other_allele)
    
    if(len_effect_allele==1&len_other_allele==1)
    {
      if((queried_effect_allele==effect_allele & queried_other_allele==other_allele)|
         (queried_effect_allele==other_allele & queried_other_allele==effect_allele))
      {
        dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                     queried_eaf,
                                     1-queried_eaf)
        dat_i$eaf <- dat_i$eaf.exposure 
        dat_i$reliability.exposure <- "high"
      }
      else
      {
        r_queried_effect_allele <- snp_reverse_base(queried_effect_allele)
        r_queried_other_allele <- snp_reverse_base(queried_other_allele)
        if((r_queried_effect_allele==effect_allele & r_queried_other_allele==other_allele)|
           (r_queried_effect_allele==other_allele & r_queried_other_allele==effect_allele))
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == r_queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "high"
        }
        else
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "low"
        }
      }
    }
    
    else
    {
      # To identify the potential DEL/ INS
      short_allele <- ifelse(len_effect_allele==1,
                             effect_allele,
                             other_allele)
      short_allele_eaf <- ifelse(short_allele == queried_effect_allele, 
                                 queried_eaf, 
                                 1-queried_eaf)
      dat_i$eaf.exposure <- ifelse(effect_allele == short_allele,
                                   short_allele_eaf,
                                   1-short_allele_eaf)
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "low"
    }
    
    dat_i[name_output]
  })
  
  return(do.call(rbind, res_tab))
}
nrow(exp_data)
exp_data<-snp_add_eaf(exp_data)
exp_data[1:187,]<-snp_add_eaf(exp_data[1:187,])
exp_data[189:399,]<-snp_add_eaf(exp_data[189:399,])
exp_data[188,]$eaf.exposure<-0.05

exp_data1<-exp_data[,c("SNP","exposure","pval.exposure","chr.exposure","pos.exposure",
                       "effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure",
                       "se.exposure","samplesize.exposure")]
write.xlsx(exp_data1, file = "D:/yanwork/IIT-I-231222663-MR/data/IV_ALL.xlsx", row.names = F)

table(exp_data$exposure)
#################结局###################
out_Hypertrophic <- extract_outcome_data(snps = exp_data$SNP, outcomes = "finn-b-L12_HYPETROPHICSCAR")
out_Hypertrophic$outcome<-"Hypertrophic scar"
out_Hypertrophic$samplesize.outcome<-766+207482


# 2、计算F值和R^2值----

#公式：R² = 2*Beta^2*EAF*(1-EAF)/(2*Beta^2*EAF*(1-EAF)+ 2*SE^2*SampleSize*EAF*(1-EAF))

exp_data$r2 <- (2 * exp_data$beta.exposure^2 * exp_data$eaf.exposure * (1 - exp_data$eaf.exposure)) / 
  (2 * exp_data$beta.exposure^2 * exp_data$eaf.exposure * (1 - exp_data$eaf.exposure) + 
     2 * exp_data$samplesize.exposure * exp_data$eaf.exposure * (1 - exp_data$eaf.exposure) * exp_data$se.exposure^2)
r2 <- sum(exp_data$r2)

#公式：F = R²*(SampleSize-2)/(1-R²)

exp_data$F <- exp_data$r2 * (exp_data$samplesize.exposure - 2) / (1 - exp_data$r2)
exp_data_meanF <- mean(exp_data$F)
quantile(exp_data$F)

# 3、harmonise----
nrow(exp_data)
har_Hypertrophic <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_Hypertrophic)

har_Hypertrophic2<-har_Hypertrophic[har_Hypertrophic$mr_keep==TRUE,c("SNP",	"chr","pos.exposure","effect_allele.exposure","other_allele.exposure",
                                                   "eaf.exposure","exposure","beta.exposure","se.exposure","pval.exposure","r2","F","outcome","beta.outcome","se.outcome",
                                                   "pval.outcome")]
colnames(har_Hypertrophic2)<-c("SNP","Chr","Pos","EA","OA","EAF","exp","beta","se","pval","r2","F","outcome","beta1","se1","pval1")

write.csv(har_Hypertrophic2, file = "D:/yanwork/IIT-I-231222663-MR/result/IV_out_Hypertrophic.csv", row.names = F)

find_exp<-function(exposure,singleexp,singloutcome,harmonise){
  a<-list()
  a$`SNP summary`<-cbind(sprintf("%.0f",nrow(exposure[exposure$exposure==singleexp,])),sprintf("%.2f",mean(exposure[exposure$exposure==singleexp,]$F)),
                         sprintf("%.2f",min(exposure[exposure$exposure==singleexp,]$F)),sprintf("%.2f",max(exposure[exposure$exposure==singleexp,]$F)))
  colnames(a$`SNP summary`)<-c("The number of all SNP","Mean F","Mini F","Max F")
  a$`The number of all unmatched SNP`<-length(exposure[exposure$exposure==singleexp,][!(exposure[exposure$exposure==singleexp,"SNP"] %in%
                                                                                          harmonise[harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"SNP"]),"SNP"])+
    nrow(harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,])
  a$`All unmatched SNP`<-c(unlist(exposure[exposure$exposure==singleexp,][!(exposure[exposure$exposure==singleexp,"SNP"] %in%harmonise[harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"SNP"]),"SNP"]),
                           unlist(harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"SNP"]))
  a$`SNP not found in summary`<-exposure[exposure$exposure==singleexp,][!(exposure[exposure$exposure==singleexp,"SNP"] %in%harmonise[harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"SNP"]),"SNP"]
  a$`Proxy SNP`<-cbind(harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"target_snp.outcome"],
                       harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"proxy_snp.outcome"],
                       round(harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,"r2"],8))
  colnames(a$`Proxy SNP`)<-c("target_snp.outcome","proxy_snp.outcome","r2")
  a$`The number of proxy SNP`<-nrow(harmonise[!is.na(harmonise$`proxy.outcome`)&harmonise$exposure==singleexp&harmonise$outcome==singloutcome,])
  a$`MR keep FALSE`<-harmonise[harmonise$exposure==singleexp&harmonise$outcome==singloutcome&harmonise$mr_keep==FALSE,"SNP"]
  return(a)
}


table(exp_data$exposure)
find_exp(exp_data,"TRAIL","Hypertrophic scar",har_Hypertrophic)

# 4、MR分析----
mr_method_list()
res <- mr(dat = har_Hypertrophic, 
          method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode","mr_wald_ratio",'mr_ivw_mre'))
res_or <- generate_odds_ratios(mr_res = res)
res_or$`OR (95% CI)`<-paste0(sprintf("%.3f",res_or$or)," (",sprintf("%.3f",res_or$or_lci95)," - ",sprintf("%.3f",res_or$or_uci95),")")
res_or$P<-sprintf("%.3f",res_or$pval)
res_or$`FDR adjusted P`<-sprintf("%.3f",p.adjust(res_or[,"pval"],method = "BH"))
write.xlsx(res_or, file = "D:/yanwork/IIT-I-231222663-MR/result/OR_Hypertrophic scar.xlsx", row.names = F)

res_or[res_or$P<0.05,]

p_scatter <- mr_scatter_plot(mr_results = res_or, dat = har_Hypertrophic)
p_scatter[8]
unique(res_or$exposure)
plot_name<-c("CTACK","Bngf","VEGF","MIF","TRAIL","TNFb","TNFa","SCGFb","SCF","IL-16",
             "RANTES","PDGFbb","MIP1b","MIP1a","MIG","MCSF","MCP3","MCP1","IL-12p70","IP10",
             "IL-18","IL-17","IL-13","IL-10","IL-8","IL-6","IL-1ra","IL-1b","HGF","IL-9",
             "IL-7","IL-5","IL-4","IL-2ra","IL-2","IFNg","GROa","GCSF","FGFBasic","Eotaxin")

for (i in c(1:40)){
  jpeg(filename = paste0("D:/yanwork/IIT-I-231222663-MR/result/Figure1_",
                         plot_name[i],"_scatter.jpeg"),units = "in",width = 10,height = 7,res = 600)
  print(p_scatter[i])
  dev.off()
}

# 5、异质性检验----
res_hetergeneity <- mr_heterogeneity(dat = har_Hypertrophic,
                                     method_list=c("mr_egger_regression",'mr_ivw_mre'))
res_hetergeneity[res_hetergeneity$method=="Inverse variance weighted (multiplicative random effects)","FDR adjusted P"]<-
  sprintf("%.2f",p.adjust(res_hetergeneity[res_hetergeneity$method=="Inverse variance weighted (multiplicative random effects)","Q_pval"],method = "BH"))
res_hetergeneity

write.xlsx(res_hetergeneity, file = "D:/yanwork/IIT-I-231222663-MR/result/hetergeneity.xlsx", row.names = F)
# 6、水平多效性分析----
res_pleiotropy <- mr_pleiotropy_test(dat = har_Hypertrophic)
res_pleiotropy$`FDR adjusted P`<-sprintf("%.2f",p.adjust(res_pleiotropy[,"pval"],method = "BH"))
write.xlsx(res_pleiotropy, file = "D:/yanwork/IIT-I-231222663-MR/result/pleiotropy.xlsx", row.names = F)

presso<-TwoSampleMR::run_mr_presso(har_Hypertrophic)
presso1<-function(presso,i){
  presso1<-presso[[i]]$`Main MR results`
  presso1$OR <- round(exp(presso1$`Causal Estimate`), 2)
  presso1$LOR <- round(exp(presso1$`Causal Estimate` - 1.96*presso1$Sd),2)
  presso1$UOR <- round(exp(presso1$`Causal Estimate` + 1.96*presso1$Sd),2)
  presso1$`OR (95% CI)`<-paste0(sprintf("%0.2f",presso1$OR)," (",sprintf("%0.2f",presso1$LOR)," - ",sprintf("%0.2f",presso1$UOR),")")
  presso1$P<-sprintf("%0.2f",presso1$`P-value`)
  presso1$`Global Test Pvalue`<-presso[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  presso1$`Distortion Test Pvalue`<-presso[[i]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
  presso1$`Distortion Test Outliers Indices`<-presso[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  return(presso1)
}
#[1] "MCP1"        "IL-18"       "IL-10"       "IFNg"        "IL-4"        "IL-1ra"      "IL-6"        "IL-12p70"    "VEGF"        "IP10"       
#[11] "TNFa"        "MIP1a"       "MCP3"        "TNFb"        "HGF"         "Eotaxin"     "GROa/CXCL1"  "MIG"         "RANTES/CCL5" "IL-2ra"     
#[21] "PDGFbb"      "MIP1b"       "SCF"         "MIF"         "CTACK/CCL27" "GCSF"        "IL-16"       "MCSF"        "IL-8"        "IL-17"      
#[31] "IL-5"        "FGFBasic"    "SCGFb"       "IL-7"        "IL-13"       "IL-9"        "IL-2"        "TRAIL"       "IL-1b"       "Bngf" 
presso1(presso,5)

presso1<-function(presso,number,exposure){
  out<-list()
  for (i in c(1:number)) {
    presso1<-presso[[i]]$`Main MR results`
    presso1$OR <- sprintf("%0.2f",exp(presso1$`Causal Estimate`))
    presso1$LOR <- sprintf("%0.2f",exp(presso1$`Causal Estimate` - 1.96*presso1$Sd))
    presso1$UOR <- sprintf("%0.2f",exp(presso1$`Causal Estimate` + 1.96*presso1$Sd))
    presso1$`OR (95% CI)`<-paste0(presso1$OR," (",presso1$LOR," - ",presso1$UOR,")")
    presso1$P<-sprintf("%0.2f",presso1$`P-value`)
    presso1$`Global Test Pvalue`<-presso[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
    presso1$Exposure<-exposure[i]
    out$P<-rbind(out$P,presso1)
    out$`Main MR results`$`Global Test`$Pvalue<-rbind(out$`Main MR results`$`Global Test`$Pvalue,
                                                      presso[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue)
    row.names(out$`Main MR results`$`Global Test`$Pvalue)[i]<-exposure[i]

  }
  return(out)
}
presso_all<-presso1(presso,40,c("MCP1","IL-18","IL-10","IFNg","IL-4","IL-1ra","IL-6","IL-12p70","VEGF","IP10","TNFa",
                                "MIP1a","MCP3","TNFb","HGF","Eotaxin","GROa/CXCL1","MIG","RANTES/CCL5","IL-2ra","PDGFbb",
                                "MIP1b","SCF","MIF","CTACK/CCL27","GCSF","IL-16","MCSF","IL-8","IL-17","IL-5",
                                "FGFBasic","SCGFb","IL-7","IL-13","IL-9","IL-2","TRAIL","IL-1b","Bngf"))

write.xlsx(presso_all$P, file = "D:/yanwork/IIT-I-231222663-MR/result/presso.xlsx",row.names = FALSE)

presso[23]

# 、harmonise----remove outlier
nrow(exp_data)
exp_data[exp_data$exposure=="SCF","SNP"][2]
exp_data[exp_data$SNP=="rs13412535",]
har_Hypertrophic[har_Hypertrophic$SNP=="rs13412535",]

mr_method_list()
res_remove <- mr(dat = har_Hypertrophic[-231,], 
          method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode","mr_wald_ratio",'mr_ivw_mre'))
res_or_remove <- generate_odds_ratios(mr_res = res_remove)
res_or_remove$`OR (95% CI)`<-paste0(sprintf("%.3f",res_or_remove$or)," (",sprintf("%.3f",res_or_remove$or_lci95)," - ",sprintf("%.3f",res_or_remove$or_uci95),")")
res_or_remove$P<-sprintf("%.3f",res_or_remove$pval)
res_or_remove$`FDR adjusted P`<-sprintf("%.3f",p.adjust(res_or_remove[,"pval"],method = "BH"))
write.xlsx(res_or_remove, file = "D:/yanwork/IIT-I-231222663-MR/result/OR_Hypertrophic scar_remove outlier.xlsx", row.names = F)

p_scatter_remove <- mr_scatter_plot(mr_results = res_or_remove, dat = har_Hypertrophic[-231,])
for (i in c(1:40)){
  jpeg(filename = paste0("D:/yanwork/IIT-I-231222663-MR/result/Figure1_",
                         plot_name[i],"_scatter remove outlier.jpeg"),units = "in",width = 10,height = 7,res = 600)
  print(p_scatter_remove[i])
  dev.off()
}
# 7、单个SNP效应分析----
res_single <- mr_singlesnp(dat = har_Hypertrophic[-231,])
## Forest plot
p_forest <- mr_forest_plot(singlesnp_results = res_single)
for (i in c(1:40)){
  jpeg(filename = paste0("D:/yanwork/IIT-I-231222663-MR/result/Figure2_",plot_name[i],"_Forest.jpeg"),
       units = "in",width = 10,height = 7,res = 600)
  print(p_forest[i])
  dev.off()
}

## Funnel plot
p_funnel <- mr_funnel_plot(singlesnp_results = res_single)
for (i in c(1:40)){
  jpeg(filename = paste0("D:/yanwork/IIT-I-231222663-MR/result/Figure3_",plot_name[i],"_funnel.jpeg"),
       units = "in",width = 10,height = 7,res = 600)
  print(p_funnel[i])
  dev.off()
}

# 8、Leave-one-out分析----
res_loo <- mr_leaveoneout(dat = har_Hypertrophic[-231,])
p_loo <- mr_leaveoneout_plot(leaveoneout_results = res_loo)

for (i in c(1:40)){
  jpeg(filename = paste0("D:/yanwork/IIT-I-231222663-MR/result/Figure4_",plot_name[i],"_leaveoneout.jpeg"),
       units = "in",width = 10,height = 7,res = 600)
  print(p_loo[i])
  dev.off()
}

save.image(file ="D:/yanwork/IIT-I-231222663-MR/result/Hypertrophic scar.Rdata")
