library(tidyverse)
library(ggpubr)
library(mgcv)
library(ggcorrplot)
library(gtools)

#change these as needed
path_to_wd<-"~/Microbiologie/RichardsLab/Sars-CoV-2/Somalogic/cytokine.epi.project"
path_to_data<-"soma_shuffled.txt"

#name of analysis, either "A2" or "A3"
analysis<-"A2"

setwd(path_to_wd)
if(!dir.exists(paste0(analysis))){
  dir.create(analysis)
}

#days of symptom onset cutoff.
time_cutoff<-14

# the following loads the data (soma), here it is assumed to be a tab delimited file, but adjust as needed.
# the data has the following columns:
# SubjectID
# age
# sex
# days_symptom
# A2 (coded as 0 for control and 1 for case)
# A3 (coded as 0 for control and 1 for case)
# then all 5284ish proteins from somalogic as a distinct column, with the column name being the somamer ID
# for example my first protein is CRYBB2.10000.28
# the order of columns does not matter, as long as the proteins are columns 7 to the last

soma_raw<-read_tsv(path_to_data)%>%
  filter(days_symptom<=time_cutoff) %>%
  mutate(A3=factor(A3,
                   labels=c("Control",
                            "Case"))) %>%
  mutate(A2=factor(A2,
                   labels=c("Control",
                            "Case"))) %>%
  mutate(sex=factor(sex))
soma_raw[,7:ncol(soma_raw)]<-scale(log(soma_raw[,7:ncol(soma_raw)]))

#need to take log, scale, and remove outliers
soma_outlier<-as.matrix(soma_raw[,7:ncol(soma_raw)])
soma_outlier[which(soma_outlier>3 | soma_outlier<(-3))]<-NA
soma_outlier[which(soma_outlier<=3 & soma_outlier>=(-3))]<-1

soma<-soma_raw
soma[,7:ncol(soma)]<-soma_raw[,7:ncol(soma_raw)]*soma_outlier

#the variables below are the proteins we'll look at, classified in general groups: CCL, CXCL, IL, soluble IL receptors, interferons, and others
colsoma<-data.frame(cols=colnames(soma))

#interleukins
interleukins<-colsoma %>% 
  filter(startsWith(cols, "IL")) %>%
  filter(!grepl("R",cols)) %>%
  filter(!grepl("ILF2",cols)) %>%
  filter(!grepl("ILF3",cols)) %>%
  filter(!grepl("ILK",cols)) %>%
  filter(!grepl("IL18BP",cols)) %>%
  filter(!grepl("IL1F5",cols)) %>%
  filter(!grepl("IL6ST",cols)) %>%
  filter(!grepl("IL12B.13733.5", cols)) %>%
  mutate(names=cols) %>%
  separate(names, into="names", sep="[.]", extra="drop")

interleukins$names[which(interleukins$cols=="IL12B.IL23A.10365.132")]<-"IL23"
interleukins$names[which(interleukins$cols=="IL12A.IL12B.10367.62")]<-"IL12"

interleukins <- interleukins %>% 
  mutate(names=factor(names, levels=mixedsort(names))) %>%
  arrange(names)

#soluble interleukin receptors
silr<-colsoma %>% 
  filter(startsWith(cols, "IL")) %>%
  filter(grepl("R",cols)) %>%
  mutate(names=cols) %>%
  separate(names, into="names", sep="[.]", extra="drop")

silr$names[which(silr$cols=="IL15RA.14054.17")]<-"IL15RA.soma1"
silr$names[which(silr$cols=="IL15RA.3445.53")]<-"IL15RA.soma2"
silr$names[which(silr$cols=="IL10RA.8104.21")]<-"IL10RA.soma1"
silr$names[which(silr$cols=="IL10RA.10344.334")]<-"IL10RA.soma2"

silr <- silr %>% 
  mutate(names=factor(names, levels=mixedsort(names))) %>%
  arrange(names)

#cxcl
cxcl<-colsoma %>% 
  filter(startsWith(cols, "CXCL")) %>%
  filter(!grepl("CXCL12.9278.9", cols)) %>%
  mutate(names=cols) %>%
  separate(names, into="names", sep="[.]", extra="drop") 

cxcl$names[which(cxcl$cols=="	CXCL9.11593.21")]<-"CXCL9.soma1"
cxcl$names[which(cxcl$cols=="CXCL9.9188.119")]<-"CXCL9.soma2"

cxcl <- cxcl %>%
  mutate(names=factor(names, levels=mixedsort(names)))%>%
  arrange(names)

#ccl
ccl<-colsoma %>% 
  filter(startsWith(cols, "CCL")) %>%
  filter(!grepl("CCL23.3028.36", cols)) %>%
  mutate(names=cols) %>%
  separate(names, into="names", sep="[.]", extra="drop") %>%
  filter(!grepl("CCL3L1", cols)) %>%
  filter(!grepl("CCL4L1", cols)) %>%
  mutate(names=factor(names, levels=mixedsort(names))) %>%
  arrange(names)

#interferons
interferons<-colsoma %>% 
  filter(startsWith(cols, "IFN")) %>%
  filter(!grepl("R",cols)) %>%
  mutate(names=cols) %>%
  separate(names, into="names", sep="[.]", extra="drop") %>%
  mutate(names=factor(names, levels=mixedsort(names))) %>%
  arrange(names)

#others
#may have to double check these names. I expect the first part (prior to first period) to be the same, but there's a chance that the numbers that follow might differ.
others<-data.frame(cols=c("CSF3.4840.73", 
                          "CSF2.4697.59", 
                          "CSF1.3738.54", 
                          "MIF.8221.19", 
                          "TNF.5936.53", 
                          "LTA.4703.87",
                          "TLR1.11149.3",
                          "TLR1.16324.38",
                          "TLR2.3835.11",
                          "TLR3.16918.198",
                          "TLR4.11101.18",
                          "TLR5.18935.14",
                          "IGHA1.IGHA2.11089.7",
                          "IGHD.IGK.IGL.4916.2",
                          "IGHE.IGK.IGL.4135.84",
                          "IGHG1.IGHG2.IGHG3.IGHG4.IGK.IGL.2744.57",
                          "IGHM.IGJ.IGK.IGL.3069.52"),
                   names=c("G-CSF", 
                           "GM-CSF", 
                           "M-CSF", 
                           "MIF", 
                           "TNF-\u03b1", 
                           "LT-\u03b1 / TNF-\u03b2",
                           "TLR1.soma1",
                           "TLR1.soma2",
                           "TLR2",
                           "TLR3",
                           "TLR4",
                           "TLR5",
                           "IgA",
                           "IgD",
                           "IgE",
                           "IgG",
                           "IgM")) %>%
  mutate(names=factor(names, levels=names))

######functions######
gam_fxn<-function(protein, name, k=15, scale_x=FALSE){
  loess_plot<-soma %>%
    ggplot(aes_string(x="days_symptom", 
                      y=protein,
                      group=analysis,
                      colour=analysis)) +
    geom_smooth() +
    xlab("Days since symptoms onset")+
    ylab(name)+
    scale_x_continuous(breaks=seq(0,15,5))+
    scale_color_manual(values=c("dodgerblue4", "red3"))+
    theme_bw()
  
  gam_full<-soma %>%
    gam(as.formula(paste(protein,"~s(days_symptom, by=",analysis,", bs=\"tp\", k=k)+
        s(days_symptom, by=sex, bs=\"tp\",k=k)+
        age+",
                         analysis,"+
        sex")),
        data=.,
        method = "REML")
  
  gam_null_outcome<-soma %>%
    gam(as.formula(paste(protein,"~s(days_symptom, by=sex, bs=\"tp\",k=k)+
        age+
        sex")),
        data=.,
        method = "REML")
  
  gam_null_sex<-soma %>%
    gam(as.formula(paste(protein,"~s(days_symptom, by=",analysis,", bs=\"tp\",k=k)+
        age+",
                         analysis)),
        data=.,
        method = "REML")
  
  gam_plot_M<-predict_graph(data=test_data_M, model=gam_full, name=name, scale_x=scale_x)
  gam_plot_F<-predict_graph(data=test_data_F, model=gam_full, name=name, scale_x=scale_x)
  
  gam_anova_outcome<-anova(gam_null_outcome, gam_full, test="Chisq")
  gam_anova_sex<-anova(gam_null_sex, gam_full, test="Chisq")
  
  return(list(loess_plot, gam_full, gam_plot_M, gam_plot_F, gam_anova_outcome, gam_anova_sex))
}


predict_graph<-function(data, model, name, scale_x=FALSE){
  fits<-data %>%
    predict(model, newdata=., type="response", se=TRUE)
  
  predicts <- data.frame(data, fits) %>% 
    mutate(lower = fit - 1.96*se.fit,
           upper = fit + 1.96*se.fit)
  if(scale_x==TRUE){
    result_plot <- ggplot(aes(x=days_symptom,
                              y=fit,
                              group=.data[[analysis]],
                              colour=.data[[analysis]]), 
                          data=predicts) +
      geom_ribbon(aes(ymin = lower, ymax=upper, group=.data[[analysis]], colour=.data[[analysis]], fill = .data[[analysis]]),
                  alpha=0.4) +
      geom_line(color="black")+    
      xlab("Days since symptoms onset")+
      ylab(name)+
      scale_x_continuous(breaks=seq(0,15,3))+
      scale_color_manual(values=c("dodgerblue4", "red3"))+
      scale_fill_manual(values=c("dodgerblue4", "red3"))+
      #ylim(c(-3,3))+
      theme_bw()
    
    return(result_plot)
  }
  else {
    result_plot <- ggplot(aes(x=days_symptom,
                              y=fit,
                              group=.data[[analysis]],
                              colour=.data[[analysis]]), 
                          data=predicts) +
      geom_ribbon(aes(ymin = lower, ymax=upper, group=.data[[analysis]], colour=.data[[analysis]], fill = .data[[analysis]]),
                  alpha=0.4) +
      geom_line(color="black")+    
      ylab(name)+
      scale_x_continuous(breaks=seq(0,15,3))+
      scale_color_manual(values=c("dodgerblue4", "red3"))+
      scale_fill_manual(values=c("dodgerblue4", "red3"))+
      #ylim(c(-3,3))+
      theme_bw()+
      theme(axis.title.x=element_blank())
    
    return(result_plot)
  }
}


#######correlation plots between cytokines######
method<-"spearman"

#corr_plot_cases
corr_cases <- soma %>%
  filter(days_symptom<=time_cutoff) %>%
  filter(.data[[analysis]]=="Case") %>%
  dplyr::select(c(interleukins$cols,
                  ccl$cols,
                  cxcl$cols,
                  interferons$cols,
                  silr$cols,
                  others$cols)) %>%
  stats::cor(x=.,use="pairwise.complete.obs",
             method=method)

corr_p_cases<- soma %>%
  filter(days_symptom<=time_cutoff) %>%
  filter(.data[[analysis]]=="Case") %>%
  dplyr::select(c(interleukins$cols,
                  ccl$cols,
                  cxcl$cols,
                  interferons$cols,
                  silr$cols,
                  others$cols)) %>%
  cor_pmat( method = method, exact=FALSE)

colnames(corr_cases) <- c(as.character(interleukins$names),
                          as.character(ccl$names),
                          as.character(cxcl$names),
                          as.character(interferons$names),
                          as.character(silr$names),
                          as.character(others$names))
rownames(corr_cases) <- colnames(corr_cases)

corr_plot_cases<-ggcorrplot(corr_cases, 
                            hc.order = TRUE, 
                            ggtheme = ggplot2::theme_bw(), 
                            colors= c("dodgerblue4","white","red3"), 
                            tl.srt = 90,
                            #type="upper",
                            p.mat = corr_p_cases,
                            insig = "blank") 
saveRDS(corr_plot_cases,file=paste0(analysis, "/",analysis,"_correlation_cytokines_cases_edgar.rds"))

order_cluster_case<-ggplot_build(corr_plot_cases)$layout$panel_scales_x[[1]]$range$range

#corr_plot_controls

proteins_df<-rbind(interferons,interleukins,silr,ccl,cxcl,others) %>%
  mutate(names=factor(names, levels=order_cluster_case)) %>%
  arrange(names, order_cluster_case)

corr_controls <- soma %>%
  filter(days_symptom<=time_cutoff) %>%
  filter(.data[[analysis]]=="Control") %>%
  dplyr::select(c(interleukins$cols,
                  ccl$cols,
                  cxcl$cols,
                  interferons$cols,
                  silr$cols,
                  others$cols)) %>%
  stats::cor(x=.,use="pairwise.complete.obs",
             method=c(method))

corr_controls<-corr_controls[proteins_df$cols,proteins_df$cols]

colnames(corr_controls) <- proteins_df$names
rownames(corr_controls) <- proteins_df$names

corr_p_controls<- soma %>%
  filter(days_symptom<=time_cutoff) %>%
  filter(.data[[analysis]]=="Control") %>%
  dplyr::select(c(interleukins$cols,
                  ccl$cols,
                  cxcl$cols,
                  interferons$cols,
                  silr$cols,
                  others$cols)) %>%
  cor_pmat( method = method, exact=FALSE)

corr_p_controls<-corr_p_controls[proteins_df$cols,proteins_df$cols]

colnames(corr_p_controls) <- proteins_df$names
rownames(corr_p_controls) <- proteins_df$names


corr_plot_controls<-ggcorrplot(corr_controls,
                               ggtheme = ggplot2::theme_bw(), 
                               colors= c("dodgerblue4","white","red3"), 
                               tl.srt = 90,
                               #type="lower",
                               p.mat = corr_p_controls,
                               insig = "blank")
#corr_plot_controls
saveRDS(corr_plot_controls,file=paste0(analysis, "/",analysis,"_correlation_cytokines_controls_edgar.rds"))

#######test data for prediction######
mean_age<-mean(soma$age)
test_data_M<-data.frame(days_symptom=rep(seq(0,15,by=0.1),2)) %>%
  mutate(!!analysis:=factor(c(rep("Control", nrow(.)/2),rep("Case", nrow(.)/2)), levels=c("Control", "Case"))) %>%
  mutate(sex=factor(rep("M", nrow(.)))) %>%
  mutate(age=rep(mean_age, nrow(.)))

test_data_F<-data.frame(days_symptom=rep(seq(0,15,by=0.1),2)) %>%
  mutate(!!analysis:=factor(c(rep("Control", nrow(.)/2),rep("Case", nrow(.)/2)), levels=c("Control", "Case"))) %>%
  mutate(sex=factor(rep("F", nrow(.)))) %>%
  mutate(age=rep(mean_age, nrow(.)))

##### go by families ###### 
n_il<-nrow(interleukins)
n_cxcl<-nrow(cxcl)
n_ccl<-nrow(ccl)
n_ifn<-nrow(interferons)
n_others<-nrow(others)
n_silr<-nrow(silr)

######interleukins######
res_il<-list()
for(i in 1:n_il){
  res_il[[i]]<-gam_fxn(protein=interleukins$cols[i],
                       name=interleukins$names[i])
}
saveRDS(res_il,file=paste0(analysis, "/",analysis,"_interleukins_edgar.rds.gz", compress = TRUE))

# #if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# il_65F<-list()
# for(i in 1:n_il){
#   il_65F[[i]]<-res_il[[i]][[3]]
# }
# adj_plot_male<-ggarrange(plotlist=il_65F,
#                          labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_il],
#                          ncol = 7,
#                          nrow = 6,
#                          common.legend = TRUE,
#                          font.label = list(size = 12),
#                          legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred IL Levels, X y.o. F"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_male


######cxcl#####
res_cxcl<-list()
for(i in 1:n_cxcl){
  res_cxcl[[i]]<-gam_fxn(protein=cxcl$cols[i],
                         name=cxcl$names[i])
}

saveRDS(res_cxcl,file=paste0(analysis, "/",analysis,"_cxcl_edgar.rds.gz", compress = TRUE))

# #if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# cxcl_65F<-list()
# for(i in 1:n_il){
#   cxcl_65F[[i]]<-res_il[[i]][[3]]
# }
# 
# adj_plot_male<-ggarrange(plotlist=cxcl_65F,
#                          labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_il],
#                          ncol = 4, 
#                          nrow = 4,
#                          common.legend = TRUE, 
#                          font.label = list(size = 12),
#                          legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred cxcl Levels, X y.o. F"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_male

######silr#####
res_silr<-list()
for(i in 1:n_silr){
  res_silr[[i]]<-gam_fxn(protein=silr$cols[i],
                         name=silr$names[i])
}
saveRDS(res_silr,file=paste0(analysis, "/",analysis,"_silr_edgar.rds.gz", compress = TRUE))

# #if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# silr_65F<-list()
# for(i in 1:n_silr){
#   silr_65F[[i]]<-res_silr[[i]][[4]]
# }
# 
# adj_plot_male<-ggarrange(plotlist=silr_65F,
#                          labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_il],
#                          ncol = 7, 
#                          nrow = 6,
#                          common.legend = TRUE, 
#                          font.label = list(size = 12),
#                          legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred soluble interleukin receptor Levels, X y.o. F"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_male

#######ccl######
res_ccl<-list()
for(i in 1:n_ccl){
  res_ccl[[i]]<-gam_fxn(protein=ccl$cols[i],
                        name=ccl$names[i])
}
saveRDS(res_ccl,file=paste0(analysis, "/",analysis,"_ccl_edgar.rds.gz", compress = TRUE))

#if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# ccl_65F<-list()
# for(i in 1:n_ccl){
#   ccl_65F[[i]]<-res_ccl[[i]][[4]]
# }
# 
# adj_plot_female<-ggarrange(plotlist=ccl_65F,
#                            labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_ccl],
#                            ncol = 5, 
#                            nrow = 5,
#                            common.legend = TRUE, 
#                            font.label = list(size = 12),
#                            legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred CCL Levels, X y.o. M"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_female


#######interferons######
res_ifn<-list()
for(i in 1:n_ifn){
  res_ifn[[i]]<-gam_fxn(protein=interferons$cols[i],
                        name=interferons$names[i])
}
saveRDS(res_ifn,file=paste0(analysis, "/",analysis,"_ifn_edgar.rds.gz", compress = TRUE))

# #if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# ifn_65F<-list()
# for(i in 1:n_ifn){
#   ifn_65F[[i]]<-res_ifn[[i]][[4]]
# }
# 
# adj_plot_female<-ggarrange(plotlist=ifn_65F,
#                            labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_ifn],
#                            ncol = 5, 
#                            nrow = 4,
#                            common.legend = TRUE, 
#                            font.label = list(size = 12),
#                            legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred Interferon Levels, X y.o. M"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_female

#####others######
res_others<-list()
for(i in 1:n_others){
  res_others[[i]]<-gam_fxn(protein=others$cols[i],
                           name=others$names[i])
}
saveRDS(res_others,file=paste0(analysis, "/",analysis,"_others_edgar.rds.gz", compress = TRUE))


# #if you want to see the resulting plot, simply uncomment below. Case should be red, and control should be blue.
# others_65F<-list()
# for(i in 1:n_others){
#   others_65F[[i]]<-res_others[[i]][[4]]
# }
# 
# adj_plot_female<-ggarrange(plotlist=others_65F,
#                            labels = c(LETTERS,paste0(LETTERS,LETTERS))[1:n_others],
#                            ncol = 5, 
#                            nrow = 4,
#                            common.legend = TRUE, 
#                            font.label = list(size = 12),
#                            legend = "bottom") %>%
#   annotate_figure(.,
#                   fig.lab = paste0(analysis, ", Inferred Others Levels, X y.o. M"),
#                   fig.lab.face = "bold",
#                   fig.lab.pos="bottom.left")
# adj_plot_female


