#H0: There is no difference between the two groups
#H1: There is a survival differential between the two groups
KM <- function(){
        #install.packages(c("survival", "survminer"))
        library("survival") # computing survival analyses
        library("survminer") #for visualizing survival analysis results

        data <- read.table("TMMnorm.level.cutpointr.out",header=T, stringsAsFactors=F, sep='\t')
        levels(factor(data$cenos))      #"0" "1"
        
        fit <- survfit(Surv(os , cenos2) ~ group, data = data, type = "kaplan-meier")
        survp <- ggsurvplot(fit, data=data, palette = c("blue","red"), pval=T, risk.table=T, conf.int=F, break.time.by= 5, legend.title = "Expression level", risk.table.fontsize = 7, surv.plot.height = 0.45, legend="right", censor = TRUE, pval.coord = c(0, 0.03), tables.height = 0.2)
        ggsave(file = "ggsurv.pdf", print(survp))
        
        #Log Rank Test in R-Survival Curve Comparison
        survdiff(formula = Surv(os , cenos) ~ group, data = data))
        }
