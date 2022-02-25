        group<-c("+","-","+","-")
        os<-c("CR","CR","PD","PD")
        count<-c(6,15,26,103) 
   
        dat<-data.frame(group,os,count)
        u<-xtabs(count~group+os,data=dat)

        chisq = chisq.test(u)
        fisher = fisher.test(u)
        chisq$p.value
        fisher$p.value
