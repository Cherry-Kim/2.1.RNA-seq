glm <- function(){
        data <- read.table("CRPD_glm_input.txt",header=T, stringsAsFactors=F)
        #levels(factor(data$response))
#best_response   response        sample gene1     gene2
#A      0       SAM1 13.61937        13.37922        
#B      1       SAM0 13.17767        12.04530        

        result <- glm(response ~ gene1, family=binomial(), data=data)
}

multinorm <- function(){
        library(nnet)
        test <- multinom(response ~ gene1 + gene2, data = data)
        summary(test)
        #pvalue calculate
        z = summary(test)$coefficients/summary(test)$standard.errors
        p_gene1_gene2 = (1-pnorm(abs(z), 0, 1))*2
}
