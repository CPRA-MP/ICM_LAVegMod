#############################################################
# model.R 
#
#############################################################
    remove(list=ls())
    assign("last.warning", NULL, envir = baseenv())

    stop_quietly <- function() {
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
    }

#############################################################
# Load libraries
#
#############################################################
library(tidyverse)
library(bbmle)

#############################################################
# Load the data and get it ready
#
#############################################################

data        = read.csv('sal_vs_depth_by_site_year.csv')
data$depthm = data$depth..ft. * 0.3048
data        = subset(data, !is.na(salinity) & !is.na(depthm) )

#############################################################
# 
#
#############################################################
negLogLikelihood = function(beta0, beta1, beta2, beta3, beta4, data )
{
    mu     = beta0 + beta1*data$salinity
    sigma  = beta2 + beta3*exp(beta4*data$salinity)
    L      = dnorm( data$depthm, mean = mu, sd = sigma )
    logLik = - (L %>% log() %>% sum())
    return(logLik)
}

estStart = list( beta0 = 0, beta1 = 0, beta2 = 0.02, beta3 = 0.2, beta4 = -0.2)
model1 = mle2( negLogLikelihood, estStart, fixed=list(beta0 = 0, beta1=0), data = list(data = data) )
model2 = mle2( negLogLikelihood, estStart, data = list(data = data) )

#############################################################
# Report results
#
#############################################################
estML1 = coef(model1) %>% as.list()
estML2 = coef(model2) %>% as.list()

cat('model.R: Msg: Initial Guess at model parameters (shown in blue)\n')
with( estStart, {
    cat('model.R: Msg: beta0 = ', beta0, '\n')
    cat('model.R: Msg: beta1 = ', beta1, '\n')
    cat('model.R: Msg: beta2 = ', beta2, '\n')
    cat('model.R: Msg: beta3 = ', beta3, '\n')
    cat('model.R: Msg: beta4 = ', beta4, '\n')
    cat('model.R: Msg: Likelihood = ',
        negLogLikelihood(beta0, beta1, beta2, beta3, beta4, data),'\n' )
} )
cat('\n')

cat('model.R: Msg: ML model1 (beta0, beta1 fixed) parameters (shown in green)\n')
with( estML1, {
    cat('model.R: Msg: beta0 = ', beta0, '\n')
    cat('model.R: Msg: beta1 = ', beta1, '\n')
    cat('model.R: Msg: beta2 = ', beta2, '\n')
    cat('model.R: Msg: beta3 = ', beta3, '\n')
    cat('model.R: Msg: beta4 = ', beta4, '\n')
    cat('model.R: Msg: Likelihood = ',
        negLogLikelihood(beta0, beta1, beta2, beta3, beta4, data),'\n' )
})
cat('\n')

cat('model.R: Msg: ML model2 parameters (shown in green)\n')
with( estML2, {
    cat('model.R: Msg: beta0 = ', beta0, '\n')
    cat('model.R: Msg: beta1 = ', beta1, '\n')
    cat('model.R: Msg: beta2 = ', beta2, '\n')
    cat('model.R: Msg: beta3 = ', beta3, '\n')
    cat('model.R: Msg: beta4 = ', beta4, '\n')
    cat('model.R: Msg: Likelihood = ',
        negLogLikelihood(beta0, beta1, beta2, beta3, beta4, data),'\n' )
})
cat('\n')


#############################################################
# Visualize the data and the fit
#
#############################################################
fitStart = with( estStart, {
    tibble( sal = data$salinity %>% unique() %>% sort(), 
                    mean = beta0 + beta1 * sal,
                    sd   = beta2 + beta3*exp(beta4 * sal),
                    uc   = mean + sd * 1.96,
                    lc   = mean - sd * 1.96 )
})

fitML1 = with ( estML1, {
    tibble( sal = data$salinity %>% unique() %>% sort(), 
                    mean = beta0 + beta1 * sal,
                    sd   = beta2 + beta3*exp(beta4 * sal),
                    uc   = mean + sd * 1.96,
                    lc   = mean - sd * 1.96 )
})

fitML2 = with ( estML2, {
    tibble( sal = data$salinity %>% unique() %>% sort(), 
                    mean = beta0 + beta1 * sal,
                    sd   = beta2 + beta3*exp(beta4 * sal),
                    uc   = mean + sd * 1.96,
                    lc   = mean - sd * 1.96 )
})
    
pl = ggplot(data, aes(x = salinity, y = depthm) ) +
        geom_point() +
        #geom_line(data = fitStart, aes(sal, mean), col='blue') +
        #geom_line(data = fitStart, aes(sal, uc), col='cyan') +
        #geom_line(data = fitStart, aes(sal, lc), col='cyan') +
        #geom_line(data = fitML1, aes(sal, mean), col='dark greed') +
        geom_line(data = fitML1, aes(sal, uc), col='green', size=2) +
        geom_line(data = fitML1, aes(sal, lc), col='green', size=2) +
        geom_line(data = fitML2, aes(sal, uc), col='blue', size=2) +
        geom_line(data = fitML2, aes(sal, lc), col='blue', size=2) +
        theme_classic() +
        xlab('Salinity (ppt)') +
        ylab('Depth (m)')
plot(pl)

#############################################################
# Check results
#
#############################################################
data$mean = with( estML1, { beta0 + beta1*data$salinity   })
data$sd   = with( estML1, { beta2 + beta3*exp( beta4*data$salinity)   })
data$uc   = data$mean + data$sd * 1.96
data$lc   = data$mean - data$sd * 1.96
data$inCI = data$lc < data$depthm & data$depthm < data$uc

cat('model.R: Msg: % of data within CI = ',sum(data$inCI)/nrow(data) * 100,'% \n')







