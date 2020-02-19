
###################################################
### Shiny app for visualisation of surveillance for preventing sustained transmission of nCov-2019
###################################################


library(shiny)
library(ggplot2)
library(gridExtra, quietly = T)
library(grid)
library(rstan, quietly = T)
library(rootSolve)
library(loo)




###############################
#  Utility functions & paramaters
###############################

R.SARS<-3
gamma.SARS<-1/3.8
beta<- R.SARS*gamma.SARS

p.exp<-function(gamma, rho, m){
  pmax(0, 1-(1/(beta*(1-rho)/gamma))^m)
}

# k<-param1[1]; r<-param2[1]; rho<-0.2; m<-10
p.gam<-function(k, r, rho, m){
  theta<-1/r
  R<-beta*(1-rho)*k*theta
  fun<-function(p) (1-p)*(1+p*R/k)^k-1
  p.sol<-uniroot.all(fun, c(0, 1))
  p<-max(p.sol)
  if(!is.numeric(p)) p<-0
  q<-1-p
  return(pmax(0,1-q^m))
}

##########################
#  shiny interface
##########################

ui = fluidPage(

	#The title
	titlePanel("Influence of surveillance for preventing sustained transmission in new locations"),
	
	
	#The sidebar for parameter input
	sidebarPanel(
	  p("Based on the analysis at R.N. Thompson (2020)"),
	tags$a(href="https://www.mdpi.com/2077-0383/9/2/498", "doi: 10.3390/jcm9020498",  target="_blank"),
	br(),
	br(),
	  p("Input data file: it should be a column with calendar days from symptom onset to hospitalisation."),
	  fileInput("file", "Choose CSV File",
	                       accept = c(
	                         "text/csv",
	                         "text/comma-separated-values,text/plain",
	                         ".csv")
	),
	
	sliderInput('max.imp.cases', 
                    label = "Number of imported cases",
                    min = 1, max = 100, step = 1, value = 10),
    sliderInput('red.transm', 
                    label = "Reduction in transmission period",
                    min = 0, max = 1, step = 0.01, value = 0.1),
  	 selectInput('distribution', "Distribution of time from symptoms onset to hospitalisation:", choices = c("Exponential", "Gamma"),  selected = "Exponential")
 	),
	#Main panel for figures and equations
	mainPanel(     
	# Output     
		plotOutput('plot'))
) #End of ui()


server = function(input, output) {
 
   readdatafile <- reactive({
    inFile <- input$file
    if (is.null(inFile))  return(NULL)
    
    tbl <- read.csv(inFile$datapath)
    
    return(tbl)
  })
  
   output$plot <- renderPlot({
   	inFile <- input$file
    rho<-input$red.transm
    m<-input$max.imp.cases
    
    if (is.null(inFile))  {
      x <- read.csv("earlyCases.csv")
     } else {
		  x <- readdatafile()
     }
    df <- data.frame(id=1:nrow(x),days=x[,1])
    df$upper <- ifelse(df$days==0,0,df$days+1)
    df$lower <- ifelse(df$days==0,0,df$days-1)
    
    fig1<-ggplot(data = df, aes(x=days, color='red')) +
      geom_histogram( binwidth  = 1, fill="white", show.legend = FALSE, size=1.1) +
    ggtitle('Observed time from symptoms onset to hospitalisation')  + xlab('Days')+ ylab("Counts")
    
    days<-seq(0.1, 15, 0.1)
    
  
   if(input$distribution=='Exponential'){
     
     fit1 <- stan(file = "exp_fit.stan",
                  data = list(N=nrow(df),low=df$lower,up=df$upper,lam_mean=mean(df$days)))
     
     res <- extract(fit1)
     param<-res$lambda
     
     ind<-sample(1:length(param), 50, replace=F)
     dots <- data.frame(x=days, y=dexp(days, param[ind[1]]))
       
     fig2 <-ggplot(data = dots, aes(x=x,y=y))+ geom_line(color= 'blue', size=0.5)  +
       labs(x = 'Days', y='PDF') + ggtitle('Fitted time from symptoms onset to hospitalisation') 
     
     for(i in 2:length(ind)){
       dots <- data.frame(x=days, y=dexp(days, param[ind[i]]))
       fig2<- fig2 +  geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5)
     }
     
     prob0<-p.exp(param, 0, m)
     prob1<-p.exp(param, rho, m)
     lm0<-as.numeric(quantile(prob0,c(0.05, 0.95)))
     lm1<-as.numeric(quantile(prob1,c(0.05, 0.95)))
     
     ans<-data.frame(cases=c("No intensification", "With intensification"), prob=c(mean(prob0), mean(prob1)), low=c(lm0[1], lm1[1]), upp=c(lm0[2], lm1[2]))
     ans$cases<-as.factor(ans$cases)
     
     log_lik1 <- extract_log_lik(fit1, merge_chains = FALSE)
     rel_n_eff <- relative_eff(exp(log_lik1))
     loo.sum<-loo(log_lik1, r_eff = rel_n_eff, cores = 2)
     
     fig3 <- ggplot(ans, aes(x=cases, y=prob)) +    geom_bar(color="black", fill="blue", stat = "identity")  + labs(x="Surveilance", y = "Prob(Sustained transmission)") + theme_bw() +ggtitle("Surveilance intensification")  +ylim(c(0,1)) # +  geom_errorbar(aes(ymin=low, ymax=upp), width=.2,  position=position_dodge(.9)) + theme(text = element_text(size=12))
     fig4<-tableGrob(loo.sum$estimates)
     }  else {
     
     fit1 <- stan(file = "gam_fit.stan",
                  data = list(N=nrow(df),low=df$lower,up=df$upper))
     
     res <- extract(fit1)
     
     param1<-res$alpha 
     param2<-res$beta
     
     ind<-sample(1:length(param1), 50, replace=F)
     dots <- data.frame(x=days, y=dgamma(days, shape=param1[ind[1]], rate=param2[ind[1]]))
     
     fig2 <-ggplot(data = dots, aes(x=x,y=y))+ geom_line(color= 'blue', size=0.5)  +
       labs(x = 'Days', y='PDF') + ggtitle('Fitted infection period') 
     
     for(i in 2:length(ind)){
       dots <- data.frame(x=days, y=dgamma(days, shape= param1[ind[i]], rate=param2[ind[i]]))
       fig2<- fig2 +  geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5)
     }
     

     prob0<-sapply(1:length(param1), function(a) p.gam(param1[a], param2[a], 0, m))
     prob1<-sapply(1:length(param1), function(a) p.gam(param1[a], param2[a], rho, m))
     lm0<-as.numeric(quantile(prob0,c(0.05, 0.95)))
     lm1<-as.numeric(quantile(prob1,c(0.05, 0.95)))
     
     ans<-data.frame(cases=c("No intensification", "With intensification"), prob=c(mean(prob0), mean(prob1)), low=c(lm0[1], lm1[1]), upp=c(lm0[2], lm1[2]))
     ans$cases<-as.factor(ans$cases)
     
     log_lik1 <- extract_log_lik(fit1, merge_chains = FALSE)
     rel_n_eff <- relative_eff(exp(log_lik1))
     loo.sum<-loo(log_lik1, r_eff = rel_n_eff, cores = 2)
     
     
     fig3 <- ggplot(ans, aes(x=cases, y=prob)) +    geom_bar(color="black", fill="blue", stat = "identity")  + labs(x="Surveilance", y = "Prob(Sustained transmission)") + theme_bw() +ggtitle("Surveilance intensification")  +ylim(c(0,1)) # +  geom_errorbar(aes(ymin=low, ymax=upp), width=.2,  position=position_dodge(.9)) + theme(text = element_text(size=12))
     fig4<-tableGrob(loo.sum$estimates)
   }
    grid.arrange(fig1, fig4, fig2 ,fig3, nrow=2, widths=c(2,1.5))
  })
   

 } #End of server()

shinyApp(ui, server)
