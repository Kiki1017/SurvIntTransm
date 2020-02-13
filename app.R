
###################################################
### Shiny app for visualisation of surveillance for preventing sustained transmission of nCov-2019
###################################################

# setwd("~/Epidemiology/nConv2019/Shiny app")
# install.packages("shiny",dependencies=TRUE, repos='http://cran.rstudio.com/')

library(shiny)
library(ggplot2)
library(gridExtra, quietly = T)
library(rstan)



##########################
#  Utility functions
##########################


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
                    label = "Reduction in transmission",
                    min = 0, max = 1, step = 0.001, value = 0.1),
  	 selectInput('distribution', "Distribution:", choices = c("Exponential", "Gamma"),  selected = "Exponential")
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
   	beta<-3/3.8
    rho<-input$red.transm
    m<-input$max.imp.cases
    if (is.null(inFile))  {
      x <- read.csv("earlyCases.csv")
     } else {
		  x <- readdatafile()
     }
    df <- data.frame(id=1:nrow(x),days=x$X6)
    df$upper <- ifelse(df$days==0,0,df$days+1)
    df$lower <- ifelse(df$days==0,0,df$days-1)
    
   if(input$distribution=='Exponential'){
     
     fit1 <- stan(file = "exp_fit.stan",
                  data = list(N=nrow(df),low=df$lower,up=df$upper,lam_mean=mean(df$days)))
     
     res <- extract(fit1)
     param<-res$lambda
     
     days<-seq(0, 15, 0.1)
     ind<-sample(1:length(param), 50, replace=F)
     dots <- data.frame(x=days, y=dexp(days, param[ind[1]]))
     
     fig<-ggplot(data = df, aes(x=days, color='red')) +
       geom_histogram(aes(y=..density..), binwidth  = 1, fill="white", show.legend = FALSE, size=1.1) +
       geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5)  +
       labs(title= 'My histogram', x = 'Infection period (days)', y='') +ggtitle("Fitted time from symptom onset to hospitalisation")
     
     for(i in 2:length(ind)){
       dots <- data.frame(x=days, y=dexp(days, param[ind[i]]))
       fig<- fig +  geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5)
     }
     
     fig
     
   }  else {
     
     fit1 <- stan(file = "gam_fit.stan",
                  data = list(N=nrow(df),low=df$lower,up=df$upper))
     
     res <- extract(fit1)
     summary(fit1)
     
     param1<-res$alpha 
     param2<-res$beta
     
     days<-seq(0, 15, 0.1)
     
     ind<-sample(1:length(param1), 50, replace=F)
     dots <- data.frame(x=days, y=dgamma(days, shape=param1[ind[1]], rate=param2[ind[1]]))
     
     fig<-ggplot(data = df, aes(x=days, color='red')) +
       geom_histogram(aes(y=..density..), binwidth  = 1, fill="white", show.legend = FALSE, size=1.1) +
       geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5) +
       labs(title= 'My histogram', x = 'Infection period (days)', y='') +ggtitle("Fitted time from symptom onset to hospitalisation")
     
     for(i in 2:length(ind)){
       dots <- data.frame(x=days, y=dgamma(days, shape= param1[ind[i]], rate=param2[ind[i]]))
       fig<- fig +  geom_line(data = dots, aes(x=x,y=y), color= 'blue', size=0.5)
     }
     
     fig
     
  	}
  })

 } #End of server()

shinyApp(ui, server)
