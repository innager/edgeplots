library(shiny)
source("helperEdge.R")

filenames <- c("ParrAlg4aGam3.RData", "ParrAlg4aGam50.RData",
               "ParrAlg4aLNorm04.RData", "ParrAlg4aLNorm05.RData",
               "ParrAlg4aMixNorm.RData", "ParrAlg4aNorm.RData")

require(RColorBrewer)
col.choices <- c(brewer.pal(5, "Set1")[2:5], brewer.pal(4, "Dark2")[4])
cols <- col.choices[rep(c(1, 5, 2:4), 2)]
cols_transp <- adjustcolor(cols, alpha.f = 0.1)
cols_border <- adjustcolor(cols, alpha.f = 0.3)

n <- 10
#xalpha <- 0.05
m <- 2^(0:20)
lim <- 20

xlabs <- parse(text = paste("2^", 0:20, sep = ""))
ordstr <- c("1st", "2nd", "3rd", "4th", "5th")

#----------------------------------------------------------------------#

ui <- fluidPage(
  titlePanel("Simulations with Edgeworth expansions, Algorithm 4a"),
  fluidRow(
    column(4, 
      wellPanel(     
      selectInput(inputId = "dist", label = "Select a distribution", 
                  choices = c("Gamma"          = "Gam", 
                              "Log-Normal"     = "LNorm", 
                              "Mixture Normal" = "MixNorm", 
                              "Normal"         = "Norm")),
      uiOutput(outputId = "ui"),
      tags$br(),
      actionButton(inputId = "go", label = "Update distribution")
      )),
    
    column(8,
      wellPanel(  
      fluidRow(
        tags$p("Following tail diagnostic, the method chooses a specified order
                 of Edgeworth expansion or t-distribution", style = "color:SteelBlue;")
      ), 
      fluidRow(
        column(6, offset = 3,
           selectInput(inputId = "order", label = "Specify an order for Edgeworth expansion",
                      choices = c(2:5, "summary for orders 2 - 5")))
      ),
#      tags$br(),
      fluidRow(
        column(12, 
#      tags$br(),
#      tags$p("some text"),
#      tags$hr(),
          sliderInput(inputId = "yscale", label = "Choose y-scale", 
                       min = 0, max = 5, value = 3, step = 0.01, width = "100%"))
      )
    ))
  ),
  
  plotOutput(outputId = "plot1", height = 400)
  
)

#----------------------------------------------------------------------#

server <- function(input, output) {
  
  rvals <- reactiveValues(
    file  = filenames[1],
    mainr = bquote("Gamma" ~ (shape == 3) ~ ", " ~ n == 10)
  )
  
  observeEvent(input$go, {  
    param <- switch(input$dist, 
                    Gam     = input$gampar,
                    LNorm   = input$lnormpar,
                    MixNorm = "",
                    Norm    = "aNorm")
    rvals$file <- filenames[grepl(param, filenames) & grepl(input$dist, filenames)]
    rvals$mainr <- switch(input$dist,
                          Gam   = bquote("Gamma" ~ (shape == .(param)) ~ ", " ~ n == 10),
                          LNorm = bquote("Log-Normal" ~ (sigma^2 == .(as.numeric(param)/10)) 
                                           ~ ", " ~ n == 10),
                          MixNorm = bquote("Mixture Normal, " ~ n == 10),
                          Norm    = bquote("Normal Distribution, " ~ n == 10)
    )                          
  })

  output$ui <- renderUI({
    switch(input$dist,
           "Gam"   = selectInput(inputId = "gampar", label = "shape parameter",
                                 choices = c(3, 50)),
           "LNorm" = selectInput(inputId = "lnormpar", label = "log variance parameter",
                                 choices = c("0.4" = "04",
                                             "0.5" = "05")),
           "MixNorm" = selectInput(inputId = "mixpar", label = "parameters", 
                                   choices = "0.5 N(-1.5, 0.25) + 0.5 N(1.5, 4)"),
           "Norm"    = selectInput(inputId = "normpar", label = "parameters",
                                   choices = "")
    )       
  })
  
  output$plot1 <- renderPlot({
    load(rvals$file)
    ylim <- c(0, input$yscale)
    mainr <- rvals$mainr
    
    if (input$order == "summary for orders 2 - 5") {
      mainl <- bquote("Edgeworth 2 - 5 Orders with Tail Diagnostic" ~ "")
      legendstr <- paste("Edgeworth", ordstr[2:5], "order")
      
      par(mfrow = c(1, 2))
      par(mar = c(5, 4, 4, 0) + 0.1)	    
      plot(NULL, xlim = rev(range(m)), ylim = ylim, log = "x", main = mainl,
           xaxt = "n", xlab = "number of tests", ylab = "error rates")
      axis(1, at = m, labels = xlabs)
      lines(rev(m), rev(pertl), col = cols[1], lwd = 2)
      
      for (term in 1:4) {
        ord <- term + 1
        Psummary <- sdLog(Parrs[[term]])
        #	Psummary <- sdReg(Parrs[[term]])
        Pmean <- Psummary$Pmean
        Plow  <- Psummary$Plow
        Phi   <- Psummary$Phi
        lines(rev(m), rev(Pmean[, 1]), col = cols[ord], lwd = 2)
        polygon(x = c(m, rev(m)), y = c(Plow[, 1], rev(Phi[, 1])),
                col = cols_transp[ord], border = cols_border[ord])
      }
      
      par(mar = c(5, 2, 4, 2) + 0.1)	
      plot(NULL, xlim = range(m), ylim = ylim, log = "x", main = mainr,
           xaxt = "n", yaxt = "n", xlab = "number of tests", ylab = "")
      axis(1, at = m, labels = xlabs) 
      axis(2, labels = FALSE)    
      lines(m, pertr, col = cols[1], lwd = 2)
      
      for (term in 1:4) {
        ord <- term + 1
        Psummary <- sdLog(Parrs[[term]])
        #	Psummary <- sdReg(Parrs[[term]])
        Pmean <- Psummary$Pmean
        Plow  <- Psummary$Plow
        Phi   <- Psummary$Phi
        lines(m, Pmean[, 2], col = cols[ord], lwd = 2)
        polygon(x = c(m, rev(m)), y = c(Plow[, 2], rev(Phi[, 2])),
                col = cols_transp[ord], border = cols_border[ord])
      }	 
      legend(2^6, ylim[2], legend = c("t-distribution", legendstr), seg.len = 3,
             lwd = 10, col = c(NA, cols_transp[2:5]), bty = "n")
      legend(2^6, ylim[2], legend = rep("", 5), bty = "n", seg.len = 3,
             lwd = 2, col = cols)  
      
    } else {
  
    ord <- as.numeric(input$order)
    term <- ord - 1
    
    legendstr <- paste("Edgeworth", ordstr[ord], "order")
    mainl <- bquote("Edgeworth" ~ .(ordstr[ord]) ~ "Order with Tail Diagnostic")
    
    Psummary <- sdLog(Parrs[[term]])
    Pmean <- Psummary$Pmean
    Plow  <- Psummary$Plow
    Phi   <- Psummary$Phi
    
    par(mfrow = c(1, 2))
    par(mar = c(5, 4, 4, 0) + 0.1)	    
    plot(NULL, xlim = rev(range(m)), ylim = ylim, log = "x", main = mainl,
         xaxt = "n", xlab = "number of tests", ylab = "error rates")
    axis(1, at = m, labels = xlabs)
    lines(rev(m), rev(pertl), col = cols[1], lwd = 2)
    lines(rev(m), rev(Pmean[, 1]), col = cols[ord], lwd = 2)
    polygon(x = c(m, rev(m)), y = c(Plow[, 1], rev(Phi[, 1])),
            col = cols_transp[ord], border = cols_border[ord])
    
    par(mar = c(5, 2, 4, 2) + 0.1)	
    plot(NULL, xlim = range(m), ylim = ylim, log = "x", main = mainr,
         xaxt = "n", yaxt = "n", xlab = "number of tests", ylab = "")
    axis(1, at = m, labels = xlabs) 
    axis(2, labels = FALSE)    
    lines(m, pertr, col = cols[1], lwd = 2)
    lines(m, Pmean[, 2], col = cols[ord], lwd = 2)
    polygon(x = c(m, rev(m)), y = c(Plow[, 2], rev(Phi[, 2])),
            col = cols_transp[ord], border = cols_border[ord])
    legend(2^6, ylim[2], legend = c("t-distribution", legendstr), seg.len = 3,
           lwd = 10, col = c(NA, cols_transp[ord]), bty = "n")
    legend(2^6, ylim[2], legend = rep("", 2), bty = "n", seg.len = 3,
           lwd = 2, col = cols[c(1, ord)])
    }
  })
  
}

#----------------------------------------------------------------------#

shinyApp(ui = ui, server = server)




