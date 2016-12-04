library(shiny)
source("helperEdge.R")

filenames <- c("ParrayGam3lr.RData", "ParrayGam50lr.RData",
               "ParrayLNorm04lr.RData", "ParrayLNorm05lr.RData",
               "ParrayMixNormlr.RData", "ParrayNormlr.RData")

require(RColorBrewer)
col.choices <- c(brewer.pal(5, "Set1")[2:5], brewer.pal(4, "Dark2")[4])
cols <- col.choices[rep(c(1, 5, 2:4), 2)]
cols_transp <- adjustcolor(cols, alpha.f = 0.1)
cols_border <- adjustcolor(cols, alpha.f = 0.3)

n <- 10
alpha <- 0.05
m <- 2^(0:20)
lim <- 20

xlabs <- parse(text = paste("2^", 0:20, sep = ""))
ordstr <- c("1st", "2nd", "3rd", "4th", "5th")

#----------------------------------------------------------------------#

ui <- fluidPage(
  titlePanel("Simulations with empirical Edgeworth expansions"),
  fluidRow(
    column(4, 
      wellPanel(     
      selectInput(inputId = "dist", label = "Select a distribution", 
                  choices = c("Gamma"          = "Gam", 
                              "Log-Normal"     = "LNorm", 
                              "Mixture Normal" = "MixNorm", 
                              "Normal"         = "Norm")),
      uiOutput(outputId = "uipar"),
      tags$br(),
      actionButton(inputId = "go", label = "Update distribution")
      )),
    
    column(8,
      wellPanel(  
      fluidRow(
        tags$p("Compare error rates for different orders of approximation 
               (truncated Edgeworth expansions).", 
               style = "color:SteelBlue;")
      ), 
      fluidRow(
        column(6, offset = 3,
          uiOutput(outputId = "uicond"))       
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
    maind = bquote("Gamma" ~ (shape == 3) ~ ", " ~ n == 10),  # dist name for plot
    conds = uniq.cond
  )
  
  observeEvent(input$go, {  
    param <- switch(input$dist, 
                    Gam     = input$gampar,
                    LNorm   = input$lnormpar,
                    MixNorm = "",
                    Norm    = "yNorm")
    rvals$file <- filenames[grepl(param, filenames) & grepl(input$dist, filenames)]
    rvals$maind <- switch(input$dist,
                          Gam   = bquote("Gamma" ~ (shape == .(param)) ~ ", " ~ n == 10),
                          LNorm = bquote("Log-Normal" ~ (sigma^2 == .(as.numeric(param)/10)) 
                                           ~ ", " ~ n == 10),
                          MixNorm = bquote("Mixture Normal, " ~ n == 10),
                          Norm    = bquote("Normal Distribution, " ~ n == 10))
    load(rvals$file) 
    rvals$conds <- getUniqCond(Parray)
  })

  output$uipar <- renderUI({
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
  
  output$uicond <- renderUI({
    selectInput(inputId = "cond", label = "Tail diagnostic combination",
                choices = rvals$conds)
  })
  
  output$plot1 <- renderPlot({
    load(rvals$file)
    S <- dim(Parray)[3]
    combs <- getCond(Parray)
    combNA <- combs$comb
    tbl.cond <- combs$tbl
    cond <- input$cond
    ylim <- c(0, input$yscale)
    mainl <- bquote(atop("combination: left tail" ~ .(substr(cond, 1, 4)),
                         .(rvals$maind)))
    mainr <- bquote(atop("combination: right tail" ~ .(substr(cond, 6, 9)),
                         "proportion of samples " ~ .(tbl.cond[cond]/S)))
    
    Pcond <- Parray[, , combNA[, cond]]
    #	Psummary <- sdReg(Pcond)
    Psummary <- sdLog(Pcond)
    Pmean <- Psummary$Pmean
    Plow  <- Psummary$Plow
    Phi   <- Psummary$Phi
    
    par(mfrow = c(1, 2))
    par(mar = c(5, 4, 4, 0) + 0.1)	    
    plot(NULL, xlim = rev(range(m)), ylim = ylim, log = "x", main = mainl,
         xaxt = "n", xlab = "number of tests", ylab = "error rates")
    axis(1, at = m, labels = xlabs)
    lines(rev(m), rev(Pmean[, 1]), col = cols[1], lwd = 2)
    for (i in 2:5) {
      lines(rev(m), rev(Pmean[, i]), col = cols[i], lwd = 2)
      polygon(x = c(m, rev(m)), y = c(Plow[, i], rev(Phi[, i])),
              col = cols_transp[i], border = cols_border[i])
    }
    
    par(mar = c(5, 2, 4, 2) + 0.1)	
    plot(NULL, xlim = range(m), ylim = ylim, log = "x", main = mainr,
         xaxt = "n", yaxt = "n", xlab = "number of tests", ylab = "")
    axis(1, at = m, labels = xlabs) 
    axis(2, labels = FALSE)    
    lines(m, Pmean[, 6], col = cols[6], lwd = 2)
    for (i in 7:10) {
      lines(m, Pmean[, i], col = cols[i], lwd = 2)
      polygon(x = c(m, rev(m)), y = c(Plow[, i], rev(Phi[, i])),
              col = cols_transp[i], border = cols_border[i])
    }
    legend("topleft", legend = paste(0:4, "-term", sep = ""), seg.len = 3,
           lwd = 10, col = c(NA, cols_transp[-1]), bty = "n")
    legend("topleft", legend = rep("", 5), bty = "n", seg.len = 3,
           lwd = 2, col = cols)
  })
  
}

#----------------------------------------------------------------------#

shinyApp(ui = ui, server = server)




