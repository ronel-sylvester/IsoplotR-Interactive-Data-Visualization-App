library(shiny)
library(shinyFiles)
library(ggplot2)  # for the diamonds dataset
library(DT)
library(IsoplotRgui)
library(IsoplotR)
library(reshape2)

ui <- fluidPage(
  title = "Examples of DataTables",
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.dataset === "mtcars"',
        helpText("Click the column header to sort a column.")
      ),
      conditionalPanel(
        'input.dataset === "iris"',
        helpText("Display 5 records by default.")
      ),
      shinyDirButton('directory2', 'Folder select', 'Please select a folder')
    ),
    mainPanel(
      textOutput("text"),
      verbatimTextOutput('dblclickIndex'),
      tabsetPanel(
        id = 'dataset',
        tabPanel("Directory", DT::dataTableOutput("mytable1")),
        tabPanel("Data", radioButtons("disp", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "head"),
                 tableOutput('contents'),
        ),
        tabPanel("Brush & Concordia", checkboxGroupInput("cnt", "Select:",
                                    choices = c("Pb206", "Pb207", "Th232", "U238", 
                                                "Final206_238", "Final207_206", 
                                                "Pb204_CPS", "U238_CPS", "Raw_U_Th_Ratio", 
                                                "err238_206", "err207_206"),
                                    selected = "Final206_238", inline = TRUE),
                 splitLayout(cellWidths = c("50%", "50%"),
                             plotOutput("plot1", brush = "plot1_brush"),
                             plotOutput("plot2")
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  volumes <- getVolumes()
  transfer <- NULL
  t238 <- 0.003   #measurement time in sec
  t206 <- 0.010
  t207 <- 0.010
  observe({
    shinyDirChoose(input, "directory2", roots=volumes, session=session)
    
    if(!is.null(input$directory2)){
      path1 <- parseDirPath(volumes, input$directory2)
      file_info <- data.frame(file.info(list.files(path1, full.names=TRUE)))
      file_info <- cbind(dir_paths = rownames(file_info), file_info)  
      file_names <- data.frame(list.files(path1))
      colnames(file_names)[1] <- "file_names"
      diamonds2 <- cbind(file_names, file_info, row.names = NULL)
      diamonds2[] <- lapply(diamonds2, as.character)
      
      
      
      output$mytable1 <- renderDataTable({
        DT_Out <- datatable(diamonds2
                            ,rownames = F
                            ,callback = JS("                                                
                                          table.on('dblclick', 'tr', function() {
                                              var tr = $(this);
                                              var info_out = table.row( this ).data();
                                              Shiny.onInputChange('dblclickIndexJS', info_out[1]);
                                          });"
                            )
        )
        return(DT_Out)
      })

        
      output$dblclickIndex <- renderText({
        UI_out <- input$dblclickIndexJS
        return(paste("Pathway:", UI_out))
      })
      
      
  
      
      output$contents <- renderTable({
          
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
        
        UI_out <- input$dblclickIndexJS
        uploaded  <- read.delim(UI_out, header=FALSE)
        names(uploaded) <- as.character(unlist(uploaded[1,]))
        uploaded <- uploaded[-1,]
        uploaded[] <- lapply(uploaded, function(x) as.numeric(as.character(x)))
        cts238 <- uploaded$U238*t238   # counts per sec U238 * count time of U238  = total counts U238 
        cts206 <- uploaded$Pb206*t206
        cts207 <- uploaded$Pb207*t207
        rsd206 <- sqrt(cts206)/cts206      # counting statistics uncertainties as relative standard deviation.
        rsd238 <- sqrt(cts238)/cts238  
        rsd207 <- sqrt(cts207)/cts207
        uploaded$Final206_238 = 1/uploaded$Final206_238
        err238_206 = sqrt((rsd238)^2+(rsd206)^2) * uploaded$Final206_238
        err207_206 = sqrt((rsd207)^2+(rsd206)^2) * uploaded$Final207_206
        filtered = as.data.frame(cbind(uploaded, err238_206, err207_206))
        
        if(input$disp == "head") {
          return(head(uploaded))
        }
        else {
          return(uploaded)
        }
      })
      
      output$plot1 = renderPlot({
        
        UI_out <- input$dblclickIndexJS
        uploaded  <- read.delim(UI_out, header=FALSE)
        names(uploaded) <- as.character(unlist(uploaded[1,]))
        uploaded <- uploaded[-1,]
        uploaded[] <- lapply(uploaded, function(x) as.numeric(as.character(x)))
        cts238 <- uploaded$U238*t238   # counts per sec U238 * count time of U238  = total counts U238 
        cts206 <- uploaded$Pb206*t206
        cts207 <- uploaded$Pb207*t207
        rsd206 <- sqrt(cts206)/cts206      # counting statistics uncertainties as relative standard deviation.
        rsd238 <- sqrt(cts238)/cts238  
        rsd207 <- sqrt(cts207)/cts207
        uploaded$Final206_238 = 1/uploaded$Final206_238
        err238_206 = sqrt((rsd238)^2+(rsd206)^2) * uploaded$Final206_238
        err207_206 = sqrt((rsd207)^2+(rsd206)^2) * uploaded$Final207_206
        filtered = as.data.frame(cbind(uploaded, err238_206, err207_206))
        
        plot.data <- melt(filtered, id.vars = 'ElapsedTime_s')
        plot.data <- plot.data[plot.data$variable %in% input$cnt, ]
        ggplot(plot.data) +
          geom_point(mapping = aes(x = ElapsedTime_s, y = value, colour = variable)) + scale_y_log10()+
          geom_line(mapping = aes(x = ElapsedTime_s, y = value, colour = variable)) + scale_y_log10()+
          labs (x = "Time(sec)", title = "17YN04A-2") + 
          scale_colour_discrete(guide=FALSE)
      })
      
      output$plot2 = renderPlot({
        
        UI_out <- input$dblclickIndexJS
        uploaded  <- read.delim(UI_out, header=FALSE)
        names(uploaded) <- as.character(unlist(uploaded[1,]))
        uploaded <- uploaded[-1,]
        uploaded[] <- lapply(uploaded, function(x) as.numeric(as.character(x)))
        cts238 <- uploaded$U238*t238   # counts per sec U238 * count time of U238  = total counts U238 
        cts206 <- uploaded$Pb206*t206
        cts207 <- uploaded$Pb207*t207
        rsd206 <- sqrt(cts206)/cts206      # counting statistics uncertainties as relative standard deviation.
        rsd238 <- sqrt(cts238)/cts238  
        rsd207 <- sqrt(cts207)/cts207
        uploaded$Final206_238 = 1/uploaded$Final206_238
        err238_206 = sqrt((rsd238)^2+(rsd206)^2) * uploaded$Final206_238
        err207_206 = sqrt((rsd207)^2+(rsd206)^2) * uploaded$Final207_206
        filtered = as.data.frame(cbind(uploaded, err238_206, err207_206))
        
        transfer = brushedPoints(filtered, input$plot1_brush, xvar = "ElapsedTime_s", yvar = "Final206_238")
        transfer = as.data.frame(subset(transfer, select = c(Final206_238, err238_206, Final207_206, err207_206)))
        transfer = as.data.frame(transfer[ -c(1:14), ])
        concordia(read.data(transfer,method = 'U-Pb', format = 2))  
      })
    }
  })
  
  ## sorted columns are colored now because CSS are attached to them
  #output$mytable2 <- DT::renderDataTable({
  #  DT::datatable(mtcars, options = list(orderClasses = TRUE))
  #})
  
  # customize the length drop-down menu; display 5 rows per page by default
  #output$mytable3 <- DT::renderDataTable({
  #  DT::datatable(iris, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  #})
  
}

shinyApp(ui, server)
#install.packages("dplyr")
