library(shiny)
library(factoextra)

# Define UI for application that draws a dendrogram
ui <- fluidPage(

    # Application title
    titlePanel("Clustering by genomic profiles"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("cex",
                        "Label size:",
                        min = 0.1,
                        max = 1,
                        value = 0.6)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("dendrogram")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$dendrogram <- renderPlot({
        # generate bins based on input$bins from ui.R
        # x    <- faithful[, 2]
        # bins <- seq(min(x), max(x), length.out = input$bins + 1)
        # # draw the histogram with the specified number of bins
        # hist(x, breaks = bins, col = 'darkgray', border = 'white',
        #      xlab = 'Waiting time to next eruption (in mins)',
        #      main = 'Histogram of waiting times')
        
        # loading file ------------------------------------------------------------

        sample_table <-
          read.delim(
            "data/genotypes_exp_val_genes.tsv",
            header = T,
            stringsAsFactors = F,
            sep = "\t",
            row.names = "v"
          )
        # base R
        sample_table <- sample_table[, colSums(is.na(sample_table)) == 0]
        # transform to numeric
        sample_df <-
          as.data.frame(sapply(sample_table, as.numeric), 
                        row.names = rownames(sample_table))
        # extract controls
        sample_df_2 = sample_df[grep("CONTROL",
                                     colnames(sample_df),
                                     invert = TRUE ,
                                     value = T)]


        
        # jaccard index -----------------------------------------------------------
        
        # do all possible pair combinations
        combinations = combn(colnames(sample_df_2), 2)
        (combinations)[0:6]
        
        # generate similarity matrix (samples vs samples)
        sim_mat = matrix(
          NA,
          ncol = ncol(sample_df_2),
          nrow = ncol(sample_df_2),
          dimnames = list(colnames(sample_df_2), colnames(sample_df_2))
        )
        head(sim_mat)
        
        for (i in 1:ncol(combinations)) {
          sim_mat[combinations[2, i], combinations[1, i]] =
            (sum(sample_df_2[, combinations[1, i]] != 0  &
                   sample_df_2[, combinations[1, i]] ==
                   sample_df_2[, combinations[2, i]])) /
            sum(sample_df_2[, combinations[2, i]] != 0 |
                  sample_df_2[, combinations[1, i]] != 0)
          sim_mat[combinations[1, i], combinations[2, i]] =
            (sum(sample_df_2[, combinations[1, i]] != 0  &
                   sample_df_2[, combinations[1, i]] ==
                   sample_df_2[, combinations[2, i]])) /
            sum(sample_df_2[, combinations[2, i]] != 0 |
                  sample_df_2[, combinations[1, i]] != 0)
        }
        diag(sim_mat) = 1
        sim_mat[is.nan(sim_mat) | is.na(sim_mat)] <- 0
        dis_mat = 1 - sim_mat
        

        # hierarchical clustering -------------------------------------------------

        index = "ball"
        res <-
          NbClust(
            data = dis_mat_2,
            diss = as.dist(dis_mat_2),
            distance = NULL,
            method = "ward.D2",
            index = index,
            min.nc = 3
          )
        fit = hclust(d = as.dist(dis_mat_2), method = "ward.D2")
        
        # plot dendrogram
        fviz_dend(
          x = as.dendrogram(fit),
          k = res$Best.nc[1][[1]],
          k_colors =rainbow(length(1:(res$Best.nc[1][[1]]))),
          rect = TRUE,
          rect_fill = TRUE,
          cex = input$cex,
          color_labels_by_k = F,
          labels_track_height = 0.2
        )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
