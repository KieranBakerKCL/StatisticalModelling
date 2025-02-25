library(shiny)

ui <- fluidPage(
  titlePanel("QQ Plot & Distribution Comparison"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Choose Distribution:",
                  choices = list("Normal" = "normal",
                                 "Exponential" = "exponential",
                                 "Student's t" = "student",
                                 "Uniform" = "uniform",
                                 "Poisson" = "poisson",
                                 "Binomial" = "binomial")),
      checkboxInput("center_std", "Center & Standardise Distribution", value = FALSE),
      conditionalPanel(
        condition = "input.center_std == true",
        helpText("Note: The distribution has been centered and standardised.")
      ),
      sliderInput("n", "Sample Size:", min = 10, max = 1000, value = 100),
      sliderInput("quantile", "Select Quantile (0 to 1):", 
                  min = 0, max = 1, value = 0.5, step = 0.01),
      # Distribution-specific parameters:
      conditionalPanel(
        condition = "input.dist == 'exponential'",
        numericInput("rate", "Rate:", value = 1)
      ),
      conditionalPanel(
        condition = "input.dist == 'student'",
        numericInput("df", "Degrees of Freedom:", value = 5, min = 1)
      ),
      conditionalPanel(
        condition = "input.dist == 'uniform'",
        numericInput("min", "Minimum:", value = 0),
        numericInput("max", "Maximum:", value = 1)
      ),
      conditionalPanel(
        condition = "input.dist == 'poisson'",
        numericInput("lambda", "Lambda:", value = 3, min = 0)
      ),
      conditionalPanel(
        condition = "input.dist == 'binomial'",
        numericInput("size", "Number of Trials:", value = 10, min = 1),
        numericInput("prob", "Probability of Success:", 
                     value = 0.5, min = 0, max = 1, step = 0.1)
      ),
      wellPanel(
        h4("Area under Density Curves up to Selected Quantile"),
        htmlOutput("areaInfo")
      )
    ),
    mainPanel(
      plotOutput("qqPlot"),
      br(),
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  
  # Reactive function to return the theoretical parameters for the selected distribution
  getParams <- reactive({
    switch(input$dist,
           normal = list(mu = 0, sigma = 1),
           exponential = list(mu = 1 / input$rate, sigma = 1 / input$rate),
           student = list(mu = 0, sigma = if (input$df > 2) sqrt(input$df/(input$df - 2)) else 1),
           uniform = list(mu = (input$min + input$max) / 2, 
                          sigma = (input$max - input$min) / sqrt(12)),
           poisson = list(mu = input$lambda, sigma = sqrt(input$lambda)),
           binomial = list(mu = input$size * input$prob, 
                           sigma = sqrt(input$size * input$prob * (1 - input$prob)))
    )
  })
  
  # Render the QQ plot with updated axis limits and dashed lines to the axes.
  output$qqPlot <- renderPlot({
    set.seed(123)
    n <- input$n
    params <- getParams()
    data <- switch(input$dist,
                   normal = rnorm(n),
                   exponential = rexp(n, rate = input$rate),
                   student = rt(n, df = input$df),
                   uniform = runif(n, min = input$min, max = input$max),
                   poisson = rpois(n, lambda = input$lambda),
                   binomial = rbinom(n, size = input$size, prob = input$prob)
    )
    if (input$center_std) {
      data <- (data - params$mu) / params$sigma
    }
    
    # Determine limits based on the 0.01 and 0.99 quantiles for theoretical and sample data.
    theo_lower <- qnorm(0.01)
    theo_upper <- qnorm(0.99)
    sample_lower <- as.numeric(quantile(data, 0.01))
    sample_upper <- as.numeric(quantile(data, 0.99))
    
    xlims <- c(theo_lower, theo_upper)
    ylims <- c(min(theo_lower, sample_lower), max(theo_upper, sample_upper))
    
    qqnorm(data,
           xlim = xlims,
           ylim = ylims,
           main = paste("QQ Plot:", input$dist, if (input$center_std) "(Standardised)" else "", "vs Normal"),
           xlab = "Theoretical Quantiles (Normal)",
           ylab = "Sample Quantiles")
    
    if (input$center_std) {
      abline(0, 1, col = "red", lwd = 2)
    } else {
      qqline(data, col = "red", lwd = 2)
    }
    
    # Compute the quantile point for the sample and its corresponding normal quantile.
    sample_q <- as.numeric(quantile(data, input$quantile))
    normal_q <- qnorm(input$quantile)
    points(normal_q, sample_q, col = "blue", pch = 19, cex = 2)
    
    # Dashed lines from the quantile dot to the axes.
    segments(x0 = normal_q, y0 = sample_q, x1 = normal_q, y1 = 0, col = "darkgray", lty = 2)
    segments(x0 = normal_q, y0 = sample_q, x1 = 0, y1 = sample_q, col = "darkgray", lty = 2)
    
    # Solid black lines for the x-axis and y-axis.
    abline(h = 0, col = "black", lwd = 2)
    abline(v = 0, col = "black", lwd = 2)
  })
  
  # Render area information as bold, colored text.
  output$areaInfo <- renderUI({
    q <- input$quantile
    params <- getParams()
    if (input$dist %in% c("normal", "exponential", "student", "uniform")) {
      if (input$center_std) {
        q_selected_std <- (switch(input$dist,
                                  normal = qnorm(q, 0, 1),
                                  exponential = qexp(q, rate = input$rate),
                                  student = qt(q, df = input$df),
                                  uniform = qunif(q, min = input$min, max = input$max)) - params$mu) / params$sigma
        area_selected <- pnorm(q_selected_std)
      } else {
        area_selected <- switch(input$dist,
                                normal = pnorm(qnorm(q, 0, 1), 0, 1),
                                exponential = pexp(qexp(q, rate = input$rate), rate = input$rate),
                                student = pt(qt(q, df = input$df), df = input$df),
                                uniform = punif(qunif(q, min = input$min, max = input$max), min = input$min, max = input$max))
      }
      area_normal <- pnorm(qnorm(q))
      
      HTML(paste0("<b style='color:blue;'>Selected Distribution Area: ", round(area_selected, 3), "</b><br>",
                  "<b style='color:red;'>Standard Normal Area: ", round(area_normal, 3), "</b>"))
    } else {
      if (input$dist == "poisson") {
        q_val <- qpois(q, lambda = input$lambda)
        area_selected <- ppois(q_val, lambda = input$lambda)
      } else if (input$dist == "binomial") {
        q_val <- qbinom(q, size = input$size, prob = input$prob)
        area_selected <- pbinom(q_val, size = input$size, prob = input$prob)
      }
      area_normal <- pnorm(qnorm(q))
      
      HTML(paste0("<b style='color:blue;'>Selected Distribution Area: ", round(area_selected, 3), "</b><br>",
                  "<b style='color:red;'>Standard Normal Area: ", round(area_normal, 3), "</b>"))
    }
  })
  
  # Render the density (or PMF) plot with updated x and y limits.
  output$distPlot <- renderPlot({
    params <- getParams()
    q <- input$quantile
    
    # For continuous distributions, define helper functions.
    densFunc <- function(x) {
      switch(input$dist,
             normal = dnorm(x, mean = 0, sd = 1),
             exponential = dexp(x, rate = input$rate),
             student = dt(x, df = input$df),
             uniform = dunif(x, min = input$min, max = input$max)
      )
    }
    quantFunc <- function(p) {
      switch(input$dist,
             normal = qnorm(p, mean = 0, sd = 1),
             exponential = qexp(p, rate = input$rate),
             student = qt(p, df = input$df),
             uniform = qunif(p, min = input$min, max = input$max)
      )
    }
    
    if (input$dist %in% c("poisson", "binomial")) {
      # Discrete distributions.
      if (input$dist == "poisson") {
        x_vals <- 0:qpois(0.999, lambda = input$lambda)
        probs <- dpois(x_vals, lambda = input$lambda)
        q_val <- qpois(q, lambda = input$lambda)
        lower_bound <- min(qpois(0.01, lambda = input$lambda), floor(qnorm(0.01)))
        upper_bound <- max(ceiling(qpois(0.999, lambda = input$lambda)), ceiling(qnorm(0.999)))
      } else {
        x_vals <- 0:input$size
        probs <- dbinom(x_vals, size = input$size, prob = input$prob)
        q_val <- qbinom(q, size = input$size, prob = input$prob)
        lower_bound <- min(qbinom(0.01, size = input$size, prob = input$prob), floor(qnorm(0.01)))
        upper_bound <- max(input$size, ceiling(qnorm(0.999)))
      }
      
      if (input$center_std) {
        x_vals <- (x_vals - params$mu) / params$sigma
        q_val <- (q_val - params$mu) / params$sigma
        lower_bound <- min((if(input$dist == "binomial") (qbinom(0.01, size = input$size, prob = input$prob) - params$mu)/params$sigma 
                            else (qpois(0.01, lambda = input$lambda) - params$mu)/params$sigma), qnorm(0.01))
      }
      ylim_max <- max(probs)
      plot(x_vals, probs, type = "h", lwd = 2, col = "blue",
           main = if(input$center_std) 
             "Discrete Distribution (Standardised) vs Standard Normal" 
           else "Discrete Distribution vs Standard Normal",
           xlab = if(input$center_std) "Standardised x" else "x", ylab = "Probability",
           xlim = c(lower_bound, upper_bound), ylim = c(0, ylim_max))
      points(x_vals, probs, pch = 16, col = "blue")
      for (i in seq_along(x_vals)){
        if (x_vals[i] <= q_val) {
          rect(xleft = x_vals[i] - 0.3, ybottom = 0, 
               xright = x_vals[i] + 0.3, ytop = probs[i],
               col = rgb(0.53, 0.81, 0.98, 0.5), border = NA)
        }
      }
      # Overlay standard normal density.
      x_norm <- seq(floor(qnorm(0.01)), ceiling(qnorm(0.999)), length.out = 500)
      y_norm <- dnorm(x_norm)
      ylim_max <- max(ylim_max, max(y_norm))
      lines(x_norm, y_norm, col = "red", lwd = 2)
      x_shade_norm <- seq(floor(qnorm(0.01)), qnorm(q), length.out = 500)
      y_shade_norm <- dnorm(x_shade_norm)
      polygon(c(x_shade_norm[1], x_shade_norm, tail(x_shade_norm, 1)),
              c(0, y_shade_norm, 0),
              col = rgb(1, 0.75, 0.8, 0.5), border = NA)
      abline(v = q_val, col = "blue", lwd = 2, lty = 2)
      abline(v = qnorm(q), col = "red", lwd = 2, lty = 2)
      legend("topright", legend = c("Selected Distribution", "Standard Normal"),
             col = c("blue", "red"), lwd = 2)
      
    } else {
      # Continuous distributions.
      if (input$center_std) {
        orig_lower <- switch(input$dist,
                             normal = qnorm(0.01, 0, 1),
                             exponential = qexp(0.01, rate = input$rate),
                             student = qt(0.01, df = input$df),
                             uniform = qunif(0.01, min = input$min, max = input$max))
        orig_upper <- switch(input$dist,
                             normal = qnorm(0.999, 0, 1),
                             exponential = qexp(0.999, rate = input$rate),
                             student = qt(0.999, df = input$df),
                             uniform = qunif(0.999, min = input$min, max = input$max))
        lower_sel_std <- (orig_lower - params$mu) / params$sigma
        upper_sel_std <- (orig_upper - params$mu) / params$sigma
        lower_norm <- qnorm(0.01)
        upper_norm <- qnorm(0.999)
        lower_bound <- min(lower_sel_std, lower_norm)
        upper_bound <- max(upper_sel_std, upper_norm)
        
        y_seq <- seq(lower_bound, upper_bound, length.out = 500)
        density_selected <- params$sigma * switch(input$dist,
                                                  normal = dnorm(params$mu + params$sigma * y_seq, 0, 1),
                                                  exponential = dexp(params$mu + params$sigma * y_seq, rate = input$rate),
                                                  student = dt(params$mu + params$sigma * y_seq, df = input$df),
                                                  uniform = dunif(params$mu + params$sigma * y_seq, min = input$min, max = input$max))
        density_normal <- dnorm(y_seq)
        ylim_max <- max(density_selected, density_normal)
        q_orig <- quantFunc(q)
        q_selected_std <- (q_orig - params$mu) / params$sigma
        q_normal <- qnorm(q)
        plot(y_seq, density_selected, type = "l", col = "blue", lwd = 2,
             xlim = c(lower_bound, upper_bound), ylim = c(0, ylim_max),
             main = paste("Densities:", input$dist, "(Standardised) vs Standard Normal"),
             xlab = "Standardised x", ylab = "Density")
        idx <- which(y_seq <= q_selected_std)
        polygon(c(y_seq[idx][1], y_seq[idx], tail(y_seq[idx], 1)),
                c(0, density_selected[idx], 0),
                col = rgb(0.53, 0.81, 0.98, 0.5), border = NA)
        abline(v = q_selected_std, col = "blue", lwd = 2, lty = 2)
        lines(y_seq, density_normal, col = "red", lwd = 2)
        idx_norm <- which(y_seq <= q_normal)
        polygon(c(y_seq[idx_norm][1], y_seq[idx_norm], tail(y_seq[idx_norm], 1)),
                c(0, density_normal[idx_norm], 0),
                col = rgb(1, 0.75, 0.8, 0.5), border = NA)
        abline(v = q_normal, col = "red", lwd = 2, lty = 2)
        legend("topright", legend = c("Selected Distribution", "Standard Normal"),
               col = c("blue", "red"), lwd = 2)
        
      } else {
        orig_lower <- switch(input$dist,
                             normal = qnorm(0.01, 0, 1),
                             exponential = qexp(0.01, rate = input$rate),
                             student = qt(0.01, df = input$df),
                             uniform = qunif(0.01, min = input$min, max = input$max))
        orig_upper <- switch(input$dist,
                             normal = qnorm(0.999, 0, 1),
                             exponential = qexp(0.999, rate = input$rate),
                             student = qt(0.999, df = input$df),
                             uniform = qunif(0.999, min = input$min, max = input$max))
        lower_bound <- min(qnorm(0.01), orig_lower)
        upper_bound <- max(qnorm(0.999), orig_upper)
        q_selected <- quantFunc(q)
        x_sel <- switch(input$dist,
                        normal = seq(lower_bound, upper_bound, length.out = 500),
                        exponential = seq(lower_bound, upper_bound, length.out = 500),
                        student = seq(lower_bound, upper_bound, length.out = 500),
                        uniform = seq(lower_bound, upper_bound, length.out = 500))
        y_sel <- switch(input$dist,
                        normal = dnorm(x_sel, 0, 1),
                        exponential = dexp(x_sel, rate = input$rate),
                        student = dt(x_sel, df = input$df),
                        uniform = dunif(x_sel, min = input$min, max = input$max))
        x_norm <- seq(lower_bound, upper_bound, length.out = 500)
        y_norm <- dnorm(x_norm)
        ylim_max <- max(max(y_sel), max(y_norm))
        plot(x_sel, y_sel, type = "l", col = "blue", lwd = 2,
             xlim = c(lower_bound, upper_bound), ylim = c(0, ylim_max),
             main = "Densities: Selected Distribution (blue) & Standard Normal (red)",
             xlab = "x", ylab = "Density")
        idx_sel <- which(x_sel <= q_selected)
        polygon(c(x_sel[idx_sel][1], x_sel[idx_sel], tail(x_sel[idx_sel], 1)),
                c(0, y_sel[idx_sel], 0),
                col = rgb(0.53, 0.81, 0.98, 0.5), border = NA)
        abline(v = q_selected, col = "blue", lwd = 2, lty = 2)
        lines(x_norm, y_norm, col = "red", lwd = 2)
        polygon(c(x_norm[1], x_norm, tail(x_norm, 1)),
                c(0, y_norm, 0),
                col = rgb(1, 0.75, 0.8, 0.5), border = NA)
        abline(v = qnorm(q), col = "red", lwd = 2, lty = 2)
        legend("topright", legend = c("Selected Distribution", "Standard Normal"),
               col = c("blue", "red"), lwd = 2)
      }
    }
  })
}

shinyApp(ui = ui, server = server)
