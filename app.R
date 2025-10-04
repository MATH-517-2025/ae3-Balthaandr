library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# ---------------------------
# Helper functions (self-contained)
# ---------------------------

make_blocks <- function(X, Y, N){
  n <- length(X)
  ord <- order(X)
  Xs <- X[ord]; Ys <- Y[ord]
  sizes <- rep(floor(n / N), N)
  remainder <- n - sum(sizes)
  if(remainder > 0) sizes[1:remainder] <- sizes[1:remainder] + 1
  blocks <- vector("list", N)
  start <- 1
  for(j in seq_len(N)){
    idx_local <- start:(start + sizes[j] - 1)
    idx_global <- ord[idx_local]
    blocks[[j]] <- list(indices = idx_global, X = Xs[idx_local], Y = Ys[idx_local])
    start <- start + sizes[j]
  }
  blocks
}

fit_block_poly4 <- function(block){
  Xb <- block$X; Yb <- block$Y
  df <- data.frame(x = Xb, y = Yb)
  fit <- tryCatch(lm(y ~ x + I(x^2) + I(x^3) + I(x^4), data = df), error = function(e) NULL)
  if(is.null(fit)){
    a0 <- mean(Yb)
    fitted <- rep(a0, length(Yb))
    secder <- rep(0, length(Yb))
    rss <- sum((Yb - fitted)^2)
    dfres <- length(Yb) - 1
    return(list(coef = c(a0 = a0, a1 = 0, a2 = 0, a3 = 0, a4 = 0),
                fitted = fitted, secder = secder, rss = rss, df.resid = dfres))
  }
  coefj <- coef(fit)
  a0 <- ifelse(is.na(coefj["(Intercept)"]), 0, coefj["(Intercept)"])
  a1 <- ifelse(is.na(coefj["x"]), 0, coefj["x"])
  a2 <- ifelse(is.na(coefj["I(x^2)"]), 0, coefj["I(x^2)"])
  a3 <- ifelse(is.na(coefj["I(x^3)"]), 0, coefj["I(x^3)"])
  a4 <- ifelse(is.na(coefj["I(x^4)"]), 0, coefj["I(x^4)"])
  fitted <- a0 + a1 * Xb + a2 * Xb^2 + a3 * Xb^3 + a4 * Xb^4
  secder <- 2 * a2 + 6 * a3 * Xb + 12 * a4 * Xb^2
  rss <- sum((Yb - fitted)^2)
  dfres <- df.residual(fit)
  list(coef = c(a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4),
       fitted = fitted, secder = secder, rss = rss, df.resid = dfres)
}

estimate_theta_sigma <- function(X, Y, N){
  n <- length(X)
  blocks <- make_blocks(X, Y, N)
  fits <- lapply(blocks, fit_block_poly4)
  fitted_all <- numeric(n)
  secder_all <- numeric(n)
  total_rss <- 0
  total_df_resid <- 0
  for(j in seq_along(blocks)){
    idx <- blocks[[j]]$indices
    fitted_all[idx] <- fits[[j]]$fitted
    secder_all[idx] <- fits[[j]]$secder
    total_rss <- total_rss + fits[[j]]$rss
    total_df_resid <- total_df_resid + fits[[j]]$df.resid
  }
  theta22_hat <- mean(secder_all^2)
  denom <- n - 5 * N
  sigma2_hat <- if(denom <= 0) NA_real_ else total_rss / denom
  list(theta22 = theta22_hat, sigma2 = sigma2_hat, fitted = fitted_all, secder = secder_all, rss = total_rss)
}

compute_h_AMISE <- function(n, sigma2_hat, theta22_hat, support_length){
  if(is.na(sigma2_hat) || is.na(theta22_hat) || theta22_hat <= 0) return(NA_real_)
  factor <- 35 * sigma2_hat * support_length / theta22_hat
  h <- n^(-1/5) * factor^(1/5)
  as.numeric(h)
}

rss_for_blocks <- function(blocks){
  rss <- 0
  for(j in seq_along(blocks)) rss <- rss + sum((blocks[[j]]$Y - blocks[[j]]$fitted)^2)
  rss
}

compute_Cp_table <- function(X, Y, N_candidates = NULL){
  n <- length(X)
  Nmax <- max(min(floor(n / 20), 5), 1)
  if(is.null(N_candidates)) N_candidates <- seq_len(Nmax)
  blocks_max <- make_blocks(X, Y, Nmax)
  fits_max <- lapply(blocks_max, fit_block_poly4)
  for(j in seq_along(blocks_max)) blocks_max[[j]]$fitted <- fits_max[[j]]$fitted
  rss_max <- rss_for_blocks(blocks_max)
  denom <- n - 5 * Nmax
  if(denom <= 0) stop("Degrees of freedom <= 0 for Nmax; increase n or reduce Nmax.")
  cp_vals <- sapply(N_candidates, function(N){
    blocksN <- make_blocks(X, Y, N)
    fitsN <- lapply(blocksN, fit_block_poly4)
    for(j in seq_along(blocksN)) blocksN[[j]]$fitted <- fitsN[[j]]$fitted
    rssN <- rss_for_blocks(blocksN)
    rssN / (rss_max / denom) - (n - 10 * N)
  })
  data.frame(N = N_candidates, Cp = cp_vals)
}

choose_optimal_N <- function(X, Y, N_candidates = NULL){
  cp_tab <- compute_Cp_table(X, Y, N_candidates)
  N_opt <- cp_tab$N[which.min(cp_tab$Cp)]
  list(N_opt = N_opt, cp_table = cp_tab)
}

generate_data <- function(n, alpha = 2, beta = 5, sigma2 = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  X <- rbeta(n, shape1 = alpha, shape2 = beta)
  mfun <- function(x) sin((x^3 + 0.1)^(-1))
  Y <- mfun(X) + rnorm(n, mean = 0, sd = sqrt(sigma2))
  list(X = X, Y = Y, mtrue = mfun(X))
}

parse_num_list <- function(txt){
  txt <- gsub("\\s+", "", txt)
  if(txt == "") return(NULL)
  parts <- strsplit(txt, ",", fixed = TRUE)[[1]]
  as.numeric(parts)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Diagnostics for h_AMISE (single-dataset)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode", "Mode", choices = c("Single run" = "single", "Grid-run" = "grid"), selected = "single"),
      hr(),
      conditionalPanel(
        "input.mode == 'single'",
        sliderInput("n_single", "Sample size n", min = 50, max = 5000, value = 500, step = 10),
        sliderInput("N_single", "Number of blocks N (single run primary)", min = 1, max = 200, value = 5, step = 1),
        sliderInput("N_range_max", "Evaluate h for N = 1.. (select upper bound)", min = 2, max = 100, value = 20, step = 1),
        sliderInput("alpha_single", "Beta alpha", min = 0.1, max = 10, value = 2, step = 0.1),
        sliderInput("beta_single", "Beta beta", min = 0.1, max = 10, value = 5, step = 0.1),
        sliderInput("sigma2_single", "Noise variance sigma^2", min = 0, max = 5, value = 1, step = 0.1),
        checkboxInput("autoN_single", "Auto-select N via Mallow's Cp (for single dataset)", value = FALSE),
        checkboxInput("flag_low_counts", "Flag blocks with low counts (<10)", value = TRUE),
        selectizeInput("Ns_for_m2", "Select up to 4 N values to compare m''(x) across x", choices = c(2,5,10,20), selected = c(2,5), multiple = TRUE),
        actionButton("run_single", "Run single dataset")
      ),
      conditionalPanel(
        "input.mode == 'grid'",
        textInput("n_grid", "n values (comma-separated)", value = "200,500,1000"),
        textInput("N_grid", "N values (comma-separated)", value = "2,5,10,20"),
        textInput("alpha_grid", "alpha values (comma-separated)", value = "0.5,1,2,5,10"),
        textInput("beta_grid", "beta values (comma-separated)", value = "0.5,1,2,5,10"),
        numericInput("sigma2_grid", "Noise variance sigma^2", value = 1, step = 0.1),
        actionButton("run_grid", "Run grid (single dataset per combo)"),
        hr(),
        h4("Cp vs n"),
        textInput("nseq_cp", "n values for Cp vs n (comma-separated)", value = "100,200,400,800,1600"),
        actionButton("compute_cp_vs_n", "Compute Cp-selected N for n sequence")
      ),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Single run",
                 verbatimTextOutput("single_info"),
                 fluidRow(column(6, plotOutput("plotXdens")), column(6, plotOutput("h_vs_N_plot"))),
                 hr(),
                 fluidRow(column(12, plotOutput("theta_sigma_vs_N"))),
                 hr(),
                 fluidRow(column(6, plotOutput("m2_blockwise_plot")), column(6, plotOutput("theta_sigma_by_block"))),
                 hr(),
                 plotOutput("localFitPlot")
        ),
        tabPanel("Grid results",
                 verbatimTextOutput("grid_info"),
                 fluidRow(column(6, plotOutput("grid_h_vs_N")), column(6, plotOutput("grid_h_vs_n"))),
                 hr(),
                 plotOutput("grid_heatmap"),
                 hr(),
                 dataTableOutput("grid_table")
        ),
        tabPanel("Cp vs n",
                 verbatimTextOutput("cpn_info"),
                 plotOutput("cp_vs_n_plot"),
                 plotOutput("h_vs_n_for_fixed_N_plot")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session){
  
  # ---------------- Single run ----------------
  single_res <- eventReactive(input$run_single, {
    n <- input$n_single; alpha <- input$alpha_single; beta <- input$beta_single; sigma2 <- input$sigma2_single
    dat <- generate_data(n, alpha, beta, sigma2, seed = as.integer(Sys.time()) %% 1e6)
    
    Nmax_eval <- max(1, as.integer(input$N_range_max))
    Ns_eval <- seq_len(Nmax_eval)
    hvals <- numeric(length(Ns_eval))
    theta_vals <- numeric(length(Ns_eval))
    sigma_vals <- numeric(length(Ns_eval))
    block_counts_list <- vector("list", length(Ns_eval))
    fits_list <- vector("list", length(Ns_eval))
    
    for(i in seq_along(Ns_eval)){
      Ni <- Ns_eval[i]
      est <- estimate_theta_sigma(dat$X, dat$Y, Ni)
      supp_len <- diff(range(dat$X))
      hvals[i] <- compute_h_AMISE(n = n, sigma2_hat = est$sigma2, theta22_hat = est$theta22, support_length = supp_len)
      theta_vals[i] <- est$theta22
      sigma_vals[i] <- est$sigma2
      blocks <- make_blocks(dat$X, dat$Y, Ni)
      fits <- lapply(blocks, fit_block_poly4)
      counts <- sapply(blocks, function(b) length(b$X))
      block_counts_list[[i]] <- counts
      fits_list[[i]] <- list(blocks = blocks, fits = fits)
    }
    
    cp_tab <- NULL; chosenN <- NULL
    if(isTRUE(input$autoN_single)){
      cp_info <- tryCatch(choose_optimal_N(dat$X, dat$Y), error = function(e) NULL)
      if(!is.null(cp_info)){
        chosenN <- cp_info$N_opt
        cp_tab <- cp_info$cp_table
      }
    }
    
    list(dat = dat, Ns_eval = Ns_eval, hvals = hvals, theta_vals = theta_vals, sigma_vals = sigma_vals,
         block_counts = block_counts_list, fits_list = fits_list, cp_tab = cp_tab, chosenN = chosenN)
  }, ignoreNULL = FALSE)
  
  output$single_info <- renderPrint({
    res <- single_res()
    if(is.null(res)) return("Press 'Run single dataset'")
    cat("Single dataset summary:\n")
    cat(sprintf("n = %d, alpha = %g, beta = %g\n", input$n_single, input$alpha_single, input$beta_single))
    if(!is.null(res$chosenN)) cat(sprintf("Cp selected N = %d\n", res$chosenN))
    cat("h range (min, median, max): ", paste(signif(range(res$hvals, na.rm=TRUE),3), collapse = ", "), "\n")
  })
  
  output$plotXdens <- renderPlot({
    res <- single_res()
    if(is.null(res)) return(NULL)
    ggplot(data.frame(x = res$dat$X), aes(x)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "grey80", color = "black") +
      geom_density() + labs(title = "Density of X", x = "X", y = "density")
  })
  
  output$h_vs_N_plot <- renderPlot({
    res <- single_res()
    if(is.null(res)) return(NULL)
    df <- data.frame(N = res$Ns_eval, h = res$hvals, theta = res$theta_vals, sigma2 = res$sigma_vals)
    p <- ggplot(df, aes(x = N, y = h)) + geom_point() + geom_line() +
      labs(title = "Estimated ĥ vs N", x = "N (blocks)", y = "ĥ")
    if(!is.null(res$chosenN)) p <- p + geom_vline(xintercept = res$chosenN, linetype = "dashed", color = "red")
    p + theme_minimal()
  })
  
  output$theta_sigma_vs_N <- renderPlot({
    res <- single_res()
    if(is.null(res)) return(NULL)
    df <- data.frame(N = res$Ns_eval, theta = res$theta_vals, sigma2 = res$sigma_vals) %>%
      pivot_longer(cols = c("theta","sigma2"), names_to = "stat", values_to = "value")
    ggplot(df, aes(x = N, y = value, color = stat)) +
      geom_point() + geom_line() +
      scale_y_log10() +
      labs(title = expression(paste("log-scale: ", hat(theta)[22], " and ", hat(sigma)^2, " vs N")),
           x = "N", y = "value (log10)") +
      theme_minimal()
  })
  
  # ---------------- m'' blockwise plot ----------------
  output$m2_blockwise_plot <- renderPlot({
    res <- single_res()
    if (is.null(res)) return(NULL)
    
    Ns_chosen <- as.integer(input$Ns_for_m2)
    if (length(Ns_chosen) == 0) {
      plot.new()
      text(0.5, 0.5, "Select some N values to compare m''(x) blockwise")
      return()
    }
    
    df_list <- list()
    for (Ni in Ns_chosen) {
      idx <- which(res$Ns_eval == Ni)
      if (length(idx) == 0) next
      info <- res$fits_list[[idx]]
      blocks <- info$blocks
      fits <- info$fits
      
      for (bi in seq_along(blocks)) {
        Xb <- blocks[[bi]]$X
        sec <- fits[[bi]]$secder
        if (length(Xb) == 0) next
        df_list[[length(df_list) + 1]] <- data.frame(
          x = Xb,
          secder = sec,
          N = paste0("N=", Ni),
          block = bi,
          stringsAsFactors = FALSE
        )
      }
    }
    
    dfall <- dplyr::bind_rows(df_list)
    if (nrow(dfall) == 0) {
      plot.new()
      text(0.5, 0.5, "No block data")
      return()
    }
    
    qhi <- quantile(dfall$secder, probs = 0.995, na.rm = TRUE)
    qlo <- quantile(dfall$secder, probs = 0.005, na.rm = TRUE)
    dfall$secder_trim <- pmin(pmax(dfall$secder, qlo), qhi)
    
    ggplot(dfall, aes(x = x, y = secder_trim, color = N)) +
      geom_point(alpha = 0.5, size = 1) +
      labs(
        title = expression(paste("Blockwise estimates of ", m^{''}, "(x)")),
        x = "x",
        y = expression(m^{''}(x)),
        color = "N"
      ) +
      theme_minimal()
  })
  
  # ---------------- Local fits ----------------
  output$localFitPlot <- renderPlot({
    res <- single_res()
    if (is.null(res)) return(NULL)
    
    Nuse <- input$N_single
    idx <- which(res$Ns_eval == Nuse)
    if (length(idx) == 0) return(NULL)
    
    info <- res$fits_list[[idx]]
    blocks <- info$blocks
    fits <- info$fits
    
    fitcurves <- do.call(rbind, lapply(seq_along(blocks), function(j) {
      Xgrid <- seq(min(blocks[[j]]$X), max(blocks[[j]]$X), length.out = 80)
      coefj <- fits[[j]]$coef
      ygrid <- coefj["a0"] + coefj["a1"] * Xgrid +
        coefj["a2"] * Xgrid^2 + coefj["a3"] * Xgrid^3 + coefj["a4"] * Xgrid^4
      data.frame(X = Xgrid, Y = ygrid, block = factor(j))
    }))
    
    df_points <- data.frame(X = res$dat$X, Y = res$dat$Y)
    
    ggplot() +
      geom_point(data = df_points, aes(x = X, y = Y), alpha = 0.6) +
      geom_line(data = fitcurves, aes(x = X, y = Y, color = block), size = 1) +
      labs(title = paste("Block polynomial fits (N =", Nuse, ")")) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  # ---------------- Grid (single run per combo) ----------------
  grid_res <- eventReactive(input$run_grid, {
    n_vals <- parse_num_list(input$n_grid); N_vals <- parse_num_list(input$N_grid)
    a_vals <- parse_num_list(input$alpha_grid); b_vals <- parse_num_list(input$beta_grid)
    sigma2 <- input$sigma2_grid
    if(is.null(n_vals) || is.null(N_vals) || is.null(a_vals) || is.null(b_vals)){
      showModal(modalDialog("Please supply comma-separated numeric lists for n,N,alpha,beta", easyClose = TRUE)); return(NULL)
    }
    combos <- expand.grid(n = n_vals, N = N_vals, alpha = a_vals, beta = b_vals, stringsAsFactors = FALSE)
    out <- list()
    for(i in seq_len(nrow(combos))){
      row <- combos[i, ]
      dat <- generate_data(row$n, row$alpha, row$beta, sigma2, seed = as.integer(Sys.time()) + i)
      est <- estimate_theta_sigma(dat$X, dat$Y, row$N)
      supp_len <- diff(range(dat$X))
      hhat <- compute_h_AMISE(n = row$n, sigma2_hat = est$sigma2, theta22_hat = est$theta22, support_length = supp_len)
      cp_info <- tryCatch({ choose_optimal_N(dat$X, dat$Y) }, error = function(e) NULL)
      out[[i]] <- data.frame(n = row$n, N = row$N, alpha = row$alpha, beta = row$beta,
                             hhat = hhat, theta22 = est$theta22, sigma2_hat = est$sigma2,
                             N_cp = if(!is.null(cp_info)) cp_info$N_opt else NA)
    }
    do.call(rbind, out)
  }, ignoreNULL = FALSE)
  
  output$grid_info <- renderPrint({
    df <- grid_res()
    if(is.null(df)) return("Run grid to compute.")
    cat("Grid complete. Example summary:\n")
    print(df %>% group_by(n) %>% summarise(mean_h = mean(hhat, na.rm=TRUE)))
  })
  
  output$grid_h_vs_N <- renderPlot({
    df <- grid_res(); if(is.null(df)) return(NULL)
    ggplot(df, aes(x = factor(N), y = hhat, color = factor(n), group = interaction(n,alpha,beta))) +
      geom_point() + geom_line(aes(group = interaction(alpha,beta,n)), alpha = 0.7) +
      facet_grid(alpha ~ beta, labeller = label_both) + theme_minimal() +
      labs(title = "h_hat vs N (single dataset per combo)", x = "N", y = "h_hat", color = "n")
  })
  
  output$grid_h_vs_n <- renderPlot({
    df <- grid_res(); if(is.null(df)) return(NULL)
    ggplot(df, aes(x = factor(n), y = hhat, color = factor(N), group = interaction(N,alpha,beta))) +
      geom_point() + geom_line(aes(group = interaction(N,alpha,beta)), alpha = 0.7) +
      facet_grid(alpha ~ beta, labeller = label_both) + theme_minimal() +
      labs(title = "h_hat vs n (single dataset per combo)", x = "n", y = "h_hat", color = "N")
  })
  
  output$grid_heatmap <- renderPlot({
    df <- grid_res(); if(is.null(df)) return(NULL)
    heat <- df %>% group_by(alpha, beta) %>% summarise(hhat_mean = mean(hhat, na.rm = TRUE))
    ggplot(heat, aes(x = factor(alpha), y = factor(beta), fill = hhat_mean)) +
      geom_tile() + scale_fill_viridis_c(na.value = "grey50") + theme_minimal() +
      labs(title = "Heatmap: mean h_hat over (alpha,beta)", x = "alpha", y = "beta", fill = "mean h")
  })
  
  output$grid_table <- renderDataTable({ grid_res() }, options = list(pageLength = 10, scrollX = TRUE))
  
  # ---------------- Cp vs n ----------------
  cp_vs_n_res <- eventReactive(input$compute_cp_vs_n, {
    nseq <- parse_num_list(input$nseq_cp)
    if(is.null(nseq)) { showModal(modalDialog("Provide n values")); return(NULL) }
    out <- data.frame(n = integer(), N_cp = integer(), stringsAsFactors = FALSE)
    for(i in seq_along(nseq)){
      n <- as.integer(nseq[i])
      dat <- generate_data(n, input$alpha_single, input$beta_single, input$sigma2_single, seed = as.integer(Sys.time()) + i)
      cp_info <- tryCatch({ choose_optimal_N(dat$X, dat$Y) }, error = function(e) NULL)
      Ncp <- if(!is.null(cp_info)) cp_info$N_opt else NA
      out <- rbind(out, data.frame(n = n, N_cp = Ncp))
    }
    out
  })
  
  output$cpn_info <- renderPrint({
    res <- cp_vs_n_res()
    if(is.null(res)) return("Press 'Compute Cp-selected N for n sequence'")
    print(res)
  })
  
  output$cp_vs_n_plot <- renderPlot({
    res <- cp_vs_n_res(); if(is.null(res)) return(NULL)
    ggplot(res, aes(x = n, y = N_cp)) + geom_point() + geom_line() +
      labs(title = "Cp-selected N vs n (single datasets generated)", x = "n", y = "N_cp") + theme_minimal()
  })
  
  # ---------------- h vs n for fixed N ----------------
  output$h_vs_n_for_fixed_N_plot <- renderPlot({
    nseq <- parse_num_list(input$nseq_cp)
    if(is.null(nseq)) nseq <- c(100,200,400,800)
    Nfixed <- parse_num_list(input$N_grid)
    if(is.null(Nfixed)) Nfixed <- c(5,10)
    out <- data.frame()
    for(N in Nfixed){
      for(n in nseq){
        dat <- generate_data(n, input$alpha_single, input$beta_single, input$sigma2_single, seed = as.integer(Sys.time()) + n + N)
        est <- estimate_theta_sigma(dat$X, dat$Y, N)
        supp_len <- diff(range(dat$X))
        hhat <- compute_h_AMISE(n = n, sigma2_hat = est$sigma2, theta22_hat = est$theta22, support_length = supp_len)
        out <- rbind(out, data.frame(n = n, N = N, hhat = hhat))
      }
    }
    ggplot(out, aes(x = factor(n), y = hhat, group = factor(N), color = factor(N))) +
      geom_point() + geom_line() +
      labs(title = "h_hat vs n for fixed N values", x = "n", y = "h_hat") +
      theme_minimal()
  })
  
}

shinyApp(ui, server)
