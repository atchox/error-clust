library(ggplot2)
library(patchwork)
library(parallel)

coverage_helper <- function(tstats, n, df, alpha) {
  t_crit <- qt(1 - alpha / 2, df = df)
  covered <- abs(tstats) <= t_crit
  p_hat <- mean(covered)

  mc_se <- sqrt(p_hat * (1 - p_hat) / R)
  z <- qnorm(1 - alpha / 2)
  lower <- max(0, p_hat - z * mc_se)
  upper <- min(1, p_hat + z * mc_se)

  list(
    coverage_df = data.frame(
      n = n,
      coverage = p_hat,
      mc_se = mc_se,
      lower = lower,
      upper = upper
    ),
    t_df = data.frame(
      n = n,
      tstat = tstats
    )
  )
}

covplot_helper <- function(cvg_df, alpha) {
  ggplot(cvg_df, aes(x = n, y = coverage)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "red") +
    labs(
      x = "sample size",
      y = "coverage probability"
    )
}

ttplot_helper <- function(t_df, p) {
  n_vals <- sort(unique(t_df$n))

  theory_df <- do.call(
    rbind,
    lapply(n_vals, function(n) {
      x <- seq(-4, 4, length.out = 400)
      data.frame(
        n = n,
        x = x,
        density = dt(x, df = n - p)
      )
    })
  )

  ggplot() +
    geom_density(
      data = t_df,
      aes(x = tstat, fill = "Empirical"),
      alpha = 0.3,
      color = "steelblue"
    ) +
    geom_line(
      data = theory_df,
      aes(x = x, y = density, color = "Theoretical"),
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    scale_fill_manual(values = "steelblue") +
    scale_color_manual(values = "black") +
    guides(
      fill = guide_legend(title = NULL),
      color = guide_legend(title = NULL)
    ) +
    facet_wrap(~n, ncol = 2) +
    labs(x = "t-statistic", y = "density") +
    theme(legend.position = "top")
}

qqplot_helper <- function(t_df, p) {
  qq_list <- list()
  n_vals <- sort(unique(t_df$n))

  for (i in seq_along(n_vals)) {
    n <- n_vals[i]
    q_emp <- sort(t_df$tstat[t_df$n == n])
    R <- length(q_emp)
    p <- ppoints(R)

    q_theory <- qt(p, df = n - p)

    qq_list[[i]] <- data.frame(
      n = n,
      q_theory = q_theory,
      q_emp = q_emp
    )
  }

  qq_df <- do.call(rbind, qq_list)

  ggplot(qq_df, aes(x = q_theory, y = q_emp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.3, size = 0.6) +
    facet_wrap(~n, ncol = 3) +
    labs(
      x = "theoretical t-quantiles",
      y = "empirical t-quantiles",
    )
}

alpha <- 0.05

n_vals <- seq(5, 100, 5)
R <- 10000
beta1 <- 1
beta2 <- 2
sigma <- 5
p <- 2

res_list <- mclapply(
  n_vals,
  function(n) {
    set.seed(42 + n)
    df <- n - p

    X1 <- rnorm(n)
    X2 <- rnorm(n)

    tstats <- sapply(seq_len(R), function(r) {
      eps <- rnorm(n, sd = sigma)
      Y <- beta1 * X1 + beta2 * X2 + eps

      fit <- lm(Y ~ X1 + X2 - 1)
      summ <- summary(fit)

      beta1_hat <- coef(fit)["X1"]
      se_beta1_hat <- summ$coefficients["X1", "Std. Error"]
      (beta1_hat - beta1) / se_beta1_hat
    })

    coverage_helper(tstats, n, df, alpha)
  },
  mc.cores = detectCores() - 1
)

coverage_df <- do.call(rbind, lapply(res_list, \(x) x$coverage_df))
coverage_df$n <- as.factor(coverage_df$n)
t_df <- do.call(rbind, lapply(res_list, \(x) x$t_df))

covp <- covplot_helper(coverage_df, alpha)
ttp <- ttplot_helper(t_df[t_df$n %in% c(10, 20, 40, 60, 80, 100), ], p)
qqp <- qqplot_helper(t_df[t_df$n %in% c(10, 20, 40, 60, 80, 100), ], p)

fig <- ((covp / qqp) | ttp) +
  plot_layout(
    widths = c(1.4, 1),
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

ggsave(
  "../figures/iid.pdf",
  fig,
  width = 10,
  height = 8,
  dpi = 300
)


ng_vals <- c(2, 3, 4, 5, 6, 7, 8)
G <- 15
R <- 10000
beta1 <- 1
beta2 <- 2
lambda <- 5
base_sigma <- 0.1
p <- 2

res_list <- mclapply(
  ng_vals,
  function(ng) {
    set.seed(42 + ng)
    n <- ng * G
    df <- n - p

    X1 <- rnorm(n)
    X2 <- rnorm(n)

    tstats <- sapply(seq_len(R), function(r) {
      eps <- unlist(lapply(seq_len(G), function(x) {
        lambda * rep(rnorm(1), ng)
      }))
      eps <- eps + rnorm(n, sd = base_sigma)

      Y <- beta1 * X1 + beta2 * X2 + eps

      fit <- lm(Y ~ -1 + X1 + X2)
      summ <- summary(fit)

      beta1_hat <- coef(fit)["X1"]
      se_beta1_hat <- summ$coefficients["X1", "Std. Error"]
      (beta1_hat - beta1) / se_beta1_hat
    })

    coverage_helper(tstats, n, df, alpha)
  },
  mc.cores = detectCores() - 1
)

coverage_df <- do.call(rbind, lapply(res_list, \(x) x$coverage_df))
coverage_df$n <- as.factor(coverage_df$n)
t_df <- do.call(rbind, lapply(res_list, \(x) x$t_df))

covp <- covplot_helper(coverage_df, alpha)
ttp <- ttplot_helper(t_df[t_df$n %in% c(30, 45, 60, 75, 90, 105), ], p)
qqp <- qqplot_helper(t_df[t_df$n %in% c(30, 45, 60, 75, 90, 105), ], p)

fig <- ((covp / qqp) | ttp) +
  plot_layout(
    widths = c(1.4, 1),
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

ggsave(
  "../figures/cluster.pdf",
  fig,
  width = 10,
  height = 8,
  dpi = 300
)
