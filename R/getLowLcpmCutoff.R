#' Function to empirically determine a log2 CPM cutoff based on ERCC RNA spike-in
#'
#' @description This function uses spike-in ERCC data, known control RNA probes,
#'    and paired samples to fit a 3rd order polynomial to determine an expression
#'    cutoff that meets the specified correlation between expected and observed
#'    fold changes.  The \code{obs} data frame used as input for the observed
#'    expression of the 92 ERCC RNA spike-ins and stores the coverage-normalized
#'    read log2 counts per million (LCPM) that mapped to the respective ERCC
#'    sequences.  Typically, prior to LCPM calculation, the read count data is
#'    normalized for any systematic differences in read coverage between samples,
#'    for example, by using the TMM normalization method as implemented in
#'    the \code{edgeR} package.
#'
#'    For each bootstrap replicate, the paired samples are subsampled with
#'    replacement.  The mean LCPM of each ERCC transcript is determined by
#'    first calculating the average LCPM value for each paired sample, and
#'    then taking the mean of those averages. The ERCC transcripts are sorted
#'    based on these means, and are then grouped into \code{n.bins} ERCC bins.
#'    Next, the Spearman correlation metric is used to calculate the association
#'    between the empirical and expected log fold change (LFC) of the ERCCs in
#'    each bin for each sample.
#'    Additionally, the average LCPM for the ERCCs in each bin are calculated
#'    for each sample. This leads to a pair of values - the average LCPM and the
#'    association value - for each sample and each ERCC bin.  Outliers within
#'    each ERCC bin are identified and removed based on >1.5 IQR.
#'    A 3rd order polynomial is fit with the explanatory variable being the
#'    average LCPM and the response variable being the Spearman correlation
#'    value between expected and observed log2 fold changes.
#'    The fitted curve is used to identify the average LCPM value with a Spearman
#'    correlation of \code{cor.value}. The results are output as an "empLCM"
#'    object as described below.  The \code{\link{summary.empLCPM}} function can
#'    be used to extract a summary of the results, and the
#'    \code{\link{plot.empLCPM}} function to plot the results for visualization.
#'
#' @param obs A data frame  of observed spike-in ERCC data.  Each row is an ERCC
#'     transcript, and each column is a sample.  Data are read
#'     coverage-normalized log2 counts per million (LCPM).
#' @param exp A data frame of expected ERCC Mix 1 and Mix 2 ratios with a column
#'    titled `expected_lfc_ratio` containing the expected log2 fold-change
#'     ratios. This data can be obtained from 'ERCC Controls Analysis' manual
#'     located on Thermo Fisher's ERCC RNA Spike-In Mix product
#'     [page](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).
#'     The 'exp_input' data frame mirrors the fields shown in the ERCC manual.
#'     For the LCPM cutoff calculation, the last column containing the log2
#'     expected fold change ratios are used.  Ensure that this column is titled
#'     "expected_lfc_ratio". See the example code below for formatting the data.     #
#' @param pairs A 2-column data frame where each row indicates a sample pair
#'     with the first column indicating the sample that received ERCC spike-ins
#'     from Mix 1 and the second column indicating the sample receiving Mix 2.
#' @param n.bins Integer.  The number of abundance bins to create.  Default is 7.
#' @param rep Integer.  The number of bootstrap replicates.  Default is 1000.
#' @param ci Numeric.  The confidence interval.  Default is 0.95.
#' @param cor.value Numeric.  The desired Spearman correlation between the
#'     empirical log2 fold change across the ERCC transcripts.  Default is 0.9.
#' @param remove.outliers If TRUE (default) outliers are identified as exceeding
#'     1.5 IQR, and are removed prior to fitting the polynomial. Set to FALSE
#'     to keep all points.
#' @param seed Integer.  The reproducibility seed.  Default is 20220719.
#'
#' @return An "empLCPM" object is returned, which contains the following named
#' elements:
#' \tabular{ll}{
#'    \code{cutoff} \tab a vector containing 3 values: the threshold value,
#'    upper confidence interval, \cr \tab and the lower confidence interval value. \cr
#'    \code{args} \tab a key: value list of arguments that were provided. \cr
#'    \code{res} \tab a list containing the main results and other
#'    information from the input. \cr \tab The \code{\link{summary.empLCPM}}
#'    function should be used to extract a summary table. \cr
#'  }
#' @seealso \code{\link{summary.empLCPM}}, \code{\link{plot.empLCPM}},
#' \code{\link{print.empLCPM}}
#' @export
#' @import stats
#' @examples
#' library(CpmERCCutoff)
#' ##############################
#' # Load and wrangle input data:
#' ##############################
#' # Load observed read counts
#' data("obs_input")
#'
#' # Set ERCC Ids to rownames
#' rownames(obs_input) = obs_input$X
#'
#' # Load expected ERCC data:
#' data("exp_input")
#'
#' # Order rows by ERCC ID.
#' exp_input = exp_input[order(exp_input$ercc_id), ]
#' rownames(exp_input) = exp_input$ercc_id
#'
#' # Load metadata:
#' data("mta_dta")
#'
#' # Pair samples that received ERCC Mix 1 with samples that received ERCC Mix 2.
#' # The resulting 2-column data frame is used for the 'pairs' argument.
#' # Note: the code here will depend on the details of the given experiment. In
#' #       this example, the post-vaccination samples (which received Mix 2)
#' #       for each subject are paired to their pre-vaccination samples (which
#' #       received Mix 1).
#' pairs_input = cbind(
#'   mta_dta[mta_dta$spike == 2, 'samid'],
#'   mta_dta[match(mta_dta[mta_dta$spike == 2, 'subid'],
#'                 mta_dta[mta_dta$spike == 1,'subid']), 'samid'])
#' # Put Mix 1 in the first column and Mix 2 in the second.
#' pairs_input = pairs_input[, c(2, 1)]
#'
#' ###############################
#' # Run getLowLcpmCutoff Function:
#' ###############################'
#' # Note: Here we use `rep = 10` for only 10 bootstrap replicates
#' #       to decrease the run time for this example; a lager number
#' #       should be used in practice (default = 1000).
#' res = getLowLcpmCutoff(obs = obs_input,
#'                        exp = exp_input,
#'                        pairs = pairs_input,
#'                        n.bins = 7,
#'                        rep = 10,
#'                        cor.value = 0.9,
#'			                  remove.outliers = TRUE,
#'                        seed = 20220719)
#'
#' # Print a short summary of the results:
#' res
#'
#' # Extract a summary table of the results:
#' summary(res)
#'
#' # Create a plot of the results:
#' plot(x = res,
#'      main = "Determination of Empirical Minimum Expression Cutoffs using ERCCs",
#'      col.trend = "blue",
#'      col.outlier = c("black", "red"))
#'
getLowLcpmCutoff = function(obs,
                            exp,
                            pairs,
                            n.bins    = 7,
                            rep       = 1000,
                            ci        = 0.95,
                            cor.value = 0.9,
			                      remove.outliers = TRUE,
                            seed      = 20220719) {

  if(!("expected_lfc_ratio" %in% colnames(exp))) {
    stop("`exp` must contain a column named 'expected_lfc_ratio' containing ",
         "the expected log2 fold change between ERCC Mix 1 and Mix 2.")
  }
  if(!all(rownames(obs) == rownames(exp))) {
    if(length(setdiff(rownames(obs), rownames(exp)) != 0)) {
      stop("Row names of `obs` and `exp` do not match.")
    }
    warning("Row names of `obs` and `exp` are not aligned. Aligning...")
    obs = obs[match(rownames(obs), rownames(exp)), ]
  }

  # Variable to hold bootstrap replicates.
  sampled = NULL

  # For each bootstrap replicate plus one for the actual:
  for (z in c(1:(rep + 1))) {

    if (z == 1) {
      # Use the original sample.
      pairs.sampled = pairs
    } else {
      # Set the seed for reproducibility of each bootstrap sample.
      if(!is.null(seed)) set.seed(seed + z - 1)
      # Subsample observed pairs with replacement for bootstrap replicate.
      pairs.sampled = pairs[sample(1:nrow(pairs), nrow(pairs), replace = T), ]
    }

    # Calculate average ERCC expression across pairs.
    avg.lcpm.mtx = (obs[, pairs.sampled[, 1]] + obs[, pairs.sampled[, 2]]) / 2
    # Create observed differential expression data set across pairs.
    # The first column in pairs is assumed to contain the ERCC Mix 1, while the
    # second column contains Mix 2.
    # Calculate the conc1 - conc2 log fold change for each paired sample:
    obs.mtx = obs[, pairs.sampled[, 1]] - obs[, pairs.sampled[, 2]]
    # Created expected differential expression data set across pairs.
    exp.mtx = as.data.frame(matrix(rep(exp$expected_lfc_ratio, ncol(obs.mtx)),
                                   ncol = ncol(obs.mtx)))

    ######################################################
    #
    # GROUP 92 ERCCS INTO SEPARATE BINS ORDERED BY AVERAGE
    # OF THE AVERAGE PRE AND POSTVAC LCPMS
    #
    ######################################################

    # Order ERCC spike-ins based on average of the average sample values.
    ord.idx = order(rowMeans(avg.lcpm.mtx), decreasing = F)

    # Create ERCC bins (split the 92 ERCCs in equal parts).
    bin.list = split(ord.idx, cut(1:92, n.bins))
    names(bin.list) = lapply(bin.list, function(x) {
      round(mean(unlist(avg.lcpm.mtx[x, ])), 2)
    })

    # Calculate association metrics for each pair over the subset of ERCCs
    # in each bin.
    cor.list = lapply(bin.list, function(x) {
        sapply(1:ncol(obs.mtx), function(y) {
            tryCatch(cor(obs.mtx[x, y], exp.mtx[x, y], method = "spearman"),
                     warning = function(w) 0)
        })
    })

    # Calculate average abundance for each sample based on ERCCs per bin.
    avg.abundance = lapply(bin.list, function(x){
      sapply(1:ncol(obs.mtx), function(y){ mean(avg.lcpm.mtx[x, y]) })
    })

    # Identify any outliers in each bin.
    if(!remove.outliers) {
      keep.bool = lapply(1:length(bin.list), function(b) rep(T, length(cor.list[[b]])))
    } else {
      keep.bool = lapply(1:length(bin.list), function(b) {
        # Note: if outlier criterion is too stringent, the bootstrap CI widens.
        Q = quantile(cor.list[[b]], c(0.25, 0.5, 0.75))
        iqr = Q[3] - Q[1]
        return(cor.list[[b]] <= Q[3] + 1.5 * iqr &
                 cor.list[[b]] >= Q[1] - 1.5 * iqr)
      })
    }
    y.clean = unlist(lapply(1:length(cor.list),
                            function(x){ cor.list[[x]][keep.bool[[x]]] }))
    x.clean = unlist(lapply(1:length(avg.abundance),
                            function(x){ avg.abundance[[x]][keep.bool[[x]]] }))

    # Fit 3rd order polynomial.
    lm.res = lm(y.clean ~ poly(x.clean, 3, raw = TRUE))

    # Get the coefficients from the fitted model.
    coef = round(lm.res$coefficients, 6)
    pol.f = function(x) coef[1] + coef[2] * x + coef[3] * x^2 + coef[4] * x^3

    # Find the threshold that provides correlation value >= 0.9.
    x.thrsh = tryCatch({
      uniroot(function(x){ pol.f(x) - cor.value },
              c(min(x.clean), max(x.clean)))$root
    }, error = function(e) {
      NA
    })

    # Save results from the original dataset.
    if (z == 1) {
      actual = x.thrsh
      avg.abundance.res = avg.abundance
      cor.list.res =  cor.list
      x.prd = seq(min(x.clean), max(x.clean), by = 0.01)
      y.prd = predict(lm.res, newdata = data.frame(x.clean = x.prd))
      # Save (x, y) points used to fit the cubic.
      x.plot = unlist(lapply(1:length(avg.abundance), function(x){ avg.abundance[[x]]}))
      y.plot = unlist(lapply(1:length(cor.list), function(x){ cor.list[[x]]}))
      # Save the outlier indicator for the points.
      outlier =  unlist(lapply(1:length(keep.bool), function(x){ !keep.bool[[x]]}))
    # Save results from the bootstrap replicates.
    } else {
      sampled = c(sampled, x.thrsh)
    }
  }

  # Calculate the quantiles of the bootstrap sampled LCPM thresholds.
  conf.res = quantile(sampled, (1 - ci) / 2 * c(1, -1) + c(0, 1), na.rm = TRUE)
  n.NA = sum(is.na(sampled)) # Number of bootstrap samples producing NA threshold.

  # Create results list and return.
  res = list(
    cutoff = list(threshold_value = signif(actual,4),
                  lower_ci = signif(conf.res[1],4),
                  upper_ci = signif(conf.res[2],4)),
    args   = list(n.bins    = n.bins,
                  rep       = rep,
                  ci        = ci,
                  cor.value = cor.value,
                  remove.outliers = remove.outliers,
                  seed      = seed),
    res    = list(boot.res  = sampled,
                  n.NA      = n.NA,
                  n.pairs   = nrow(pairs),
                  lm.res    = lm.res,
                  x.plot    = x.plot,
                  y.plot    = y.plot,
                  outlier   = outlier,
                  cor.list  = cor.list.res,
                  avg.abundance = avg.abundance.res))

  class(res) = "empLCPM"

  return(res)
}

#' Summarizing the empirically derived LCPM cutoff
#'
#' Summary method for class "empLCPM".
#'
#' @param object A 'empLCPM' object from \code{\link{getLowLcpmCutoff}}.
#' @param ... Additional arguments are ignored.
#' @return The function \code{summary.empLCPM} extracts the key results and
#' returns a string of the summary.
#' @export
summary.empLCPM = function(object, ...) {
  s = data.frame(n.pairs = object$res$n.pairs,
                 n.bins = object$args$n.bins,
                 cor.val = object$args$cor.value,
                 n.boot = object$args$rep,
                 conf.level = object$args$ci,
                 lcpm.cutoff = object$cutoff$threshold_value,
                 lcpm.lower = object$cutoff$lower_ci,
                 lcpm.upper = object$cutoff$upper_ci)
  rownames(s) = NULL
  return(s)
}

#' Print the empirically derived LCPM cutoff results
#'
#' Print method for class "empLCPM".
#'
#' @param x A 'empLCPM' object from \code{\link{getLowLcpmCutoff}}.
#' @param ... Additional arguments are ignored.
#' @return The function \code{print.empLCPM} extracts the key results and
#' returns a string of the summary.
#' @export
print.empLCPM = function(x, ...) {
  s = summary(x)
  cat(paste0("LCPM Cutoff Value: ", s$lcpm.cutoff,
             " (", round(s$conf.level * 100, 1), "% bootstrap CI: ",
             x$cutoff$lower_ci, ', ',x$cutoff$upper_ci, ")\n"))
}

#' Plot the empirically derived LCPM cutoff
#'
#' Plot method for class "empLCPM" that plots the 3rd order polynomial fit for
#' the Spearman correlation between the expected and observed ERCC fold changes
#' across log2 CPM abundance bins and highlights the empirically derived LCPM
#' cutoff.
#'
#' @param x A 'empLCPM' object from \code{\link{getLowLcpmCutoff}}.
#' @param main An (optional) title for the plot.
#' @param col.trend Color used for the 3rd order polynomial fit.
#' @param col.outlier A vector specifying the default color for the points
#' in the scatterplot, and the color for the outlier points. The default is
#' to color all non-outlier points black and the outliers red.
#' @param ... Additional arguments are ignored.
#' @return The function \code{plot.empLCPM} creates a plot that summarizes 3rd
#' order polynomial fit for the Spearman correlation between the expected and
#' observed ERCC fold changes across log2 CPM abundance bins.  Vertical lines
#' are used to indicate the optimal cutoff value, along with the lower and upper
#'95\% bootstrap confidence interval highlighted by dashed vertical lines.
#' @export
plot.empLCPM = function(x, main = "",
                        col.trend = "purple",
                        col.outlier = c("black", "red"),
                        ...) {
  if (x$res$n.NA / x$args$rep > 0.05) {
    # This may primarily happen if the fitted cubic function is shifted
    # too low, with no LCPM cutoff resulting in correlation of 0.9.
    warning(x$res$n.NA, " of ", x$args$rep, "bootstrap cutoffs resulted in NA.")
  }

  # Plot all points - color code outliers differently.
  plot(x = unlist(x$res$x.plot),
       y = unlist(x$res$y.plot),
       cex      = 0.7,
       pch      = 20,
       cex.axis = 0.7,
       xlab     = '',
       ylab     = '',
       main     = main,
       col      = col.outlier[x$res$outlier + 1])
  # Add x-axis and y-axis labels.
  mtext('Spearman Correlation',
        side = 2,
        line = 2.6,
        cex  = 1)
  mtext(expression('Log'[2] * ' CPM'),
        side = 1,
        line = 2.6,
        cex  = 1)

  ## Add predicted values to plot.
  x.prd = seq(min(x$res$x.plot), max(x$res$x.plot), by = 0.01)
  y.prd = predict(x$res$lm.res, newdata = data.frame(x.clean = x.prd))
  lines(x.prd,
        y.prd,
        pch = 20,
        col = col.trend,
        lwd = 2)


  ## Add actual min cutoff and 95% CI to plot.
  abline(v   = x$cutoff$threshold_value,
         col = 'darkblue',
         lwd = 2)
  abline(v   = x$cutoff$lower_ci,
         col = 'lightblue',
         lty = 2,
         lwd = 2)
  abline(v   = x$cutoff$upper_ci,
         col = 'lightblue',
         lty = 2,
         lwd = 2)

  ## Add cutoffs.
  legend(par('usr')[2] - ((par('usr')[2] - par('usr')[1]) * .4),
         par('usr')[3] + ((par('usr')[4] - par('usr')[3]) * 0.15),
         title  = expression('Miniumum Log'[2] * ' CPM Cut Off [95% CI]'),
         legend = paste(round(x$cutoff$threshold_value, 2),
                        ' [', round(x$cutoff$lower_ci, 2), ', ',
                        round(x$cutoff$upper_ci, 2), ']'),
         bty    = 'n',
         cex    = 0.75)
  legend('topleft',
         legend = paste0('Bin Size =', x$args$n.bins),
         bty    = 'n',
         cex    = 0.75)

  legend(par('usr')[2] - ((par('usr')[2] - par('usr')[1]) * .3),
         par('usr')[3] + ((par('usr')[4] - par('usr')[3]) * 0.45),
         legend = c('Cubic polynomial Fit',
                    'Cut Off Point',
                    'Outliers Excluded From\nCurve Fitting'),
         pch = c(NA, NA, 16),
         lty = c(1, 2, NA),
         lwd = c(2, 2, NA),
         col = c(col.trend, 'darkblue', col.outlier[2]),
         cex = 0.7,
         bty = 'n')
}
