#' A data frame of observed spike in ERCC normalized LCPM data
#'
#' A data frame of observed gene expression results for ERCC RNA that was spiked
#' into samples.  Data are read coverage-normalized log2 counts per million (LCPM).
#'
#' @format A data frame with 92 rows and 50 Columns: Each row is an ERCC
#'     transcript, and each column is a sample.  Data are read
#'     coverage-normalized LCPM.
#' \describe{
#'   \item{X}{This first column is ERCC spike in mRNA Ids (ERCC-00002 -- ERCC-00171)}
#'   \item{SAM15}{Sample 15 containing normalized log2 counts per million}
#'   \item{SAM36}{Sample 36 containing normalized log2 counts per million}
#'   \item{SAM19}{Sample 19 containing normalized log2 counts per million}
#'   \item{SAM53}{Sample 53 containing normalized log2 counts per million}
#'   \item{SAM42}{Sample 42 containing normalized log2 counts per million}
#'   \item{SAM32}{Sample 32 containing normalized log2 counts per million}
#'   \item{SAM18}{Sample 18 containing normalized log2 counts per million}
#'   \item{SAM48}{Sample 48 containing normalized log2 counts per million}
#'   \item{SAM26}{Sample 26 containing normalized log2 counts per million}
#'   \item{SAM37}{Sample 37 containing normalized log2 counts per million}
#'   \item{SAM38}{Sample 38 containing normalized log2 counts per million}
#'   \item{SAM29}{Sample 29 containing normalized log2 counts per million}
#'   \item{SAM17}{Sample 17 containing normalized log2 counts per million}
#'   \item{SAM41}{Sample 41 containing normalized log2 counts per million}
#'   \item{SAM09}{Sample 09 containing normalized log2 counts per million}
#'   \item{SAM07}{Sample 07 containing normalized log2 counts per million}
#'   \item{SAM14}{Sample 14 containing normalized log2 counts per million}
#'   \item{SAM02}{Sample 02 containing normalized log2 counts per million}
#'   \item{SAM05}{Sample 05 containing normalized log2 counts per million}
#'   \item{SAM25}{Sample 25 containing normalized log2 counts per million}
#'   \item{SAM08}{Sample 08 containing normalized log2 counts per million}
#'   \item{SAM28}{Sample 28 containing normalized log2 counts per million}
#'   \item{SAM44}{Sample 44 containing normalized log2 counts per million}
#'   \item{SAM04}{Sample 04 containing normalized log2 counts per million}
#'   \item{SAM10}{Sample 10 containing normalized log2 counts per million}
#'   \item{SAM31}{Sample 31 containing normalized log2 counts per million}
#'   \item{SAM21}{Sample 21 containing normalized log2 counts per million}
#'   \item{SAM20}{Sample 20 containing normalized log2 counts per million}
#'   \item{SAM52}{Sample 52 containing normalized log2 counts per million}
#'   \item{SAM46}{Sample 46 containing normalized log2 counts per million}
#'   \item{SAM01}{Sample 01 containing normalized log2 counts per million}
#'   \item{SAM13}{Sample 13 containing normalized log2 counts per million}
#'   \item{SAM39}{Sample 39 containing normalized log2 counts per million}
#'   \item{SAM49}{Sample 49 containing normalized log2 counts per million}
#'   \item{SAM30}{Sample 30 containing normalized log2 counts per million}
#'   \item{SAM50}{Sample 50 containing normalized log2 counts per million}
#'   \item{SAM11}{Sample 11 containing normalized log2 counts per million}
#'   \item{SAM35}{Sample 35 containing normalized log2 counts per million}
#'   \item{SAM06}{Sample 06 containing normalized log2 counts per million}
#'   \item{SAM27}{Sample 27 containing normalized log2 counts per million}
#'   \item{SAM33}{Sample 33 containing normalized log2 counts per million}
#'   \item{SAM22}{Sample 22 containing normalized log2 counts per million}
#'   \item{SAM24}{Sample 24 containing normalized log2 counts per million}
#'   \item{SAM16}{Sample 16 containing normalized log2 counts per million}
#'   \item{SAM34}{Sample 34 containing normalized log2 counts per million}
#'   \item{SAM03}{Sample 03 containing normalized log2 counts per million}
#'   \item{SAM47}{Sample 47 containing normalized log2 counts per million}
#'   \item{SAM40}{Sample 40 containing normalized log2 counts per million}
#'   \item{SAM23}{Sample 23 containing normalized log2 counts per million}
#' }
#'
"obs_input"


#' A data frame of expected ERCC1 and ERCC2 ratios
#'
#' A data frame that contains the expected spike in ERCC data.  This data can be obtained
#' from 'ERCC Controls Analysis' manual located on Thermo Fisher's ERCC RNA
#' Spike-In Mix product
#' [page](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).
#' The 'exp_input' data frame mirrors the fields shown in the ERCC manual.
#' For the LCPM cutoff calculation, the last column of the log2 expected fold
#' change ratios are used.  Ensure that this column is titled
#' "expected_lfc_ratio". See the example code below for formatting the data.
#'
#' @format A data frame with 92 rows and 4 columns: Each row represents an ERCC
#'     transcript.  Columns are described below:
#' \describe{
#'   \item{ercc_id}{ERCC spike-in mRNA Ids (ERCC-00002 -- ERCC-00171)}
#'   \item{subgroup}{ERCC subgroups (A -- D)}
#'   \item{ercc1_conc}{ERCC1 concentration (0.014 -- 30,000)}
#'   \item{ercc2_conc}{ERCC2 concentration (0.007 -- 30,000)}
#'   \item{expected_fc_ratio}{Expected fold change ratio (.5 -- 4)}
#'   \item{expected_lfc_ratio}{Expected log2 fold change ratio (-1 -- 2)}
#'  }
#' @examples
#' # Order rows by ERCC ID and assign to rownames.
#' exp_input = exp_input[order(exp_input$ercc_id), ]
#' rownames(exp_input) = exp_input$ercc_id
#'
"exp_input"

#' A data frame containing sample-level ERCC meta data
#'
#' A data frame containing sample-level ERCC meta data.  In this experiment,
#' subjects had repeated measures and the baseline samples were spiked with
#' ERCC1 and four post-vaccination samples were each spiked with ERCC2.   This
#' meta data is used to create the 2-column data frame of sample pairs where the
#' first column will contain sample IDs that received ERCC1, and the second
#' column will contain sample IDs that received ERCC2.  The pairs data frame is
#' used as input to the `pairs` argument in the `getLowCpmCutOff` function.  Use
#' `?getLowCpmcutOff` to review example code regarding the formatting of the
#' `pairs` data frame.
#'
#' @format A data frame with 49 rows and 4 columns: Each row is a sample, and each
#' column contains metadata such as subject IDs, spike in control type, and
#' treatment groups.  For this study, data was collected at various time points,
#' however, under different experiment conditions, the 'day' column can be
#' represented as treatment group.
#' \describe{
#'   \item{samid}{Sample IDs (SAM01 -- SAM49)}
#'   \item{subid}{Subject/Participant ID (SUB01 -- SUB10)}
#'   \item{day}{collection day (0 -- 14)}
#'   \item{spike}{ERCC spike in control Mix (1 or 2)}
#'  }
#'
"mta_dta"
