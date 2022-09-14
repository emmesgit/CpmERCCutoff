## LOAD DATA:
library(CpmERCCutoff)

#Load observed read counts
data("obs_input")
rownames(obs_input) = obs_input$X

#Load expected ERCC data:
data("exp_input")
exp_input <- exp_input[order(exp_input$ercc_id), ]
rownames(exp_input) = exp_input$ercc_id

# #Calculate log2 fold change ratio
# exp_input$expected_lfc_ratio =
#   log2(exp_input$ercc1_conc) - log2(exp_input$ercc2_conc)

#Load metadata:
data("mta_dta")
#Pair samples that received spike 1 with samples that received spike 2.
#The resulting 2-column data frame is used for the 'pairs' argument.
pairs_input = cbind(
  mta_dta[mta_dta$spike == 2, 'samid'],
  mta_dta[match(mta_dta[mta_dta$spike == 2, 'subid'],
                mta_dta[mta_dta$spike == 1,'subid']), 'samid'])
pairs_input = pairs_input[, c(2, 1)] # Put Mix 1 in first column, Mix 2 in second.


##TESTING:

test_that(
  "Return object is 'empLCPM' class and of length 3.", {
  res <- getLowLcpmCutoff(obs = obs_input,
                         exp = exp_input,
                         pairs = pairs_input,
                         rep=10)
  expect_s3_class(getLowLcpmCutoff(obs_input, exp_input, pairs_input, rep=10),
              "empLCPM")
  expect_length(res,3)
})

test_that("Reproducible Result is as expected.", {
  res <- getLowLcpmCutoff(obs = obs_input,
                         exp = exp_input,
                         pairs = pairs_input,
                         rep = 10,
                         seed = 20220718)
  expect_equal(res$cutoff$threshold_value, 2.987)
})
