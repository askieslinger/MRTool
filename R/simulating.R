#' Simulate U,X and Y and run MR on them
#'
#' Simulates confounder U, exposure X and outcome Y based on genotype data and betacoefficients.
#' For examples see vignette browseVignettes('MRTool')
#'
#' @param SNP_data
#' @param SNP
#' @param my_Parameters
#' @param reverse=FALSE TRUE if MR should be conducted with X as putative outcome
#' and Y as putative exposure
#'
#' @return list(reg_res,sim_data,mr_res)
#' @export
#'
#' @import data.table
sim_and_mr <- function(SNP_data, SNP,  my_Parameters,reverse=FALSE) {
  # create data for example scenario
  my_MR_data <- SimulateMRData(
    SNP = SNP,
    SNP_data =SNP_data,
    Parameters = my_Parameters
  )
  setkey(my_MR_data$results,Scenario,term)
  if(reverse){
    my_MR_data$results[beta=='G_X',beta:='gxtmp']
    my_MR_data$results[beta=='G_Y',beta:='G_X']
    my_MR_data$results[beta=='gxtmp',beta:='G_Y']
  }
  # create the input object
  # create DT with only results from X and Y
  res_xy <- rbind(my_MR_data$results[.(unique(Scenario),c('U'))],my_MR_data$results[.(unique(Scenario),c('X'))])
  res_snps <- rbindlist(lapply(paste0('`',SNP,'`'), FUN=function(snp) my_MR_data$results[.(unique(Scenario),snp)]))
  setkey(res_snps,term)
  res<-lapply(
    res_snps[, unique(term)],
    FUN = function(snp) {
      mr_data_snp <- rbind(res_xy, res_snps[.(snp)])
      mr_objects <- list()
      for (scen in my_Parameters$Scenario) {
        my_MR_object <-
          CreateMRInputObject(MR_data = list(results=mr_data_snp[Scenario == scen],data=my_MR_data$data),
                              MR_Scenario = scen,correlation=F)
        # use the MR package
        MR_main_object <- MendelianRandomization::mr_allmethods(my_MR_object, method = "ivw")
        setnames(MR_main_object$Values, 4, '95% CI lower')
        setnames(MR_main_object$Values, 5, '95% CI upper')
        mr_objects[[scen]] = MR_main_object
      }
      mr_objects
    }
  )
  names(res) <-  res_snps[, unique(term)]
  if(reverse){
    my_MR_data$results[beta=='G_X',beta:='gxtmp']
    my_MR_data$results[beta=='G_Y',beta:='G_X']
    my_MR_data$results[beta=='gxtmp',beta:='G_Y']
  }
  list(
    reg_res = my_MR_data$results,
    sim_data = my_MR_data$data,
    mr_res = res
  )

}


#' Simulate MR Data
#'
#' For examples see vignette browseVignettes('MRTool')
#' @param SNP_data
#' @param SNP
#' @param Parameters
#'
#' @return
#' @export
#'
#' @importFrom foreach %do%
#' @import data.table
SimulateMRData <- function(SNP_data = SNP_data,
                           # The genetic data
                           SNP = NULL,
                           Parameters = NULL # Parameters generated with the SetMRParams function
)
{
  # set the weight to the error in the data generation

  if (is.null(Parameters) ||
      is.na(Parameters)) {
    stop(
      "Please create an parameter object via the SetMRParams() function and supply it to the function as 'Parameters =' or select >2 SNP for the analysis."
    )
  }

  All_Scenarios <-
    foreach::foreach(myScenario = unique(Parameters$Scenario)) %do% {

      params <-
        Parameters[Scenario == myScenario, ]
      simsd_X <- params[, unique(simsd_X)]
      simsd_Y <- params[, unique(simsd_Y)]
      simsd_U <- params[, unique(simsd_U)]



      # Start with the data of each SNP
      G <- SNP_data[, .SD, .SDcols = SNP]


      # Change order for later multiplication of the proper G_U with each matching SNP
      if (attributes(Parameters)$Vary_between_SNP == T) {
        setcolorder(x = G, neworder = SNP)
      }
      if (length(params$G_X) > 1) {
        setcolorder(x = G, neworder = SNP)
      }
      ### data creation: U ---------------------------------------------------
      #==================#
      # Data creation: U #
      #==================#

      # Generate the data for U, X and Y based on the parameters provided

      if (length(params$G_U) > 1) {
        U <-
          as.data.table(as.matrix(G) %*% params$G_U + simsd_U * rnorm(n = nrow(G)))
      } else{
        U <- rowSums(params$G_U * G) + simsd_U * rnorm(n = nrow(G))
      }
      # }

      ### data creation: x --------------------------------------------------------------------
      ###

      each_G_X <- params$G_X
      if (length(each_G_X) > 1) {
        X <-
          as.data.table(as.matrix(G) %*% each_G_X  + unique(params$U_X) * U + simsd_X * rnorm(n = nrow(G)))
      } else{
        X <-
          rowSums(each_G_X * G) + params$U_X * U + simsd_X * rnorm(n = nrow(G)) # old error mod: (1 - params$G_X - params$U_X)
      }

      ###
      ### data creation: Y -------------------------------------------------------------------------------------
      #==================#
      # Data creation: Y #
      #==================#

      if (length(params$G_Y) > 1) {
        Y <-
          as.data.table(
            unique(params$U_Y) * U + unique(params$X_Y) * X + (as.matrix(G) %*% params$G_Y) + simsd_Y * rnorm(n = nrow(G))
          )
      } else{
        Y <-
          params$U_Y * U + params$X_Y * X + rowSums(params$G_Y * G) + params$simsd_Y * rnorm(n = nrow(G))
      }

      # Create Data.Table of vector U and X
      U_X_Y <- data.table(U = U,
                          X = X,
                          Y = Y)
      setnames(U_X_Y, old = names(U_X_Y), new = c('U', 'X', 'Y'))

      # merge the Data
      mydata <-
        cbind(G, U_X_Y)

      ###===========================###
      # 190201 Correct Implementation #
      ###===========================###

      # Create a readable vector with all my SNP -> lm needs `` to distinguish between : and other parts
      # allSNP <- paste0("`", SNP[SNP %in% params$SNP], "`")
      mySNP <- SNP

      # loop through each SNP and fit the required models and extract the desired statistics
      all_betas <-
        lapply(mySNP, function(my.snp) {
          # debug
          # my.snp <- mySNP[1]

          # necessary for model formula (doesn't like ":" etc)
          my.snp.raw <- my.snp
          my.snp <- paste0("`", my.snp, "`")

          # create formulas X ~ SNP & Y ~ SNP
          G_X_lm_formula <-
            as.formula(paste0("X ~ ", my.snp))
          G_Y_lm_formula <-
            as.formula(paste0("Y ~ ", my.snp))
          G_U_lm_formula <-
            as.formula(paste0("U ~ ", my.snp))

          # fit the lm for
          G_X_lm <-
            summary(lm(formula = G_X_lm_formula, data = mydata))
          G_Y_lm <-
            summary(lm(formula = G_Y_lm_formula, data = mydata))
          G_U_lm <-
            summary(lm(formula = G_U_lm_formula, data = mydata))
          #print(paste('fstat: ',G_X_lm$fstatistic[1]))
          # get r2
          G_X_res <-
            as.data.table(broom::tidy(G_X_lm))
          G_Y_res <-
            as.data.table(broom::tidy(G_Y_lm))
          G_U_res <-
            as.data.table(broom::tidy(G_U_lm))

          # remove intercept
          G_X_res <-
            G_X_res[term == (my.snp), ]
          G_Y_res <-
            G_Y_res[term == (my.snp), ]
          G_U_res <-
            G_U_res[term == (my.snp), ]

          # add r2, n and identifier for which beta
          G_X_res[, `:=`(
            adjRsquared = G_X_lm$adj.r.squared,
            fstat = G_X_lm$fstatistic[1],
            n = length(residuals(G_X_lm)),
            beta = "G_X"
          )]
          G_Y_res[, `:=`(
            adjRsquared = G_Y_lm$adj.r.squared,
            fstat = G_Y_lm$fstatistic[1],
            n = length(residuals(G_Y_lm)),
            beta = "G_Y"
          )]
          G_U_res[, `:=`(
            adjRsquared = G_U_lm$adj.r.squared,
            fstat = G_U_lm$fstatistic[1],
            n = length(residuals(G_U_lm)),
            beta = "G_U"
          )]

          # assign the specific effects to each row
          if (nrow(params) > 1) {
            G_X_res[, assigned.effect := params[SNP == (my.snp.raw), G_X]]
            G_Y_res[, assigned.effect := params[SNP == (my.snp.raw), G_Y]]
            G_U_res[, assigned.effect := params[SNP == (my.snp.raw), G_U]]
          } else if (nrow(params) == 1) {
            G_X_res[, assigned.effect := params[, G_X]]
            G_Y_res[, assigned.effect := params[, G_Y]]
            G_U_res[, assigned.effect := params[, G_U]]
          } else {
            stop("Something is wrong! Please check the parameter file.")
          }

          G_X_res[, `:=`(sim_sd = simsd_X)]
          G_Y_res[, `:=`(sim_sd = simsd_Y)]
          G_U_res[, `:=`(sim_sd = simsd_U)]

          # consolidate and return
          res <-
            rbindlist(list(G_X_res, G_Y_res, G_U_res))
          return(res)
        })

      # consolidate in one data table
      results <- rbindlist(all_betas)

      # Also get via lm and incorporate it into results (but not in the G-loop)
      U_X_lm_formula <- as.formula("X ~ U")
      U_Y_lm_formula <- as.formula("Y ~ U")
      X_Y_lm_formula <- as.formula("Y ~ X  ")
      Y_X_lm_formula <- as.formula("X ~ Y ")

      # fit the lm for
      U_X_lm <-
        summary(lm(formula = U_X_lm_formula, data = mydata))
      U_Y_lm <-
        summary(lm(formula = U_Y_lm_formula, data = mydata))
      X_Y_lm <-
        summary(lm(formula = X_Y_lm_formula, data = mydata))
      Y_X_lm <-
        summary(lm(formula = Y_X_lm_formula, data = mydata))

      # get r2
      U_X_res <-
        as.data.table(broom::tidy(U_X_lm))
      U_Y_res <-
        as.data.table(broom::tidy(U_Y_lm))
      X_Y_res <-
        as.data.table(broom::tidy(X_Y_lm))
      Y_X_res <-
        as.data.table(broom::tidy(Y_X_lm))

      # remove intercept/just keep the effects of U on X/Y and X on Y
      U_X_res <- U_X_res[term == "U", ]
      U_Y_res <- U_Y_res[term == "U", ]
      X_Y_res <- X_Y_res[term == "X", ]
      Y_X_res <- Y_X_res[term == "Y", ]

      # add r2, n and identifier for which beta
      U_X_res[, `:=`(
        adjRsquared = U_X_lm$adj.r.squared,
        fstat = U_X_lm$fstatistic[1],
        n = length(residuals(U_X_lm)),
        beta = "U_X",
        assigned.effect = unique(params$U_X),
        sim_sd = simsd_U
      )]
      U_Y_res[, `:=`(
        adjRsquared = U_Y_lm$adj.r.squared,
        fstat = U_Y_lm$fstatistic[1],
        n = length(residuals(U_Y_lm)),
        beta = "U_Y",
        assigned.effect = unique(params$U_Y),
        sim_sd = simsd_U
      )]
      X_Y_res[, `:=`(
        adjRsquared = X_Y_lm$adj.r.squared,
        fstat = X_Y_lm$fstatistic[1],
        n = length(residuals(X_Y_lm)),
        beta = "X_Y",
        assigned.effect = unique(params$X_Y),
        sim_sd = simsd_U
      )]
      Y_X_res[, `:=`(
        adjRsquared = Y_X_lm$adj.r.squared,
        fstat = Y_X_lm$fstatistic[1],
        n = length(residuals(Y_X_lm)),
        beta = "Y_X",
        assigned.effect = 0,
        sim_sd = simsd_U
      )]

      # consolidate
      results.rest <-
        rbindlist(list(U_X_res, U_Y_res, X_Y_res,Y_X_res))
      # results.rest[ , Rep := (myRep)]

      all.results <-
        rbindlist(list(results, results.rest))
      all.results[, Scenario := myScenario]
      mydata[, Scenario := myScenario]

      # consolidate data and results in one list
      combined.output <-
        list(all.results, mydata)

      # name the list elements according to their scenario
      names(combined.output) <-
        c("results", "data")
      #print(paste('consolidation finished: ',myScenario))


      # return
      return(combined.output)
    }

  ###
  # Loops have ended here, now the output needs to be consolidated
  ###
  #print('left foreach')
  # extract all model results
  res.results <-
    lapply(All_Scenarios, function(my.scenario) {
      # my.scenario <- All_Scenarios[[1]]
      my.scenario$results
    })
  all.results <- rbindlist(res.results)

  # extract all model data
  res.data <-
    lapply(All_Scenarios, function(my.scenario) {
      # my.scenario <- All_Scenarios[[1]]
      my.scenario$data
    })
  all.data <- rbindlist(res.data)
  setkey(all.data, Scenario)
  # Output
  return(list(results = all.results,
              data = all.data))
}


#' Create MRInput Object
#'
#' For examples see vignette browseVignettes('MRTool')
#' @param MR_data
#' @param MR_Scenario
#' @param correlation
#'
#' @return
#' @export
#'
CreateMRInputObject <- function(MR_data = my_MR_data,
                                MR_Scenario = "V1",
                                correlation = T) {
  # Remind the user to input data if missing
  if (any(is.na(MR_data) | is.na(MR_Scenario)) == T) {
    return(print(
      "Please complete or correct data input. Something seems to be missing!"
    ))
    stop()
  }

  # What are the names of the SNP I am using as IVs
  mySNP <- names(MR_data$data[Scenario == (MR_Scenario)])
  mySNP <- mySNP[!mySNP %in% c("U", "X", "Y", "Scenario")]

  # First get the results and the data for the scenario I want to look at
  myData <- MR_data$data[Scenario == (MR_Scenario)]
  myResults <- MR_data$results[Scenario == (MR_Scenario)]

  # Create the Correlation Matrix for the IVs
  # First create the Covariance Matrix
  C <- cov(myData[, .SD, .SDcols  = mySNP], use = "na.or.complete")

  # Then pull the standard deviations from the Covariance Matrix
  S <- diag(C) ^ (-1 / 2)
  if (length(S) > 1) {
    S <- diag(S)
  }
  # Calc the Correlation Matrix
  CM <- S %*% C %*% S


  # Create the input object used by MendelianRandomization package
  if (correlation) {
    MRInputObject <-
      MendelianRandomization::mr_input(
        bx = myResults[Scenario == (MR_Scenario) & beta == "G_X", estimate],
        bxse = myResults[Scenario == (MR_Scenario) &
                           beta == "G_X", std.error],
        by = myResults[Scenario == (MR_Scenario) &
                         beta == "G_Y", estimate],
        byse = myResults[Scenario == (MR_Scenario) &
                           beta == "G_Y", std.error],
        correlation = CM,
        # exposure = "Simulated Exposure",
        # outcome = "Simulated Outcome",
        snps = myResults[Scenario == (MR_Scenario) &
                           beta == "G_X", term]
      )
  } else{
    MRInputObject <-
      MendelianRandomization::mr_input(
        bx = myResults[Scenario == (MR_Scenario) & beta == "G_X", estimate],
        bxse = myResults[Scenario == (MR_Scenario) &
                           beta == "G_X", std.error],
        by = myResults[Scenario == (MR_Scenario) &
                         beta == "G_Y", estimate],
        byse = myResults[Scenario == (MR_Scenario) &
                           beta == "G_Y", std.error],
        # exposure = "Simulated Exposure",
        # outcome = "Simulated Outcome",
        snps = myResults[Scenario == (MR_Scenario) &
                           beta == "G_X", term]
      )
  }
  return(MRInputObject)
}
