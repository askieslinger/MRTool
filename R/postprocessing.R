#' Transform result of iterations of simulation and MR
#'
#' Transforms results of x iterations with each n scenario simulations into median per scenario
#' for relevant variables.
#' For examples see vignette browseVignettes('MRTool')
#'
#' @param results list of results of single iteration. result of single iteration is list(reg_res,sim_data,mr_res)
#' @param SNP names of snps used
#' @param iterations integer how many times the simulation was repeated
#' @param rev true if mr was executed in the anticausal direction. false for causal direction
#' @param pve_grid grid of percentages of variance explained as output of cal_betas
#'
#' @return list(res_list,sim_data). res_list
#'     is list of regression results and mr results. sim_data is result of simulation.
#'     'power' is proportion of p-value <0.05 of MR estimate per scenario over all iterations.
#'     caus_eff ist Lineare Regression zwischen X und Y in die Richtung der MR, also estimate_X_Y
#'     für die kausale Richtung und estimate_Y_X für die antikausale MR.
#'     Ergebnisse der MR sind 'Estimate', 'Std Error', 'P-Value', '95% CI upper'
#'     und '95% CI lower'.
#'     Im datatable davor sind Ergebnisse der LR und Parameter der Simulation.
#'     Am Ende des datatable sind die pves.
#'
#' @export
#'
transform_results <-
  function(results, SNP, iterations, pve_grid, rev = FALSE) {
    ###---- processing of regression results-----------------------------------------------------------------------
    # bind regression results of every iteration into one DT
    reg_res <- rbindlist(lapply(results, function(x)
      x$reg_res))
    # dcast into a wide DT (instead of long DT), take median of values for each scenario
    a <- data.table::dcast(
      rbindlist(lapply(SNP, function(x)
        reg_res[term == paste0('`', x, '`')])),
      Scenario + term + n ~ beta,
      value.var = c(
        'assigned.effect',
        'estimate',
        'std.error',
        'statistic',
        'p.value',
        'adjRsquared',
        'fstat',
        'sim_sd'
      ),
      fun = median
    )
    b <- data.table::dcast(
      rbindlist(lapply(c('U', 'X', 'Y'), function(x)
        reg_res[term == x])),
      Scenario + n ~ beta,
      value.var = c(
        'assigned.effect',
        'estimate',
        'std.error',
        'statistic',
        'p.value',
        'adjRsquared',
        'fstat'
      ),
      fun = median
    )
    reg_res_mean <- merge(a, b)
    ###------------processing of mr results------------------------------------
    # combine the mr results per snp, per scenario, per iteration into one data.table
    mrs <- rbindlist(lapply(
      results,
      FUN = function(iter)
        rbindlist(lapply(
          iter$mr_res,
          FUN = function(snp)
            rbindlist(lapply(
              snp,
              FUN = function(scen)
                scen$Values
            ), idcol = 'Scenario')
        ), idcol = 'snp')
    ), idcol = 'iteration')
    setDT(mrs)

    # data.table with median results of ivw method
    mrs_ivw_median <-
      mrs[Method == 'IVW', lapply(.SD, median), by = c('Scenario', 'snp'), .SDcols = c('Estimate',
                                                                                       'Std Error',
                                                                                       'P-value',
                                                                                       '95% CI upper',
                                                                                       '95% CI lower')]

    setkey(mrs_ivw_median, Scenario, snp)
    # calculate power
    ans <-
      mrs[Method == 'IVW', .(power = sum(`P-value` <= 0.05) / iterations), by = c('Scenario', 'snp')]
    setkey(ans, Scenario, snp)
    mrs_ivw_median <- merge(mrs_ivw_median, ans)

    ###------------------- processing sim_data -------------------------------------------------------------
    # check variance of X per Scenario
    rbindlist(lapply(
      results,
      FUN = function(x)
        x$sim_data[, var(X), by = Scenario]
    ))[, mean(V1), by = Scenario]
    # combine sim data for every iteration into one data.table
    sim_data <- rbindlist(lapply(results, function(x)
      x$sim_data), idcol = 'iteration')


    ###-----------Merging results of regression and mr for easy plotting ---------------------------------------------
    #rename term into snp_id
    reg_res_mean[, `:=`(snp_id = term)]
    reg_res_mean[, `:=`(term = NULL)]
    setkey(reg_res_mean, Scenario, snp_id)

    if (rev) {
      reg_res_mean[, caus_eff := estimate_Y_X]
    } else{
      reg_res_mean[, caus_eff := estimate_X_Y]
    }

    #split mr results into two DT by snp as IV, combine with regression results
    res_list <-
      lapply(
        split(mrs_ivw_median, by = 'snp'),
        FUN = function(snp_res) {
          reg_res_mean[snp_res, on = 'Scenario']
        }
      )
    #combine that with pves as assigned. two data.tables ,one for each snp as IV
    res_list <-
      lapply(
        res_list,
        FUN = function(snp_res)
          snp_res[pve_grid, on = c('Scenario', 'snp_id')]
      )

    # return
    list(res_list = res_list, sim_data = sim_data)
  }
