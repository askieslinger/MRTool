#' calculate correlation for a datatable of snps
#'
#' returns 2 snps each that have lowest,
#' positive high, positive low, negative high, negative low correlation. uses pearson correlation
#' For examples see vignette browseVignettes('MRTool')
#'
#' @param genedose_test snp datatable. first two columns identifiers
#' @param method used to compute correlation coefficient. 'pearson' (default), 'kendall' or 'spearman'
#' @return
#' @export
#'
get_corrs <- function(genedose_test,method='pearson') {

  snp_names <- names(genedose_test)[-(1:2)]
  cor_m <- cor(genedose_test[, ..snp_names],method=method)
  cor_m[cor_m == 1] = NA
  snp_cor <-
    data.table(cor = c('pos_high', 'pos_low', 'neg_high', 'neg_low', 'lowest'))

  snp_cor=tryCatch({
  inds <- arrayInd(which.min(abs(cor_m)), .dim = dim(cor_m))
  snp_cor[cor == 'lowest', `:=`(snp1 = rownames(cor_m)[inds[, 1]],
                                snp2 = colnames(cor_m)[inds[, 2]])]
}, error=function(e){
  print('Es wurden keine SNPs mit einer negativ-hohen Korrelation gefunden.')
  })


  snp_cor=tryCatch({
  inds <-
    arrayInd(sample(which(cor_m > -0.31 &
                            cor_m < (-0.29)), 1), .dim = dim(cor_m))
  snp_cor[cor == 'neg_high', `:=`(snp1 = rownames(cor_m)[inds[, 1]],
                                  snp2 = colnames(cor_m)[inds[, 2]])]
  }, error=function(e){
    print('Es wurden keine SNPs mit einer negativ-hohen Korrelation gefunden.')
    snp_cor=snp_cor[cor!='neg_high']
    return(snp_cor)}
  )

  snp_cor=tryCatch({
  inds <-
    arrayInd(sample(which(cor_m > 0.29 & cor_m < 0.31), 1), .dim = dim(cor_m))
  snp_cor[cor == 'pos_high', `:=`(snp1 = rownames(cor_m)[inds[, 1]],
                                  snp2 = colnames(cor_m)[inds[, 2]])]
  }, error=function(e){
    print('Es wurden keine SNPs mit einer positiv-hohen Korrelation gefunden.')
    snp_cor=snp_cor[cor!='pos_high']
    return(snp_cor)}
  )

  snp_cor=tryCatch({
  inds <-
    arrayInd(sample(which(cor_m > 0.09 & cor_m < 0.11), 1), .dim = dim(cor_m))
  snp_cor[cor == 'pos_low', `:=`(snp1 = rownames(cor_m)[inds[, 1]],
                                 snp2 = colnames(cor_m)[inds[, 2]])]
  }, error=function(e){
    print('Es wurden keine SNPs mit einer positiv-niedrigen Korrelation gefunden.')
    snp_cor=snp_cor[cor!='pos_low']
    return(snp_cor)}
  )

  snp_cor=tryCatch({
  inds <-
    arrayInd(sample(which(cor_m > -0.11 &
                            cor_m < (-0.09)), 1), .dim = dim(cor_m))
  snp_cor[cor == 'neg_low', `:=`(snp1 = rownames(cor_m)[inds[, 1]],
                                 snp2 = colnames(cor_m)[inds[, 2]])]
  }, error=function(e){
    print('Es wurden keine SNPs mit einer negativ-niedrigen Korrelation gefunden.')
    snp_cor=snp_cor[cor!='neg_low']
    return(snp_cor)}
  )

  if(nrow(snp_cor==1)){
    snp_cor[, cor_val := cor_m[snp1, snp2]]
  }else{
  snp_cor[, cor_val := diag(cor_m[snp1, snp2])]
  }

  snp_cor
}


#' calculates betacoefficients from explained variances
#'
#' For examples see vignette browseVignettes('MRTool')
#'
#' @param pves a named list of parameters and percentages of variance explained for each one.
#'     If multiple snps, each param is a list of snp-vectors, each vector the scenarios for this snp
#' @param maf a named list of maf for every snp
#' @param gy_direct should there be a direct effect of G on Y. default true. does not need to be changed.
#' @return a list of two data.tables. a data.table containing all valid combinations
#'     of coefficients, one row is one snp in one scenario. A data.table with the adjusted combinations
#'     of percentages of variance explained
#' @export
#'
#'
#' @import data.table
#'
cal_betas <- function(pves, maf, gy_direct = T) {

  if((sum(unlist(pves)>1)>0)|| (sum(unlist(pves)<0)>0)){
    stop('All percentages of variance explained must be between 0 and 1!')
  }

  #unterscheidung 1 oder mehrere snps
  if (length(maf) == 1) {
    #grid mit den kombinationen
    pve_grid <- data.table(expand.grid(pves))
    pve_grid[, Scenario := paste0("V", 1:.N)]
    setkey(pve_grid, Scenario)
    # varianz von genotyp
    var_g <- 2 * maf * (1 - maf)
    # neues grid mit den koeffs
    beta_grid <- data.table(Scenario = pve_grid$Scenario)
    setkey(beta_grid, Scenario)
    # begin calculating coefficients
    beta_grid[, `:=`(
      G_U = sqrt(pve_grid$G_U / var_g),
      simsd_U = sqrt(1 - pve_grid$G_U)
    )]
    if(sum(is.nan(unlist(beta_grid)))>0){
      stop('percentages of variance explained for GU is impossible')
    }
    beta_grid[,`:=`(
      U_X = sqrt(pve_grid$U_X/(simsd_U^2))
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for UX is impossible')
    }

    beta_grid[, `:=`(
      G_X = sqrt(pve_grid$G_X / var_g) - U_X * G_U,
      simsd_X = sqrt(1 - pve_grid$G_X - simsd_U ^ 2 * U_X ^ 2)
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for GX and or UX is impossible')
    }
    beta_grid[,`:=`(
      X_Y = sqrt(pve_grid$X_Y/(simsd_X^2))
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for XY is impossible')
    }
    beta_grid[,`:=`(
      U_Y = pmax(0,sqrt(pve_grid$U_Y/(simsd_U^2)) - X_Y * U_X)
    )]

    pve_grid[,U_Y:=beta_grid[,(U_Y+X_Y*U_X)^2*simsd_U^2]]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for UY is impossible')
    }
    pve_grid[,snp_id:=names(maf)]
    pve_grid[,snp_id:=lapply(.SD,FUN=function(name)paste0('`',name,'`')),.SDcols='snp_id']

    # unterscheidung ob direkter effekt von g auf y vorhanden
    # nicht mehr nötig
    if (gy_direct == F) {
      beta_grid[, `:=`(
        G_Y = 0,
        simsd_Y = sqrt(
          1 - (X_Y * (G_X + U_X * G_U) + U_Y * G_U) ^ 2 * var_g - X_Y ^ 2 * simsd_X ^
            2 - pve_grid$U_Y * simsd_U ^ 2
        )
      )]
    } else{

      beta_grid[,  `:=`(
        G_Y =pmax(0,sqrt(pve_grid$G_Y / var_g) - (X_Y * (G_X + U_X * G_U) + U_Y * G_U))
      )]
      beta_grid[,  `:=`(simsd_Y=sqrt(
        1 - (X_Y * (G_X + U_X * G_U) + U_Y * G_U+G_Y) ^ 2 * var_g - X_Y ^ 2 * simsd_X ^
          2 - pve_grid$U_Y * simsd_U ^ 2
      )
      )]
    }
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for GY and or XY and or UY is impossible')
    }

    setorder(x = beta_grid, G_U, G_Y, X_Y, G_X)

    setattr(beta_grid, "Vary_between_SNP", F)
    setattr(beta_grid, "G_X_randomization", F)
    setattr(beta_grid, "G_X_between", NA)
    list(beta_grid,pve_grid)
  } else{
    ### mehr als 1 SNP ######################################

    no_scens <- unique(lengths(pves$G_X))
    if (length(no_scens)>1 || length(unique(lengths(pves$G_Y))) >1 || length(unique(lengths(pves$G_U))) >1){
      stop('number of scenarios per snp not identical. \nPlease input a list(snp1=c(scen1,scen2,..),snp2=(scen1,scen2,..),...) for every parameter G_*')
    }
    no_scens <- unlist(no_scens)
    if(!isTRUE(all.equal(length(pves$G_X),length(pves$G_U))) || !isTRUE(all.equal(length(pves$G_X),length(pves$G_Y)))){
      stop('number of snps in G_X,G_U,G_Y not identical. \nPlease input a list(snp1=c(scen1,scen2,..),snp2=(scen1,scen2,..),...)for every parameter G_*')
    }
    if(!isTRUE(all.equal(unique(lengths(pves$G_Y)),no_scens) )|| !isTRUE(all.equal(unique(lengths(pves$G_U)),no_scens))){
      stop('number of scenarios per snp not identical. \nPlease input a list(snp1=c(scen1,scen2,..),snp2=(scen1,scen2,..),...)for every parameter G_*')
    }
    # kleine dt mit snps als columns
    gx <- as.data.table(pves$G_X)
    setnames(gx, names(maf))
    gy <- as.data.table(pves$G_Y)
    setnames(gy, names(maf))
    gu <- as.data.table(pves$G_U)
    setnames(gu, names(maf))

    # pro snp, mache grid mit kombis der ux,uy,xy, g variablen in fester kombi drangehaengt
    tmp <- lapply(names(maf), function(snp) {
      x <-
        data.table(expand.grid(c(pves[4:6], list(G_X =unlist( gx[, ..snp])
        )
        )
        )
        )[, Scenario := paste0('V', 1:.N)]
      y <-
        data.table(expand.grid(c(pves[4:6], list(G_Y = unlist(gy[, ..snp])))))[, Scenario := paste0('V', 1:.N)]
      u <-
        data.table(expand.grid(c(pves[4:6], list(G_U = unlist(gu[, ..snp]))
        )
        )
        )[, Scenario := paste0('V', 1:.N)]
      x[,`:=`(G_Y=y[,G_Y],G_U=u[,G_U])]
    })

    names(tmp)<-names(maf)

    # fuege dts pro snp untereinander zusammen
    pve_grid <- rbindlist(tmp, idcol = 'snp_id')
    setkey(pve_grid, Scenario)
    # fuege varianz pro snp hinzu
    var_g <- 2 * maf * (1 - maf)
    varg_dt <-as.data.table(var_g,keep.rownames = 'snp_id')
    setkey(varg_dt,snp_id)
    # neues grid fuer koeffizienten
    beta_grid <-
      data.table(Scenario = pve_grid$Scenario, snp_id = pve_grid$snp_id)
    setkeyv(beta_grid,c('Scenario','snp_id'))
    beta_grid <- beta_grid[varg_dt,on='snp_id']
    setkey(beta_grid, Scenario)
    beta_grid[, `:=`(
      G_U = sqrt(pve_grid$G_U / var_g)
    )]

    beta_grid <- beta_grid[pve_grid[,.(simsd_U=sqrt(1-sum(G_U))),by=Scenario]]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for GU is impossible')
    }
    beta_grid[,`:=`(
      U_X = sqrt(pve_grid$U_X/(simsd_U^2))
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for UX and or GX is impossible')
    }

    beta_grid <- beta_grid[pve_grid[,.(simsd_X=sum(G_X)),by=Scenario]]

    beta_grid[, `:=`(
      G_X = sqrt(pve_grid$G_X / var_g) - U_X * G_U,
      simsd_X =sqrt(1-simsd_X -simsd_U^2*U_X^2)
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for GX and or UX is impossible')
    }
    beta_grid[,`:=`(
      X_Y = sqrt(pve_grid$X_Y/(simsd_X^2))
    )]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for XY is impossible')
    }

    beta_grid[,`:=`(
      U_Y = pmax(0,sqrt(pve_grid$U_Y/(simsd_U^2)) - X_Y * U_X)
    )]
    pve_grid[,U_Y:=beta_grid[,(U_Y+X_Y*U_X)^2*simsd_U^2]]
    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for UY and or XY is impossible')
    }

    # unterscheidung, ob direkter effect g auf y vorhanden
    if (gy_direct == F) {
      beta_grid[,simsd_Y:=sum((X_Y * (G_X + U_X * G_U) + U_Y * G_U) ^ 2 * var_g),by=Scenario]
      beta_grid[, `:=`(
        G_Y = 0,
        simsd_Y = sqrt(
          1 - simsd_Y - X_Y ^ 2 * simsd_X ^
            2 - pve_grid$U_Y * simsd_U ^ 2
        )
      )]
    } else{
      beta_grid[,  `:=`(
        G_Y =pmax(0,sqrt(pve_grid$G_Y / var_g) - (X_Y * (G_X + U_X * G_U) + U_Y * G_U))
      )]
      beta_grid <- beta_grid[,simsd_Y:=sum((X_Y * (G_X + U_X * G_U) + U_Y * G_U+G_Y) ^ 2 * var_g),by=Scenario]
      beta_grid[,simsd_Y := sqrt(1-simsd_Y - X_Y ^ 2 * simsd_X ^ 2 - pve_grid$U_Y * simsd_U ^ 2 )]
    }

    if(sum(sapply(beta_grid,is.nan))>0){
      stop('percentages of variance explained for GY and or XY and or UY is impossible')
    }
    setorder(x = beta_grid, G_U, G_Y, X_Y, G_X)
    setkeyv(beta_grid,c('Scenario','snp_id'))
    # legacy
    setattr(beta_grid, "Vary_between_SNP", F)
    setattr(beta_grid, "G_X_randomization", F)
    setattr(beta_grid, "G_X_between", NA)
    pve_grid[,snp_id:=lapply(.SD,FUN=function(name)paste0('`',name,'`')),.SDcols='snp_id']

    # return
    list(beta_grid,pve_grid)
  }
}






#' Set MR Parameters (betacoefficients)
#'
#' Create a testplan for different parameter combinations
#' or create a custom plan for different IV parameters. If Vary_between_SNP, The Parameters set
#' for G_U, G_X and G_Y will not create individual scenarios but one, where the different
#' parameters are distributed to the individual SNP, cycled through all give SNP, based on the
#' number of parameters provided by the user, accepting only one parameter (the first)
#' for U_Y, U_X and X_Y
#' For examples see vignette browseVignettes('MRTool')
#'
#'
#' @param Vary_between_SNP
#' @param G_X_randomization
#' @param SNP
#' @param G_U
#' @param G_Y
#' @param X_Y
#' @param G_X
#' @param U_X
#' @param U_Y
#' @param sim_sd
#' @param maf
#'
#'
#' @return
#' @export
#'
SetMRParams <- function(Vary_between_SNP = F,
                        G_X_randomization = F,
                        SNP = NULL,
                        G_U = c(0, 0.1),
                        G_Y = c(0, 0.1),
                        X_Y = c(0, 0.1),
                        G_X = c(0.03, 0.1),
                        U_X = c(0, 0.1),
                        U_Y = c(0, 0.1),
                        sim_sd = 1,
                        maf) {
  if (is.null(SNP)) {
    stop("Please select >0 SNP to analyze!")
  }

  # Error can't be 0, which will cause singularity issues
  if (length(sim_sd) == 1 && sim_sd == 0) {
    print(
      "Sorry, an error of 0 will cause singularity issues later on and thus causes errors. Please choose an error weight >0"
    )
    stop()
  }


  # means create different scenarios
  if (Vary_between_SNP == F & G_X_randomization == T) {
    # set boundaries for a runif creation of G_X parameter for each SNP

    # Create the Testplan
    testplan2 <-
      data.table(expand.grid(
        G_U = G_U,
        # MR CONDITION (must be 0)
        G_Y = G_Y,
        # MR CONDITION (must be 0)
        X_Y = X_Y,
        # causal effect
        G_X = paste0(G_X[1], " < G_X < ", G_X[2]),
        # MR CONDITION
        U_X = U_X,
        # not measured in Burgess
        U_Y = U_Y
      )) # not measured in Burgess

    setorder(x = testplan2, G_U, G_Y, X_Y, G_X)
    testplan2[, Scenario := paste0("V", 1:.N)]

    # make the Vary_between_SNP info accessible accessible
    attr(testplan2, "Vary_between_SNP") <- F
    attr(testplan2, "G_X_randomization") <- T

    # make the G_X_between values accessible
    attr(testplan2, "G_X_between") <- c(G_X[1], G_X[2])
    #print('vary=F,G_X_rando=T')
    #print(testplan2[])
    return(testplan2[])

    # create testplan with single set G_X parameter for each SNP
  } else if (Vary_between_SNP == F & G_X_randomization == F) {
    # Create the Testplan
    testplan2 <-
      data.table(expand.grid(
        G_U = G_U,
        # MR CONDITION (must be 0)
        G_Y = G_Y,
        # MR CONDITION (must be 0)
        X_Y = X_Y,
        # causal effect
        G_X = G_X,
        # MR CONDITION
        U_X = U_X,
        # not measured in Burgess
        U_Y = U_Y,
        sim_sd = sim_sd
      )) # not measured in Burgess
    setorder(x = testplan2, G_U, G_Y, X_Y, G_X)
    testplan2[, Scenario := paste0("V", 1:.N)]
    attr(testplan2, "Vary_between_SNP") <- F
    attr(testplan2, "G_X_randomization") <- F
    attr(testplan2, "G_X_between") <- NA
    #print('vary=F,G_X_rando=F')
    #print(testplan2[])
    return(testplan2[])

    # tell to set parameter properly END OF THIS SUB-LOOP
  } else if (Vary_between_SNP == T & G_X_randomization == T) {
    # Create only one scenario with variation between SNP
    # create a differing scenario with random G_X effect set between 2 boundaries
    testplan <- data.table(
      SNP = SNP,
      G_U = G_U,
      G_Y = G_Y,
      X_Y = X_Y[1],
      G_X = paste0(G_X[1], " < G_X < ", G_X[2]),
      U_X = U_X[1],
      U_Y = U_Y[1]
    )
    testplan[, Scenario := "V1"]
    attr(testplan, "Vary_between_SNP") <- T
    attr(testplan, "G_X_between") <- c(G_X[1], G_X[2])
    attr(testplan, "G_X_randomization") <- T

    #print('vary=T,G_X_rando=T')
    #print(testplan[])

    return(testplan[])
  } else if (Vary_between_SNP == T & G_X_randomization == F) {
    # # create part of test plan that is G specific =======================
    # testplan <- expand.grid(SNP = SNP,
    #                         G_U = G_U,
    #                         G_Y = G_Y,
    #                         G_X = G_X,
    #                         stringsAsFactors = F)
    # testplan <- testplan[rep(row.names(testplan), length(SNP)),]
    #
    # SNP <- rep(SNP, each = length(SNP))
    #

    testplan <- data.frame(
      SNP = SNP,
      G_U = G_U,
      G_Y = G_Y,
      G_X = G_X,
      stringsAsFactors = F
    )

    # create the remaining test plan
    non.g <- expand.grid(
      X_Y = X_Y,
      U_X = U_X,
      U_Y = U_Y,
      stringsAsFactors = F
    )

    # replicate for each SNP per scenario
    non.g.extend <- non.g[rep(row.names(non.g),
                              each = length(SNP)), ]

    testplan2 <- cbind(testplan, non.g.extend)
    setDT(testplan2)
    scenarios <- paste0("V", rep(1:nrow(non.g),
                                 each = length(SNP)))
    testplan2[, Scenario := (scenarios)]

    # add attributes for easier handling
    attr(testplan2, "Vary_between_SNP") <- T
    attr(testplan2, "G_X_randomization") <- F
    attr(testplan2, "G_X_between") <- NA

    #print('vary=T,G_X_rando=F')
    #print(testplan2[])

    return(testplan2[])

  } else{
    print(
      "Please set a valid (TRUE or FALSE) 'G_X_randomization' and 'Vary_between_SNP parameter"
    )
  }
}
