#===============================================================================
#
# MAST w/ RANDOM EFFECT
#
#===============================================================================
# See code provided by Ricardo Albanus (https://github.com/rdalbanus):
# https://github.com/satijalab/seurat/issues/3712#issuecomment-1379578940

# We want a different behavior from MASTDETest, where the latent variable is
# used as a random effect in the model
rem_MASTDETest <- function(
    data.use,
    cells.1,
    cells.2,
    latent.vars = NULL,
    verbose = TRUE,
    # New option - random effect variable (should be included in latent.vars)
    re.var = NULL,
    ...
) {
  # Check for MAST
  if (!PackageCheck('MAST', error = FALSE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  
  require(MAST)
  require(glue)
  
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  
  # This is the main change in the code - we want ~ ( 1 | re.var) in the formula:
  # ~ condition + lat.vars + (1 | re.var)
  if (!is.null(re.var)) {
    if (!re.var %in% latent.vars.names) {
      stop("Random effect variable should be included in latent variables!")
    }
    latent.vars.names <- latent.vars.names[!latent.vars.names %in% re.var]
    fmla <- as.formula(
      object = paste0(
        " ~ ", paste(latent.vars.names, collapse = "+"), glue("+ (1|{re.var})")
      )
    )
    # print(fmla)
    # We need glmer to make this work
    method <-  "glmer" # trying to troubleshoot this - it can clash with the already existing method var in the function call    
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, method = "glmer", ...)
  } else {
    # Original code
    fmla <- as.formula(
      object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
    )
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  }
  summaryCond <- MAST::summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable
  
  # The output format is slightly different, so we need adapt the code
  if(!is.null(re.var)) {
    p_val <- summaryDt[summaryDt$"component" == "H", 4]$`Pr(>Chisq)`
    genes.return <- summaryDt[summaryDt$"component" == "H", 1]$primerid
  } else {
    p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
    genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  }
  
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

# This is the zlm function with a minor change to fix a buggy variable being passed possibly due to namespace clash
my_zlm <- function (formula, sca, method = "bayesglm", silent = TRUE, ebayes = TRUE, 
                    ebayesControl = NULL, force = FALSE, hook = NULL, parallel = TRUE, 
                    LMlike, onlyCoef = FALSE, exprs_values = MAST::assay_idx(sca)$aidx, 
                    ...) 
{
  dotsdata = list(...)$data
  if (!is.null(dotsdata)) {
    if (!missing(sca)) 
      stop("Cannot provide both `sca` and `data`")
    sca = dotsdata
  }
  if (!inherits(sca, "SingleCellAssay")) {
    if (inherits(sca, "data.frame")) {
      if (!is.null(dotsdata)) {
        return(.zlm(formula, method = method, silent = silent, 
                    ...))
      }
      else {
        return(.zlm(formula, data = sca, method = method, 
                    silent = silent, ...))
      }
    }
    else {
      stop("`sca` must inherit from `data.frame` or `SingleCellAssay`")
    }
  }
  if (missing(LMlike)) {
    method <- match.arg(method, MAST:::methodDict[, keyword])
    method <- MAST:::methodDict[keyword == method, lmMethod]
    if (!is(sca, "SingleCellAssay")) 
      stop("'sca' must be (or inherit) 'SingleCellAssay'")
    if (!is(formula, "formula")) 
      stop("'formula' must be class 'formula'")
    Formula <- MAST:::removeResponse(formula)
    priorVar <- 1
    priorDOF <- 0
    if (ebayes) {
      # if (!methodDict[lmMethod == method, implementsEbayes]) 
      if (!all(MAST:::methodDict[lmMethod == method, implementsEbayes]))
        stop("Method", method, " does not implement empirical bayes variance shrinkage.")
      ebparm <- MAST:::ebayes(t(SummarizedExperiment:::assay(sca, exprs_values)), ebayesControl, 
                              model.matrix(Formula, SummarizedExperiment::colData(sca)))
      priorVar <- ebparm[["v"]]
      priorDOF <- ebparm[["df"]]
      stopifnot(all(!is.na(ebparm)))
    }
    obj <- new_with_repaired_slots(classname = method, design = SummarizedExperiment::colData(sca), 
                                   formula = Formula, priorVar = priorVar, priorDOF = priorDOF, 
                                   extra = list(...))
  }
  else {
    if (!missing(formula)) 
      warning("Ignoring formula and using model defined in 'objLMLike'")
    if (!inherits(LMlike, "LMlike")) 
      stop("'LMlike' must inherit from class 'LMlike'")
    obj <- LMlike
  }
  ee <- t(SummarizedExperiment:::assay(sca, exprs_values))
  genes <- colnames(ee)
  ng <- length(genes)
  MM <- MAST:::model.matrix(obj)
  coefNames <- colnames(MM)
  listEE <- setNames(seq_len(ng), genes)
  obj <- MAST:::fit(obj, ee[, 1], silent = silent)
  nerror <- totalerr <- 0
  pb = progress::progress_bar$new(total = ng, format = " Completed [:bar] :percent with :err failures")
  .fitGeneSet <- function(idx) {
    hookOut <- NULL
    tt <- try({
      obj <- MAST:::fit(obj, response = ee[, idx], silent = silent, 
                        quick = TRUE)
      if (!is.null(hook)) 
        hookOut <- hook(obj)
      nerror <- 0
    })
    if (is(tt, "try-error")) {
      obj@fitC <- obj@fitD <- NULL
      obj@fitted <- c(C = FALSE, D = FALSE)
      nerror <- nerror + 1
      totalerr = totalerr + 1
      if (nerror > 5 & !force) {
        stop("We seem to be having a lot of problems here...are your tests specified correctly?  \n If you're sure, set force=TRUE.", 
             tt)
      }
    }
    pb$tick(tokens = list(err = totalerr))
    if (onlyCoef) 
      return(cbind(C = coef(obj, "C"), D = coef(obj, "D")))
    summaries <- MAST:::summarize(obj)
    structure(summaries, hookOut = hookOut)
  }
  if (!parallel || getOption("mc.cores", 1L) == 1) {
    listOfSummaries <- lapply(listEE, .fitGeneSet)
  }
  else {
    listOfSummaries <- parallel::mclapply(listEE, .fitGeneSet, 
                                          mc.preschedule = TRUE, mc.silent = silent)
  }
  if (onlyCoef) {
    out <- do.call(abind, c(listOfSummaries, rev.along = 0))
    return(aperm(out, c(3, 1, 2)))
  }
  cls <- sapply(listOfSummaries, function(x) class(x))
  complain <- if (force) 
    warning
  else stop
  if (mean(cls == "try-error") > 0.5) 
    complain("Lots of errors here..something is amiss.")
  hookOut <- NULL
  if (!is.null(hook)) 
    hookOut <- lapply(listOfSummaries, attr, which = "hookOut")
  message("\nDone!")
  summaries <- MAST:::collectSummaries(listOfSummaries)
  summaries[["LMlike"]] <- obj
  summaries[["sca"]] <- sca
  summaries[["priorVar"]] <- obj@priorVar
  summaries[["priorDOF"]] <- obj@priorDOF
  summaries[["hookOut"]] <- hookOut
  summaries[["exprs_values"]] <- exprs_values
  summaries[["Class"]] <- "ZlmFit"
  zfit <- do.call(new, as.list(summaries))
  zfit
}

# We will replace the original function with ours inside the Seurat namespace
assignInNamespace('MASTDETest', rem_MASTDETest, asNamespace("Seurat"))
# getFromNamespace("MASTDETest", "Seurat")

assignInNamespace('zlm', my_zlm, asNamespace("MAST"))
# getFromNamespace("zlm", "MAST")


## FROM MAST
# https://github.com/RGLab/MAST/blob/devel/R/helper-methods.R
new_with_repaired_slots = function(classname, ..., extra){
  #... were actually named in caller, so assumed good.
  # extra are dot args, and might have been passed in erroneously by a caller.
  bad_slots = setdiff(names(extra), slotNames(classname))
  if(length(bad_slots) > 0){
    warning(sprintf("Dropping illegal slot(s) %s for class %s.  
                    This likely indicates a bug in an upstream package.", 
                    paste0(bad_slots, collapse = ', '), classname))
  }
  do.call(new, c(Class = classname, list(...), extra[setdiff(names(extra), bad_slots)]))
}

