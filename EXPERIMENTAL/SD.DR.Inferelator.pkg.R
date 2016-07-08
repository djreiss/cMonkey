DATE <-
"Mon Jan 31 10:43:00 2011"
VERSION <-
"0.0.6"
.onLoad <-
function( libname, pkgname ) { ##.onAttach
    cat( "Loading ", pkgname, " version ", VERSION, " (", DATE, ")\n", sep="" )
    cat( "Copyright (C) David J Reiss, Institute for Systems Biology; dreiss.isb@gmail.com.org.\n" )
    cat( "Modified by Sam Danziger, Institute for Systems Biology; sdanziger@systemsbiology.org.\n" )
    cat( "http://baliga.systemsbiology.net/cmonkey\n" )
    cat( "\nNOTE that this package is still sloppy in that it relies upon some global variables:\n" )
    cat( "\nThis version will accept priors\n" )
    cat( "'predictor.mats', 'envMap', 'colMap', and optionally 'predictors'.\n" )
  }

combine.symbol <-
"~~"
cv.glmnet <-
function (x, y, lambda, K = 10, cv.reps = 10, trace = FALSE, 
    plot.it = TRUE, se = TRUE, weights = NA, ...) 
{
    all.folds <- do.call("c", lapply(1:cv.reps, function(i) cv.folds(length(y), 
        K)))
    apply.func <- get.apply.func()
    if (is.na(weights[1])) 
        weights <- rep(1, length(y))
    residmat <- do.call(cbind, apply.func(1:length(all.folds), 
        function(i) {
            if (trace) 
                cat("CV Fold", i, "\n")
            omit <- all.folds[[i]]
            fit <- my.glmnet(x[-omit, ], y[-omit], lambda = lambda, 
                weights = weights[-omit], ...)
            fit <- predict.elnet(fit, x[omit, , drop = FALSE], 
                ...)
            if (length(omit) == 1) {
                fit <- matrix(fit, nrow = 1)
            }
            apply((y[omit] - fit)^2, 2, mean)
        }))
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(cv = cv, cv.error = cv.error, lambda = lambda, 
        fraction = log(lambda))
    if (plot.it) {
        plot(lambda, cv, type = "b", ylim = range(cv, cv + cv.error, 
            cv - cv.error), ...)
        if (se) 
            error.bars(lambda, cv + cv.error, cv - cv.error, 
                width = 1/length(lambda))
    }
    invisible(object)
}
fill.in.time.series <-
function (cols, col.map, fill.all.ts.frac = 1.25, remove.all.ts.frac = 0.1, 
    fill.ts.gap.size = 1, ...) 
{
    out.cols <- cols
    all.cols <- rownames(col.map)
    firsts <- which(col.map$is1stLast == "f")
    lasts <- which(col.map$is1stLast == "l")
    ts <- lapply(1:length(firsts), function(i) all.cols[firsts[i]:lasts[i]])
    cols.frac <- length(cols)/length(all.cols)
    for (j in 1:length(ts)) {
        frac <- mean(ts[[j]] %in% cols)
        if (frac >= cols.frac * fill.all.ts.frac) {
            new.cols <- ts[[j]][!(ts[[j]] %in% cols)]
            out.cols <- unique(c(out.cols, new.cols))
        }
        else if (frac < cols.frac * remove.all.ts.frac) {
            remove.cols <- ts[[j]][ts[[j]] %in% cols]
            out.cols <- out.cols[!(out.cols %in% remove.cols)]
        }
        else {
            for (q in 1:fill.ts.gap.size) {
                is.in <- which(ts[[j]] %in% out.cols)
                expand <- unique(sort(c(is.in, is.in - 1, is.in + 
                  1)))
                expand <- expand[expand > 0 & expand <= length(ts[[j]])]
                out.cols <- unique(c(out.cols, ts[[j]][expand]))
            }
        }
    }
    out.cols
}
filter.by.aic <-
function (mean.profile, predictor.matrix, top.aic.to.keep, r.min.cutoff = 0.2, 
    force.pos.neg = NA, ...) 
{
    cors <- apply(predictor.matrix, 1, cor, mean.profile, use = "pairwise")
    paic <- rep(Inf, nrow(predictor.matrix))
    apply.func <- get.apply.func()
    tmp <- unlist(apply.func(which(abs(cors) >= r.min.cutoff), 
        function(j) AIC(lm(mean.profile ~ as.numeric(predictor.matrix[j, 
            ])))))
    paic[which(abs(cors) >= r.min.cutoff)] <- tmp
    names(paic) <- rownames(predictor.matrix)
    if (is.na(force.pos.neg)) {
        best.preclusts <- names(sort(paic))[1:min(length(paic), 
            top.aic.to.keep)]
    }
    else {
        best.preclusts <- c(names(sort(cors, decreasing = T))[1:min(length(cors), 
            round(top.aic.to.keep * force.pos.neg))], names(sort(cors))[1:min(length(cors), 
            round(top.aic.to.keep * (1 - force.pos.neg)))])
    }
    result <- predictor.matrix[best.preclusts, , drop = F]
    return(result)
}
get.apply.func <-
function (plot = F) 
if (multicore:::isChild() || plot || (exists("DEBUG") && DEBUG)) lapply else mclapply
get.boot.coef.quantiles <-
function (coef.obj, use.grep = F) 
{
    coef.quantiles <- NULL
    n.boot <- length(coef.obj$coeffs.boot)
    if (n.boot > 1) {
        tmp <- unlist(coef.obj$coeffs.boot)
        tmp2 <- table(names(tmp))
        coef.quantiles <- t(sapply(names(tmp2), function(i) {
            if (!use.grep) 
                tmp3 <- tmp[names(tmp) == i]
            else tmp3 <- tmp[grepl(i, names(tmp))]
            tmp3 <- c(tmp3, rep(0, max(0, n.boot - length(tmp3))))
            c(n = sum(names(tmp) == i)/n.boot, quantile(abs(tmp3), 
                prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)))
        }))
        coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
            1, function(i) all(i[-1] == 0)), ]
        print(coef.quantiles, digits = 3)
    }
    invisible(coef.quantiles)
}
get.cluster.predictors <-
function (cluster.rows, cluster.profile, predictor.mat, predictor.mat.ands, 
    tf.groups, env.names, r.cutoff = 1.0, aic.filter.cutoff = 25, 
    quiet = F, force.include = env.names, ...) 
{
    tf.groups <- tf.groups[names(tf.groups) %in% rownames(predictor.mat)]
    force.include <- force.include[force.include %in% names(tf.groups)]
    predictor.mat.env <- predictor.mat[force.include, , drop = F]
    possibly.regulates <- numeric()
    is.in <- sapply(tf.groups, function(j) mean(j %in% cluster.rows, 
        na.rm = T))
    high.cor <- apply(predictor.mat[names(tf.groups), ], 1, cor, 
        cluster.profile, use = "pairwise")
    high.cor[is.na(high.cor)] <- 0
    if (any(is.na(high.cor))) {
        tmp <- which(is.na(high.cor))
        is.in[tmp] <- 1
        high.cor[tmp] <- 1
    }
    if (any(is.in > 0 | high.cor > r.cutoff)) {
        possibly.regulates <- c(is.in[is.in > 0], high.cor[is.in <= 
            0 & high.cor > r.cutoff])
        is.in <- which(is.in > 0 | high.cor > r.cutoff)
        predictor.mat <- predictor.mat[-is.in, ]
        is.in <- unique(sort(unlist(lapply(names(is.in), function(j) grep(j, 
            rownames(predictor.mat.ands))))))
        predictor.mat.ands <- predictor.mat.ands[-is.in, ]
    }
    if (!quiet && length(possibly.regulates) > 0) 
        cat("POSSIBLY REGULATES (but removed from regression):", 
            paste(unique(names(possibly.regulates)), sep = ", "), 
            "\n")
    if (!is.na(aic.filter.cutoff) && aic.filter.cutoff != 0) {
        best.singletons <- filter.by.aic(cluster.profile, predictor.matrix = predictor.mat, 
            top.aic.to.keep = aic.filter.cutoff, ...)
        best.combos <- NULL
        if (!is.null(predictor.mat.ands)) {
            best.combos <- filter.by.aic(cluster.profile, predictor.matrix = predictor.mat.ands, 
                top.aic.to.keep = nrow(best.singletons) * 20, 
                ...)
            if (!is.na(r.cutoff) && r.cutoff < 1) {
                tmp.cor <- t(cor(t(best.singletons), t(best.combos)))
                to.elim <- which(tmp.cor > r.cutoff - 0.2, arr = T)
                best.combos <- best.combos[!rownames(best.combos) %in% 
                  rownames(to.elim), , drop = F]
            }
            if (nrow(best.combos) > aic.filter.cutoff) 
                best.combos <- best.combos[1:aic.filter.cutoff, 
                  ]
            if (nrow(best.combos) > nrow(best.singletons) + length(env.names)) 
                best.combos <- best.combos[1:(nrow(best.singletons) + 
                  length(env.names)), ]
        }
        predictor.mat <- unique(rbind(predictor.mat[rownames(best.singletons), 
            , drop = F], predictor.mat.env, predictor.mat.ands[rownames(best.combos), 
            , drop = F]))
    }
    else {
        predictor.mat <- unique(rbind(predictor.mat, predictor.mat.env, 
            predictor.mat.ands))
    }
    list(possibly.regulates = possibly.regulates, predictors = rownames(predictor.mat))
}
get.input.matrix <-
function (profile, predictor.mat, conds.use = "ALL", col.map = NULL, 
    tau = 10, ratio.cutoff = 3, quiet = F) 
{
    out.tmp <- profile
    in.tmp <- predictor.mat
    if (!is.null(col.map) && !is.na(tau) && tau > 0) {
        if (!quiet) 
            cat("Time series data supplied, converting predictors to difference equation... \n")
        if (!quiet) 
            cat("Tau =", tau, "\n")
        conds <- colnames(predictor.mat)
        cm <- col.map[conds, ]
        good.i <- ((cm$isTs == TRUE) & (cm$is1stLast %in% c("m", 
            "l"))) | (cm$isTs == FALSE & cm$is1stLast == "e")
        curs <- as.character(cm$condName[good.i])
        prevs <- as.character(cm$prevCol[good.i])
        prevs[is.na(prevs)] <- curs[is.na(prevs)]
        del.ts <- as.numeric(as.character(cm$delta.t[good.i]))
        del.ts[del.ts < 1] <- 1
        tmp <- curs %in% names(out.tmp) & prevs %in% names(out.tmp) & 
            prevs %in% colnames(predictor.mat)
        prevs <- prevs[tmp]
        curs <- curs[tmp]
        del.ts <- del.ts[tmp]
        out.tmp <- ((tau/del.ts) * (out.tmp[curs] - out.tmp[prevs])) + 
            out.tmp[prevs]
        in.tmp <- predictor.mat[, prevs]
        colnames(in.tmp) <- names(out.tmp)
    }
    else {
        if (!quiet) 
            cat("Time series data NOT supplied, using Tau = 0.\n")
    }
    out.tmp[out.tmp > ratio.cutoff] <- ratio.cutoff
    out.tmp[out.tmp < -ratio.cutoff] <- -ratio.cutoff
    in.tmp[in.tmp > ratio.cutoff] <- ratio.cutoff
    in.tmp[in.tmp < -ratio.cutoff] <- -ratio.cutoff
    if (conds.use[1] == "ALL") 
        conds.use <- names(out.tmp)
    df.tmp <- t(in.tmp[, names(out.tmp) %in% conds.use])
    df.tmp <- df.tmp[, !is.na(apply(df.tmp, 2, var, use = "pair")) & 
        apply(df.tmp, 2, var, use = "pair") > 0.01]
    output <- as.numeric(out.tmp[names(out.tmp) %in% conds.use])
    names(output) <- names(out.tmp)[names(out.tmp) %in% conds.use]
    df.tmp[is.na(df.tmp)] <- 0
    list(inp = df.tmp, outp = output)
}
#' Get the predictor matrix for inferelator
#' 
#' @param  predictors  The genes to use as predictors
#' @param  data  The input matrix: rows are gene names, cols are experiments
#' @param  gene.prefix  A prefix for gene names (DEFAULT: "VNG")
#' @param  preclust.k  The number of clusters to precluster the predictors into (DEFAULT: 30)
#' @param  funcs  Used for producing "and" matrix. For "make.combined.predictors" (DEFAULT: "min")
#' @param  quiet  Set to TRUE to supress output (DEFAULT: FALSE)
#' @param ...   
#' @return  list(predictor.mat, predictor.mat.ands, genetic.names, tf.groups, env.names)
#' @usage predictors <- get.predictor.matrices(predictors, data, gene.prefix = "VNG", preclust.k = 30, funcs = "min", quiet = F)
#' @export
get.predictor.matrices <-
function (predictors, data, gene.prefix = "VNG", preclust.k = 30, 
    funcs = "min", quiet = F, ...) 
{
    if (!quiet) 
        cat("Computing predictor matrices...\n")
    predictors <- predictors[predictors %in% rownames(data)]
    predictors.genetic <- grep(paste("^", gene.prefix, sep = ""), 
        predictors, value = T)
    env.names <- setdiff(predictors, predictors.genetic)
    if (!is.na(preclust.k) && preclust.k != 0 && preclust.k < 
        length(predictors)) {
        tmp <- preclust.tfs.kmeans(data = data, tfs = predictors.genetic, 
            clust.count = preclust.k, ...)
        predictor.mat <- tmp$result
        tf.groups <- tmp$tf.groups
        rm(tmp)
    }
    else {
        predictor.mat <- data[predictors, ]
        rownames(predictor.mat) <- predictors
        tf.groups <- as.list(predictors)
        names(tf.groups) <- rownames(predictor.mat)
    }
    if (length(env.names) > 0 && any(!env.names %in% rownames(predictor.mat))) 
        predictor.mat <- rbind(predictor.mat, data[env.names[!env.names %in% 
            rownames(predictor.mat)], ])
    tmp <- unique(predictor.mat)
    if (any(!rownames(tmp) %in% rownames(predictor.mat))) 
        cat("Predictor", rownames(tmp)[!rownames(tmp) %in% rownames(predictor.mat)], 
            "is not unique. Removing.\n")
    predictor.mat <- tmp
    rm(tmp)
    predictor.mat.ands <- NULL
    if (!is.na(funcs) && length(funcs) > 0) {
        if (!quiet) 
            cat("Computing combined predictor matrix...\n")
        predictor.mat.ands <- make.combined.predictors(predictor.mat, 
            funcs = funcs, ...)
    }
    list(predictor.mat = predictor.mat, predictor.mat.ands = unique(predictor.mat.ands), 
        genetic.names = unique(predictors.genetic), tf.groups = tf.groups, 
        env.names = unique(env.names))
}
inferelate.one.cluster <-
function (cluster, predictors, data, col.map = NULL, conds.use = c("clust", 
    "ALL")[1], quiet = F, plot = T, shrink.opt = c("glmnet", 
    "lars")[1], predictor.mats = NULL, weighted = T, exp.weights = NULL,...) 
{
    if (!exists("predictor.mats") || is.null(predictor.mats)) 
        predictor.mats <<- get.predictor.matrices(predictors, 
            data, quiet = quiet, ...)
    cluster.rows <- cluster$rows
    cluster.conds <- colnames(data)
    if (conds.use == "clust") {
        if (exists("col.map") && !is.null(col.map)) 
            cluster.conds <- fill.in.time.series(cluster$cols, 
                col.map, fill.all.ts.frac = 1.25, remove.all.ts.frac = 0.1, 
                fill.ts.gap.size = 1)
        else cluster.conds <- cluster$cols
    }
    cluster.conds <- cluster.conds[cluster.conds %in% colnames(data)]
    cluster.profile <- apply(data[cluster.rows, , drop = F], 
        2, mean, na.rm = T)
    cluster.weights <- NA
    if (weighted) {
        cluster.vars <- apply(data[cluster.rows, , drop = F], 
            2, var, na.rm = T)
        cluster.vars <- cluster.vars/(abs(cluster.profile) + 
            0.05)
        cluster.vars[is.na(cluster.vars) | cluster.vars == 0] <- 1
        cluster.weights <- 1/cluster.vars

    }

    #SD 10/11/10  take into account any input weights from exp.weights
    if (! is.null(exp.weights) ) {
        if(!is.na(cluster.weights)[1]) { 
        	cluster.weights<-exp.weights*cluster.weights 
        } else {
        	cluster.weights<-exp.weights
        }
    }

    
    tmp <- get.cluster.predictors(cluster.rows, cluster.profile[cluster.conds], 
        predictor.mats$predictor.mat[, cluster.conds], predictor.mats$predictor.mat.ands[, 
            cluster.conds], predictor.mats$tf.groups, predictor.mats$env.names, 
        quiet = quiet, ...)
    #SD 9/28/10  Add in any value that has a prior > .5
    pGTE.5<-names(cluster$priors)[cluster$priors>=.5]
    possibleNames<-c(rownames(predictor.mats$predictor.mat.ands),rownames(predictor.mats$predictor.mat))
    tmp$predictors<-c(tmp$predictors,pGTE.5[pGTE.5 %in% possibleNames])
       
    #SD 10/11/10 Remove any values from "possibly regulates" that has been added because of a prior
    tmp$possibly.regulates<-tmp$possibly.regulates[!(tmp$possibly.regulates %in% tmp$predictors)]
        
    possibly.regulates <- tmp$possibly.regulates
    predictor.mat <- rbind(predictor.mats$predictor.mat, predictor.mats$predictor.mat.ands)
    predictor.mat <- predictor.mat[rownames(predictor.mat) %in% 
        tmp$predictors, , drop = F]
    predictor.mat <- mean.variance.normalize(predictor.mat, filter = NA)
    if (!quiet) 
        cat("Inferelating on biclust #", cluster$k, "using", 
            length(cluster.conds), "conditions and", nrow(predictor.mat), 
            "predictors:\n", paste(rownames(predictor.mat), sep = ", "), 
            "\n")
    coeffs<-NULL #SD 02-11-2011, If shrinkage fails for any reason, don't have any coeffs
    if (shrink.opt == "glmnet") {
    	#SD 9/28/10 to include possible priors
    	penalties<-NULL
    	if ( ! is.null(cluster$priors) ) {penalties<-1-cluster$priors+.5}
    	
        try(coeffs <- inferelator.enet(cluster.profile, predictor.mat, 
            cluster.conds, col.map = col.map, quiet = quiet, 
            weights = cluster.weights, penalties = penalties, ...))
    }
    else if (shrink.opt == "lars") {
        try(coeffs <- inferelator(cluster.profile, predictor.mat, 
            cluster.conds, col.map = col.map, quiet = quiet, 
            ...))
    }
    singleresult <- coeffs$coeffs
    coeffs.boot <- coeffs$coeffs.boot
    coef.quantiles <- coeffs$coef.quantiles
    all.inputs <- coeffs$all.inputs
    coeffs$coeffs <- coeffs$coeffs.boot <- coeffs$coef.quantiles <- coeffs$all.inputs <- NULL
    coeffs$main <- paste("Bicluster", cluster$k, cluster$nrows, 
        "genes")
    if (length(singleresult) > 0) {
        conds <- c(cluster.conds, colnames(data)[!colnames(data) %in% 
            cluster.conds])
        coeffs$predictor.mat <- predictor.mat[names(singleresult), 
            conds, drop = F]
        coeffs$colors <- c("red", ifelse(singleresult > 0, "#ffaaaa", 
            "#aaffaa"))
        return(list(k = cluster$k, coeffs = singleresult, possibly.regulates = possibly.regulates, 
            cluster.conds = cluster.conds, coeffs.boot = coeffs.boot, 
            coef.quantiles = coef.quantiles, all.inputs = all.inputs, 
            plot.info = coeffs))
    }
    else {
        return(list(k = cluster$k, coeffs = numeric(), possibly.regulates = possibly.regulates, 
            cluster.conds = cluster.conds, coeffs.boot = coeffs.boot, 
            coef.quantiles = coef.quantiles, all.inputs = all.inputs, 
            plot.info = coeffs))
    }
}
inferelator <-
function (profile, predictor.mat, conds.use, col.map = NULL, 
    tau = 10, ratio.cutoff = 3, coef.cutoff = 0.02, cv.k = 10, 
    cv.choose = "min+2se", n.boot = 1, boot.opt = c("resample", 
        "cv"), rescale.coeffs = T, quiet = T, max.coeffs = NA, 
    min.coeffs = NA, ...) 
{
    if (cv.choose == "min") 
        cv.choose <- "min+0se"
    tmp <- get.input.matrix(profile, predictor.mat, conds.use, 
        col.map = col.map, tau = tau, ratio.cutoff = ratio.cutoff, 
        quiet = quiet)
    df.tmp <- tmp$inp
    output <- tmp$outp
    rm(tmp)
    apply.func <- get.apply.func()
    out.coe <- apply.func(1:n.boot, function(boot) {
        cols <- 1:length(output)
        if (boot > 1 && boot.opt == "resample") 
            cols <- sample(cols, replace = T)
        lars.obj <- try(lars(df.tmp[cols, ], output[cols], type = "lasso", 
            trace = F), silent = quiet)
        if (class(lars.obj) == "try-error") {
            tries <- 1
            while (tries <= 20 && class(lars.obj) == "try-error") {
                lars.obj <- try(lars(df.tmp[cols, ], output[cols], 
                  type = "lasso", trace = F), silent = quiet)
                tries <- tries + 1
            }
        }
        if (class(lars.obj) == "try-error") 
            return(numeric())
        cv.lars.obj <- try(cv.lars(df.tmp[cols, ], output[cols], 
            K = cv.k, type = "lasso", plot.it = F, trace = F), 
            silent = quiet)
        if (class(cv.lars.obj) == "try-error") {
            tries <- 1
            while (tries <= 20 && class(cv.lars.obj) == "try-error") {
                cv.lars.obj <- try(cv.lars(df.tmp[cols, ], output[cols], 
                  K = cv.k, type = "lasso", plot.it = F, trace = F), 
                  silent = quiet)
                tries <- tries + 1
            }
        }
        if (class(cv.lars.obj) == "try-error") 
            return(numeric())
        min.i <- which.min(cv.lars.obj$cv)
        min.err <- cv.lars.obj$cv.error[min.i]
        if (grepl("+", cv.choose[1], fixed = T)) {
            se <- as.numeric(gsub("min+", "", gsub("se", "", 
                cv.choose[1])))
            best.s <- min(which(cv.lars.obj$cv <= min(cv.lars.obj$cv) + 
                se * min.err))
        }
        else best.s <- which.min(cv.lars.obj$cv)
        orig.coeffs <- coeffs <- coef.lars(lars.obj, s = cv.lars.obj$fraction[best.s], 
            mode = "fraction")
        sorted <- names(sort(abs(coeffs), decreasing = T))
        coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
        if (!is.na(min.coeffs) && length(coeffs) < min.coeffs) 
            coeffs <- orig.coeffs[sorted[1:min.coeffs]]
        if (!is.na(max.coeffs) && length(coeffs) > max.coeffs) 
            coeffs <- orig.coeffs[sorted[1:max.coeffs]]
        if (!quiet) 
            cat(boot, cv.choose[1], min.i, min.err, best.s, cv.lars.obj$cv[best.s], 
                cv.lars.obj$fraction[best.s], length(coeffs), 
                "\n")
        if (rescale.coeffs && length(coeffs) > 0) {
            ins <- df.tmp[, names(coeffs), drop = F]
            coeffs.s <- coef(lm(output[cols] ~ ins[cols, ] - 
                1))
            names(coeffs.s) <- names(coeffs)
            coeffs <- coeffs.s[abs(coeffs.s) >= coef.cutoff]
        }
        if (boot == 1) {
            out <- list(coeffs = coeffs, lars.obj = lars.obj, 
                cv.lars.obj = cv.lars.obj, best.s = best.s, se = se, 
                min.err = min.err)
            return(out)
        }
        coeffs
    })
    lars.obj <- out.coe[[1]]$lars.obj
    cv.lars.obj <- out.coe[[1]]$cv.lars.obj
    best.s <- out.coe[[1]]$best.s
    se <- out.coe[[1]]$se
    min.err <- out.coe[[1]]$min.err
    out.coe[[1]] <- out.coe[[1]]$coeffs
    coeffs <- out.coe[[1]]
    coeffs <- coeffs[order(abs(coeffs), decreasing = T)]
    coef.quantiles <- NULL
    if (n.boot > 1) {
        tmp <- unlist(out.coe)
        tmp2 <- table(names(tmp))
        coef.quantiles <- t(sapply(names(tmp2), function(i) {
            tmp3 <- tmp[names(tmp) == i]
            tmp3 <- c(tmp3, rep(0, n.boot - length(tmp3)))
            c(n = sum(names(tmp) == i)/n.boot, quantile(abs(tmp3), 
                prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) * 
                sign(mean(tmp3[tmp3 != 0], na.rm = T)))
        }))
        coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
            1, function(i) all(i[-1] == 0)), ]
        if (!quiet) 
            print(coef.quantiles, digits = 3)
    }
    if (!quiet) {
        cat("NONZERO COEFFS:\n")
        print(coeffs)
    }
    return(list(coeffs = coeffs, coeffs.boot = out.coe, coef.quantiles = coef.quantiles, 
        lars.obj = lars.obj, cv.lars.obj = cv.lars.obj, cv.choose = cv.choose, 
        best.s = best.s, se = se, min.err = min.err, all.inputs = rownames(df.tmp)))
}
inferelator.enet <-
function (profile, predictor.mat, conds.use, col.map = NULL, 
    tau = 10, ratio.cutoff = 3, coef.cutoff = 0.02, cv.k = 10, 
    cv.choose = "min+2se", n.boot = 1, boot.opt = c("resample", 
        "cv"), rescale.coeffs = T, quiet = T, alpha = 0.9, weights = NA, 
    penalties = NA, max.coeffs = NA, min.coeffs = NA, ...) 
{
    if (cv.choose == "min") 
        cv.choose <- "min+0se"
    tmp <- get.input.matrix(profile, predictor.mat, conds.use, 
        col.map = col.map, tau = tau, ratio.cutoff = ratio.cutoff, 
        quiet = quiet)

    #Filter out NA experiments to avoid shrinkage errors
    #  SD, 9/22/10
    df.tmp <- tmp$inp[!is.na(tmp$outp),]
    output <- tmp$outp[!is.na(tmp$outp)]
    
    rm(tmp)
    in.penalties <- rep(1, ncol(df.tmp))
    names(in.penalties) <- colnames(df.tmp)
    if (!is.na(penalties) && any(names(penalties) %in% names(in.penalties))) {
        penalties <- penalties[names(penalties) %in% names(in.penalties)]
        in.penalties[names(penalties)] <- penalties
    }
    if (is.na(weights[1])) 
        weights <- rep(1, length(output))
    else weights <- weights[names(output)]
    names(weights) <- names(output)
    apply.func <- get.apply.func()
    if (!quiet) 
        cat("Alpha =", alpha, "\n")
    out.coe <- apply.func(1:n.boot, function(boot) {
        cols <- 1:length(output)
        if (boot > 1 && boot.opt == "resample") 
            cols <- sample(cols, replace = T)
        glmnet.obj <- my.glmnet(df.tmp[cols, ], output[cols], 
            penalty.factor = in.penalties, weights = weights[cols], 
            alpha = if (alpha == "cv.choose") 
                0
            else alpha, ...)
        if ("try-error" %in% class(glmnet.obj)) {
            tries <- 1
            while (tries <= 20 && "try-error" %in% class(glmnet.obj)) {
                glmnet.obj <- try(my.glmnet(df.tmp[cols, ], output[cols], 
                  penalty.factor = in.penalties[cols], weights = weights[cols], 
                  alpha = if (alpha == "cv.choose") 
                    0
                  else alpha, ...), silent = quiet, priors)
                tries <- tries + 1
            }
        }
        if ("try-error" %in% class(glmnet.obj)) 
            return(numeric())
        cv.glmnet.obj <- try(cv.glmnet(df.tmp[cols, ], output[cols], 
            lambda = glmnet.obj$lambda, K = cv.k, trace = F, 
            penalty = in.penalties, weights = weights[cols], 
            alpha = if (alpha == "cv.choose") 
                1
            else alpha, plot.it = F), silent = quiet)
        if ("try-error" %in% class(cv.glmnet.obj)) {
            tries <- 1
            while (tries <= 20 && "try-error" %in% class(cv.glmnet.obj)) {
                cv.glmnet.obj <- try(cv.glmnet(df.tmp[cols, ], 
                  output[cols], lambda = glmnet.obj$lambda, K = cv.k, 
                  trace = F, penalty = in.penalties, weights = weights[cols], 
                  alpha = if (alpha == "cv.choose") 
                    1
                  else alpha, plot.it = F), silent = quiet)
                tries <- tries + 1
            }
        }
        if ("try-error" %in% class(cv.glmnet.obj)) 
            return(numeric())
        cv.glmnet.obj$alpha <- alpha
        min.i <- which.min(cv.glmnet.obj$cv)
        min.err <- cv.glmnet.obj$cv.error[min.i]
        se <- 1
        if (grepl("+", cv.choose[1], fixed = T)) {
            se <- as.numeric(gsub("min+", "", gsub("se", "", 
                cv.choose[1])))
            best.s <- min(which(cv.glmnet.obj$cv <= min(cv.glmnet.obj$cv) + 
                se * min.err))
        }
        else best.s <- which.min(cv.glmnet.obj$cv)
        coeffs.tmp <- as.matrix(coef(glmnet.obj, s = glmnet.obj$lambda[best.s]))
        coeffs <- coeffs.tmp[coeffs.tmp != 0, ]
        names(coeffs) <- rownames(coeffs.tmp)[coeffs.tmp != 0]
        orig.coeffs <- coeffs <- coeffs[names(coeffs) != "(Intercept)"]
        sorted <- names(sort(abs(coeffs), decreasing = T))
        coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
        if (!is.na(min.coeffs) && length(coeffs) < min.coeffs && 
            length(orig.coeffs) >= min.coeffs) 
            coeffs <- orig.coeffs[sorted[1:min.coeffs]]
        if (!is.na(max.coeffs) && length(coeffs) > max.coeffs) 
            coeffs <- orig.coeffs[sorted[1:max.coeffs]]
        if (!quiet) 
            cat(boot, cv.choose[1], min.i, min.err, best.s, cv.glmnet.obj$cv[best.s], 
                glmnet.obj$lambda[best.s], length(coeffs), "\n")
        if (rescale.coeffs && length(coeffs) > 0) {
            ins <- df.tmp[, names(coeffs), drop = F]
            glmnet.obj2 <- my.glmnet(ins[cols, , drop = F], output[cols], 
                alpha = alpha, penalty.factor = in.penalties, weights = weights[cols], 
                ...)
            coeffs.s <- t(as.matrix(coef(glmnet.obj2)))
            coeffs <- coeffs.s[nrow(coeffs.s), ]
            coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
            coeffs <- coeffs[names(coeffs) != "(Intercept)"]
        }
        if (boot == 1) {
            out <- list(coeffs = coeffs, lars.obj = glmnet.obj, 
                cv.lars.obj = cv.glmnet.obj, best.s = best.s, 
                se = se, min.err = min.err)
            return(out)
        }
        coeffs
    })
    lars.obj <- out.coe[[1]]$lars.obj
    cv.lars.obj <- out.coe[[1]]$cv.lars.obj
    best.s <- out.coe[[1]]$best.s
    se <- out.coe[[1]]$se
    min.err <- out.coe[[1]]$min.err
    out.coe[[1]] <- out.coe[[1]]$coeffs
    coeffs <- out.coe[[1]]
    coeffs <- coeffs[order(abs(coeffs), decreasing = T)]
    coef.quantiles <- NULL
    if (n.boot > 1) {
        tmp <- unlist(out.coe)
        tmp2 <- sort(table(names(tmp)), decreasing = T)
        coef.quantiles <- t(sapply(names(tmp2), function(i) {
            tmp3 <- tmp[names(tmp) == i]
            tmp3 <- c(tmp3, rep(0, n.boot - length(tmp3)))
            c(n = sum(names(tmp) == i)/n.boot, quantile(abs(tmp3), 
                prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) * 
                sign(mean(tmp3[tmp3 != 0], na.rm = T)))
        }))
        coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
            1, function(i) all(i[-1] == 0)), ]
        if (!quiet) 
            print(coef.quantiles, digits = 3)
    }
    if (!quiet) {
        cat("NONZERO COEFFS:\n")
        print(coeffs)
    }
    return(list(coeffs = coeffs, coeffs.boot = out.coe, coef.quantiles = coef.quantiles, 
        lars.obj = lars.obj, cv.lars.obj = cv.lars.obj, cv.choose = cv.choose, 
        best.s = best.s, se = se, min.err = min.err, all.inputs = rownames(df.tmp)))
}
load.egrin.data <-
function (path = ".", ...) 
{
    load(paste(path, "data_orig_EGRIN/egrin_newcode_workspace.RData", 
        sep = "/"))
    load(paste(path, "data_orig_EGRIN/env_map_egrin.RData", sep = "/"))
    relevant.env <- c("oxygen", "illumination", "Fe", "Cu", "Co", 
        "Mn", "Zn", "Ni", "gamma", "uv")
    env.map <- env.map[relevant.env, ]
    load(paste(path, "data_orig_EGRIN/col_map_egrin.RData", sep = "/"))
    for (i in 1:(length(colMap) - 1)) {
        col.map <- rbind(col.map, as.data.frame(colMap[[i]]))
    }
    rownames(col.map) <- names(colMap)[1:(length(colMap) - 1)]
    colnames(col.map)[colnames(col.map) == "del.t"] <- "delta.t"
    pc <- as.character(col.map$prevCol)
    pc[is.na(pc)] <- as.character(col.map$condName[is.na(pc)])
    col.map$prevCol <- as.factor(pc)
    col.map$delta.t[is.na(col.map$delta.t)] <- 9999
    predictors <- c(readLines(paste(path, "data/halo/halo_tfs.txt", 
        sep = "/")), rownames(env.map))
    data <- rbind(ratios.egrin, env.map)
    load(paste(path, "data_orig_EGRIN/egrin_coeffs.RData", sep = "/"))
    invisible(list(col.map = col.map, env.map = env.map, predictors = predictors, 
        data = data, clusterStack.egrin = clusterStack.egrin))
}
make.combined.predictors <-
function (predictor.mat, predictors = rownames(predictor.mat), 
    funcs = "min", r.filter = 0.8, ...) 
{
    if (is.null(funcs) || is.na(funcs)) 
        return(predictor.mat)
    result <- NULL
    tmp <- t(combn(predictors, 2))
    tmp <- tmp[tmp[, 1] != tmp[, 2], ]
    apply.func <- get.apply.func()
    for (func in funcs) {
        tmp2 <- do.call(rbind, apply.func(1:nrow(tmp), function(i) apply(predictor.mat[tmp[i, 
            ], ], 2, func, na.rm = T)))
        tmp2[is.infinite(tmp2)] <- NA
        rownames(tmp2) <- paste(tmp[, 1], tmp[, 2], rep(func, 
            nrow(tmp)), sep = combine.symbol)
        result <- rbind(result, tmp2)
        rm(tmp2)
    }
    cat("Combined", funcs, "predictor matrix is", nrow(result), 
        "x", ncol(result), "\n")
    out <- NULL
    if (!is.na(r.filter) && r.filter > 0 && r.filter < 1) {
        apply.func <- get.apply.func()
        all.cors <- cor(t(predictor.mat), t(result), use = "pairwise")
        tmp <- apply.func(1:nrow(result), function(i) {
            nm <- strsplit(rownames(result)[i], combine.symbol, 
                fixed = T)[[1]]
            tmp.out <- NULL
            ttmp <- all.cors[nm[1:2], i]
            ttmp[is.na(ttmp)] <- 0
            if (!any(ttmp > r.filter)) 
                tmp.out <- rownames(result)[i]
            tmp.out
        })
        tmp <- do.call("c", tmp)
        out <- result[tmp, ]
        cat("Filtered for cor <=", r.filter, ", combined predictor matrix is now ", 
            nrow(out), "x", ncol(out), "\n")
    }
    if(r.filter >= 1) { out<- result }
    attr(out, "r.filter") <- r.filter
    return(out)
}

#' Turn a matrix into a matrix of Z-scores.
#' Means and SDs are calculated for each row (gene)
#' 
#' @param input.matrix  The input matrix
#' @param filter The minimum allowable fraction of non-NA data points (DEFAULT: 0.04)
#' @param byCols Set to true to calculate Z-scores for columns (i.e. experiments).  May be better for predictors matrix  (DEFAULT: FALSE)
#' @return  A matrix of Z scores
#' @usage data <- mean.variance.normalize(input.matrix, filter = 0.04, byCols = F)
#' @export
mean.variance.normalize <-
function (input.matrix, filter = 0.04, byCols = F) 
{
    if (!is.na(filter)) {
        which.good <- which(apply(input.matrix, 1, function(i) mean(!is.na(i), 
            na.rm = T)) >= filter)
        input.matrix <- input.matrix[which.good, ]
    }
    
    if (byCols) { 
    	input.matrix <- apply(input.matrix, 2, function(x) { (x-mean(x,na.rm=T))/(sd(x,na.rm=T)) })
    } else { 
	means <- apply(input.matrix, 1, mean, na.rm = T)
    	sds <- apply(input.matrix, 1, sd, na.rm = T)
	sds[sds == 0 | is.na(sds)] <- 1
	input.matrix <- apply(input.matrix, 2, "-", means)
	input.matrix <- apply(input.matrix, 2, "/", sds)
    }
    
    return(input.matrix)
}
my.glmnet <-
function (x, y, family = c("gaussian", "binomial", "poisson", 
    "multinomial", "cox"), weights, offset = NULL, alpha = 1, 
    nlambda = 100, lambda.min = 1e-06, lambda = NULL, standardize = TRUE, 
    thresh = 1e-04, dfmax = ncol(x) + 1, pmax = min(dfmax * 1.2, 
        ncol(x)), exclude, penalty.factor = rep(1, ncol(x)), 
    maxit = 100, HessianExact = FALSE, type = c("covariance", 
        "naive"), ...)  {
        #browser()
  glmnet(x, y, family, weights, offset, alpha, nlambda, lambda.min, 
      lambda, standardize, thresh, dfmax, pmax, exclude, penalty.factor, 
      maxit, HessianExact, type)
}
nwInf.package <-
function (install = T, update.web = F, check = F, version = "0.0.5") 
{
    #source.files <- c("runnit.R", "inferelator.R", "inferelator_enet.R","predictelator.R", "write.inf.network.R")
    if (exists("envMap") && !is.null(envMap)) {
        cat("Packaging Halo data...\n")
        halo <- list(ratios = ratios, envMap = envMap, colMap = colMap, 
            tfs = halo_tfs)
        halo <<- halo
    }
    onLoad <- function(libname, pkgname) {
        cat("Loading ", pkgname, " version ", VERSION, " (", 
            DATE, ")\n", sep = "")
        cat("Copyright (C) David J Reiss, Institute for Systems Biology; dreiss.isb@gmail.com.org.\n")
        cat("http://baliga.systemsbiology.net/cmonkey\n")
        cat("\nNOTE that this package is still sloppy in that it relies upon some global variables:\n")
        cat("'predictor.mats', 'envMap', 'colMap', and optionally 'predictors'.\n")
    }
    source("~/scratch/halo/generic_scripts/construct.package.R")
    construct.package("cMonkeyNwInf", version = version, source.files = source.files, 
        functions.visible = c("runnit.wrapper.halo","runnit.wrapper","runnit.wrapper.sce",
            "plot.coeff.obj","write.cytoscape.files"), data = if (exists("halo") && 
            !is.null(halo)) 
            list(halo = "halo")
        else NULL, required = c("lars", "glmnet", "multicore", 
            "Matrix"), suggested = c("cMonkey", "foreach", "doMC", 
            "ff", "igraph"), short.desc = "Inferelator-like network inference on cMonkey biclusters", 
        long.desc = "Inferelator-like network inference on cMonkey biclusters", 
        onLoad = onLoad)
    if (install) 
        system(sprintf("R CMD INSTALL lib/DR.Inferelator_%s.tar.gz", 
            version))
    if (check) {
        cwd <- setwd("lib")
        system(sprintf("R CMD CHECK DR.Inferelator_%s.tar.gz", 
            version))
        setwd(cwd)
    }
    if (update.web) {
        system(sprintf("scp lib/DR.Inferelator_%s.tar.gz bragi:/local/apache2/htdocs/cmonkey/DR.Inferelator_latest.tar.gz", 
            version))
    }
}
plot.cluster.coeffs <-
function (coefs, scale = 1, cex = 0.5, ...) 
{
    require(igraph)
    network <- data.frame()
    comb.cnt <- 1
    node.types <- character()
    for (coe in coefs) {
        if (length(coe$coeffs) <= 0) {
            network <- rbind(network, data.frame(n1 = sprintf("bic%s", 
                coe$k), n2 = sprintf("bic%s", coe$k), weight = NA, 
                mode = "-"))
        }
        else {
            for (i in 1:length(coe$coeffs)) {
                n <- strsplit(names(coe$coeffs)[i], combine.symbol, 
                  fixed = T)[[1]]
                if (length(n) == 1) {
                  network <- rbind(network, data.frame(n1 = n, 
                    n2 = sprintf("bic%s", coe$k), weight = coe$coeffs[i], 
                    mode = ">"))
                }
                else {
                  n2 <- paste("AND", comb.cnt, sep = "")
                  network <- rbind(network, data.frame(n1 = n2, 
                    n2 = sprintf("bic%s", coe$k), weight = coe$coeffs[i], 
                    mode = ">"))
                  network <- rbind(network, data.frame(n1 = n[1], 
                    n2 = n2, weight = 0, mode = "-"))
                  network <- rbind(network, data.frame(n1 = n[2], 
                    n2 = n2, weight = 0, mode = "-"))
                  comb.cnt <- comb.cnt + 1
                }
            }
        }
        if (!is.null(coe$possibly.regulates) && length(coe$possibly.regulates) > 
            0) {
            for (i in 1:length(coe$possibly.regulates)) {
                network <- rbind(network, data.frame(n1 = names(coe$possibly.regulates)[i], 
                  n2 = sprintf("bic%s", coe$k), weight = 0, mode = "*"))
            }
        }
    }
    gr <- graph.edgelist(as.matrix(network[, 1:2]), directed = T)
    gr.layout <- layout.fruchterman.reingold.grid(gr, niter = 3000 * 
        length(coefs)^2, coolexp = 0.5, ...)
    gr.layout <- layout.norm(gr.layout, -1, 1, -1, 1)
    node.names <- get.vertex.attribute(gr, "name")
    node.sizes <- rep(15, length(node.names))
    names(node.sizes) <- node.names
    node.sizes[grepl("^bic", node.names)] <- 25
    node.sizes[grepl("^AND", node.names)] <- 10
    node.sizes <- node.sizes * scale/length(coefs)
    node.colors <- rep("lightgreen", length(node.names))
    names(node.colors) <- node.names
    node.colors[grepl("^bic", node.names)] <- "lightblue"
    node.colors[grepl("^AND", node.names)] <- "gray"
    node.frame.colors <- rep("black", length(node.names))
    names(node.frame.colors) <- node.names
    node.frame.colors[!node.names %in% names(predictor.mats$tf.groups)] <- "red"
    node.frame.colors[grepl("^bic", node.names)] <- "blue"
    node.frame.colors[grepl("^AND", node.names)] <- "gray"
    node.shapes <- rep("circle", length(node.names))
    names(node.shapes) <- node.names
    node.shapes[grepl("^bic", node.names)] <- "square"
    node.names[grepl("^AND", node.names)] <- ""
    node.names <- gsub("TFGROUP", "tf", node.names)
    edge.colors <- ifelse(is.na(network$weight), "white", ifelse(network$weight > 
        0, "red", ifelse(network$weight < 0, "green", "blue")))
    edge.colors[as.character(network$mode) == "*"] <- "black"
    edge.widths <- abs(network$weight) * 6 + 0.25
    edge.widths[is.na(edge.widths)] <- 0.25
    edge.widths[as.character(network$mode) == "*"] <- 0.25
    tmp <- as.character(network$mode)
    tmp[tmp == "*"] <- "-"
    network.mode <- as.factor(tmp)
    plot(gr, layout = gr.layout, axes = F, margin = 0, rescale = F, 
        vertex.label = node.names, vertex.size = node.sizes, 
        vertex.color = node.colors, vertex.shape = node.shapes, 
        edge.arrow.size = 0.5, vertex.frame.color = node.frame.colors, 
        vertex.label.cex = cex, edge.color = edge.colors, edge.width = edge.widths, 
        edge.arrow.mode = as.character(network$mode))
    invisible(cbind(network, edge.colors, edge.widths))
}
plot.coeff.obj <-
function (coeffs, do.scattersmooth = F, ...) 
{
    layout(matrix(c(1, 1, 1, 2, 2, 2, 4, 4, 1, 1, 1, 2, 2, 2, 
        4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5, 
        3, 3, 3, 3, 5, 5, 5, 5), nrow = 5, ncol = 8, byrow = T))
    my.plotCVLars <- function(cv.lars.object, se = TRUE, ...) {
        attach(cv.lars.object)
        plot(fraction, cv, type = "b", ylim = range(cv, cv + 
            cv.error, cv - cv.error), ...)
        if (se) 
            error.bars(fraction, cv + cv.error, cv - cv.error, 
                width = 1/length(fraction))
        detach(cv.lars.object)
        invisible()
    }
    if (!is.null(coeffs$n.boot)) 
        n.boot <- 1
    else if (!is.null(coeffs[[1]]$n.boot)) {
        n.boot <- coeffs[[1]]$n.boot
        coeb <- coeffs
        coeffs <- coeffs[[1]]
    }
    pi <- coeffs$plot.info
    require(lars)
    require(glmnet)
    if ("glmnet" %in% class(pi$lars.obj)) 
        plot(pi$lars.obj, "lambda")
    else plot(pi$lars.obj)
    lines(rep(pi$cv.lars.obj$fraction[pi$best.s], 2), c(-999, 
        999), col = 2, lty = 2, lwd = 3)
    my.plotCVLars(pi$cv.lars.obj, se = TRUE, main = class(pi$lars.obj)[1])
    if ("glmnet" %in% class(pi$lars.obj)) 
        legend("bottomleft", pi$cv.choose)
    else legend("topright", pi$cv.choose)
    lines(rep(pi$cv.lars.obj$fraction[pi$best.s], 2), c(-999, 
        999), col = 2, lty = 2, lwd = 3)
    if (grepl("+", pi$cv.choose, fixed = T)) {
        lines(rep(pi$cv.lars.obj$fraction[which.min(pi$cv.lars.obj$cv)], 
            2), rep(min(pi$cv.lars.obj$cv), 2) + c(0, pi$se * 
            pi$min.err), col = 2, lty = 2, lwd = 1)
        lines(c(pi$cv.lars.obj$fraction[pi$best.s], pi$cv.lars.obj$fraction[which.min(pi$cv.lars.obj$cv)]), 
            rep(min(pi$cv.lars.obj$cv), 2) + pi$se * pi$min.err, 
            col = 2, lty = 2, lwd = 1)
    }
    if (length(coeffs$coeffs) > 0) {
        matplot(t(rbind(coeffs$observed[pi$clust.conds.plot], 
            pi$predictor.mat[, pi$clust.conds.plot])), col = pi$colors, 
            ylab = "Normalized expression", xlab = "Conditions", 
            type = "l", main = pi$main)
        legend("bottomright", c("biclust", names(coeffs$coeffs)), 
            col = pi$colors, lty = 1, cex = 0.5)
        lines(pi$cluster.profile, col = "red")
    }
    else {
        plot(coeffs$observed[pi$clust.conds.plot], col = "red", 
            ylab = "Normalized expression", xlab = "Conditions", 
            type = "l", main = pi$main)
        legend("bottomright", "biclust", col = "red", lty = 1, 
            cex = 0.5)
    }
    if (!is.null(coeffs$pred.ts) && nrow(coeffs$pred.ts) > 1) {
        matlines(t(apply(coeffs$pred.ss[, pi$clust.conds.plot], 
            2, quantile, prob = c(0.1, 0.9))), col = rep("lightblue", 
            2), lty = 1, lwd = 3)
        matlines(t(apply(coeffs$pred.ts[, pi$clust.conds.plot], 
            2, quantile, prob = c(0.1, 0.9))), col = rep("gray", 
            2), lty = 1, lwd = 3)
    }
    if (n.boot > 1) {
        coeb <- coeb[sapply(coeb, length) > 0]
        pred.ts <- t(sapply(coeb, "[[", "pred.ts"))
        tmp <- t(apply(pred.ts, 2, quantile, prob = c(0.05, 0.5, 
            0.95)))
        rownames(tmp) <- colnames(coeb[[1]]$pred.ts)
        matlines(tmp[pi$clust.conds.plot, ], typ = "l", lty = 1, 
            col = c("gray", "red", "gray"), lwd = 3)
    }
    lines(coeffs$pred.ss[1, pi$clust.conds.plot], col = "blue")
    lines(coeffs$pred.ts[1, pi$clust.conds.plot], col = "black")
    lines(rep(pi$n.conds, 2), c(-999, 999), col = "gray", lty = 2, 
        lwd = 3)
    legend("bottomleft", c("pred.ss", "pred.ts"), col = c("blue", 
        "black"), lty = 1, cex = 0.5)
    lines(coeffs$observed[pi$clust.conds.plot], col = "red")
    out.net <- plot.cluster.coeffs(list(coeffs))
    if (do.scattersmooth && !is.null(coeffs$pred.ts)) {
        source("~/scratch/halo/generic_scripts/scattersmooth.R")
        scattersmooth(coeffs$observed[coeffs$cluster.conds][!is.na(coeffs$pred.ts[1, 
            coeffs$cluster.conds])], coeffs$pred.ts[1, coeffs$cluster.conds][!is.na(coeffs$pred.ts[1, 
            coeffs$cluster.conds])])
    }
    invisible(out.net)
}
preclust.tfs.kmeans <-
function (data, tfs, clust.count, n.iter = 200, n.start = 25, 
    seed = 31337, r.cutoff = 0.85, ...) 
{
    if (!is.na(seed)) 
        set.seed(seed)
    tf.matrix <- data[tfs[tfs %in% rownames(data)], ]
    
    #Set any NA values to 0
    tf.matrix[is.na(tf.matrix)]<-0
    
    data.c <- kmeans(tf.matrix, clust.count, iter.max = n.iter, 
        nstart = n.start)
    result <- data.c$centers
    rownames(result) <- paste("TFGROUP", 1:clust.count, sep = "")
    tf.groups <- lapply(1:length(data.c$size), function(i) names(which(data.c$cluster == 
        i)))
    names(tf.groups) <- rownames(result)
    cat("Preclustered with k-means, predictor matrix is", nrow(result), 
        "x", ncol(result), "\n")
    tmp <- apply(result, 1, function(i) apply(data[tfs, ], 1, 
        cor, i))
    if (any(tmp > r.cutoff)) {
        high.cors <- apply(tmp, 2, function(i) which(i > r.cutoff))
        to.be.added <- lapply(names(tf.groups), function(i) names(high.cors[[i]])[!names(high.cors[[i]]) %in% 
            tf.groups[[i]]])
        for (i in 1:length(tf.groups)) tf.groups[[i]] <- unique(c(tf.groups[[i]], 
            to.be.added[[i]]))
        result <- t(sapply(tf.groups, function(i) apply(data[i, 
            , drop = F], 2, mean, na.rm = T)))
    }
    for (i in 1:length(tf.groups)) if (length(tf.groups[[i]]) == 
        1) 
        names(tf.groups)[i] <- rownames(result)[i] <- tf.groups[[i]][1]
    return(list(result = result, tf.groups = tf.groups))
}
predictelate <-
function (cluster.rows, coeffs, ratios, predictor.mats = NULL, 
    tf.groups = NULL, col.map = NULL, tau = 10, max.coeffs = length(coeffs), 
    ...) 
{
    if (length(coeffs) <= 0) {
        out <- ratios[1, ] * 0
        return(out)
    }
    if (max.coeffs < length(coeffs)) 
        coeffs <- sort(coeffs, decreasing = T)[1:max.coeffs]
    coeff.names <- unique(unlist(strsplit(names(coeffs), combine.symbol, 
        fixed = T)))
    if (is.null(predictor.mats) && !is.null(tf.groups) && any(coeff.names %in% 
        names(tf.groups))) {
        tfgroup.ratios <- t(sapply(tf.groups[which(names(tf.groups) %in% 
            coeff.names)], function(i) apply(ratios[i, , drop = F], 
            2, mean)))
        rownames(tfgroup.ratios) <- names(tf.groups)[names(tf.groups) %in% 
            coeff.names]
        ratios <- rbind(ratios, tfgroup.ratios)
    }
    else {
        ratios <- rbind(ratios, predictor.mats$predictor.mat[, 
            colnames(ratios)], predictor.mats$predictor.mat.ands[, 
            colnames(ratios)])
    }
    out.ss <- 0
    for (j in 1:length(coeffs)) {
        if (coeffs[j] == 0) 
            next
        nodes <- unlist(strsplit(names(coeffs)[j], "~~", fixed=T))
        if (length(nodes) == 1) {
            if (!nodes[1] %in% rownames(ratios)) 
                next
            tmp <- ratios[nodes, ]
        }
        else if (length(nodes) == 3) {
            if (names(coeffs)[j] %in% rownames(ratios)) {
                tmp <- ratios[names(coeffs)[j], ]
            }
            else {
                if (!all(nodes[1:2] %in% rownames(ratios))) 
                  next
                tmp <- apply(ratios[c(nodes[1], nodes[2]), ], 
                  2, FUN = nodes[3], na.rm = T)
            }
        }
        tmp[is.na(tmp)] <- 0
        out.ss <- out.ss + tmp * coeffs[j]
    }
    out <- out.ss
    if (!is.null(col.map) && !is.na(tau) && tau > 0) {
        conds <- colnames(ratios)
        prevs <- as.character(col.map[conds, "prevCol"])
        del.ts <- as.numeric(as.character(col.map[conds, "delta.t"]))
        del.ts[del.ts < 1] <- 1
        cluster.prof <- apply(ratios[cluster.rows, , drop = F], 
            2, mean, na.rm = T)
        tmp1 <- cluster.prof[prevs]
        tmp1[is.na(tmp1)] <- 0
        tmp2 <- out.ss[conds]
        tmp2[is.na(tmp2)] <- 0
        out.ts <- (tau * tmp1 + del.ts * tmp2)/(tau + del.ts)
        names(out.ts) <- conds
        out <- out.ts
    }
    out[out > 3] <- 3
    out[out < -3] <- -3
    out
}

#' Run DR.Inferelator on some data to infer regulation of biclusters
#'   Inferelator uses regression to predict which genes cause a clusters expression level to change
#'   d(bicluster expression level)/dt = B1*G1 + B2*G2 + ... Bn*Gn
#'   Where Gx is gene expression level and Bx is a coefficient to be determined by Elastic Net shrinkage
#'
#' @param ks  Which clusters in clusterStack to infer regulation for
#' @param data  A gene expression matrix with rows as genes expression levels and cols as experiments
#' @param col.map  An Inferelator dataframe containing experimental times and other information  (http://err.bio.nyu.edu/inferelator/)
#' @param predictors  The genes to be considered as predictors
#' @usage clusterStack  A list of biclusters from cMonkey.  Includes $priors, prior probability of regulation
#' @param tau  The supposed time delay between TF mRNA expression and target mRNA expression (DEFAULT: 10)
#' @param pred.data  Data used only for making predictors, e.g. raw matrix data (DEFAULT: NULL, i.e. use "data")
#' @param filter.pred.by.col Normalization algorithm for pred.data. One of FALSE, TRUE, or NULL (DEFAULT: FALSE, i.e. filter by row, same as for data)
#' @param plot Set to TRUE to plot out shrinkage (DEFAULT: T)
#' @param coeffs  UNKNOWN, SD (DEFAULT: NULL)
#' @param tf.groups  UNKNOWN, SD (DEFAULT: 72)
#' @param n.boot  UNKNOWN, SD (DEFAULT: 1)
#' @param boot.opt  Bootstrapping algorithm one of "resample.lars", "resample.rows", "resample", or "lars" (DEFAULT: "resample.lars")
#' @param exp.weights  Used to weight certain experiments (i.e. cols in data) (DEFAULT: NULL)
#' @param ...  
#' @return Coefficients for putative regulators
#' @usage coeffs<-runnit(ks, data, col.map, predictors, clusterStack, tau = 10, pred.data = NULL, filter.pred.by.col = F, plot = T, coeffs = NULL, tf.groups = 72, n.boot = 1, boot.opt = c("resample.lars", "resample.rows", "resample", "lars")[1], exp.weights=NULL, ...)
#' @export
runnit <-
function (ks, data, col.map, predictors, clusterStack, tau = 10, pred.data = NULL, filter.pred.by.col = F,
    plot = T, coeffs = NULL, tf.groups = 72, n.boot = 1, boot.opt = c("resample.lars", 
        "resample.rows", "resample", "lars")[1], exp.weights=NULL, ...) 
{
    if (is.null(pred.data)) { pred.data <- data }    
    data <- mean.variance.normalize(data, filter = 0.04)

    if (! is.null(filter.pred.by.col)) {
        pred.data <- mean.variance.normalize(pred.data, filter = 0.04, byCols = filter.pred.by.col)
    }

    predictors <- predictors[predictors %in% rownames(pred.data)]
    if (!exists("predictor.mats") || ((is.na(tf.groups) || tf.groups == 
        0 || tf.groups >= length(predictors)) && length(predictor.mats$tf.groups) != 
        length(predictors)) || (!is.na(tf.groups) && tf.groups != 
        0 && tf.groups < length(predictors) && length(predictor.mats$tf.groups) != 
        tf.groups)) {
        predictor.mats <<- get.predictor.matrices(predictors, 
            pred.data, preclust.k = tf.groups, ...)
    }
    n.boot.lars <- 1
    boot.opt.lars <- "resample"
    if (n.boot > 1 && boot.opt %in% c("resample.lars", "lars")) {
        n.boot.lars <- n.boot
        n.boot <- 1
        if (boot.opt == "lars") 
            boot.opt.lars <- "cv"
    }
    apply.func <- get.apply.func(plot)
    if (n.boot > 1) 
        apply.func <- lapply
    out <- apply.func(ks, function(i) {
        cluster <- clusterStack[[i]]
        k <- cluster$k
        apply.func <- get.apply.func()
        if (n.boot == 1) 
            apply.func <- lapply
        out.k <- apply.func(1:n.boot, function(boot) {
            cat("***      BICLUSTER:", k, boot, "\n")
            clust <- cluster
            if (boot > 1) {
                if (boot.opt %in% c("resample", "resample.rows")) 
                  clust$rows <- sample(clust$rows, replace = T)
                if (boot.opt == "resample") 
                  clust$cols <- sample(colnames(data), length(clust$cols), 
                    replace = F)
            }
            
            #SD 1/10/11 Remove any genes that appear in the clusterStack but aren't in the data.
            #  This may occur during filtering if using only a subset of the data
  	    missingBool <- !clust$rows %in% rownames(data)
	    if ( any (missingBool) ) {
	    	warning(paste('Gene(s):',clust$rows[missingBool],'are missing for the data matrix.'))
	    	clust$rows <- clust$rows[!missingBool]
	    }		

            coeffs <- inferelate.one.cluster(clust, predictors, 
                data, predictor.mats = predictor.mats, tau = tau, 
                col.map = col.map, n.boot = n.boot.lars, boot.opt = boot.opt.lars, 
                quiet = n.boot > 1, exp.weights=exp.weights, ...)
            clust.conds <- sort(coeffs$cluster.conds)
            
            observed <- apply(data[clust$rows, , drop = F], 2, 
                mean, na.rm = T)
            apply.func <- get.apply.func()
            pred.ss <- do.call(rbind, apply.func(coeffs$coeffs.boot, 
                function(b) predictelate(clust$rows, b, data, 
                  predictor.mats = predictor.mats, tau = tau, 
                  ...)))
            pred.ts <- do.call(rbind, apply.func(coeffs$coeffs.boot, 
                function(b) predictelate(clust$rows, b, data, 
                  predictor.mats = predictor.mats, tau = tau, 
                  col.map = col.map, ...)))
            if (is.null(pred.ss)) 
                pred.ss <- t(observed * 0)
            if (is.null(pred.ts)) 
                pred.ts <- t(observed * 0)
            if ("weighted" %in% names(list(...)) && list(...)$weighted == 
                TRUE) {
                vars <- apply(data[clust$rows, , drop = F], 2, 
                  var, na.rm = T)
                vars <- vars/(abs(observed) + 0.05)
                vars[is.na(vars) | vars == 0] <- 1
                weights <- 1/vars
                weights <- weights/sum(weights) * length(weights)
            }
            else {
                weights <- rep(1, ncol(data))
                names(weights) <- colnames(data)
            }
            rmsd.ss <- sqrt(weighted.mean((pred.ss[nrow(pred.ss), 
                ] - observed)[clust.conds]^2, weights[clust.conds], 
                na.rm = T))
            rmsd.ts <- sqrt(weighted.mean((pred.ts[nrow(pred.ts), 
                ] - observed)[clust.conds]^2, weights[clust.conds], 
                na.rm = T))
            not.clust.conds <- colnames(data)[!colnames(data) %in% 
                clust.conds]
            rmsd.ts.out <- sqrt(weighted.mean((pred.ts[nrow(pred.ts), 
                ] - observed)[not.clust.conds]^2, weights[not.clust.conds], 
                na.rm = T))
            coeffs$plot.info$main <- paste("Bicluster", cluster$k, 
                cluster$nrows, "genes")
            coeffs$plot.info$clust.conds.plot <- c(clust.conds, 
                sort(colnames(data)[!colnames(data) %in% clust.conds]))
            coeffs$plot.info$n.conds <- length(clust.conds)
            if (n.boot <= 1) 
                cat(k, tau, rmsd.ss, rmsd.ts, rmsd.ts.out, "\n")
            coeffs$pred.ss <- pred.ss
            coeffs$pred.ts <- pred.ts
            coeffs$rmsd <- c(ss = rmsd.ss, ts = rmsd.ts, ts.out = rmsd.ts.out)
            coeffs$observed <- observed
            coeffs$n.boot <- n.boot
            coeffs$boot.opt <- boot.opt
            attr(coeffs, "class") <- "coeff.obj"
            if (boot > 1) 
                coeffs$plot.info <- NULL
            coeffs
        })
        names(out.k) <- paste(k, 1:n.boot, sep = ".")
        if (n.boot > 1) {
            cc.tmp <- out.k
            nb <- max(n.boot, n.boot.lars)
            cc.tmp <- cc.tmp[sapply(cc.tmp, length) > 0]
            cc <- lapply(cc.tmp, "[[", "coeffs")
            tmp <- cc
            names(tmp) <- NULL
            tmp <- unlist(tmp)
            tmp2 <- sort(table(names(tmp)), decreasing = T)
            coef.quantiles <- t(sapply(names(tmp2), function(i) {
                tmp3 <- tmp[names(tmp) == i]
                tmp3 <- c(tmp3, rep(0, nb - length(tmp3)))
                c(n = sum(names(tmp3) == i)/nb, quantile(abs(tmp3), 
                  prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) * 
                  sign(mean(tmp3[tmp3 != 0], na.rm = T)))
            }))
            coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
                1, function(i) all(i[-1] == 0)), ]
            out.k[[1]]$coef.quantiles <- coef.quantiles
        }
        else if (n.boot.lars > 1) {
            print(out.k[[1]]$coef.quantiles, digits = 3)
        }
        if (plot) {
            try(plot.coeff.obj(out.k, ...))
        }
        attr(out.k, "class") <- "coeff.obj"
        out.k
    })
    out <- do.call(c, out)
    attr(out, "CALL") <- match.call(expand.dots = T)
    attr(out, "class") <- "coeff.obj"
    invisible(out)
}
runnit.egrin.data <-
function (ks = 1:300, tau = 10, plot = T, coeffs = NULL, tf.groups = 72, 
    n.boot = 1, boot.opt = c("resample.lars", "resample.rows", 
        "resample", "lars")[1], ...) 
{
    if (!"egrin.data" %in% searchpaths()) {
        egrin.data <- load.egrin.data(...)
        attach(egrin.data)
    }
    out <- runnit(ks, data, col.map, predictors, clusterStack.egrin, 
        tau = tau, plot = plot, coeffs = coeffs, tf.groups = tf.groups, 
        n.boot = n.boot, boot.opt = boot.opt, ...)
    detach(egrin.data)
    invisible(out)
}
runnit.wrapper.halo <-
function (f, ks = "all", ...) 
{
    if (is.character(f) && file.exists(f) && (!exists("e") || 
        e$tmp.file != f)) {
        load(f, envir = .GlobalEnv)
        print(f)
        assign("tmp.file", f, env = e)
    }
    else if (is.environment(f)) {
        e <- f
        rm(f)
    }
    ratios <- e$get.cluster.matrix()
    if (ks[1] == "all") 
        ks <- 1:e$k.clust
    if (!exists("envMap")) 
        envMap <- NULL
    if (!exists("colMap")) 
        colMap <- NULL
    if (!exists("predictors")) 
        predictors <- readLines("data/halo/halo_tfs.txt")
    data <- ratios
    if (!is.null(envMap)) {
        envMap <- envMap[, !is.na(apply(envMap, 2, var, use = "pair")) & 
            apply(envMap, 2, var, use = "pair") > 0.01, drop = F]
        envMap <- envMap[rownames(envMap) %in% colnames(ratios), 
            , drop = F]
        ratios <- ratios[, colnames(ratios) %in% rownames(envMap), 
            drop = F]
        data <- rbind(ratios, t(as.matrix(envMap)))
        predictors <- c(predictors, colnames(envMap))
    }
    if (!is.null(predictors)) 
        predictors <- predictors[predictors %in% rownames(data)]
    out <- runnit(ks, data, colMap, predictors, e$clusterStack, 
        gene.prefix = e$genome.info$gene.prefix, ...)
    invisible(out)
}

#SD 10/11/10 Run for Yeast with best options
runnit.wrapper.sce <-
function(e) 
{
    out<-runnit.wrapper(e,cv.choose="min+4se", tf.groups=999, alpha=0.8, 
		n.boot=1, tau=10, r.cutoff=2, r.filter=1.0, weighted=T, 
		aic.filter=25, plot=F)
    coeffs<-out
    save(coeffs,e,file=paste("coeffs.e.",Sys.Date(),sep=""))
    invisible(out)
}

#SD 10/11/10 Should work for Halo or Yeast
runnit.wrapper <-
function (f, ks = "all", filter.pred.by.col = T, predictors = predictors, ...) 
{
    if (is.character(f) && file.exists(f) && (!exists("e") || 
        e$tmp.file != f)) {
        load(f, envir = .GlobalEnv)
        print(f)
        assign("tmp.file", f, env = e)
    }
    else if (is.environment(f)) {
        e <- f
        rm(f)
    }
    exp.weights<-NULL #Used to weight based on conditions
    row.weights<-e$row.weights
    ratios <- get.cluster.matrix(e$ratios)
    pred.data <- ratios.raw <- get.cluster.matrix(e$ratios.raw)
    if (ks[1] == "all") 
        ks <- 1:e$k.clust
    if (!exists("envMap")) 
        envMap <- NULL
    if (!exists("colMap")) 
        colMap <- NULL
    if (!exists("predictors")) 
        predictors <- rownames(ratios.raw)
    if (!is.null(envMap)) {
        envMap <- envMap[, !is.na(apply(envMap, 2, var, use = "pair")) & 
            apply(envMap, 2, var, use = "pair") > 0.01, drop = F]
        envMap <- envMap[rownames(envMap) %in% colnames(ratios), 
            , drop = F]
        ratios <- ratios[, colnames(ratios) %in% rownames(envMap), 
            drop = F]
        ratios.raw <- ratios.raw[, colnames(ratios) %in% rownames(envMap), 
            drop = F]
        pred.data <- rbind(ratios.raw, t(as.matrix(envMap)))
        ratios <- rbind(ratios, t(as.matrix(envMap)))
        predictors <- c(predictors, colnames(envMap))
        
        #SD 10/11/10  Construct the "weights" matrix for priors
        if (!is.null(row.weights)) {
        	exp.weights<-rep(0,ncol(pred.data))
        	orderIdx<-match(colnames(pred.data),rownames(envMap))
        	for (cond in names(row.weights)) {
        		exp.weights<-exp.weights+row.weights[cond]*envMap[orderIdx,cond]/max(envMap[orderIdx,cond])
        	}
        	names(exp.weights)<-colnames(pred.data)
        	
        	#SD 02/10/11  Remove anything from predictors that has zero weight
        	if (any(exp.weights == 0))  {
			keepBool <- colnames(ratios) %in% names(exp.weights[exp.weights != 0])
			ratios <- ratios [,keepBool]
			keepBool <- colnames(pred.data) %in% names(exp.weights[exp.weights != 0])
			pred.data <- pred.data [,keepBool]
			exp.weights <- exp.weights[exp.weights != 0]
			ratios <- ratios[,match(colnames(pred.data),colnames(ratios))]
			colMap <- colMap[rownames(colMap) %in% names(exp.weights),]
			ts.inds <- unique(colMap$ts.ind)
			new.ts.inds <- rep(0,length(colMap$ts.ind))
			for (i in 1:length(ts.inds)) {
				tsBool <- colMap$ts.ind == ts.inds[i]
				new.ts.inds[tsBool] <- i
				minExps <- rownames(colMap[tsBool & colMap$time==min(colMap[tsBool,'time']),])
				minExp <- rownames(colMap[minExps,])[!colMap[minExps,'prevCol'] %in% rownames(colMap)]
				colMap[minExp,'is1stLast'] <- 'f'
				colMap[minExp,'prevCol'] <- NA
				maxExps <- rownames(colMap[tsBool & colMap$time==max(colMap[tsBool,'time']),])
				maxExp <- rownames(colMap[maxExps,])[!rownames(colMap[maxExps,]) %in% colMap[maxExps,'prevCol']]
				colMap[maxExp,'is1stLast'] <- 'l'
			}
			colMap$ts.ind <- new.ts.inds
			colMap$numTS <- i
        	} #if (any(exp.weights == 0))
        } #if (!is.null(row.weights))
    }
    #Remove any row differences between predictor and data matrices
    pred.data <- pred.data[rownames(pred.data) %in% rownames(ratios),]
    ratios <- ratios[rownames(ratios) %in% rownames(pred.data),]
    #newClusterStack <- e$clusterStack 
    newClusterStack <- list()
 
    for ( i in 1:length(e$clusterStack) ) { 
    	newClust <- e$clusterStack[[i]]
    	if (any(newClust$rows %in% rownames(ratios)) & any(newClust$cols %in% colnames(ratios))) {		
    		newClust$rows <- newClust$rows[newClust$rows %in% rownames(ratios)]
    		newClust$cols <- newClust$cols[newClust$cols %in% colnames(ratios)]
    		newClusterStack[[length(newClusterStack)+1]] <- newClust
    	} 
    }
    if (!is.null(predictors)) 
        predictors <- predictors[predictors %in% rownames(pred.data)]
    out <- runnit(ks=ks, data=ratios, col.map=colMap, predictors=predictors, clusterStack=newClusterStack, 
        pred.data=pred.data, gene.prefix = e$genome.info$gene.prefix, exp.weights=exp.weights, filter.pred.by.col=filter.pred.by.col, ...)
    invisible(out)
}
write.cytoscape.files <-
function (inf.result, clusterStack, sif.filename) 
{
    out <- unlist(inf.result)
    write("weight (java.lang.Double)", "weights.eda")
    gatecount = 1
    for (j in 1:length(out)) {
        nodes <- strsplit(names(out)[j], "\\.")[[1]]
        if (length(nodes) == 2) {
            write(paste(nodes[2], "activates", nodes[1]), sif.filename, 
                append = T)
            write(paste(nodes[2], "(activates)", nodes[1], "=", 
                out[j]), "weights.eda", append = T)
        }
        else if (length(nodes) == 4) {
            write(paste(nodes[2], "combines", paste("AND-", gatecount, 
                sep = "")), sif.filename, append = T)
            write(paste(nodes[3], "combines", paste("AND-", gatecount, 
                sep = "")), sif.filename, append = T)
            write(paste(paste("AND-", gatecount, sep = ""), "activates", 
                nodes[1]), sif.filename, append = T)
            write(paste(paste("AND-", gatecount, sep = ""), "(activates)", 
                nodes[1], "=", out[j]), "weights.eda", append = T)
            write(paste(paste("AND-", gatecount, sep = ""), "=", 
                "(logicGate)"), "types.noa", append = T)
            gatecount = gatecount + 1
        }
        else if (length(nodes) == 1) {
            write(nodes[1], sif.filename, append = T)
        }
    }
    write("clusterGenes", "clusterGenes.noa")
    write("clusterConditions", "clusterConditions.noa")
    write("clusterGeneCount", "clusterGeneCount.noa")
    write("clusterConditionCount", "clusterConditionCount.noa")
    write("clusterMotifPValues", "clusterMotifPValues.noa")
    write("clusterMotifs", "clusterMotifs.noa")
    write("clusterResiduals", "clusterResiduals.noa")
    for (i in 1:length(clusterStack)) {
        write(paste(paste(i, " = ", "(", sep = ""), paste(clusterStack[[i]]$rows, 
            collapse = "::"), ")", sep = ""), "clusterGenes.noa", 
            append = T)
        write(paste(paste(i, " = ", "(", sep = ""), paste(clusterStack[[i]]$cols, 
            collapse = "::"), ")", sep = ""), "clusterConditions.noa", 
            append = T)
        write(paste(i, " = ", clusterStack[[i]]$nrows, sep = ""), 
            "clusterGeneCount.noa", append = T)
        write(paste(i, " = ", clusterStack[[i]]$ncols, sep = ""), 
            "clusterConditionCount.noa", append = T)
        write(paste(i, " = ", clusterStack[[i]]$e.val, sep = ""), 
            "clusterMotifPValues.noa", append = T)
        write(paste(i, " = ", clusterStack[[i]]$resid, sep = ""), 
            "clusterResiduals.noa", append = T)
        write(paste(i, " = (cluster)", sep = ""), "types.noa", 
            append = T)
        if (length(clusterStack.redox[[i]]$motif.out$pssms) > 
            0) {
            write(paste(i, " = (", gsub(" ", "::", paste(lapply(clusterStack.redox[[i]]$motif.out$pssm, 
                pssm.to.string), collapse = " ")), ")", sep = ""), 
                "clusterMotifs.noa", append = T)
        }
        else {
            write(paste(i, "=", "()"), "clusterMotifs.noa", append = T)
        }
    }
}

#' Convert a ratios list into a matrix
#'  Copied from cMonkey source on 1/28/11 and added ratios input
#'
#' @param ratios  the list of ratios matrices
#' @param rows  the rows to select (DEFAULT: NULL, i.e. all)
#' @param cols  the cols to select (DEFAULT: NULL, i.e. all)
#' @param matrices  the submatrices to pull from (DEFAULT: names( ratios ), i.e. all)
#' @usage ratios<-get.cluster.matrix(ratios, rows=NULL, cols=NULL, matrices=names( ratios ))
#' @return A vector of similar elements to remove
get.cluster.matrix <- function(ratios, rows=NULL, cols=NULL, matrices=names( ratios ) ) {
  if (! is.null(ratios) ) {
	if ( is.null( attr( ratios, "rnames" ) ) ) attr( ratios, "rnames" ) <- unique(unlist(sapply(ratios,rownames))) #Match list up with cMonkey ratios list
	if ( is.null( rows ) ) rows <- attr( ratios, "rnames" )
	if ( is.null( attr( ratios, "cnames" ) ) ) attr( ratios, "cnames" ) <- unique(unlist(sapply(ratios,colnames)))
	if ( is.null( cols ) ) cols <- attr( ratios, "cnames" )
	cols.b <- attr( ratios, "cnames" )[ attr( ratios, "cnames" ) %in% cols ]

	rats <- matrix( NA, nrow=length( rows ), ncol=length( cols.b ) )
	rownames( rats ) <- rows; colnames( rats ) <- cols.b
	cnames <- character()
	for ( n in matrices ) { ## names( ratios ) ) {
	  r.tmp <- ratios[[ n ]][ rows[ rows %in% rownames( ratios[[ n ]] ) ],
				   cols.b[ cols.b %in% colnames( ratios[[ n ]] ) ], drop=F ]
	  if ( is.null( r.tmp ) || all( is.na( r.tmp ) ) ) next
	  if ( is.vector( r.tmp ) ) { r.tmp <- t( r.tmp ); rownames( r.tmp ) <- rows }
	  cnames <- c( cnames, colnames( r.tmp ) )
	  rats[ rownames( r.tmp ), colnames( r.tmp ) ] <- r.tmp; rm( r.tmp )
	}
  	out <- rats[ ,colnames( rats ) %in% cnames, drop=F ]
  } else {
  	out <- NULL
  }
  return(out)
}
