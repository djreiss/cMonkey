DATE <-
"Fri Nov 15 11:58:34 2013"
VERSION <-
"4.9.11"
.onLoad <-
function( libname, pkgname ) { ##.onAttach
    packageStartupMessage( "Loading ", pkgname, " version ", VERSION, " (", DATE, ")" )
    packageStartupMessage( "Copyright (C) David J Reiss, Institute for Systems Biology; dreiss@systemsbiology.org." )
    packageStartupMessage( "https://github.com/dreiss-isb/cmonkey" )
    if ( grepl( "beta", VERSION ) ) return()
    vers <- try( readLines( "http://rawgithub.com/dreiss-isb/cmonkey/master/VERSION" ), silent=T )
    if ( class( vers ) != "try-error" ) {
      vers <- gsub( " ", "", vers )
      if ( vers != VERSION ) packageStartupMessage( "\nYou are not using the most current version of cMonkey.\nPlease consider upgrading to v", vers, " via:\n\n> install.packages('devtools', dep=T)\n> require(devtools)\n> install_github('cmonkey', 'dreiss-isb', subdir='cMonkey')" )
      else packageStartupMessage( "Congratulations! You are using the latest version of cMonkey.\n" )
    } else {
      packageStartupMessage( "Could not check to see if you are using the latest version of cMonkey." )
    }
  }
adjust.all.clusters <-
function (env, ks = 1:env$k.clust, force.motif = T, ...) 
{
    old.stats <- env$stats
    tmp <- env$row.col.membership.from.clusterStack(env$clusterStack)
    row.membership <- tmp$r
    col.membership <- tmp$c
    rm(tmp)
    gc()
    mc <- env$get.parallel(length(ks))
    new.rm <- mc$apply(ks, function(k) env$adjust.clust(k, row.membership, 
        ...)$r)
    rm <- do.call(cbind, new.rm)
    for (i in 1:nrow(rm)) {
        tmp <- unique(rm[i, rm[i, ] != 0])
        rm[i, ] <- c(tmp, rep(0, ncol(rm) - length(tmp)))
    }
    rm <- rm[, apply(rm, 2, sum) != 0, drop = F]
    colnames(rm) <- NULL
    env$clusterStack <- lapply(1:env$k.clust, function(k) list(rows = rownames(which(rm == 
        k, arr = T)), cols = env$clusterStack[[k]]$cols))
    env$clusterStack <- env$get.clusterStack(ks = 1:k.clust)
    env$post.adjust <- FALSE
    env$cmonkey.one.iter(env, dont.update = T, force.row = T, 
        force.col = T, force.motif = if (force.motif & !no.genome.info) 
            "run.meme", force.net = T)
    print(rbind(OLD = old.stats[nrow(old.stats), ], NEW = env$stats[nrow(env$stats), 
        ]))
    invisible(env)
}
adjust.clust <-
function (k, row.memb = get("row.membership"), expand.only = T, 
    limit = 100, scores = "r.scores", quant.cutoff = 0.33, force.expand = 0) 
{
    if (scores == "rr.scores" || scores == "r.scores") {
        tmp <- get.combined.scores(quant = F)
        r.scores <- tmp$r
        if (scores == "rr.scores") {
            scores <- get.density.scores(ks = 1:k.clust)$r
            scores <- 1 - scores[, ]
        }
        else {
            scores <- r.scores
        }
        rm(r.scores, tmp)
        gc()
    }
    else {
        scores <- get(scores)
    }
    get.rows2 <- function(k, rm) rownames(which(rm == k, arr = T))
    scores <- scores[, ]
    old.rows <- get.rows(k)
    if (force.expand == 0) {
        wh <- names(which(scores[which(!attr(ratios, "rnames") %in% 
            old.rows), k] < quantile(scores[old.rows, k], quant.cutoff, 
            na.rm = T)))
    }
    else {
        expand.only <- TRUE
        wh <- names(sort(scores[!attr(ratios, "rnames") %in% 
            old.rows, k], decreasing = F)[1:force.expand])
    }
    if (length(wh) > limit) {
        warning("Surpassing limit.")
        return(invisible(list(r = row.memb)))
    }
    else if (length(wh) <= 0) 
        return(invisible(list(r = row.memb)))
    tries <- 0
    while (length(wh) > 0 && tries < 50) {
        wh2 <- names(which.max(scores[wh, k]))
        wh2.scores <- scores[wh2, row.memb[wh2, ]]
        wh2a <- names(which.max(scores[get.rows2(k, rm = row.memb), 
            k]))
        for (col in 1:ncol(row.memb)) if (all(row.memb[wh2, col] == 
            0)) 
            break
        if (col == ncol(row.memb) && any(row.memb[wh2, col] != 
            0)) {
            row.memb <- cbind(row.memb, rep(0, nrow(row.memb)))
            col <- col + 1
        }
        row.memb[wh2, col] <- k
        if (!expand.only) 
            row.memb[wh2a, row.memb[wh2a, ] == k] <- 0
        if (force.expand == 0) {
            wh <- names(which(scores[which(!attr(ratios, "rnames") %in% 
                get.rows2(k, rm = row.memb)), k] < quantile(scores[get.rows2(k, 
                rm = row.memb), k], quant.cutoff, na.rm = T)))
        }
        else {
            wh <- wh[!wh %in% wh2]
        }
        if (length(get.rows2(k, rm = row.memb)) > cluster.rows.allowed[2]) 
            break
        tries <- tries + 1
    }
    new.rows <- get.rows2(k, rm = row.memb)
    rm(scores)
    gc()
    if (any(!new.rows %in% old.rows) || any(!old.rows %in% new.rows)) 
        cat("ADJUSTED CLUSTER:", k, length(old.rows), length(new.rows), 
            sum(!old.rows %in% new.rows), "\n")
    row.memb <- t(apply(row.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    row.memb <- row.memb[, apply(row.memb, 2, sum) != 0, drop = F]
    colnames(row.memb) <- NULL
    invisible(list(r = row.memb))
}
cluster.pclust <-
function (k, mot.inds = "COMBINED") 
{
    inds <- mot.inds
    if (is.null(inds)) 
        return(list(p.clusts = NA, e.vals = NA))
    if (mot.inds[1] == "COMBINED") 
        inds <- names(get("mot.weights"))
    if (is.null(inds)) 
        return(list(p.clusts = NA, e.vals = NA))
    rows <- get.rows(k)
    p.clusts <- sapply(inds, function(n) {
        ms <- meme.scores[[n]][[k]]
        out <- NA
        if (length(rows) > 0 && !is.null(ms$pv.ev) && !is.null(ms$pv.ev[[1]])) {
            if ("p.value" %in% colnames(ms$pv.ev[[1]])) 
                out <- mean(log10(ms$pv.ev[[1]][rownames(ms$pv.ev[[1]]) %in% 
                  rows, "p.value"]), na.rm = T)
            else if ("pvals" %in% colnames(ms$pv.ev[[1]])) 
                out <- mean(log10(ms$pv.ev[[1]][rownames(ms$pv.ev[[1]]) %in% 
                  rows, "pvals"]), na.rm = T)
        }
        out
    })
    e.vals <- sapply(inds, function(n) {
        ms <- meme.scores[[n]][[k]]
        sapply(1:length(ms$meme.out), function(i) if (length(rows) > 
            0 && !is.null(ms$meme.out) && !is.null(ms$meme.out[[i]])) 
            ms$meme.out[[i]]$e.value
        else NA)
    })
    if (!is.matrix(e.vals)) 
        e.vals <- t(t(e.vals))
    if (mot.inds[1] == "COMBINED") {
        p.clusts <- weighted.mean(p.clusts, mot.weights[inds], 
            na.rm = T)
        e.vals <- apply(e.vals, 1, weighted.mean, mot.weights[inds], 
            na.rm = T)
    }
    else if (mot.inds[1] != "COMBINED") {
        if (length(p.clusts) < length(inds) && all(is.na(p.clusts))) 
            p.clusts <- rep(NA, length(inds))
        e.vals <- apply(e.vals, 1, function(i) if (length(i) < 
            length(inds) && all(is.na(i))) 
            rep(NA, length(inds))
        else i)
    }
    if (!is.matrix(e.vals)) 
        e.vals <- t(t(e.vals))
    if (is.matrix(e.vals) && ncol(e.vals) != length(inds)) 
        e.vals <- t(e.vals)
    if (mot.inds[1] != "COMBINED") 
        names(p.clusts) <- colnames(e.vals) <- inds
    else e.vals <- as.vector(e.vals)
    list(p.clusts = p.clusts, e.vals = e.vals)
}
cluster.resid <-
function (k, rats.inds = "COMBINED", varNorm = F, in.cols = T, 
    ...) 
{
    residual.submatrix <- function(rats, rows, cols, varNorm = F, 
        ...) {
        rows <- rows[rows %in% rownames(rats)]
        cols <- cols[cols %in% colnames(rats)]
        if (length(rows) <= 1 || length(cols) <= 1) 
            return(1)
        maxRowVar <- attr(rats, "maxRowVar")
        rats <- rats[rows, cols]
        if (is.vector(rats) || any(dim(rats) <= 1) || mean(is.na(rats)) > 
            0.95) 
            return(1)
        d.rows <- rowMeans(rats, na.rm = T)
        d.cols <- colMeans(rats, na.rm = T)
        d.all <- mean(d.rows, na.rm = T)
        rats[, ] <- rats[, ] + d.all - outer(d.rows, d.cols, 
            "+")
        average.r <- mean(abs(rats), na.rm = TRUE)
        if (varNorm && !is.null(maxRowVar)) {
            row.var <- mean(apply(rats, 1, var, use = "pairwise.complete.obs"), 
                na.rm = T)
            if (is.na(row.var) || row.var > maxRowVar) 
                row.var <- maxRowVar
            average.r <- average.r/row.var
        }
        average.r
    }
    inds <- rats.inds
    if (rats.inds[1] == "COMBINED") 
        inds <- names(get("row.weights"))
    rows <- get.rows(k)
    cols <- get.cols(k)
    resids <- sapply(ratios[inds], function(rn) {
        if (in.cols) 
            residual.submatrix(rn, rows, cols, varNorm = varNorm)
        else residual.submatrix(rn, get.rows(k), colnames(rn)[!colnames(rn) %in% 
            cols], varNorm = varNorm)
    })
    if (rats.inds[1] == "COMBINED") 
        resids <- weighted.mean(resids, row.weights[inds], na.rm = T)
    if (rats.inds[1] != "COMBINED" && length(resids) < length(inds) && 
        all(is.na(resids))) {
        resids <- rep(NA, length(inds))
        names(resids) <- inds
    }
    resids
}
cluster.summary <-
function (e.cutoff = 0.01, nrow.cutoff = 5, seq.type = names(mot.weights)[1], 
    plot = F, sort = c("score.norm", "score", "resid", "e.value1", 
        "e.value2", "nrow")) 
{
    ms <- NULL
    if (!is.null(seq.type)) 
        ms <- meme.scores[[seq.type]]
    if (is.null(ms)) 
        e.cutoff <- NA
    score <- sapply(1:k.clust, function(k) mean(row.scores[get.rows(k), 
        k], na.rm = T, trim = 0.01)) * row.scaling[iter] + if (!is.null(mot.scores)) 
        sapply(1:k.clust, function(k) mean(mot.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * mot.scaling[iter]
    else 0 + if (!is.null(net.scores)) 
        sapply(1:k.clust, function(k) mean(net.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * net.scaling[iter]
    else 0
    nrow <- sapply(1:k.clust, function(k) length(get.rows(k)))
    out <- data.frame(k = 1:k.clust, nrow = nrow, score = score, 
        resid = sapply(1:k.clust, cluster.resid, varNorm = F), 
        consensus1 = sapply(1:k.clust, function(k) if (is.null(ms) || 
            is.null(ms[[k]]$meme.out) || length(ms[[k]]) <= 3) 
            ""
        else pssm.to.string(ms[[k]]$meme.out[[1]]$pssm)), e.value1 = sapply(1:k.clust, 
            function(k) if (is.null(ms) || is.null(ms[[k]]$meme.out) || 
                length(ms[[k]]) <= 3) 
                Inf
            else ms[[k]]$meme.out[[1]]$e.value), consensus2 = sapply(1:k.clust, 
            function(k) if (is.null(ms) || is.null(ms[[k]]$meme.out) || 
                length(ms[[k]]) <= 3) 
                ""
            else if (length(ms[[k]]$meme.out) == 1) 
                ""
            else pssm.to.string(ms[[k]]$meme.out[[2]]$pssm)), 
        e.value2 = sapply(1:k.clust, function(k) if (is.null(ms) || 
            is.null(ms[[k]]$meme.out) || length(ms[[k]]) <= 3) 
            Inf
        else if (length(ms[[k]]$meme.out) <= 1) 
            Inf
        else ms[[k]]$meme.out[[2]]$e.value))
    if (all(out$consensus2 == "")) 
        out$consensus2 <- out$e.value2 <- NULL
    if (!is.na(sort[1]) && sort[1] %in% colnames(out)) 
        out <- out[order(out[[sort[1]]]), ]
    out
}
clusters.w.conds <-
function (conds, ks = 1:k.clust, p.val = F) 
{
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        cols <- get.cols(i)
        if (!p.val) 
            sum(cols %in% conds)
        else phyper(sum(cols %in% conds), length(conds), attr(ratios, 
            "ncol") - length(conds), length(cols), lower = F) * 
            length(ks)
    }))
}
clusters.w.func <-
function (func, ks = 1:k.clust, short = F, max.rows = 999, p.val = F) 
{
    if (p.val) {
        long.names <- get.long.names(attr(ratios, "rnames"), 
            short = short)
        n2 <- length(grep(func, long.names, perl = T, ignore.case = T))
    }
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        rows <- get.rows(i)
        if (length(rows) <= 1) 
            return(NA)
        rows.l <- get.long.names(rows, short = short)
        if (!p.val) {
            if (length(rows) >= max.rows) 
                NA
            else length(grep(func, rows.l, perl = T, ignore.case = T))
        }
        else {
            phyper(length(grep(func, rows.l, perl = T, ignore.case = T)), 
                n2, attr(ratios, "nrow") - n2, length(rows), 
                lower = F) * length(ks)
        }
    }))
}
clusters.w.genes <-
function (genes, ks = 1:k.clust, p.val = F) 
{
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        rows <- get.rows(i)
        if (length(rows) <= 1) 
            return(NA)
        if (!p.val) 
            sum(rows %in% genes)
        else phyper(sum(rows %in% genes), length(genes), attr(ratios, 
            "nrow") - length(genes), length(rows), lower = F) * 
            length(ks)
    }))
}
cluster.the.new.motif.clusters <-
function (c, cutoff = 0.15) 
{
    args <- attr(c, "in.args")
    motifs <- args$motifs
    if (motifs == "ALL") {
        motifs <- unlist(lapply(1:nrow(motif.widths), function(i) {
            ii <- which(motif.widths[i, ] > 0)
            if (length(ii) <= 0) 
                return(NULL)
            paste("MOT", i, ii, sep = "_")
        }))
    }
    max.motifs <- length(motifs)
    cat(length(motifs), "motifs\n")
    frac.in.coding <- coding.fracs$all.fracs[motifs]
    motifs <- motifs[!is.na(frac.in.coding) & frac.in.coding < 
        coding.fracs$mean.fracs - 0.01]
    cat(length(motifs), "motifs\n")
    dirname <- sprintf("filehash/new_motif_shadows_%s_%d", as.character(args$p.cutoff), 
        length(motifs))
    system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
        dirname, args$distance.weight.cutoff))
    system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
        dirname, args$distance.weight.cutoff))
    tmp <- fread("new.mcltemp", sep = " ")
    require(Matrix)
    m <- Matrix(0, nrow = max.motifs, ncol = max.motifs)
    m[as.matrix(tmp[, V1, V2])] <- tmp[, V3]
    m[as.matrix(tmp[, V2, V1])] <- tmp[, V3]
    rm(tmp)
    gc()
    ov.means <- Matrix(data = 0, nrow = length(c), ncol = length(c), 
        sparse = T)
    rnames <- 1:length(motifs)
    names(rnames) <- motifs
    for (i in 1:(attr(c, "mc.length") - 1)) {
        m.i <- rnames[c[[i]]]
        m.i <- m.i[!is.na(m.i)]
        if (length(m.i) <= 0) 
            next
        x.tmp <- m[m.i, ]
        for (j in (i + 1):length(c)) {
            m.j <- rnames[c[[j]]]
            m.j <- m.j[!is.na(m.j)]
            if (length(m.j) <= 0) 
                next
            xx.tmp <- x.tmp[, m.j]
            ov.means[i, j] <- mean(xx.tmp != 0)
        }
        cat("Similar:", i, max(ov.means[i, ]), which(ov.means[i, 
            ] > cutoff), "\n")
    }
    ov.means <- as.matrix(ov.means)
    ov.means[lower.tri(ov.means)] <- t(ov.means)[lower.tri(t(ov.means))]
    mc.length <- attr(c, "mc.length")
    if (!any(ov.means > cutoff)) 
        cutoff <- cutoff - 0.05
    ind <- 1
    motif.cluster.clusters <- rep(0, ncol(ov.means))
    for (i in 1:mc.length) {
        tt <- which(ov.means[i, ] > cutoff)
        ttmp2 <- tt
        ttmp2 <- ttmp2[ttmp2 != i]
        if (length(ttmp2) <= 0) 
            next
        best.match <- ttmp2[which.max(ov.means[i, ttmp2])]
        if (motif.cluster.clusters[i] != 0) {
            motif.cluster.clusters[tt] <- motif.cluster.clusters[i]
        }
        else if (best.match <= mc.length && motif.cluster.clusters[best.match] != 
            0) {
            motif.cluster.clusters[tt] <- motif.cluster.clusters[best.match]
        }
        else {
            tt <- c(i, tt)
            tt <- tt[motif.cluster.clusters[tt] == 0]
            motif.cluster.clusters[tt] <- ind
            ind <- ind + 1
        }
        cat(i, tt, motif.cluster.clusters[i], "\n")
    }
    m = c(which(motif.cluster.clusters[1:mc.length] == 0), which(motif.cluster.clusters != 
        0))
    cat("Total unique motif clusters:", length(table(motif.cluster.clusters[m])) + 
        table(motif.cluster.clusters[m])["0"] - 1, "\n")
    out2 <- new.env()
    out2$ov.means <- ov.means
    out2$motif.cluster.clusters <- motif.cluster.clusters
    unlink("new.mcltemp")
    out2
}
cmonkey <-
function (env = NULL, ...) 
{
    if (((is.null(list(...)$dont.init) || !list(...)$dont.init) && 
        (is.null(env$dont.init) || !env$dont.init) && (!exists("dont.init") || 
        !dont.init)) || is.null(env) || is.null(env$genome.info)) {
        env <- cmonkey.init(env, ...)
    }
    else {
        if (sink.number() > 0) 
            for (i in 1:sink.number()) try(sink(), silent = T)
        if (env$save.logfile != FALSE) 
            sink(env$save.logfile, split = T, append = T)
    }
    cat("\033[31mTIME STARTED:", env$time.started, "\033[0m\n")
    if ((!exists("clusterStack", envir = env) || length(env$clusterStack) < 
        env$k.clust) && exists("ratios", envir = env)) 
        env$cmonkey.re.seed(env)
    while (env$iter <= env$n.iter) {
        iter <- env$iter
        env$cmonkey.one.iter(env)
    }
    if (!is.na(env$plot.iters) && (iter %in% env$plot.iters || 
        (iter - 1) %in% env$plot.iters)) 
        try(env$plotStats(iter, plot.clust = env$favorite.cluster()))
    env$iter <- iter <- env$iter - 1
    print(env$cluster.summary())
    parent.env(env) <- globalenv()
    parent.env(env$cmonkey.params) <- env
    env$clusterStack <- env$get.clusterStack(ks = 1:env$k.clust, 
        force = T)
    print(env$cluster.summary())
    env$set.param("time.ended", date(), env$cmonkey.params)
    env$time.ended <- date()
    cat("\033[31mTIME ENDED:", env$time.ended, "\033[0m\n")
    invisible(env)
}
cmonkey.init <-
function (env = NULL, ...) 
{
    if (!exists("cmonkey.params")) {
        cmonkey.params <- new.env(hash = T)
    }
    tmp.e <- environment(cMonkey:::cmonkey)
    if (!is.null(env) && (is.list(env) || is.environment(env))) {
        for (i in names(env)) if (!i %in% names(list(...))) 
            assign(i, env[[i]])
        if (is.list(env)) 
            env <- NULL
    }
    for (i in ls(tmp.e)) {
        f2 <- NULL
        if ((!is.null(env) && exists(i, envir = env, inherit = F))) {
            f <- try(get(i, envir = env))
            f2 <- try(get(i, envir = tmp.e))
        }
        else if (exists(i, envir = .GlobalEnv, inherit = F)) {
            f <- try(get(i, envir = .GlobalEnv))
            f2 <- try(get(i, envir = tmp.e))
        }
        else if (exists(i)) {
            f <- try(get(i))
            f2 <- try(get(i, envir = tmp.e))
        }
        else {
            f <- try(get(i, envir = tmp.e))
        }
        if (class(f) == "function") 
            environment(f) <- sys.frames()[[length(sys.frames())]]
        assign(i, f)
        if (!is.null(f2) && class(f2) == "function" && object.size(f2) != 
            object.size(f)) {
            environment(f2) <- sys.frames()[[length(sys.frames())]]
            assign(paste("super", i, sep = "."), f2)
            if (!is.null(env)) {
                assign(paste("super", i, sep = "."), f2, envir = env)
                environment(env[[paste("super", i, sep = ".")]]) <- env
            }
        }
    }
    rm(f, f2, tmp.e)
    if (!is.null(env)) 
        for (i in ls(env)) assign(i, get(i, env))
    args <- mget(names(formals()), env = as.environment(-1))
    for (i in names(args)) if (!i %in% c("...", "env")) 
        set.param(i, args[[i]])
    for (i in names(list(...))) if (i != "env") 
        set.param(i, list(...)[[i]])
    rm(args)
    if (sink.number() > 0) 
        for (i in 1:sink.number()) try(sink(), silent = T)
    set.param("save.logfile", FALSE)
    if (save.logfile != FALSE) 
        sink(save.logfile, split = T, append = (exists("dont.init") && 
            dont.init) || (exists("is.inited") && !is.inited))
    if (!exists("organism")) {
        message("WARNING: No organism was set; using \"hpy\".")
        organism <- "hpy"
        Sys.sleep(3)
    }
    set.param("organism", organism)
    if (!exists("ratios") && exists("row.weights")) {
        try({
            ratios <- lapply(names(row.weights), get)
            names(ratios) <- names(row.weights)
        })
    }
    if ((exists("ratios") && !is.null(ratios))) {
        if (is.matrix(ratios) || is.data.frame(ratios)) 
            ratios <- list(ratios = load.ratios(ratios))
        else if (is.character(ratios)) 
            ratios <- lapply(ratios, load.ratios)
        else ratios <- lapply(ratios, function(r) as.matrix(load.ratios(r)))
        ratios <- ratios[sapply(ratios, function(r) all(dim(r) > 
            0))]
        attr(ratios, "rnames") <- sort(unique(unlist(lapply(ratios, 
            rownames))))
        attr(ratios, "cnames") <- sort(unique(unlist(lapply(ratios, 
            colnames))))
        attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
        attr(ratios, "ncol") <- length(attr(ratios, "cnames"))
        if (is.null(names(ratios))) {
            names(ratios) <- paste("ratios", 1:length(ratios), 
                sep = ".")
            if (exists("row.weights") && !all(names(ratios) %in% 
                names(row.weights))) 
                for (i in names(ratios)) row.weights[i] <- row.weights[1]
        }
        for (n in names(ratios)) {
            if (ncol(ratios[[n]]) > 1) {
                attr(ratios[[n]], "maxRowVar") <- mean(apply(ratios[[n]][, 
                  ], 1, var, use = "pair"), na.rm = T)
                attr(ratios[[n]], "all.colVars") <- apply(ratios[[n]][, 
                  ], 2, var, use = "pair", na.rm = T)
            }
        }
        rm(n)
    }
    if (exists("ratios") && is.null(names(ratios))) 
        names(ratios) <- paste("ratios", 1:length(ratios), sep = ".")
    if (!is.null(env) && exists("ratios")) 
        assign("ratios", ratios, envir = env)
    if (length(ratios) == 1) 
        names(ratios) <- "ratios"
    set.param("cog.org", "?")
    set.param("rsat.species", "?")
    set.param("n.iter", 2000)
    set.param("n.clust.per.row", 2)
    if (exists("ratios") && !is.null(ratios)) {
        set.param("k.clust", round(attr(ratios, "nrow") * n.clust.per.row/20))
    }
    else {
        set.param("k.clust", 100)
    }
    set.param("n.clust.per.col", if (exists("ratios") && attr(ratios, 
        "ncol") >= 60) 
        round(k.clust/2)
    else round(k.clust * 2/3))
    set.param("all.genes.even.if.not.in.ratios", FALSE)
    set.param("row.iters", seq(1, n.iter, by = 2))
    set.param("col.iters", seq(1, n.iter, by = 5))
    set.param("meme.iters", c(1, seq(100, n.iter, by = 100)))
    set.param("mot.iters", seq(1, n.iter, by = 10))
    set.param("net.iters", seq(1, n.iter, by = 7))
    set.param("row.scaling", 1)
    set.param("row.weights", c(ratios = 1))
    set.param("mot.scaling", c(rep(1e-05, 100), seq(0, 1, length = n.iter * 
        3/4)))
    set.param("mot.weights", c(`upstream meme` = 1))
    set.param("net.scaling", seq(1e-05, 0.5, length = n.iter * 
        3/4))
    set.param("net.weights", c(string = 0.5, operons = 0.5))
    set.param("grouping.weights", numeric())
    set.param("plot.iters", seq(2, n.iter, by = 25))
    set.param("post.adjust", TRUE)
    set.param("parallel.cores", TRUE)
    set.param("parallel.cores.motif", TRUE)
    set.param("max.changes", c(rows = 0.5, cols = 5))
    set.param("cluster.rows.allowed", c(3, 70))
    set.param("merge.cutoffs", c(n = 0.3, cor = 0.975))
    set.param("fuzzy.index", 0.75 * exp(-(1:n.iter)/(n.iter/4)))
    set.param("translation.tab", NULL)
    set.param("seed.method", c(rows = "kmeans", cols = "best"))
    set.param("maintain.seed", NULL)
    set.param("n.motifs", c(rep(1, n.iter/3), rep(2, n.iter/3)))
    if (file.exists("./progs") && file.exists("./progs/meme")) {
        set.param("progs.dir", "./progs/")
    }
    else if ("package:cMonkey" %in% search() && file.exists(sprintf("%s/progs/", 
        system.file(package = "cMonkey")))) {
        set.param("progs.dir", sprintf("%s/progs/", system.file(package = "cMonkey")))
    }
    else if (any(mot.scaling > 0) && (!exists("no.genome.info") || 
        !no.genome.info)) {
        message("WARNING: You do not have meme/mast/dust installed in the correct location.\nTrying to install it now.\n")
        install.binaries()
        set.param("progs.dir", sprintf("%s/progs/", system.file(package = "cMonkey")))
        if ("package:cMonkey" %in% search() && !file.exists(sprintf("%s/progs/", 
            system.file(package = "cMonkey")))) 
            message("WARNING: Could not install meme. Please see the website for installation instructions.")
    }
    set.param("meme.cmd", paste(progs.dir, "meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $compute", 
        sep = "/"))
    set.param("mast.cmd", sprintf("%s/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 999999 -mev 9999999 -mt 0.99 -seqp -remcorr", 
        progs.dir))
    set.param("dust.cmd", sprintf("%s/dust $fname", progs.dir))
    set.param("operon.shift", TRUE)
    set.param("bg.order", 3)
    set.param("recalc.bg", FALSE)
    set.param("motif.upstream.search", c(-20, 150))
    set.param("motif.upstream.scan", c(-30, 250))
    set.param("rsat.urls", c("http://rsat.ccb.sickkids.ca/", 
        "http://rsat.ulb.ac.be/rsat/", "http://embnet.ccg.unam.mx/rsa-tools"))
    set.param("stats.iters", c(1, seq(5, n.iter, by = 5)))
    set.param("cm.script.each.iter", "cm.script.each.iter.R")
    set.param("date.run", format(Sys.time(), "%y %b %d %H:%M:%S"))
    set.param("cmonkey.version", cm.version)
    set.param("session.info", unlist(list(R.version, Sys.info(), 
        Sys.getenv(), sessionInfo())), quiet = T)
    set.param("time.started", date())
    if (exists("ratios") && !is.null(ratios)) {
        set.param("cmonkey.filename", paste("cmonkey", cmonkey.version, 
            organism, paste(sapply(ratios, dim), collapse = "x"), 
            gsub(" ", "_", date.run), sep = "_"))
    }
    else {
        set.param("cmonkey.filename", paste("cmonkey", cmonkey.version, 
            organism, "0x0", gsub(" ", "_", date.run), sep = "_"))
    }
    op <- options(digits.secs = 10)
    set.param("rnd.seed", as.integer(substr(gsub("[-:. ]", "", 
        as.character(Sys.time())), 12, 20)))
    options(op)
    rm(op)
    set.seed(rnd.seed)
    set.param("big.memory", FALSE)
    set.param("big.memory.verbose", FALSE)
    if (organism == "hsa") 
        rsat.urls[1] <- rsat.urls[2]
    if (!exists("rsat.species") || rsat.species == "?" || is.na(rsat.species)) {
        err <- dlf("data/KEGG/KEGG_taxonomy.txt", "http://baliga.systemsbiology.net/cmonkey/taxonomy.txt")
        if (class(err) != "try-error") {
            tab <- read.delim("data/KEGG/KEGG_taxonomy.txt", 
                sep = "\t", comment = "#", head = F, as.is = T)
            rsat.spec <- as.character(subset(tab, V2 == organism, 
                select = "V4", drop = T))[1]
            rm(tab)
            if (any(strsplit(rsat.spec, "")[[1]] == "(")) 
                rsat.spec <- gsub("\\s\\(.*\\)", "", rsat.spec)
        }
        else {
            rsat.spec <- "?"
        }
        rsat.spec <- gsub(" ", "_", rsat.spec, fixed = T)
        kegg.spec <- rsat.spec
        if (!file.exists("data/RSAT_genomes_listing.txt")) {
            require(RCurl)
            tmp <- strsplit(getURL(paste(rsat.urls[1], "/data/genomes/", 
                sep = "")), "\n")[[1]]
            writeLines(tmp, con = "data/RSAT_genomes_listing.txt")
        }
        vals <- character()
        if (file.exists("data/RSAT_genomes_listing.txt")) {
            tmp <- readLines("data/RSAT_genomes_listing.txt")
            vals <- grep(rsat.spec, tmp, fixed = T, val = T)
        }
        if (length(vals) <= 0) {
            message("Could not find correct organism for RSAT... will try to guess...")
            max.dist <- 0.5
            vals <- rep("", 3)
            while (length(vals) > 1) {
                vals <- agrep(rsat.spec, tmp, ignore = T, max.dist = max.dist, 
                  val = T)
                max.dist <- max.dist - 0.01
                if (length(vals) <= 0) {
                  max.dist <- max.dist + 0.02
                  vals <- agrep(rsat.spec, tmp, ignore = T, max.dist = max.dist, 
                    val = T)
                  break
                }
            }
            if (length(vals) > 1) {
                rsat.spec <- sapply(strsplit(vals, "[<>/]"), 
                  "[", 8)
                message("Found ", length(rsat.spec), " matches...")
                rsat.spec <- rsat.spec[menu(rsat.spec, graphics = F, 
                  title = "Please choose one.")]
            }
            if (length(vals) == 1) {
                rsat.spec <- strsplit(vals, "[<>/]")[[1]][8]
                message("Found one match: ", rsat.spec, " ...")
                message("If this is not correct, you're not quite out of luck -- set the 'rsat.species' parameter manually.")
            }
        }
        set.param("rsat.species", rsat.spec, override = T)
        rm(tmp, rsat.spec, err, vals)
    }
    else {
        set.param("rsat.species", rsat.species)
    }
    if (!exists("taxon.id") || taxon.id == "?" || is.na(taxon.id) || 
        length(taxon.id) <= 0) {
        fname <- dlf("data/GO/proteome2taxid", "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid")
        tab <- read.delim(gzfile("data/GO/proteome2taxid"), head = F)
        taxon.id <- subset(tab, V1 == gsub("_", " ", rsat.species))$V2
        if (length(taxon.id) <= 0) 
            taxon.id <- subset(tab, grepl(gsub("_", " ", rsat.species), 
                V1))$V2[1]
        set.param("taxon.id", taxon.id, override = T)
        rm(tab, fname)
    }
    if (!exists("cog.org") || cog.org == "?" || is.na(cog.org)) {
        tmp <- strsplit(organism, "")[[1]]
        tmp[1] <- toupper(tmp[1])
        cog.o <- paste(tmp, sep = "", collapse = "")
        if (cog.o == "") 
            cog.o <- "?"
        set.param("cog.org", cog.o, override = T)
        rm(cog.o, tmp)
    }
    else {
        set.param("cog.org", cog.org)
    }
    message("Organism is ", organism, " ", cog.org, " ", rsat.species, 
        " ", taxon.id)
    genome.loc <- paste(rsat.urls[1], "/data/genomes/", rsat.species, 
        "/genome/", sep = "")
    fname <- paste("data/", rsat.species, "/organism.tab", sep = "")
    err <- dlf(fname, paste(genome.loc, "/organism.tab", sep = ""))
    org.tab <- readLines(fname)
    org.tab <- strsplit(org.tab[length(org.tab)], "\t")[[1]]
    is.eukaryotic <- any(grepl("Eukaryota", org.tab))
    cat("Is eukaryote:", is.eukaryotic, "\n")
    rm(err, org.tab, genome.loc, fname)
    if (is.eukaryotic) {
        message("Organism is a eukaryote; presuming there are no operons.")
        set.param("is.eukaryotic", TRUE, override = T)
        set.param("operon.shift", FALSE, override = T)
        set.param("discard.genome", TRUE, override = T)
        set.param("recalc.bg", FALSE, override = T)
        if ("operons" %in% names(net.weights)) {
            net.weights <- net.weights[names(net.weights) != 
                "operons"]
            set.param("net.weights", net.weights, override = T)
        }
    }
    if (get.parallel(100, verbose = T)$mc) 
        on.exit(try(kill(children(), SIGKILL)), add = T)
    on.exit({
        if (sink.number() > 0) for (i in 1:sink.number()) try(sink(), 
            silent = T)
    })
    if (sum(net.weights, na.rm = T) > 0) 
        net.weights <- net.weights/sum(net.weights, na.rm = T)
    if (sum(row.weights, na.rm = T) > 0) 
        row.weights <- row.weights/sum(row.weights, na.rm = T)
    if (sum(mot.weights, na.rm = T) > 0) 
        mot.weights <- mot.weights/sum(mot.weights, na.rm = T)
    for (i in c("n.motifs", "meme.cmd", "mast.cmd", "meme.iters", 
        "operon.shift", "bg.order", "motif.upstream.search", 
        "motif.upstream.scan")) {
        v <- get(i)
        if (all(names(mot.weights) %in% names(v))) 
            next
        if (is.vector(v) && length(v) > 1) 
            v <- list(`1` = v)
        names(v) <- names(mot.weights)[1]
        for (n in names(mot.weights)[!names(mot.weights) %in% 
            names(v)]) {
            if (is.list(v)) 
                v[[n]] <- v[[1]]
            else if (is.vector(v)) 
                v[n] <- v[1]
            names(v)[length(v)] <- names(mot.weights)[length(v)]
        }
        assign(i, v)
    }
    rm(v)
    if (!is.null(env)) 
        for (i in ls()) {
            if (i %in% c("i", "env")) 
                next
            tmp <- get(i)
            if (class(tmp) == "function") 
                environment(tmp) <- env
            assign(i, tmp, envir = env)
        }
    if (!is.na(rsat.species) && (!exists("genome.info") || genome.info$species != 
        rsat.species)) {
        cat("Initializing genome info for organism", organism, 
            "\n")
        set.param("no.genome.info", FALSE)
        genome.info <- get.genome.info()
        if (!is.null(env)) 
            assign("genome.info", genome.info, envir = env)
        if (is.na(taxon.id) || length(taxon.id) <= 0) {
            taxon.id <- genome.info$taxon.id
            set.param("taxon.id", taxon.id, override = T)
            message("Organism is ", organism, " ", cog.org, " ", 
                rsat.species, " ", taxon.id)
        }
        if (exists("ratios") && !is.null(ratios)) 
            tmp <- toupper(attr(ratios, "rnames"))
        else if (exists("genome.info") && !is.null(genome.info$feature.names)) {
            tmp <- toupper(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T))
            if (exists("ratios") && !is.null(ratios)) 
                tmp <- tmp[toupper(tmp) %in% toupper(attr(ratios, 
                  "rnames"))]
        }
        qqq <- sapply(1:4, function(nch) max(table(substr(tmp, 
            1, nch)))/length(tmp))
        nch <- 0
        if (any(qqq > 0.9)) {
            nch <- which(qqq > 0.9)
            nch <- nch[length(nch)]
        }
        else if (any(qqq > 0.6)) {
            nch <- which(qqq > 0.6)
            nch <- nch[length(nch)]
        }
        else if (any(qqq > 0.4)) {
            nch <- which(qqq > 0.4)
            nch <- nch[length(nch)]
        }
        prefix <- NA
        if (nch > 0) {
            prefix <- names(which.max(table(substr(tmp, 1, nch))))
            message("Assuming gene/probe names have common prefix '", 
                prefix, "'.")
            genome.info$gene.prefix <- prefix
        }
        else {
            message("Could not find a common gene/probe identifier prefix. This only matters if there's no expression matrix.")
            prefix <- genome.info$gene.prefix <- NA
        }
        if (TRUE) {
            tmp2 <- tmp
            if (length(unique(nchar(tmp2))) > 1) {
                nc <- max(nchar(tmp2))
                for (i in 1:length(tmp2)) {
                  tmp2[i] <- paste(tmp2[i], paste(rep(" ", nc - 
                    nchar(tmp2[i])), sep = "", collapse = ""), 
                    sep = "")
                }
            }
            tmp2 <- do.call(rbind, strsplit(tmp2, ""))
            regex <- apply(tmp2, 2, function(i) sort(unique(i)))
            for (i in 1:length(regex)) {
                ii <- as.integer(regex[[i]])
                if (!any(is.na(ii))) {
                  if (length(ii) == length(ii[1]:ii[length(ii)]) && 
                    all(sort(ii) == ii[1]:ii[length(i)])) 
                    regex[[i]] <- paste("[", paste(ii[1], ii[length(ii)], 
                      sep = "-"), "]", sep = "")
                }
                if (length(regex[[i]][regex[[i]] != " "]) > 1) 
                  regex[[i]] <- c("[", regex[[i]], "]")
                if (any(regex[[i]] == "" | regex[[i]] == " " | 
                  is.na(regex[[i]]))) 
                  regex[[i]] <- c(regex[[i]][regex[[i]] != " "], 
                    "?")
            }
            regex <- paste(unlist(lapply(regex, paste, sep = "", 
                collapse = "")), sep = "", collapse = "")
            regex <- gsub(".", "\\.", regex, fixed = T)
            message("Assuming gene/probe names have regex '", 
                regex, "'.")
            genome.info$gene.regex <- regex
            rm(tmp, tmp2)
        }
        all.names <- unique(c(as.character(genome.info$feature.names$id), 
            as.character(genome.info$feature.names$names), attr(ratios, 
                "rnames")))
        if (exists("translation.tab")) 
            all.names <- unique(c(all.names, as.character(translation.tab$V1), 
                as.character(translation.tab$V2)))
        if (!all.genes.even.if.not.in.ratios) 
            all.names <- all.names[all.names %in% attr(ratios, 
                "rnames")]
        genome.info$all.gene.names <- unique(grep(paste("^", 
            genome.info$gene.regex, sep = ""), all.names, perl = T, 
            val = T, ignore = T))
        if (length(genome.info$all.gene.names) <= length(attr(ratios, 
            "rnames"))/2) {
            genome.info$all.gene.names <- c(genome.info$all.gene.names, 
                unique(grep(paste("^", genome.info$gene.prefix, 
                  sep = ""), all.names, perl = T, val = T, ignore = T)))
        }
        if (!is.null(env)) 
            assign("genome.info", genome.info, envir = env)
        genome.info$operons <- NULL
        if ((operon.shift || "operons" %in% names(net.weights)) && 
            !no.genome.info) {
            tmp.operons <- try(get.operon.predictions("microbes.online"))
            if (class(tmp.operons) == "try-error") {
                message("Could not fetch operons file. Assuming it doesn't exist (eukaryote?)")
                set.param("is.eukaryotic", TRUE, override = T)
                set.param("operon.shift", FALSE, override = T)
                operon.shift[1:length(operon.shift)] <- FALSE
                if ("operons" %in% names(net.weights)) {
                  net.weights <- net.weights[names(net.weights) != 
                    "operons"]
                  set.param("net.weights", net.weights, override = T)
                }
            }
            else {
                genome.info$operons <- tmp.operons
            }
            rm(tmp.operons)
            if (!is.null(env)) 
                assign("genome.info", genome.info, envir = env)
        }
        if (!exists("ratios") || is.null(ratios)) {
            message("WARNING: No ratios matrix -- will generate an 'empty' one with all annotated ORFs for 'probes'.")
            if (!is.null(genome.info$gene.regex)) 
                rows <- unique(as.character(subset(genome.info$feature.names, 
                  grepl(paste("^", genome.info$gene.regex, sep = ""), 
                    names, ignore = T, perl = T), select = "names", 
                  drop = T)))
            else rows <- unique(as.character(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T)))
            ratios <- list(ratios = t(t(rep(NA, length(rows)))))
            rownames(ratios$ratios) <- rows
            attr(ratios, "rnames") <- sort(unique(rows))
            rm(rows)
            attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
            attr(ratios, "ncol") <- 1
            cat("Ratios: ", attr(ratios, "nrow"), "x", 1, "\n")
        }
        rm(nch, prefix, regex, tmp, qqq)
        if (!no.genome.info && length(mot.weights) > 0) {
            genome.info$all.upstream.seqs <- genome.info$bg.list <- list()
            genome.info$bg.fname <- character()
            for (i in names(mot.weights)) {
                cat("Pre-computing all '", i, "' seqs (", paste(motif.upstream.scan[[i]], 
                  collapse = ", "), ")...\n", sep = "")
                genome.info$all.upstream.seqs[[i]] <- get.sequences(genome.info$all.gene.names, 
                  seq.type = i, distance = motif.upstream.scan[[i]], 
                  filter = F)
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
                message(sum(!attr(ratios, "rnames") %in% names(genome.info$all.upstream.seqs[[i]])), 
                  " probes have no '", i, "' sequence.")
                if (!is.na(bg.order[i])) {
                  cat("Pre-computing '", i, "' residue bg distrib (order=", 
                    bg.order[i], ")...\n", sep = "")
                  tmp.seqs <- if (!is.null(genome.info$all.upstream.seqs[[i]])) 
                    genome.info$all.upstream.seqs[[i]]
                  else get.sequences(genome.info$all.gene.names, 
                    distance = motif.upstream.search[[i]], seq.type = i, 
                    filter = F)
                  genome.info$bg.fname[i] <- my.tempfile("meme.tmp", 
                    suf = ".bg")
                  capture.output(genome.info$bg.list[[i]] <- mkBgFile(tmp.seqs, 
                    order = bg.order[i], verbose = T, bgfname = genome.info$bg.fname[i], 
                    use.rev.comp = grepl("-revcomp", meme.cmd[i])))
                  rm(tmp.seqs)
                }
                else {
                  message("NOT USING a global sequence background distribution!")
                }
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
            }
        }
        networks <- list()
        if (!is.na(net.iters) && any(net.iters %in% 1:n.iter)) {
            if (length(grep("string", names(net.weights))) > 
                0) {
                if ("string" %in% names(net.weights)) {
                  if (exists("string.links")) {
                    string <- string.links
                  }
                  else {
                    cat("Loading STRING network.\n")
                    string <- get.STRING.links(genome.info$org.id$V1[1])
                  }
                  if (!is.null(string) && nrow(string) > 0) {
                    cat("Read in", nrow(string), "STRING edges; weight =", 
                      net.weights["string"], "\n")
                    string$combined_score <- string$combined_score/max(string$combined_score, 
                      na.rm = T) * 1000
                    string$combined_score <- 1000 * exp(string$combined_score/1000)/exp(1)
                    networks[["string"]] <- string
                  }
                  else {
                    warning("Could not load STRING network. Either", 
                      organism, "is not there or your gene names are not standard.")
                  }
                  rm(string)
                }
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if ("operons" %in% names(net.weights) && !is.null(genome.info$operons)) {
                cat("Converting operon predictions into a network...\n")
                tmp <- tapply(genome.info$operons$gene, genome.info$operons$head)
                names(tmp) <- genome.info$operons$gene
                mc <- get.parallel(length(unique(tmp)))
                out.sif <- do.call(rbind, mc$apply(unique(tmp), 
                  function(j) {
                    whch <- which(tmp == j)
                    gs <- names(whch)
                    if (length(gs) <= 1 || length(gs) > attr(ratios, 
                      "nrow")/20) 
                      return(NULL)
                    tmp.sif <- t(combn(gs, 2))
                    tmp.sif <- tmp.sif[tmp.sif[, 1] != tmp.sif[, 
                      2], , drop = F]
                    data.frame(protein1 = tmp.sif[, 1], protein2 = tmp.sif[, 
                      2], combined_score = rep(1000, nrow(tmp.sif)))
                  }))
                if (!is.null(out.sif) && nrow(out.sif) > 0) {
                  out.sif$combined_score <- rep(1000, nrow(out.sif))
                  colnames(out.sif) <- c("protein1", "protein2", 
                    "combined_score")
                  networks[["operons"]] <- out.sif
                }
                rm(tmp, mc, out.sif)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (exists("net.weights") && length(net.weights) > 
                0 && !is.null(names(net.weights))) {
                for (i in names(net.weights)) {
                  if (i %in% names(networks)) 
                    next
                  if (file.exists(i)) {
                    cat("Loading sif interactions from file:", 
                      i, "; weight =", net.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else if (exists(i) && !is.null(ncol(get(i))) && 
                    ncol(get(i)) >= 2) {
                    cat("Using network '", i, "' that exists in memory already; weight = ", 
                      net.weights[i], "\n", sep = "")
                    sif <- get(i)
                    if (ncol(sif) == 2) 
                      sif <- cbind(sif, rep(1, nrow(sif)))
                    colnames(sif) <- c("protein1", "protein2", 
                      "combined_score")
                  }
                  else {
                    next
                  }
                  networks[[basename(i)]] <- sif
                  rm(sif)
                }
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (exists("grouping.weights") && length(grouping.weights) > 
                0) {
                if (exists("net.weights")) 
                  net.weights <- c(net.weights, grouping.weights)
                else net.weights <- grouping.weights
                for (i in names(grouping.weights)) {
                  if (i %in% names(networks)) 
                    next
                  if (file.exists(i)) {
                    cat("Loading groupings from file:", i, "; weight =", 
                      grouping.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else if (exists(i) && !is.null(ncol(get(i))) && 
                    ncol(get(i)) >= 2) {
                    cat("Using groupings from '", i, "' that exists in memory already; weight = ", 
                      grouping.weights[i], "\n", sep = "")
                    sif <- get(i)
                    if (ncol(sif) == 2) 
                      sif <- cbind(sif, combined_score = rep(1, 
                        nrow(sif)))
                  }
                  colnames(sif) <- c("group", "protein", "combined_score")
                  if (sum(unique(as.character(sif$protein)) %in% 
                    attr(ratios, "rnames")) < sum(unique(as.character(sif$group)) %in% 
                    attr(ratios, "rnames"))) {
                    sif <- sif[, c(2, 1, 3)]
                    colnames(sif) <- c("group", "protein", "combined_score")
                  }
                  sif <- sif[order(sif$group), ]
                  tmp <- tapply(sif$protein, sif$group)
                  names(tmp) <- as.character(sif$protein)
                  cat("Converting", length(unique(tmp)), "groupings to a network (this may take a while for big grouping files)...")
                  mc <- get.parallel(length(unique(tmp)))
                  out.sif <- mc$apply(unique(tmp), function(j) {
                    whch <- which(tmp == j)
                    gs <- names(whch)
                    if (length(gs) <= 1 || length(gs) > attr(ratios, 
                      "nrow")/20) 
                      return(NULL)
                    ws <- sif$combined_score[whch]
                    names(ws) <- gs
                    tmp.sif <- t(combn(gs, 2))
                    tmp.sif <- tmp.sif[tmp.sif[, 1] != tmp.sif[, 
                      2], , drop = F]
                    tmp.sif <- data.frame(protein1 = tmp.sif[, 
                      1], protein2 = tmp.sif[, 2], combined_score = (ws[tmp.sif[, 
                      1]] + ws[tmp.sif[, 2]])/2)
                    rownames(tmp.sif) <- NULL
                    if (j%%100 == 0) 
                      cat(j)
                    cat(".")
                    tmp.sif
                  })
                  cat(length(unique(tmp)), "... ")
                  out.sif <- do.call(rbind, out.sif)
                  colnames(out.sif) <- c("protein1", "protein2", 
                    "combined_score")
                  networks[[basename(i)]] <- out.sif
                  cat("DONE\n")
                }
                rm(sif, tmp, out.sif, i, mc)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            for (n in names(networks)) {
                nn <- networks[[n]]
                if (nrow(nn) <= 0) {
                  message("WARNING: no edges in network", n, 
                    "... skipping.")
                  if (length(grep(n, seed.method[1])) > 0) {
                    message("ALSO, we have to change the row seeding method from", 
                      seed.method, "to 'kmeans'.")
                    seed.method["rows"] <- "kmeans"
                    set.param("seed.method", seed.method, override = T)
                  }
                  next
                }
                nodes <- unique(c(as.character(nn$protein1), 
                  as.character(nn$protein2)))
                cat(nrow(nn), "edges,", length(nodes), "nodes in network", 
                  n, "\n")
                nn <- subset(nn, as.character(protein1) != as.character(protein2))
                dupes <- duplicated(nn[, c("protein1", "protein2")])
                if (sum(dupes) > 0) {
                  cat("Merging", sum(dupes), "duplicate edges in network", 
                    n, "; this could take a while for networks with lots of nodes...\n")
                  tmp.nn <- subset(nn, dupes)
                  dupe.nodes <- unique(c(as.character(tmp.nn$protein1), 
                    as.character(tmp.nn$protein2)))
                  if (length(dupe.nodes) < 6000) {
                    tmp <- tapply(tmp.nn$combined_score, tmp.nn[, 
                      c("protein1", "protein2")], sum, na.rm = T)
                    tmp2 <- which(!is.na(tmp), arr = T)
                    nn.new <- data.frame(protein1 = rownames(tmp)[tmp2[, 
                      1]], protein2 = colnames(tmp)[tmp2[, 2]], 
                      combined_score = tmp[tmp2])
                    rm(tmp, tmp2)
                    nn <- rbind(nn.new, nn)
                    rm(nn.new)
                    nn <- nn[!duplicated(nn[, c("protein1", "protein2")]), 
                      ]
                  }
                  rm(tmp.nn, dupe.nodes)
                }
                if (exists("ratios") && !is.null(ratios) && !any(nodes %in% 
                  attr(ratios, "rnames"))) {
                  if (median(nchar(nodes)) > median(nchar(attr(ratios, 
                    "rnames"))) && any(substr(nodes, 1, median(nchar(attr(ratios, 
                    "rnames")))) %in% attr(ratios, "rnames"))) {
                    nn$protein1 <- substr(as.character(nn$protein1), 
                      1, median(nchar(attr(ratios, "rnames"))))
                    nn$protein2 <- substr(as.character(nn$protein2), 
                      1, median(nchar(attr(ratios, "rnames"))))
                    nodes <- unique(c(as.character(nn$protein1), 
                      as.character(nn$protein2)))
                  }
                  if (!is.null(genome.info$synonyms)) {
                    rr <- attr(ratios, "rnames")[!attr(ratios, 
                      "rnames") %in% nodes]
                    if (length(rr) > 0) {
                      cat("Reconciling network", n, length(rr), 
                        "node names with probe names...\n")
                      syns <- get.synonyms(rr)
                      mc <- get.parallel(length(syns))
                      is.there <- unlist(mc$apply(syns, function(i) any(i %in% 
                        nodes)))
                      syns <- syns[is.there]
                      nnc1 <- as.character(nn$protein1)
                      nnc2 <- as.character(nn$protein2)
                      nnc1.t <- !nnc1 %in% attr(ratios, "rnames")
                      nnc2.t <- !nnc2 %in% attr(ratios, "rnames")
                      mc <- get.parallel(2)
                      tmp <- mc$apply(1:2, function(ii) {
                        for (i in names(syns)) {
                          if (ii == 1) 
                            nnc1[nnc1.t & nnc1 %in% syns[[i]]] <- i
                          else nnc2[nnc2.t & nnc2 %in% syns[[i]]] <- i
                        }
                        if (ii == 1) 
                          return(nnc1)
                        else return(nnc2)
                      })
                      nnc1 <- tmp[[1]]
                      nnc2 <- tmp[[2]]
                      rm(tmp, nnc1.t, nnc2.t)
                      cat(sum(!is.there), "probes have no nodes in", 
                        n, "network (but", sum(attr(ratios, "rnames") %in% 
                          nodes, na.rm = T) + sum(is.there), 
                        "do)\n")
                      nn$protein1 <- nnc1
                      nn$protein2 <- nnc2
                      tmp <- nnc1 %in% attr(ratios, "rnames") & 
                        nnc2 %in% attr(ratios, "rnames")
                      nn <- subset(nn, tmp == TRUE)
                      rm(tmp, syns, is.there, nnc1, nnc2, nnc1.t, 
                        nnc2.t, tmp, rr, i)
                    }
                  }
                }
                else {
                  cat(sum(!attr(ratios, "rnames") %in% nodes), 
                    "probes have no nodes in", n, "network (but", 
                    sum(attr(ratios, "rnames") %in% nodes, na.rm = T), 
                    "do)\n")
                }
                ttmp <- nn[, c(2, 1, 3)]
                colnames(ttmp) <- colnames(nn)
                nn <- rbind(nn, ttmp)
                rm(ttmp)
                nn <- nn[!duplicated(nn[, c("protein1", "protein2")]), 
                  ]
                cat(n, "network filtered, symmetrized and uniquified:", 
                  nrow(nn), "edges.\n")
                networks[[n]] <- nn
                if (!is.null(env)) 
                  assign("networks", networks, envir = env)
            }
            rm(n, nn, nodes, dupes)
            if (length(networks) > 1) {
                sums <- sapply(networks, function(n) sum(n$combined_score, 
                  na.rm = T))
                ms <- max(sums[sums > 0], na.rm = T)
                if (length(sums) > 0 && !is.na(ms)) 
                  for (n in names(networks)) if (sums[n] > 0) 
                    networks[[n]]$combined_score <- networks[[n]]$combined_score/sums[n] * 
                      ms
                rm(n, sums, ms)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (!is.null(names(net.weights))) 
                names(net.weights) <- basename(names(net.weights))
        }
        if (!no.genome.info && cog.org != "" && cog.org != "?" && 
            !is.null(cog.org) && !all(is.na(plot.iters) | plot.iters == 
            0)) {
            cat("Loading COG functional codes (for plotting), org. code", 
                cog.org, ": trying NCBI whog file...\n")
            genome.info$cog.code <- get.COG.code(cog.org)
            cat(sum(!is.na(genome.info$cog.code)), "genes have a COG code (", 
                if (is.null(genome.info$cog.code)) 
                  attr(ratios, "nrow")
                else sum(is.na(genome.info$cog.code)), "do not)\n")
        }
    }
    iter <- 0
    meme.scores <- clusterStack <- list()
    for (i in names(mot.weights)) {
        meme.scores[[i]] <- list()
        meme.scores[[i]][[k.clust + 1]] <- ""
    }
    stats <- row.scores <- col.scores <- mot.scores <- net.scores <- NULL
    if (!exists("favorite.cluster")) 
        favorite.cluster <- function() min(which(sapply(1:k.clust, 
            function(k) length(get.rows(k))) > cluster.rows.allowed[1] * 
            2))
    row.scaling <- extend.vec(row.scaling)
    mot.scaling <- extend.vec(mot.scaling)
    net.scaling <- extend.vec(net.scaling)
    n.motifs <- lapply(n.motifs, extend.vec)
    fuzzy.index <- extend.vec(fuzzy.index)
    is.inited <- TRUE
    if (is.null(env)) 
        env <- new.env(hash = T, parent = globalenv())
    else parent.env(env) <- globalenv()
    parent.env(cmonkey.params) <- env
    attr(env, "class") <- c("environment", "cmonkey")
    for (i in ls()) {
        if (i %in% c("i", "env")) 
            next
        tmp <- get(i)
        if (class(tmp) == "function") 
            environment(tmp) <- env
        assign(i, tmp, envir = env)
    }
    if (exists("favorite.cluster")) 
        env$favorite.cluster <- favorite.cluster
    environment(env$favorite.cluster) <- env
    if (exists("cm.func.each.iter")) {
        env$cm.func.each.iter <- cm.func.each.iter
        environment(env$cm.func.each.iter) <- env
        try(env$cm.func.each.iter(), silent = T)
    }
    cat("INITIALIZATION IS COMPLETE.\n")
    env$iter <- env$iter + 1
    invisible(env)
}
cmonkey.one.iter <-
function (env, dont.update = F, ...) 
{
    env <- env$update.all.clusters(env, dont.update = F, ...)
    row.memb <- sapply(1:k.clust, function(k) attr(ratios, "rnames") %in% 
        get.rows(k))
    if (is.vector(row.memb)) 
        row.memb <- t(row.memb)
    rownames(row.memb) <- attr(ratios, "rnames")
    col.memb <- sapply(1:k.clust, function(k) attr(ratios, "cnames") %in% 
        get.cols(k))
    if (is.vector(col.memb)) 
        col.memb <- t(col.memb)
    rownames(col.memb) <- attr(ratios, "cnames")
    if (iter %in% stats.iters) {
        tmp.stats <- env$get.stats()
        stats <- env$stats
        if (iter > 1 && nrow(stats) > 0 && any(!colnames(tmp.stats) %in% 
            colnames(stats))) {
            for (cn in colnames(tmp.stats)[!colnames(tmp.stats) %in% 
                colnames(stats)]) stats[[cn]] <- rep(NA, nrow(stats))
        }
        env$stats <- rbind(stats, tmp.stats)
        cat(organism, as.matrix(env$stats[nrow(env$stats), ]), 
            "\n")
    }
    else {
        cat(sprintf("==> %04d %.3f %.3f %.3f\n", iter, mean(env$row.scores[, 
            ][row.memb], na.rm = T), if (!is.null(env$mot.scores)) 
            mean(env$mot.scores[, ][row.memb & env$mot.scores[, 
                ] < 0], na.rm = T, trim = 0.05)
        else NA, if (!is.null(env$net.scores)) 
            mean(env$net.scores[, ][row.memb], na.rm = T, trim = 0.05)
        else NA))
    }
    if (!is.na(plot.iters) && iter %in% plot.iters) {
        env$plotStats(iter, plot.clust = env$favorite.cluster(), 
            new.dev = T)
    }
    if (exists("cm.func.each.iter")) 
        try(cm.func.each.iter(), silent = T)
    if (any(cm.script.each.iter != "")) {
        for (f in cm.script.each.iter) {
            if (file.exists(f) && file.info(f)$size > 1) {
                tmp <- readLines(f)
                if (all(substr(tmp, 1, 1) == "#")) 
                  next
                if (tmp[1] != "## QUIET") 
                  cat("Sourcing the script '", f, "' ...\n", 
                    sep = "")
                try(source(f, echo = tmp[1] != "## QUIET", local = T), 
                  silent = T)
                rm(tmp)
            }
        }
    }
    if (get.parallel()$mc) {
        if (getDoParName() == "doMC") {
            chld <- multicore::children()
            if (length(chld) > 0) {
                try({
                  multicore::kill(chld)
                  tmp <- multicore::collect(chld)
                }, silent = T)
            }
        }
        else if (getDoParName() == "doSNOW" && "data" %in% ls(pos = foreach:::.foreachGlobals)) {
            cl <- get("data", pos = foreach:::.foreachGlobals)
            if (!is.null(data)) 
                stopCluster(cl)
        }
    }
    if (!dont.update) 
        env$iter <- env$iter + 1
    invisible(env)
}
cmonkey.re.seed <-
function (env) 
{
    if (!exists("rnd.seed", envir = env$cmonkey.params)) {
        op <- options(digits.secs = 10)
        tmp.time <- as.character(Sys.time())
        options(op)
        rm(op)
        tmp.rnd.seed <- as.integer(substr(gsub("[-:. ]", "", 
            tmp.time), 12, 20))
        cat("RESETTING RANDOM SEED: ")
        env$set.param("date.run", tmp.time, env$cmonkey.params)
        env$date.run <- env$cmonkey.params$date.run
        env$set.param("rnd.seed", tmp.rnd.seed, env$cmonkey.params)
        env$rnd.seed <- env$cmonkey.params$rnd.seed
        set.seed(env$rnd.seed)
        rm(tmp.rnd.seed)
    }
    if (!is.null(env$ratios) && attr(env$ratios, "ncol") > 1) {
        cat("Seeding all clusters using methods:", env$seed.method, 
            "\n")
        tmp <- env$seed.clusters(env$k.clust, seed.method = env$seed.method["rows"], 
            col.method = env$seed.method["cols"])
    }
    else {
        cat("Seeding all clusters using methods: rnd rnd\n")
        tmp <- env$seed.clusters(env$k.clust, seed.method = "rnd", 
            col.method = "rnd")
    }
    env$clusterStack <- lapply(1:env$k.clust, function(k) list(rows = rownames(which(tmp$rm == 
        k, arr = T)), cols = rownames(which(tmp$cm == k, arr = T))))
    attr(env$clusterStack, "iter") <- env$iter - 1
    invisible(env)
}
cm.version <-
"4.9.11"
col.let <-
c("A", "C", "G", "T")
DEBUG <-
function (...) 
{
}
dlf <-
function (f, url, msg = NULL, mode = "wb", quiet = F, ...) 
{
    err <- 0
    if (mode == "ab" || !file.exists(f) || file.info(f)$size == 
        0) {
        if (!file.exists(dirname(f))) 
            try(dir.create(dirname(f), recursive = T))
        if (!is.null(msg)) 
            cat(msg, "\n")
        err <- try(download.file(url, destfile = f, mode = mode, 
            quiet = quiet, ...))
    }
    closeAllConnections()
    err
}
ensemble.to.database <-
function (out) 
{
    dir.create("filehash/DATABASES/")
    tmp <- data.table(round(out$e$ratios$ratios * 1000)/1000)
    tmp <- data.table(gene = rownames(out$e$ratios$ratios), tmp)
    write.table(tmp, file = "filehash/DATABASES/ratios.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    rm(tmp)
    system("bzip2 -fv filehash/DATABASES/ratios.tsv", wait = F)
    tmp <- data.table(out$e$genome.info$feature.tab)
    write.table(tmp, file = "filehash/DATABASES/feature_tab.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    rm(tmp)
    system("bzip2 -fv filehash/DATABASES/feature_tab.tsv", wait = F)
    tmp <- data.table(out$e$genome.info$feature.names)
    write.table(tmp, file = "filehash/DATABASES/feature_names.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    rm(tmp)
    system("bzip2 -fv filehash/DATABASES/feature_names.tsv", 
        wait = F)
    tmp <- data.table(out$e$genome.info$operons)
    write.table(tmp, file = "filehash/DATABASES/operons.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    rm(tmp)
    system("bzip2 -fv filehash/DATABASES/operons.tsv", wait = F)
    tmp <- data.table(gene = names(out$e$genome.info$all.upstream.seqs[[1]]), 
        sequence = out$e$genome.info$all.upstream.seqs[[1]])
    write.table(tmp, file = "filehash/DATABASES/upstream_seqs.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    rm(tmp)
    system("bzip2 -fv filehash/DATABASES/upstream_seqs.tsv", 
        wait = F)
    if (!is.null(out$e$networks$operons)) {
        tmp <- data.table(out$e$networks$operons)
        write.table(tmp, file = "filehash/DATABASES/operon_network.tsv", 
            quote = F, sep = "\t", row.names = F, col.names = T, 
            append = F)
        rm(tmp)
        system("bzip2 -fv filehash/DATABASES/operon_network.tsv", 
            wait = F)
    }
    if (!is.null(out$e$networks$string)) {
        tmp <- data.table(out$e$networks$string)
        write.table(tmp, file = "filehash/DATABASES/string_network.tsv", 
            quote = F, sep = "\t", row.names = F, col.names = T, 
            append = F)
        rm(tmp)
        system("bzip2 -fv filehash/DATABASES/string_network.tsv", 
            wait = F)
    }
    bzcon1 <- "filehash/DATABASES/bicluster.tsv"
    bzcon2 <- "filehash/DATABASES/bicluster_genes.tsv"
    bzcon3 <- "filehash/DATABASES/bicluster_conds.tsv"
    bzcon4 <- "filehash/DATABASES/bicluster_motifs.tsv"
    wrote <- FALSE
    for (k in 1:out$e$k.clust) {
        if (k%%100 == 0) 
            print(k)
        tab1 <- tab2 <- tab3 <- tab4 <- NULL
        bb <- paste("BIC", k, sep = "_")
        clust <- out$get.bicluster.info(bb)[[1]]
        if (!is.null(clust)) {
            fname <- names(out$e$fnames.to.cluster[which(out$e$fnames.to.cluster == 
                k)])
            k.orig <- clust$k
            tab1 <- data.table(bic = k, nrow = clust$nrow, ncol = clust$ncol, 
                resid = clust$resid, pclust = clust$p.clust, 
                eval = min(clust$e.val, na.rm = T), k_orig = k.orig, 
                fname = fname)
            if (!is.null(tab1)) 
                write.table(tab1, bzcon1, quote = F, sep = "\t", 
                  row.names = F, col.names = !wrote, append = wrote)
        }
        genes <- clust$rows
        if (!is.null(genes)) 
            tab2 <- data.table(bic = k, gene = genes)
        if (!is.null(tab2)) 
            write.table(tab2, bzcon2, quote = F, sep = "\t", 
                row.names = F, col.names = !wrote, append = wrote)
        conds <- clust$cols
        if (!is.null(conds)) 
            tab3 <- data.table(bic = k, cond = conds)
        if (!is.null(tab3)) 
            write.table(tab3, bzcon3, quote = F, sep = "\t", 
                row.names = F, col.names = !wrote, append = wrote)
        mots <- out$get.motifs(bicluster = bb)[[1]]
        if (!is.null(mots)) 
            tab4 <- data.table(bic = k, mot = gsub("MOT_", "", 
                mots))
        if (!is.null(tab4)) 
            write.table(tab4, bzcon4, quote = F, sep = "\t", 
                row.names = F, col.names = !wrote, append = wrote)
        wrote <- TRUE
    }
    system(sprintf("bzip2 -fv -9 %s", bzcon1), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon2), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon3), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon4), wait = F)
    bzcon1 <- "filehash/DATABASES/motif.tsv"
    bzcon2 <- "filehash/DATABASES/motif_meme_posn.tsv"
    bzcon3 <- "filehash/DATABASES/motif_mast_posn.tsv"
    bzcon4 <- "filehash/DATABASES/motif_pssm.tsv"
    wrote <- FALSE
    for (k in 1:out$e$k.clust) {
        bic <- paste("BIC", k, sep = "_")
        mots <- out$get.motifs(bicluster = bic)[[1]]
        for (m in mots) {
            tab1 <- NULL
            if (k%%100 == 0) 
                print(m)
            minf <- out$get.motif.info(m)[[1]]
            mm <- gsub("MOT_", "", m)
            if (!is.null(minf)) {
                cf <- out$coding.fracs$all.fracs[m]
                tab <- data.table(mot = mm, width = minf$width, 
                  llr = minf$llr, eval = minf$e.value, sites = minf$sites, 
                  coding = cf, good = (cf < out$coding.fracs$mean.fracs - 
                    0.01))
                if (!is.null(tab)) 
                  write.table(tab, bzcon1, quote = F, sep = "\t", 
                    row.names = F, col.names = !wrote, append = wrote)
                tab <- NULL
                if (nrow(minf$posns) > 0) 
                  tab <- as.data.table(cbind(mot = mm, minf$posns))
                if (!is.null(tab)) 
                  write.table(tab, bzcon2, quote = F, sep = "\t", 
                    row.names = F, col.names = !wrote, append = wrote)
                tab <- NULL
                if (nrow(minf$mast) > 0) 
                  tab <- as.data.table(cbind(mot = mm, minf$mast))
                if (!is.null(tab)) 
                  write.table(tab, bzcon3, quote = F, sep = "\t", 
                    row.names = F, col.names = !wrote, append = wrote)
                pssm <- minf$pssm
                colnames(pssm) <- out$e$col.let
                pssm <- as.data.table(pssm)
                tab <- data.table(cbind(mot = mm, ind = 1:nrow(pssm), 
                  round(pssm * 1000)/1000))
                if (!is.null(tab)) 
                  write.table(tab, bzcon4, quote = F, sep = "\t", 
                    row.names = F, col.names = !wrote, append = wrote)
                wrote <- TRUE
            }
        }
    }
    system(sprintf("bzip2 -fv -9 %s", bzcon1), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon2), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon3), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon4), wait = F)
    bzcon1 <- "filehash/DATABASES/motif_clust.tsv"
    bzcon2 <- "filehash/DATABASES/motif_clust_combined_pssm.tsv"
    bzcon3 <- "filehash/DATABASES/motif_clust_aligned_pssms.tsv"
    wrote <- FALSE
    for (k in 1:out$mc.length) {
        if (k%%100 == 0) 
            print(k)
        mc <- paste("MOTC", k, sep = "_")
        mots <- out$get.motifs(motif.clust = mc)[[1]]
        tab <- data.table(motc = k, mot = gsub("MOT_", "", mots))
        if (!is.null(tab)) 
            write.table(tab, bzcon1, quote = F, sep = "\t", row.names = F, 
                col.names = !wrote, append = wrote)
        mc.info <- out$get.motif.cluster.info(mc)[[1]]
        pssm <- attr(mc.info, "combined.pssm")
        colnames(pssm) <- out$e$col.let
        pssm <- as.data.table(pssm)
        tab <- data.table(cbind(motc = k, ind = 1:nrow(pssm), 
            round(pssm * 1000)/1000))
        if (!is.null(tab)) 
            write.table(tab, bzcon2, quote = F, sep = "\t", row.names = F, 
                col.names = !wrote, append = wrote)
        aligned.pssms <- attr(mc.info, "aligned.pssms")
        wrote2 <- wrote
        for (m in names(aligned.pssms)) {
            pssm <- aligned.pssms[[m]]
            colnames(pssm) <- out$e$col.let
            pssm <- as.data.table(pssm)
            tab <- data.table(cbind(motc = k, motif = gsub("MOT_", 
                "", m), ind = 1:nrow(pssm), round(pssm * 1000)/1000))
            if (!is.null(tab)) 
                write.table(tab, bzcon3, quote = F, sep = "\t", 
                  row.names = F, col.names = !wrote2, append = wrote2)
            wrote2 <- TRUE
        }
        wrote <- TRUE
    }
    system(sprintf("bzip2 -fv -9 %s", bzcon1), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon2), wait = F)
    system(sprintf("bzip2 -fv -9 %s", bzcon3), wait = F)
    tmp <- out$pssm.scans
    setnames(tmp, "gene", "p.value", "posn", "mot", "bic")
    tmp$strand <- ifelse(tmp$mot > 0, "+", "-")
    tmp$mot <- abs(tmp$mot)
    write.table(tmp, file = "filehash/DATABASES/motif_scans.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    system(sprintf("bzip2 -fv -9 %s", "filehash/DATABASES/motif_scans.tsv"), 
        wait = F)
    rm(tmp)
    gc()
    load("filehash/fimo_out_1e-06.RData")
    write.table(fimo.out, file = "filehash/DATABASES/motif_scans_fimo.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    system(sprintf("bzip2 -fv -9 %s", "filehash/DATABASES/motif_scans_fimo.tsv"), 
        wait = F)
    rm(fimo.out)
    gc()
    write.table(out$meme.hits, file = "filehash/DATABASES/motif_meme_hits.tsv", 
        quote = F, sep = "\t", row.names = F, col.names = T, 
        append = F)
    system(sprintf("bzip2 -fv -9 %s", "filehash/DATABASES/motif_meme_hits.tsv"), 
        wait = F)
    genomes <- lapply(names(out$e$genome.info$genome.seqs), function(i) {
        nc <- nchar(out$e$genome.info$genome.seqs[i])
        output <- c(1, nc)
        names(output) <- c(i, "")
        output
    })
    options(cores = 4)
    bzcon1 <- "filehash/DATABASES/motif_promoter_counts.tsv"
    wrote <- FALSE
    for (mc in 1:out$mc.length) {
        for (i in genomes) {
            mots <- out$motif.clusts[[mc]]
            tmp <- out$plot.promoter.architecture(i, type = "fimo", 
                dont.plot = T, motif.filter = mots, include.bad = T)
            not.zeroes <- apply(cbind(tmp$counts[, 1], tmp$mat), 
                1, function(i) any(i != 0))
            tab <- data.table(motc = mc, loc = as.integer(names(tmp$counts[not.zeroes, 
                ])), count = tmp$counts[not.zeroes, 1], round(tmp$mat[not.zeroes, 
                ] * 1000)/1000)
            cat(mc, nrow(tab), "\n")
            if (!is.null(tab)) 
                write.table(tab, bzcon1, quote = F, sep = "\t", 
                  row.names = F, col.names = !wrote, append = wrote)
            wrote <- TRUE
        }
    }
    tmp <- out$plot.promoter.architecture(i, type = "fimo", dont.plot = T, 
        include.bad = F, count.all = T)
    not.zeroes <- apply(cbind(tmp$counts[, 1], tmp$mat), 1, function(i) any(i != 
        0))
    tab <- data.table(motc = 0, loc = as.integer(names(tmp$counts[not.zeroes, 
        1])), count = tmp$counts[not.zeroes, 1], round(tmp$mat[not.zeroes, 
        ] * 1000)/1000)
    write.table(tab, bzcon1, quote = F, sep = "\t", row.names = F, 
        col.names = F, append = T)
    system(sprintf("bzip2 -fv -9 %s", bzcon1), wait = F)
}
evaluate.all.clusterings <-
function (...) 
{
    load("filehash/fimo_out_1e-06.RData")
    load("filehash/new_motif_shadows_1e-06.RData")
    retry.new.clustering <- function(mcl.I = 1.2, mcl.pi = 2, 
        distance.weight.cutoff = 0.999999, mcl.cmd = paste("./progs/mcl-10-201/local/bin/mcl new.mcltemp -o %s --abc", 
            " -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1"), 
        meme.the.clusters = F, dirname = "filehash/new_motif_shadows_1e-06_52977") {
        outfile <- sprintf("filehash/new.mcltemp.I%s.pi%s.dcut%s", 
            gsub(".", "", sprintf("%.1f", mcl.I), fixed = T), 
            gsub(".", "", sprintf("%.1f", mcl.pi), fixed = T), 
            sprintf("%.3f", distance.weight.cutoff))
        system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
            dirname, distance.weight.cutoff))
        system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
            dirname, distance.weight.cutoff))
        mcl.cmd <- sprintf("%s  2>&1", mcl.cmd)
        mcl.cmd <- sprintf(mcl.cmd, outfile, mcl.I, mcl.pi)
        print(mcl.cmd)
        mcl.out <- system(mcl.cmd, intern = T, ignore.stderr = F)
        system(sprintf("gzip -fv %s", outfile))
        clusts <- lapply(strsplit(readLines(gzfile(sprintf("%s.gz", 
            outfile))), "\t"), as.integer)
        cat("GOT", length(clusts), "motif clusters with", length(unlist(clusts)), 
            "motifs.\n")
        clusts <- clusts[sapply(clusts, length) >= 3]
        cat("GOT", sum(sapply(clusts, length) >= 3), "motif clusters (length > 3) with", 
            length(unlist(clusts[sapply(clusts, length) >= 3])), 
            "motifs.\n")
        cat("GOT", sum(sapply(clusts, length) >= 10), "motif clusters (length > 10) with", 
            length(unlist(clusts[sapply(clusts, length) >= 10])), 
            "motifs.\n")
        attr(clusts, "mcl.cmd") <- mcl.cmd
        attr(clusts, "mcl.out") <- mcl.out
        clusts
    }
    out$retry.new.clustering <- retry.new.clustering
    environment(out$retry.new.clustering) <- out
    tmpl <- list()
    try(load("tmpl.RData"))
    sys.source("~/scratch/biclust/cmonkey-ensemble-funcs.R", 
        envir = out, chdir = T)
    for (mcl.I in rev(c(1.2, 1.5, 2, 2.5, 3, 4.5, 6))) {
        for (mcl.pi in c(1, 2, 5, 10, 20)) {
            for (dc in c(0.8, 0.9, 0.95, 0.98, 0.99, 0.999, 0.999999)) {
                if (paste(mcl.I, mcl.pi, dc) %in% names(tmpl)) 
                  next
                print(paste(mcl.I, mcl.pi, dc))
                tmp2 <- capture.output(tmp <- out$retry.new.clustering(mcl.I = mcl.I, 
                  mcl.pi = mcl.pi, distance.weight.cutoff = dc, 
                  ...))
                attr(tmp, "captured.output") <- tmp2
                tmpl[[paste(mcl.I, mcl.pi, dc)]] <- tmp
                print(paste(mcl.I, mcl.pi, dc))
                print(attr(tmp, "mcl.cmd"))
                print(attr(tmp, "mcl.out"))
                print(tmp2)
                print(date())
                save(tmpl, file = "tmpl.RData")
            }
        }
    }
    save(tmpl, file = "tmpl.RData")
    motifs <- unlist(lapply(1:nrow(out$motif.widths), function(i) {
        ii <- which(out$motif.widths[i, ] > 0)
        if (length(ii) <= 0) 
            return(NULL)
        paste("MOT", i, ii, sep = "_")
    }))
    frac.in.coding <- out$coding.fracs$all.fracs[motifs]
    motifs <- motifs[!is.na(frac.in.coding) & frac.in.coding < 
        out$coding.fracs$mean.fracs - 0.01]
    require(data.table)
    df = data.table(read.delim(bzfile(sprintf("filehash/new_motif_shadows_1e-06_%d.tsv.bz2", 
        length(motifs))), head = F))
    setkey(df, V1, V2)
    bad.inds <- integer()
    for (i in 1:length(motifs)) {
        if (!file.exists(sprintf("filehash/new_motif_shadows_1e-06_%d/%08d.tsv.bz2", 
            length(motifs), i))) {
            bad.inds = c(bad.inds, i)
            next
        }
        tmp = read.delim(bzfile(sprintf("filehash/new_motif_shadows_1e-06_%d/%08d.tsv.bz2", 
            length(motifs), i)), head = F)
        if (all(tmp$V3 == 0)) {
            bad.inds = c(bad.inds, i)
            print(i)
        }
    }
    cat(length(bad.inds), "motifs with no hits to genome...?\n")
    df = df[!(V1 %in% bad.inds) & !(V2 %in% bad.inds)]
    rm(frac.in.coding)
    gc()
    require(Matrix)
    m = Matrix(0, nrow = length(motifs), ncol = length(motifs))
    m = forceSymmetric(m)
    m[cbind(df$V1, df$V2)] <- 1 - df$V3
    m[cbind(df$V2, df$V1)] <- 1 - df$V3
    save(m, df, file = "m_and_df.RData")
    max.len <- 400
    m.inds = 1:length(motifs)
    dens.out = list()
    for (i in names(tmpl)) {
        print(i)
        if (!is.null(dens.out[[i]])) 
            next
        ii <- i
        i = tmpl[[i]]
        mc.len = sum(sapply(i, length) >= 10)
        if (mc.len > max.len) 
            mc.len = max.len
        cluster.densities = mclapply(i[1:mc.len], function(j) {
            mm = m[j, j]
            d = c(sum(mm), nrow(mm))
        })
        between.densities = mclapply(i[1:mc.len], function(j) {
            not.j = m.inds[m.inds %in% unlist(i[1:mc.len]) & 
                !m.inds %in% j]
            mm = m[j, not.j]
            d = c(sum(mm), nrow(mm), ncol(mm))
        })
        out.density = {
            non.hits = m.inds[!m.inds %in% unlist(i[1:mc.len])]
            mm = m[non.hits, non.hits]
            d = c(sum(mm), nrow(mm))
        }
        dens.out[[ii]] = list(cluster.dens = cluster.densities, 
            between.dens = between.densities, out.dens = out.density)
        gc()
    }
    save(dens.out, file = "dens.out.RData")
    tmp2 = t(sapply(names(dens.out), function(i) {
        cl.size = sapply(dens.out[[i]]$cluster.dens, function(j) j[2])
        max.ind = if (length(cl.size) > 400) 
            400
        else length(cl.size)
        in.dens = sapply(dens.out[[i]]$cluster.dens[1:max.ind], 
            function(j) j[1]/j[2]/j[2])
        betw.dens = sapply(dens.out[[i]]$between.dens[1:max.ind], 
            function(j) j[1]/j[2]/j[3])
        in.dens[in.dens == 0 & betw.dens == 0] = betw.dens[in.dens == 
            0 & betw.dens == 0] = NA
        out.dens = dens.out[[i]]$out.dens[1]/dens.out[[i]]$out.dens[1]/dens.out[[i]]$out.dens[1]
        c(sum(cl.size[1:max.ind]), max.ind, weighted.mean(in.dens, 
            cl.size[1:max.ind], na.rm = T), weighted.mean(betw.dens, 
            cl.size[1:max.ind], na.rm = T), out.dens)
    }))
    tmp <- t(sapply(strsplit(names(tmpl), " "), as.numeric))
    dc <- table(tmp[, 3])
    dc[] <- 1:length(dc)
    pi <- table(tmp[, 2])
    pi[] <- 0:(length(pi) - 1)
    mcl.I <- table(tmp[, 1])
    mcl.I[] <- 1:length(mcl.I)
    mc.len <- sapply(tmpl, function(i) sum(sapply(i, length) >= 
        10))
    par(mfrow = c(3, 3))
    plot(tmp2[, 1], log10(tmp2[, 3]), xlab = "cl.size", ylab = "log10(in.dens)", 
        cex = 0.7, col = dc[as.character(tmp[, 3])], pch = pi[as.character(tmp[, 
            2])])
    legend("bottomleft", legend = names(dc), pch = 20, cex = 0.7, 
        col = 1:length(dc), title = "distance.cut")
    plot(tmp2[, 1], log10(tmp2[, 4]), xlab = "cl.size", ylab = "log10(betw.dens)", 
        cex = 0.7, col = dc[as.character(tmp[, 3])], pch = pi[as.character(tmp[, 
            2])])
    legend("bottomleft", legend = names(pi), cex = 0.7, pch = 0:(length(pi) - 
        1), title = "mcl.PI")
    plot(tmp2[, 1], log10(tmp2[, 5]), xlab = "cl.size", ylab = "log10(out.dens)", 
        cex = 0.7, col = dc[as.character(tmp[, 3])], pch = pi[as.character(tmp[, 
            2])])
    plot(tmp2[, 2], log10(tmp2[, 3]/tmp2[, 4]), xlab = "mc.len", 
        ylab = "log10(in.dens/betw.dens)", cex = 0.7, col = dc[as.character(tmp[, 
            3])], pch = pi[as.character(tmp[, 2])])
    plot(tmp2[, 1], log10(tmp2[, 3]/tmp2[, 4]), xlab = "cl.size", 
        ylab = "log10(in.dens/betw.dens)", cex = 0.7, col = dc[as.character(tmp[, 
            3])], pch = pi[as.character(tmp[, 2])])
    plot(tmp2[, 1], log10(tmp2[, 3]/tmp2[, 4]), xlab = "cl.size", 
        ylab = "log10(in.dens/betw.dens)", cex = 0.7, col = mcl.I[as.character(tmp[, 
            1])], pch = pi[as.character(tmp[, 2])])
    legend("bottomleft", legend = names(mcl.I), pch = 20, cex = 0.7, 
        col = 1:length(mcl.I), title = "mcl.I")
    plot(mc.len, log10(tmp2[, 3]/tmp2[, 4]), xlab = "mc.len", 
        ylab = "log10(in.dens/betw.dens)", cex = 0.7, col = mcl.I[as.character(tmp[, 
            1])], pch = pi[as.character(tmp[, 2])])
    graphics.to.pdf("evaluate_all_clusterings.pdf")
    clusts <- tmpl[[which.max(log10(tmp2[, 3]/tmp2[, 4]))]]
    for (i in 1:length(clusts)) clusts[[i]] <- motifs[clusts[[i]]]
    mc.length <- sum(sapply(clusts, length) >= 10)
    for (i in 1:mc.length) {
        tmp <- out$cluster.motifs("hclust", motifs = clusts[[i]], 
            min.gene.overlap = 1, e.value.cutoff = Inf, p.value.cutoff = 1e-04, 
            resid.cutoff = Inf, n.cutoff = 1, expand = F, include.bad = T, 
            find.bad = NA, in.tt.out = NULL, k.cut = 1)
        attr(clusts[[i]], "tt.out") <- tmp$tt.out
        attr(clusts[[i]], "tt.out2") <- tmp$tt.out2[[1]]
    }
    pdf("evaluate_all_clusterings_clusters.pdf")
    par(mfrow = c(4, 4))
    for (i in 1:mc.length) {
        pssm <- attr(attr(clusts[[i]], "tt.out2"), "combined.pssm")
        e$viewPssm(pssm, main = paste(i, length(clusts[[i]])))
    }
    dev.off()
    attr(clusts, "mc.length") <- mc.length
    attr(clusts, "motifs") <- motifs
}
extend.vec <-
function (v, n = n.iter) 
{
    if (length(v) < n) 
        v <- c(v, rep(v[length(v)], n.iter - length(v)))
    v
}
filter.sequences <-
function (seqs, start.stops = NULL, seq.type = paste(c("upstream", 
    "upstream.noncod", "upstream.noncod.same.strand", "downstream", 
    "gene")[1], "meme"), distance = motif.upstream.search[[seq.type]], 
    uniquify = T, remove.repeats = T, remove.atgs = T, mask.overlapping.rgns = F, 
    verbose = F, ...) 
{
    if (uniquify) 
        seqs <- seqs[!get.dup.seqs(seqs)]
    if (remove.repeats && length(grep("NNNNNN", seqs)) <= 1) {
        if (verbose) 
            cat("Removing low-complexity regions from sequences.\n")
        seqs.new <- remove.low.complexity(seqs, seq.type = seq.type)
        if (length(seqs.new) == length(seqs)) 
            seqs <- seqs.new
        else warning("Remove low complexity failed - skipping!")
        rm(seqs.new)
    }
    if (remove.atgs && any(distance < 0)) {
        tmp <- names(seqs)
        substr(seqs, distance[2] + 1, distance[2] + 4) <- "NNNN"
        names(seqs) <- tmp
    }
    if (mask.overlapping.rgns) {
        if (is.null(start.stops)) 
            start.stops <- attr(seqs, "start.stops")
        if (!is.null(start.stops)) {
            overlaps <- apply(start.stops, 1, function(i) subset(start.stops, 
                i[4] == contig & (i[1] >= start & i[1] <= end) | 
                  (i[2] >= start & i[2] <= end)))
            overlaps <- lapply(names(overlaps), function(g) overlaps[[g]][rownames(overlaps[[g]]) != 
                g, ])
            names(overlaps) <- rownames(start.stops)
            is.overlapping <- sapply(overlaps, nrow)
            overlaps <- overlaps[is.overlapping > 0]
            for (i in names(overlaps)) {
                if (nrow(overlaps[[i]]) <= 0) 
                  next
                seq1 <- seqs[i]
                if (start.stops[i, 3] == "R") 
                  seq1 <- rev.comp(seq1)
                ss1 <- sapply(20:nchar(seq1), function(i) substr(seq1, 
                  1, i))
                ss2 <- sapply(1:(nchar(seq1) - 20), function(i) substr(seq1, 
                  i, nchar(seq1)))
                for (j in 1:nrow(overlaps[[i]])) {
                  seq2 <- seqs[rownames(overlaps[[i]])[j]]
                  if (overlaps[[i]][j, 3] == "R") 
                    seq2 <- rev.comp(seq2)
                  g1 <- sapply(sapply(ss1, grep, seq2), length)
                  rgn <- c(1, nchar(seq1))
                  if (all(g1 > 0)) {
                  }
                  else if (any(g1 > 0)) {
                    ind <- which(diff(g1) != 0)
                    rgn <- c(1, ind - 1)
                  }
                  else {
                    g2 <- sapply(sapply(ss2, grep, seq2), length)
                    if (any(g2 > 0)) {
                      ind <- which(diff(g2) != 0)
                      rgn <- c(ind + 1, nchar(seq1))
                    }
                  }
                  if (verbose) 
                    cat(sprintf("Masking region %d-%d of sequence %s (%s)\n", 
                      rgn[1], rgn[2], i, rownames(overlaps[[i]])[j]))
                  substr(seq1, rgn[1], rgn[2]) <- paste(rep("N", 
                    rgn[2] - rgn[1] + 1), collapse = "")
                  seq <- seq1
                  if (start.stops[i, 3] == "R") 
                    seq <- rev.comp(seq)
                  seqs[i] <- seq
                  other.ov <- rownames(overlaps[[i]])[j]
                  overlaps[[other.ov]] <- overlaps[[other.ov]][rownames(overlaps[[other.ov]]) != 
                    i, , drop = F]
                }
            }
        }
    }
    if (!is.null(start.stops)) 
        attr(seqs, "start.stops") <- start.stops[names(seqs), 
            , drop = F]
    seqs
}
fimo.all.motifs <-
function (p.cutoff = 1e-05) 
{
    in.args <- c(mget(names(formals()), env = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    dir.create("filehash/fimo_out")
    seqs.file <- e$my.tempfile("fimo_seqs", tmpdir = "./filehash/fimo_out")
    writeLines(paste(paste(">", names(e$genome.info$genome.seqs), 
        sep = ""), e$genome.info$genome.seqs, sep = "\n"), con = seqs.file)
    inds <- c(seq(1, e$k.clust, by = 100), e$k.clust)
    files <- mclapply(1:(length(inds) - 1), function(i) {
        mots.file <- e$all.motifs.to.mast.file(ks = inds[i]:inds[i + 
            1], seq.type = names(e$mot.weights)[1], e.value.cutoff = Inf, 
            resid.cutoff = Inf)
        out.file <- sprintf("filehash/fimo_out/fimo_out_%d_%d", 
            inds[i], inds[i + 1])
        print(out.file)
        if (file.exists(sprintf("%s.bz2", out.file))) 
            return(out.file)
        cmd <- sprintf("./progs/fimo --max-stored-scores 9999999999 --text --verbosity 2 %s %s |bzip2 -c >%s.bz2", 
            mots.file, seqs.file, out.file)
        print(cmd)
        out <- system(cmd, intern = T)
        out.file
    }, mc.preschedule = F)
    fimo.out <- system(sprintf("bunzip2 -c filehash/fimo_out/fimo_out_*.bz2 | awk '($6<=%s){print}'", 
        as.character(p.cutoff)), intern = T)
    fimo.out <- as.data.frame(do.call(rbind, lapply(fimo.out, 
        function(i) strsplit(i, "\t")[[1]])))
    colnames(fimo.out) <- strsplit(system(sprintf("bunzip2 -c %s.bz2 | head -1", 
        files[1]), intern = T), "\t")[[1]]
    fimo.out$Start <- as.integer(as.character(fimo.out$Start))
    fimo.out$Stop <- as.integer(as.character(fimo.out$Stop))
    fimo.out$`Log-odds` <- as.numeric(as.character(fimo.out$`Log-odds`))
    fimo.out$`p-value` <- as.numeric(as.character(fimo.out$`p-value`))
    tmp <- do.call(rbind, strsplit(as.character(fimo.out$Motif), 
        "_"))
    fimo.out$bic <- as.integer(tmp[, 2])
    fimo.out$mot <- as.integer(tmp[, 3])
    fimo.out$Strand <- substr(tmp[, 1], 1, 1)
    rm(tmp)
    fimo.out$Motif <- NULL
    fimo.out <- as.data.table(fimo.out)
    setkey(fimo.out, bic, mot, Seq, Start)
    attr(fimo.out, "in.args") <- in.args
    save(fimo.out, file = sprintf("filehash/fimo_out_%s.RData", 
        as.character(p.cutoff)))
    fimo.out
}
foreach.register.backend <-
function (par) 
{
    if (!require(foreach)) 
        return(NULL)
    if (par > 1 && require(doMC, quietly = T)) 
        registerDoMC(cores = par)
    else registerDoSEQ()
}
get.all.scores <-
structure(function (ks = 1:k.clust, force.row = F, force.col = F, 
    force.motif = F, force.net = F, quantile.normalize = F) 
{
    mc <- get.parallel(length(ks))
    if (is.null(row.scores)) {
        row.scores <- matrix(0, nrow = attr(ratios, "nrow"), 
            ncol = max(ks))
        rownames(row.scores) <- attr(ratios, "rnames")
        row.scores <- matrix.reference(row.scores)
    }
    if (force.row || (row.scaling[iter] > 0 && !is.na(row.iters[1]) && 
        iter %in% row.iters)) {
        if (is.null(row.scores)) {
            row.scores <- matrix(0, nrow = attr(ratios, "nrow"), 
                ncol = max(ks))
            rownames(row.scores) <- attr(ratios, "rnames")
            row.scores <- matrix.reference(row.scores)
        }
        else row.scores[, ks] <- 0
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        rownames(row.memb) <- attr(ratios, "rnames")
        for (i in names(ratios)) {
            if (row.weights[i] == 0 || is.na(row.weights[i])) 
                next
            tmp.row <- do.call(cbind, mc$apply(ks, get.row.scores, 
                ratios = ratios[[i]]))
            tmp <- is.infinite(tmp.row) | is.na(tmp.row)
            if (any(tmp)) 
                tmp.row[tmp] <- quantile(tmp.row[row.memb[rownames(tmp.row), 
                  ] & !tmp], 0.95)
            tmp <- rownames(row.scores)[rownames(row.scores) %in% 
                rownames(tmp.row)]
            row.scores[tmp, ks] <- row.scores[tmp, ks] + tmp.row[tmp, 
                ] * row.weights[i]
            rm(tmp.row, tmp)
        }
    }
    if (n.clust.per.col < k.clust && (force.col || (row.scaling[iter] > 
        0 && !is.na(col.iters[1]) && iter %in% col.iters))) {
        if (is.null(col.scores)) {
            col.scores <- matrix(0, nrow = attr(ratios, "ncol"), 
                ncol = max(ks))
            rownames(col.scores) <- attr(ratios, "cnames")
            col.scores <- matrix.reference(col.scores)
        }
        else col.scores[, ks] <- 0
        col.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "cnames") %in% get.cols(k))
        rownames(col.memb) <- attr(ratios, "cnames")
        for (i in names(row.weights)) {
            if (row.weights[i] == 0 || is.na(row.weights[i])) 
                next
            tmp.col <- do.call(cbind, mc$apply(ks, get.col.scores, 
                ratios = ratios[[i]]))
            tmp <- is.infinite(tmp.col) | is.na(tmp.col)
            if (any(tmp)) 
                tmp.col[tmp] <- quantile(tmp.col[col.memb[, ][rownames(tmp.col), 
                  ] == 1 & !tmp], 0.95)
            tmp <- rownames(col.scores)[rownames(col.scores) %in% 
                rownames(tmp.col)]
            col.scores[tmp, ks] <- col.scores[tmp, ks] + tmp.col[tmp, 
                ] * row.weights[i]
            rm(tmp.col, tmp)
        }
    }
    for (i in names(mot.weights)) {
        if (force.motif == "run.meme" || (mot.scaling[iter] > 
            0 && !is.na(meme.iters[[i]][1]) && iter %in% meme.iters[[i]] && 
            exists("genome.info") && !no.genome.info)) {
            if (mot.weights[i] == 0 || is.na(mot.weights[i])) 
                next
            tmp <- motif.all.clusters(ks, seq.type = i, verbose = T)
            meme.scores[[i]] <- tmp
        }
    }
    if (force.motif == TRUE || force.motif == "run.meme" || (mot.scaling[iter] > 
        0 && !is.na(mot.iters[1]) && iter %in% mot.iters && exists("genome.info") && 
        !no.genome.info)) {
        if (is.null(mot.scores)) {
            mot.scores <- matrix(0, nrow = attr(ratios, "nrow"), 
                ncol = max(ks))
            rownames(mot.scores) <- attr(ratios, "rnames")
            mot.scores <- matrix.reference(mot.scores)
        }
        else mot.scores[, ks] <- 0
        tmp.mots <- list()
        for (i in names(mot.weights)) {
            if (mot.weights[i] == 0 || is.na(mot.weights[i])) 
                next
            tmp.mot <- do.call(cbind, mc$apply(ks, get.motif.scores, 
                meme.scores = meme.scores, seq.type = i))
            tmp.mot[is.infinite(tmp.mot) | is.na(tmp.mot)] <- 0
            if (quantile.normalize && sum(mot.weights > 0 & !is.na(mot.weights)) > 
                1) 
                tmp.mots[[i]] <- tmp.mot
            else mot.scores[, ks] <- mot.scores[, ks] + tmp.mot[, 
                ] * mot.weights[i]
            rm(tmp.mot)
        }
        if (quantile.normalize && length(tmp.mots) > 1) {
            tmp.mots <- quantile.normalize.scores(tmp.mots, weights = mot.weights[mot.weights != 
                0])
            for (i in names(tmp.mots)) mot.scores[, ks] <- mot.scores[, 
                ks] + tmp.mots[[i]][, ] * mot.weights[i]
            rm(tmp.mots)
        }
    }
    cluster.ns <- NULL
    if (force.net || (net.scaling[iter] > 0 && !is.na(net.iters[1]) && 
        exists("genome.info") && iter %in% net.iters)) {
        if (is.null(net.scores)) {
            net.scores <- matrix(0, nrow = attr(ratios, "nrow"), 
                ncol = max(ks))
            rownames(net.scores) <- attr(ratios, "rnames")
            net.scores <- matrix.reference(net.scores)
        }
        else net.scores[, ks] <- 0
        tmp.nets <- list()
        for (i in names(networks)) {
            if (net.weights[i] == 0 || is.na(net.weights[i])) 
                next
            if (nrow(subset(networks[[i]], protein1 %in% attr(ratios, 
                "rnames") & protein2 %in% attr(ratios, "rnames"))) <= 
                0) 
                next
            tmp.net <- do.call(cbind, mc$apply(ks, get.network.scores, 
                net = networks[[i]]))
            if (all(is.na(tmp.net)) || all(is.character(tmp.net))) 
                next
            tmp.net[is.infinite(tmp.net) | is.na(tmp.net)] <- 0
            if (quantile.normalize && sum(net.weights > 0 & !is.na(net.weights)) > 
                1) 
                tmp.nets[[i]] <- tmp.net
            else net.scores[, ks] <- net.scores[, ks] + tmp.net[, 
                ] * net.weights[i]
            cluster.ns <- cbind(cluster.ns, do.call(c, mc$apply(ks, 
                function(k) mean(tmp.net[get.rows(k), k], na.rm = T, 
                  trim = 0.05))))
            colnames(cluster.ns)[ncol(cluster.ns)] <- i
            rm(tmp.net)
        }
        if (quantile.normalize && length(tmp.nets) > 1) {
            tmp.nets <- quantile.normalize.scores(tmp.nets, weights = net.weights[net.weights != 
                0])
            for (i in names(tmp.nets)) net.scores[, ks] <- net.scores[, 
                ks] + tmp.nets[[i]][, ] * net.weights[i]
            rm(tmp.nets)
        }
        cluster.ns <- cbind(cluster.ns, do.call(c, mc$apply(ks, 
            function(k) mean(net.scores[get.rows(k), k], na.rm = T, 
                trim = 0.05))))
        colnames(cluster.ns)[ncol(cluster.ns)] <- "net.scores"
    }
    list(r = row.scores, m = mot.scores, ms = meme.scores, n = net.scores, 
        c = col.scores, cns = cluster.ns)
}, version = 1)
get.clust <-
function (k, fill = T, fill.motif = T, seq.type = names(mot.weights), 
    varNorm = F, ...) 
{
    gen.clust <- function(rowNames, colNames = NA, fill = F, 
        motif = F, n.motifs = 3, ...) {
        rowNames <- rowNames[rowNames %in% attr(ratios, "rnames")]
        if (!is.null(colNames) && length(colNames) > 1 && !is.na(colNames)) 
            colNames <- colNames[colNames %in% attr(ratios, "cnames")]
        c.tmp <- list(nrows = length(rowNames), ncols = length(colNames), 
            rows = rowNames, cols = colNames, k = 999, p.clust = 1, 
            e.val = rep(999, n.motifs), resid = {
                out = rep(NA, length(row.weights))
                names(out) <- names(row.weights)
                resid = out
            })
        if (fill && c.tmp$nrows > 0 && c.tmp$ncols > 0 && !all(is.na(colNames))) 
            c.tmp$resid <- cluster.resid(k, names(row.weights), 
                varNorm = varNorm, ...)
        return(c.tmp)
    }
    cols <- get.cols(k)
    rows <- get.rows(k)
    if (length(cols) <= 0) 
        cols <- NA
    clust <- gen.clust(rows, cols, fill = fill, motif = F, n.motifs = max(unlist(n.motifs)))
    clust$k <- k
    if (fill.motif) {
        tmp <- cluster.pclust(k, seq.type)
        clust$e.val <- tmp$e.vals
        clust$p.clust <- tmp$p.clusts
    }
    clust
}
get.cluster.matrix <-
function (rows = NULL, cols = NULL, matrices = names(ratios)) 
{
    if (is.null(rows)) 
        rows <- attr(ratios, "rnames")
    if (is.null(cols)) 
        cols <- attr(ratios, "cnames")
    cols.b <- attr(ratios, "cnames")[attr(ratios, "cnames") %in% 
        cols]
    rats <- matrix(NA, nrow = length(rows), ncol = length(cols.b))
    rownames(rats) <- rows
    colnames(rats) <- cols.b
    cnames <- character()
    for (n in matrices) {
        r.tmp <- ratios[[n]][rows[rows %in% rownames(ratios[[n]])], 
            cols.b[cols.b %in% colnames(ratios[[n]])], drop = F]
        if (is.null(r.tmp) || all(is.na(r.tmp))) 
            next
        if (is.vector(r.tmp)) {
            r.tmp <- t(r.tmp)
            rownames(r.tmp) <- rows
        }
        cnames <- c(cnames, colnames(r.tmp))
        rats[rownames(r.tmp), colnames(r.tmp)] <- r.tmp
        rm(r.tmp)
    }
    rats[, colnames(rats) %in% cnames, drop = F]
}
get.clusterStack <-
function (ks = 1:k.clust, force = F, ...) 
{
    if (!force && !is.null(attr(clusterStack, "iter")) && attr(clusterStack, 
        "iter") == iter) 
        return(clusterStack)
    mc <- get.parallel(length(ks))
    clusterStack <- mc$apply(ks, get.clust, ...)
    attr(clusterStack, "iter") <- iter
    clusterStack
}
get.COG.code <-
function (org, rows = attr(ratios, "rnames")) 
{
    up.rows <- toupper(rows)
    out <- rep("-", length(rows))
    names(out) <- up.rows
    fname <- "data/COG_whog.txt"
    err <- dlf(fname, "ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog", 
        "Fetching COG codes from NCBI...")
    lines <- readLines(gzfile(fname))
    closeAllConnections()
    hits <- grep(paste(org, "\\:|COG", sep = ""), lines)
    hpy.hits <- grep(paste(org, "\\:", sep = ""), lines[hits])
    if (length(hpy.hits) <= 0) 
        return(NULL)
    genes <- gsub(paste("\\s+", org, "\\:\\s+", sep = ""), "", 
        lines[hits][hpy.hits], perl = T)
    cogs <- lines[hits][hpy.hits - 1]
    cog.codes <- sapply(strsplit(cogs, "[\\s+\\[\\]]", perl = T), 
        "[", 2)
    cog.codes <- substr(cog.codes, 1, 1)
    genes <- toupper(genes)
    mc <- get.parallel(length(genes))
    tmp <- mc$apply(1:length(genes), function(i) {
        gn <- strsplit(genes[i], " ")[[1]]
        if (length(gn) <= 0) 
            next
        gn <- gn[!is.na(gn)]
        if (!all(gn %in% up.rows)) 
            gn <- toupper(unlist(get.synonyms(gn, ignore = T)))
        if (sum(gn %in% up.rows) <= 0) 
            return(character())
        gn <- gn[gn %in% up.rows]
        out[up.rows %in% gn] <- cog.codes[i]
        out
    })
    for (t in tmp) if (length(t) > 0) 
        out[t != "-"] <- t[t != "-"]
    out[out == "-"] <- NA
    names(out) <- rows
    closeAllConnections()
    out
}
get.cols <-
function (k) 
{
    out <- clusterStack[[k]]$cols
    if (is.null(out) || is.na(out) || length(out) <= 0) 
        out <- attr(ratios, "cnames")
    out
}
get.col.scores <-
function (k, for.cols = "all", ratios = ratios[[1]], norm.method = c("mean", 
    "all.colVars", "none")[1], ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.cols[1] == "all") 
        for.cols <- colnames(ratios)
    rows <- rows[rows %in% rownames(ratios)]
    if (length(rows) <= 1) 
        return(rep(NA, length(for.cols)))
    rats <- ratios[rows, for.cols, drop = F]
    row.weights <- if (exists("get.row.weights")) 
        get.row.weights(rows, cols, ratios)
    else NA
    if (is.na(row.weights[1])) {
        rats.mn <- colMeans(rats, na.rm = T)
    }
    else {
        rats.mn <- apply(rats, 2, weighted.mean, w = row.weights[rows], 
            na.rm = T)
    }
    rats[, ] <- t(t(rats) - rats.mn)^2
    rats <- colMeans(rats, na.rm = T)
    var.norm <- 0.99
    if (norm.method == "all.colVars") {
        all.colVars <- attr(ratios, "all.colVars")
        if (!is.null(all.colVars)) 
            var.norm <- all.colVars[for.cols]
    }
    else if (norm.method == "mean") {
        var.norm <- abs(rats.mn)
    }
    rats <- rats/(var.norm + 0.01)
    rats
}
get.combined.scores <-
function (quantile.normalize = F) 
{
    r.scores <- row.scores[, ]
    r.scores <- matrix.reference(r.scores)
    if (!quantile.normalize) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        rownames(row.memb) <- attr(ratios, "rnames")
        tmp <- r.scores[, ] < -20
        r.scores[, ][tmp] <- min(r.scores[, ][!tmp], na.rm = T)
        rsm <- r.scores[, ][row.memb]
        tmp <- mad(rsm, na.rm = T)
        if (tmp != 0) 
            r.scores[, ] <- (r.scores[, ] - median(rsm, na.rm = T))/tmp
        else {
            tmp <- sd(rsm, na.rm = T)
            if (tmp != 0) 
                r.scores[, ] <- (r.scores[, ] - median(rsm, na.rm = T))/tmp
        }
        rm(tmp, rsm)
    }
    tmp <- r.scores[, ] < -20
    r.scores[, ][tmp] <- min(r.scores[, ][!tmp], na.rm = T)
    rm(tmp)
    r.scores[, ][is.infinite(r.scores[, ])] <- NA
    r.scores[, ][is.na(r.scores[, ])] <- max(r.scores[, ], na.rm = T)
    if (!quantile.normalize && !is.null(mot.scores) || !is.null(net.scores)) 
        rs.quant <- quantile(r.scores[, ], 0.01, na.rm = T)
    if (!is.null(mot.scores)) {
        m.scores <- mot.scores[, ]
    }
    else m.scores <- NULL
    if (!is.null(mot.scores) && !is.null(m.scores)) {
        tmp <- m.scores[, ] < -20
        m.scores[, ][tmp] <- min(m.scores[, ][!tmp], na.rm = T)
        rm(tmp)
        if (!quantile.normalize) {
            m.scores[, ] <- m.scores[, ] - quantile(m.scores[, 
                ], 0.99, na.rm = T)
            m.scores[, ] <- m.scores[, ]/abs(quantile(m.scores[, 
                ], 0.01, na.rm = T)) * abs(rs.quant)
        }
    }
    if (!is.null(net.scores)) {
        n.scores <- net.scores[, ]
        n.scores <- matrix.reference(n.scores)
    }
    else n.scores <- NULL
    if (!is.null(net.scores) && !is.null(n.scores)) {
        n.scores[, ] <- n.scores[, ] - quantile(n.scores[, ], 
            0.99, na.rm = T)
        if (!quantile.normalize) {
            qqq <- abs(quantile(n.scores[, ], 0.01, na.rm = T))
            if (qqq == 0) 
                qqq <- sort(n.scores[, ])[10]
            if (qqq == 0) 
                qqq <- min(n.scores[, ], na.rm = T)
            if (qqq != 0) 
                n.scores[, ] <- n.scores[, ]/qqq * abs(rs.quant)
            rm(qqq)
        }
    }
    if (!is.null(col.scores)) {
        c.scores <- col.scores[, ] * 0
        c.scores <- matrix.reference(c.scores)
        tmp <- c.scores[, ] < -20
        c.scores[, ][tmp] <- min(c.scores[, ][!tmp], na.rm = T)
        rm(tmp)
    }
    else c.scores <- NULL
    new.weights <- c(row = row.scaling[iter], mot = mot.scaling[iter], 
        net = net.scaling[iter])
    if (quantile.normalize) {
        tmp <- list(row = r.scores, mot = m.scores, net = n.scores)
        if (sum(sapply(tmp, function(i) !is.null(i))) > 1) {
            wts <- new.weights[!sapply(tmp, is.null)]
            tmp <- quantile.normalize.scores(tmp, weights = wts)
            if (!is.null(r.scores)) 
                r.scores[, ] <- tmp$row[, ]
            if (!is.null(m.scores)) 
                m.scores[, ] <- tmp$mot[, ]
            if (!is.null(n.scores)) 
                n.scores[, ] <- tmp$net[, ]
        }
        rm(tmp)
    }
    if (new.weights["row"] != 1) 
        r.scores[, ] <- r.scores[, ] * new.weights["row"]
    if (!is.null(m.scores)) {
        tmp <- !is.na(m.scores[, ])
        if (is.null(r.scores)) 
            r.scores <- m.scores[, ] * 0
        r.scores[, ][tmp] <- r.scores[, ][tmp] + m.scores[, ][tmp] * 
            new.weights["mot"]
    }
    if (!is.null(n.scores)) {
        tmp <- !is.na(n.scores[, ])
        if (is.null(r.scores)) 
            r.scores <- n.scores[, ] * 0
        r.scores[, ][tmp] <- r.scores[, ][tmp] + n.scores[, ][tmp] * 
            new.weights["net"]
    }
    r.scores <- matrix.reference(r.scores)
    c.scores <- matrix.reference(c.scores)
    invisible(list(r = r.scores, c = c.scores, scalings = new.weights))
}
get.density.scores <-
function (ks = 1:k.clust, r.scores, c.scores, plot = "none", 
    bw.scale = function(nr) exp(-nr/10) * 10) 
{
    rr <- attr(ratios, "rnames")
    rs <- r.scores
    bw.r <- max(diff(range(rs[, ], na.rm = T))/100, 0.001)
    get.rr.scores <- function(k) {
        rows <- get.rows(k)
        cols <- get.cols(k)
        rsk <- rs[, k, drop = T]
        if (length(rows) > 0 && length(cols) > 0) {
            bw <- bw.r * bw.scale(length(rows))
            d <- density(rsk[rows], na.rm = T, bw = bw, adjust = 2, 
                from = min(rsk, na.rm = T) - 1, to = max(rsk, 
                  na.rm = T) + 1, n = 256)
            p <- approx(d$x, rev(cumsum(rev(d$y))), rsk)$y
            if ("rows" %in% plot) {
                h = hist(rsk, breaks = 50, main = NULL, xlab = "Combined scores")
                tmp.scale <- round(attr(ratios, "nrow")/length(rows)/4)
                hist(rep(rsk[rows], tmp.scale), breaks = h$breaks, 
                  col = "red", border = "red", add = T)
                hist(rsk, breaks = h$breaks, add = T)
                lines(d$x, d$y/max(d$y, na.rm = T) * attr(ratios, 
                  "nrow")/50, col = "blue")
                lines(sort(rsk), p[order(rsk)]/max(p, na.rm = T) * 
                  attr(ratios, "nrow")/50, col = "green")
            }
        }
        else {
            p <- rep(1, attr(ratios, "nrow"))
        }
        return(p/sum(p, na.rm = T))
    }
    rr.scores <- NULL
    mc <- get.parallel(length(ks))
    if (!is.null(row.scores)) {
        rr.scores <- row.scores[, ] * 0
        rr.scores <- matrix.reference(rr.scores)
        rr.scores[, ] <- do.call(cbind, mc$apply(ks, get.rr.scores))
        rr.scores[, ][is.infinite(rr.scores[, ])] <- NA
    }
    cc.scores <- NULL
    if (!is.null(col.scores)) {
        cs <- c.scores
        bw.c <- max(diff(range(cs[, ], na.rm = T))/100, 0.001)
        get.cc.scores <- function(k) {
            cols <- get.cols(k)
            rows <- get.rows(k)
            csk <- cs[, k, drop = T]
            if (length(cols) > 0 && length(rows) > 0 && !all(is.na(csk[cols])) && 
                !all(is.infinite(csk[cols])) & !all(csk[cols][!is.na(csk[cols])] == 
                csk[cols[!is.na(csk[cols])][1]])) {
                d <- density(csk[cols], na.rm = T, from = min(csk, 
                  na.rm = T) - 1, to = max(csk, na.rm = T) + 
                  1, bw = bw.c, adjust = 2, n = 256)
                p <- approx(d$x, rev(cumsum(rev(d$y))), csk)$y
                if ("cols" %in% plot) {
                  h = hist(csk, breaks = 50, main = NULL, xlab = "Combined scores")
                  tmp.scale <- round(attr(ratios, "ncol")/length(cols)/4) + 
                    1
                  hist(rep(csk[cols], tmp.scale), breaks = h$breaks, 
                    col = "red", border = "red", add = T)
                  hist(csk, breaks = h$breaks, add = T)
                  lines(d$x, d$y/max(d$y, na.rm = T) * attr(ratios, 
                    "ncol")/50, col = "blue")
                  lines(sort(csk), p[order(csk)]/max(p, na.rm = T) * 
                    attr(ratios, "ncol")/50, col = "green")
                }
            }
            else {
                p <- rep(1, attr(ratios, "ncol"))
            }
            return(p/sum(p, na.rm = T))
        }
        cc.scores <- col.scores[, ] * 0
        cc.scores <- matrix.reference(cc.scores)
        if (!is.null(c.scores)) {
            cc.scores[, ] <- do.call(cbind, mc$apply(ks, get.cc.scores))
            cc.scores[, ][is.infinite(cc.scores)] <- NA
        }
    }
    invisible(list(r = rr.scores, c = cc.scores))
}
get.dup.seqs <-
function (seqs) 
{
    out <- duplicated(seqs)
    names(out) <- names(seqs)
    out
}
get.gene.coords <-
function (rows, op.shift = T, op.table = genome.info$operons, 
    ...) 
{
    if (is.null(rows)) {
        if (organism == "hal") 
            rows <- grep("^VNG", as.character(genome.info$feature.tab$canonical_Name), 
                perl = T, val = T)
        else if (organism == "sce") 
            rows <- grep(e$genome.info$gene.regex, as.character(genome.info$feature.tab$id), 
                perl = T, val = T)
        else rows <- grep("^NP_", as.character(genome.info$feature.tab$id), 
            perl = T, val = T)
    }
    rows <- unique(rows)
    syns <- get.synonyms(rows, ...)
    tab <- genome.info$feature.tab
    ids <- lapply(syns, function(s) s[s %in% tab$id])
    if (all(sapply(ids, length) < 1)) {
        warning("Could not find gene start/stop for any input genes", 
            call. = F)
        return(NULL)
    }
    if (any(sapply(ids, length) < 1)) 
        warning("Could not find gene start/stop for all input genes", 
            call. = F)
    ids <- ids[sapply(ids, length) >= 1]
    ids <- sapply(ids, "[", 1)
    ids <- data.frame(id = ids, names = names(ids))
    coos <- NULL
    if (op.shift && exists("op.table")) {
        if (attr(op.table, "source") == "rsat") {
            ops <- merge(ids, op.table, by.x = "id", by.y = "query", 
                all = F, sort = F)
            ops2 <- ops[order(ops$lead), ]
            coos <- merge(ops, tab, by.x = "lead", by.y = "name", 
                all = F, sort = F)[, c("id.x", "names", "contig", 
                "strand", "start_pos", "end_pos")]
        }
        else if (attr(op.table, "source") == "microbes.online") {
            ops <- NULL
            if (!any(ids$names %in% op.table$gene)) {
                ids2 <- lapply(syns, function(s) s[s %in% op.table$gene])
                if (all(sapply(ids2, length) < 1)) {
                  warning("Could not find operon info for any input genes", 
                    call. = F)
                }
                else {
                  if (any(sapply(ids2, length) < 1)) 
                    warning("Could not find operon info for all input genes", 
                      call. = F)
                  ids2 <- ids2[sapply(ids2, length) >= 1]
                  ids2 <- sapply(ids2, "[", 1)
                  ids2 <- data.frame(id = ids2, names = names(ids2))
                  if (any(!rows %in% ids2$names)) {
                    capture.output(prefix <- get.gene.regex(as.character(ids2$id))[1])
                    missing <- rows[!rows %in% ids2$names]
                    syns2 <- sapply(syns[missing], function(i) grep(sprintf("^%s", 
                      prefix), i, val = T)[1])
                    syns2[is.na(syns2)] <- names(syns2)[is.na(syns2)]
                    names(syns2) <- NULL
                    ids2 <- rbind(ids2, data.frame(id = syns2, 
                      names = missing))
                  }
                  ops <- merge(ids2, op.table, by.x = "id", by.y = "gene", 
                    all.x = T, incomparables = NA, sort = F)
                }
            }
            if (is.null(ops)) 
                ops <- merge(ids, op.table, by.x = "names", by.y = "gene", 
                  all.x = T, incomparables = NA, sort = F)
            if (any(is.na(ops$head))) {
                head <- as.character(ops$head)
                head[is.na(head)] <- as.character(ops$id[is.na(head)])
                ops$head <- as.factor(head)
            }
            head.genes <- unique(as.character(ops$head))
            head.genes <- head.genes[!is.na(head.genes)]
            head.syns <- get.synonyms(head.genes)
            head.ids <- lapply(head.syns, function(s) s[s %in% 
                tab$id])
            head.ids <- head.ids[sapply(head.ids, length) >= 
                1]
            head.ids <- data.frame(id = sapply(head.ids, "[", 
                1), names = names(head.ids))
            ops2 <- merge(ops, head.ids, by.x = "head", by.y = "names", 
                all.x = T, sort = F)
            coos <- merge(ops2, tab, by.x = "id.y", by.y = "id", 
                all.x = T, sort = F)[, c("id.x", "names", "contig", 
                "strand", "start_pos", "end_pos")]
        }
    }
    else {
        coos <- merge(ids, tab, by = "id", sort = F)[, c("id", 
            "names", "contig", "strand", "start_pos", "end_pos")]
    }
    colnames(coos)[1] <- "id"
    if (is.factor(coos$start_pos)) 
        coos$start_pos <- as.numeric(levels(coos$start_pos))[coos$start_pos]
    if (is.factor(coos$end_pos)) 
        coos$end_pos <- as.numeric(levels(coos$end_pos))[coos$end_pos]
    coos[!duplicated(coos[, 1:4]), ]
}
get.gene.regex <-
function (names = NULL, verbose = F) 
{
    if (!is.null(names)) 
        tmp <- names
    else {
        if (exists("ratios") && !is.null(ratios)) 
            tmp <- toupper(attr(ratios, "rnames"))
        else if (exists("genome.info") && !is.null(genome.info$feature.names)) {
            tmp <- toupper(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T))
            if (exists("ratios") && !is.null(ratios)) 
                tmp <- tmp[toupper(tmp) %in% toupper(attr(ratios, 
                  "rnames"))]
        }
    }
    qqq <- sapply(1:4, function(nch) max(table(substr(tmp, 1, 
        nch)))/length(tmp))
    nch <- 0
    if (any(qqq > 0.9)) {
        nch <- which(qqq > 0.9)
        nch <- nch[length(nch)]
    }
    else if (any(qqq > 0.6)) {
        nch <- which(qqq > 0.6)
        nch <- nch[length(nch)]
    }
    else if (any(qqq > 0.4)) {
        nch <- which(qqq > 0.4)
        nch <- nch[length(nch)]
    }
    prefix <- NA
    if (nch > 0) {
        prefix <- names(which.max(table(substr(tmp, 1, nch))))
        if (verbose) 
            message("Assuming gene/probe names have common prefix '", 
                prefix, "'.")
        genome.info$gene.prefix <- prefix
    }
    else {
        if (verbose) 
            message("Could not find a common gene/probe identifier prefix. This only matters if there's no expression matrix.")
        prefix <- genome.info$gene.prefix <- NA
    }
    tmp2 <- tmp
    if (length(unique(nchar(tmp2))) > 1) {
        nc <- max(nchar(tmp2))
        for (i in 1:length(tmp2)) tmp2[i] <- paste(tmp2[i], rep(" ", 
            nc - nchar(tmp2[i])), sep = "", collapse = "")
    }
    tmp2 <- do.call(rbind, strsplit(tmp2, ""))
    regex <- apply(tmp2, 2, function(i) sort(unique(i)))
    for (i in 1:length(regex)) {
        ii <- as.integer(regex[[i]])
        if (!any(is.na(ii))) {
            if (length(ii) == length(ii[1]:ii[length(ii)]) && 
                all(ii) == ii[1]:ii[length(i)]) 
                regex[[i]] <- paste("[", paste(ii[1], ii[length(ii)], 
                  sep = "-"), "]", sep = "")
        }
        if (length(regex[[i]][regex[[i]] != " "]) > 1) 
            regex[[i]] <- c("[", regex[[i]], "]")
        if (any(regex[[i]] == "" | regex[[i]] == " " | is.na(regex[[i]]))) 
            regex[[i]] <- c(regex[[i]][regex[[i]] != " "], "?")
    }
    regex <- paste(unlist(lapply(regex, paste, sep = "", collapse = "")), 
        sep = "", collapse = "")
    if (verbose) 
        message("Assuming gene/probe names have regex '", regex, 
            "'.")
    c(prefix, regex)
}
get.genome.info <-
function (fetch.upstream = F) 
{
    rsat.url <- rsat.urls[1]
    feature.tab <- feature.names <- genome.seqs <- operons <- org.id <- synonyms <- NULL
    genome.loc <- paste(rsat.url, "/data/genomes/", rsat.species, 
        "/genome/", sep = "")
    fname <- paste("data/", rsat.species, "/organism_names.tab", 
        sep = "")
    err <- dlf(fname, paste(genome.loc, "/organism_names.tab", 
        sep = ""))
    if (class(err) == "try-error") {
        tmp.url <- paste(rsat.url, "/data/genomes/", rsat.species, 
            "_EnsEMBL/genome/organism_names.tab", sep = "")
        err <- dlf(fname, tmp.url)
        if (class(err) != "try-error") 
            genome.loc <- paste(rsat.url, "/data/genomes/", rsat.species, 
                "_EnsEMBL/genome/", sep = "")
    }
    if (!file.exists(fname) || file.info(fname)$size <= 0) 
        stop(paste("Genome info for", rsat.species, "does not exist. Please check", 
            genome.loc, "and let me know if I am wrong"))
    nskip <- sum(substr(readLines(gzfile(fname), n = 20), 1, 
        2) == "--" | readLines(gzfile(fname), n = 20) == "")
    org.id <- read.delim(gzfile(fname), head = F, as.is = T, 
        skip = nskip)
    if (!exists("taxon.id") || is.na(taxon.id) || is.null(taxon.id)) 
        taxon.id <- org.id$V1[1]
    cat("Organism taxon id:", taxon.id, "\n")
    closeAllConnections()
    if (!no.genome.info && !all(grepl("file=", names(mot.weights)))) {
        fname <- paste("data/", rsat.species, "/feature.tab", 
            sep = "")
        use.cds <- FALSE
        err <- dlf(fname, paste(genome.loc, "feature.tab", sep = ""), 
            paste("Fetching genome annotation data from RSAT", 
                rsat.url, "..."))
        if (class(err) == "try-error") {
            err <- dlf(fname, paste(genome.loc, "cds.tab", sep = ""))
            use.cds <- TRUE
        }
        cat("Loading genome annotation data...\n")
        head <- readLines(gzfile(fname), n = 30)
        nskip <- length(grep("^--", head))
        feature.tab <- read.delim(gzfile(fname), skip = nskip, 
            head = F, comment = "", as.is = F)
        closeAllConnections()
        head <- strsplit(gsub("^-- ", "", head[grep("^-- id", 
            head, perl = T)], perl = T), "\t")[[1]]
        colnames(feature.tab) <- head[1:ncol(feature.tab)]
        fname <- paste("data/", rsat.species, "/feature_names.tab", 
            sep = "")
        err <- dlf(fname, paste(genome.loc, if (!use.cds) 
            "feature_names.tab"
        else "cds_names.tab", sep = ""))
        nskip <- sum(substr(readLines(gzfile(fname), n = 20), 
            1, 2) == "--")
        closeAllConnections()
        feature.names <- read.delim(gzfile(fname), head = F, 
            as.is = T, skip = nskip, row.names = NULL, comment = "")
        closeAllConnections()
        colnames(feature.names) <- c("id", "names", "type")
        feature.names <- unique(feature.names)
        chroms <- unique(as.character(feature.tab$contig))
        chroms <- chroms[!is.na(chroms) & chroms != ""]
        if (!is.na(mot.iters[1])) {
            genome.seqs <- list()
            for (i in chroms) {
                cat("Loading genome sequence, chromosome", i, 
                  "\n")
                fname <- paste("data/", rsat.species, "/", i, 
                  ".raw", sep = "")
                err <- dlf(fname, paste(genome.loc, i, ".raw", 
                  sep = ""))
                if (class(err) == "try-error") {
                  ii <- gsub(":", "_", i, fixed = T)
                  err <- dlf(fname, paste(genome.loc, ii, ".raw", 
                    sep = ""))
                  if (class(err) == "try-error") {
                    err <- dlf(fname, paste(genome.loc, gsub(".[0-9]$", 
                      "", i), ".raw", sep = ""))
                    if (class(err) == "try-error") 
                      cat("ERROR reading genome sequence", i, 
                        "\n")
                    else fname <- paste("data/", rsat.species, 
                      "/", gsub(".[0-9]$", "", i), ".raw", sep = "")
                  }
                  else fname <- paste("data/", rsat.species, 
                    "/", ii, ".raw", sep = "")
                }
                out <- try(readLines(gzfile(fname)), silent = T)
                if (class(out) == "try-error" || length(out) == 
                  0 || is.na(out) || out == "" || out == "CHARACTER(0)") 
                  out <- try(readLines(fname), silent = T)
                if (class(out) == "try-error" || length(out) == 
                  0 || is.na(out) || out == "" || out == "CHARACTER(0)") {
                  cat("ERROR reading genome sequence", i, "\n")
                  next
                }
                if (length(out) > 1) 
                  out <- paste(out, collapse = "", sep = "")
                out <- toupper(out)
                print(nchar(out))
                genome.seqs[[i]] <- out
            }
            if (length(genome.seqs) != length(chroms)) {
                cat("WARNING: Could not read sequence for chromosomes", 
                  chroms[!chroms %in% names(genome.seqs)], "\n")
                feature.tab <- subset(feature.tab, contig %in% 
                  names(genome.seqs))
            }
            if (length(genome.seqs) <= 0) 
                genome.seqs <- NULL
        }
    }
    synonyms <- NULL
    if (exists("synonym.thesaurus")) {
        cat("Loading synonym information from custom thesaurus...\n")
        tmp <- read.csv(gzfile(synonym.thesaurus), header = F)
        synonyms <- strsplit(toupper(as.character(tmp[, 2])), 
            ";")
        names(synonyms) <- toupper(as.character(tmp[, 1]))
        rm(tmp)
    }
    else if (exists("ratios") && !is.null(ratios)) {
        cat("Gathering all \"standard\" orf names and other synonyms for all probe names...\n")
        synonyms <- get.synonyms(attr(ratios, "rnames"), feature.names, 
            verbose = T)
    }
    if (!is.null(synonyms)) {
        cat("Mean number of synonyms per probe:", mean(sapply(synonyms, 
            length), na.rm = T), "\n")
        is.bad <- sapply(names(synonyms), function(i) length(synonyms[[i]]) == 
            0 || substr(synonyms[[i]][1], 1, 5) == "Error")
        if (sum(is.bad) > 0) {
            cat("These", sum(is.bad), "probe names have no matching ORF annotation:\n")
            print(names(which(is.bad)))
        }
        rm(is.bad)
    }
    closeAllConnections()
    invisible(list(species = rsat.species, genome.seqs = genome.seqs, 
        feature.tab = feature.tab, feature.names = feature.names, 
        org.id = org.id, taxon.id = taxon.id, synonyms = synonyms))
}
get.long.names <-
function (k, shorter = F) 
{
    if (is.numeric(k[1])) {
        rows <- get.rows(k)
    }
    else {
        rows <- k
    }
    if (is.null(genome.info$feature.tab)) {
        out <- rep("", length(rows))
        names(out) <- rows
        return(rows)
    }
    ids <- get.synonyms(rows)
    mc <- list(apply = lapply)
    if (!shorter) 
        desc <- mc$apply(ids, function(i) subset(genome.info$feature.tab, 
            id %in% i, select = c("id", "description")))
    else {
        desc <- mc$apply(ids, function(i) subset(genome.info$feature.tab, 
            id %in% i, select = c("id", "name", "description")))
        for (i in 1:length(desc)) if (length(desc[[i]]$name) > 
            0 && desc[[i]]$name %in% rows) {
            if (grepl("(", desc[[i]]$description, fixed = T)) 
                desc[[i]]$name <- strsplit(as.character(desc[[i]]$description), 
                  "[()]", perl = T)[[1]][2]
        }
    }
    out <- sapply(desc, function(i) as.character(i[1, 2]))
    out <- out[rows]
    names(out) <- rows
    out[is.na(out) | out == names(out)] <- ""
    out
}
get.mast.pvals <-
function (mast.output, in.genes = NULL) 
{
    space.pad <- function(lines, length) {
        nc <- nchar(lines)
        nc[nc >= length] <- 0
        spaces <- sapply(1:length(lines), function(i) paste(rep(" ", 
            length - nc[i]), sep = "", collapse = ""))
        paste(lines, spaces)
    }
    out <- list()
    start <- grep("SECTION III: ANNOTATED SEQUENCES", mast.output)
    if (length(start) == 0 || is.na(start)) 
        return(out)
    end <- grep("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", 
        mast.output[(start + 3):length(mast.output)]) + start + 
        3
    line.starts <- grep("LENGTH = ", mast.output[(start + 2):(start + 
        1 + end)]) + start + 1
    if (is.null(line.starts) || length(line.starts) == 0) 
        return(out)
    for (i in 1:length(line.starts)) {
        l <- line.starts[i]
        gene <- mast.output[l - 2]
        if (is.null(gene) || is.na(gene) || (!is.null(in.genes) && 
            !(gene %in% in.genes) && !(toupper(gene) %in% toupper(in.genes)))) 
            next
        l.next <- line.starts[i + 1] - 2
        if (i >= length(line.starts)) 
            l.next <- end
        if (l.next - l <= 5) 
            next
        submast <- mast.output[l:(l.next - 1)]
        l.start <- which(submast == "")[1] + 1
        if (submast[l.start] == "") 
            l.start <- l.start + 1
        q <- list()
        for (i in 1:6) q[[i]] <- space.pad(submast[seq((l.start + 
            i - 1), length(submast), by = 6)], 80)
        seq.starts <- as.integer(sapply(strsplit(q[[5]], " "), 
            "[", 1))
        char.skip <- which(strsplit(q[[5]][1], "")[[1]] %in% 
            c("G", "A", "T", "C", "N", "X"))[1]
        mots <- unlist(strsplit(gsub("[\\[\\]\\<\\>]", "", paste(substr(q[[1]], 
            char.skip, 80), collapse = ""), perl = T), "\\s+", 
            perl = T))
        mots <- as.integer(mots[!is.na(as.integer(mots))])
        mots <- mots[!is.na(mots)]
        p.vals <- strsplit(paste(substr(q[[2]], char.skip, 80), 
            collapse = ""), "\\s+")[[1]]
        p.vals <- as.numeric(p.vals[!is.na(as.numeric(p.vals))])
        posns <- integer()
        for (i in 1:length(q[[1]])) {
            posns <- c(posns, which(strsplit(substr(q[[1]][i], 
                char.skip, 80), "")[[1]] %in% c("[", "<")) + 
                seq.starts[i])
        }
        out[[gene]] <- list(pvals = p.vals, mots = mots, posns = posns)
    }
    return(out)
}
getMastPValuesAndEValues <-
function (mastOutput, get.p.values = NULL) 
{
    lines <- grep("COMBINED P-VALUE", mastOutput)
    if (length(lines) > 0) {
        splitted <- strsplit(mastOutput[lines], "[\\t\\s]+", 
            perl = T)
        out <- t(sapply(1:length(lines), function(i) {
            gene <- mastOutput[lines[i] - 2]
            splt <- splitted[[i]]
            p.val <- splt[8]
            e.val <- splt[11]
            c(gene = gene, p.value = p.val, e.value = e.val)
        }))
        out <- data.frame(gene = out[, "gene"], p.value = as.numeric(out[, 
            "p.value"]), e.value = as.numeric(out[, "e.value"]))
    }
    out2 <- data.frame()
    if (!is.null(get.p.values) && !is.na(get.p.values)) {
        tmp <- get.mast.pvals(mastOutput, in.genes = get.p.values)
        for (g in names(tmp)) {
            pv <- as.numeric(tmp[[g]]$pvals)
            pos <- as.integer(tmp[[g]]$posns)
            mots <- as.integer(tmp[[g]]$mots)
            if (!all(c(length(pv), length(pos)) == length(mots))) 
                pv <- c(pv, rep(pv[1], length(pos) - length(pv)))
            out2 <- rbind(out2, data.frame(gene = g, pvals = pv, 
                posns = pos, mots = mots))
        }
    }
    return(list(out, out2))
}
getMemeMotifInfo <-
function (memeOutput) 
{
    out <- list()
    lines <- grep("^MOTIF\\s+\\d", memeOutput, perl = T)
    if (length(lines) <= 0) 
        lines <- grep("^MOTIF\\s+", memeOutput, perl = T)
    if (length(lines) > 0) {
        pssms <- getMemeMotifPssm(memeOutput, n.motif = length(lines))
        splitted <- strsplit(memeOutput[lines], "[\\t\\s]+", 
            perl = T)
        for (i in 1:length(lines)) {
            splt <- splitted[[i]]
            motif <- as.integer(splt[2])
            width <- as.integer(splt[5])
            sites <- as.integer(splt[8])
            llr <- as.integer(splt[11])
            e.value <- as.numeric(sub("\\+", "", splt[14]))
            pssm <- pssms[[motif]]$pssm
            l2 <- grep(paste("Motif", motif, "sites sorted by position p-value"), 
                memeOutput) + 4
            l3 <- grep("--------------------------------------------------------------------------------", 
                memeOutput[(l2 + 1):length(memeOutput)])[1] + 
                l2 - 1
            posns <- do.call(rbind, strsplit(memeOutput[l2:l3], 
                "[\\t\\s]+", perl = T))[, c(1:4, 6)]
            colnames(posns) <- c("gene", "strand", "start", "p.value", 
                "site")
            posns <- data.frame(gene = posns[, "gene"], strand = posns[, 
                "strand"], start = as.integer(posns[, "start"]), 
                p.value = as.numeric(posns[, "p.value"]), site = posns[, 
                  "site"])
            out[[motif]] <- list(width = width, sites = sites, 
                llr = llr, e.value = e.value, pssm = pssm, posns = posns)
        }
    }
    out
}
getMemeMotifPssm <-
function (memeOut, n.motif = 1) 
{
    pssms <- list()
    for (i in 1:n.motif) {
        m.line1 <- grep(sprintf("Motif %d position-specific probability matrix", 
            i), memeOut)
        if (length(m.line1) > 0) {
            m.desc <- strsplit(memeOut[m.line1 + 2], " ")[[1]]
            winLen <- as.numeric(m.desc[6])
            e.val <- as.numeric(m.desc[10])
            pssm <- do.call(rbind, strsplit(memeOut[m.line1 + 
                2 + 1:winLen], "\\s+", perl = T))[, 2:5]
            pssm <- matrix(as.numeric(pssm), nrow = winLen, ncol = 4, 
                byrow = F)
            pssms[[i]] <- list(pssm = pssm, e.val = e.val)
        }
        else {
            pssms[[i]] <- list(pssm = NULL, e.val = 99999)
        }
    }
    return(pssms)
}
get.motif.scores <-
function (k, meme.scores, seq.type = "upstream meme", for.rows = "all") 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.rows[1] == "all") 
        for.rows <- attr(ratios, "rnames")
    if (length(rows) <= 1 || is.null(meme.scores[[seq.type]]$all.pv)) 
        return(rep(NA, length(for.rows)))
    m.scores <- log(meme.scores[[seq.type]]$all.pv[, k])
    m.scores <- m.scores[for.rows]
    return(m.scores)
}
get.network.scores <-
function (k, net = networks$string, for.rows = "all", p1.col = "protein1", 
    p2.col = "protein2", score.col = "combined_score", combine.func = sum) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.rows[1] == "all") 
        for.rows <- attr(ratios, "rnames")
    if (length(rows) < 1) 
        return(rep(NA, length(for.rows)))
    cons <- net[as.character(net[[p1.col]]) %in% rows, c(p2.col, 
        score.col), drop = F]
    if (is.null(cons) || nrow(cons) <= 0) 
        return(rep(NA, length(for.rows)))
    cons <- cons[as.character(cons[[p2.col]]) %in% for.rows, 
        , drop = F]
    if (is.null(cons) || nrow(cons) <= 0) 
        return(rep(NA, length(for.rows)))
    tmp <- tapply(as.numeric(cons[[score.col]]), as.character(cons[[p2.col]]), 
        combine.func, na.rm = T)/length(rows)
    scores <- rep(NA, length(for.rows))
    names(scores) <- for.rows
    scores[names(tmp)] <- tmp
    return(-log(scores + 1))
}
get.operon.predictions <-
function (fetch.predicted.operons = "microbes.online", org.id = genome.info$org.id$V1[1]) 
{
    operons <- NULL
    if (fetch.predicted.operons == "rsat") {
        rsat.url <- rsat.urls[1]
        cat("Using operon predictions from RSAT...\n")
        fname <- paste("data/", rsat.species, "/rsat_operon_predictions.html", 
            sep = "")
        err <- dlf(fname, paste(rsat.url, "/infer-operons.cgi?organism=", 
            rsat.species, "&genes=all&return_leader=on&return_operon=on&return_query=on&", 
            "output=display&dist_thr=55", sep = ""))
        operons <- readLines(gzfile(fname))
        start <- which(operons == "<INPUT type=\"hidden\" NAME=\"gene_selection\" VALUE=\"#lead\toperon\tquery") + 
            1
        end <- which(operons == "<INPUT type=\"hidden\" NAME=\"feattype\" VALUE=\"\">") - 
            2
        operons <- do.call(rbind, strsplit(operons[start:end], 
            "\t+", perl = T))
        colnames(operons) <- c("lead", "operon", "query")
        operons <- as.data.frame(operons)
    }
    else if (fetch.predicted.operons == "microbes.online") {
        cat("Using operon predictions from MicrobesOnline...\n")
        fname <- paste("data/", rsat.species, "/microbesonline_operons_gnc", 
            org.id, ".named", sep = "")
        err <- dlf(fname, paste("http://www.microbesonline.org/operons/gnc", 
            org.id, ".named", sep = ""))
        if (org.id != taxon.id && (!file.exists(fname) || file.info(fname)$size == 
            0)) {
            fname <- paste("data/", rsat.species, "/microbesonline_operons_gnc", 
                taxon.id, ".named", sep = "")
            err <- dlf(fname, paste("http://www.microbesonline.org/operons/gnc", 
                taxon.id, ".named", sep = ""))
        }
        if (file.exists(fname)) 
            cat("Succesfully fetched operon predictions. Parsing...\n")
        ops <- read.delim(gzfile(fname))
        ops2 <- ops
        on <- as.character(ops2$SysName1)[-1]
        bOp <- ops2$bOp[-1]
        on[which(on == "")] <- on[which(on == "") - 1]
        ops2$SysName1 <- c(as.character(ops2$SysName1)[1], on)
        on <- as.character(ops2$SysName2)[-1]
        on[which(on == "")] <- on[which(on == "") - 1]
        ops2$SysName2 <- c(as.character(ops2$SysName2)[1], on)
        ops2 <- subset(ops2, bOp == "TRUE")
        sn1 <- as.character(ops2$SysName1)
        sn1[sn1 == "" | is.na(sn1)] <- as.character(ops2$Name1)[sn1 == 
            "" | is.na(sn1)]
        sn2 <- as.character(ops2$SysName2)
        sn2[sn2 == "" | is.na(sn2)] <- as.character(ops2$Name2)[sn2 == 
            "" | is.na(sn2)]
        operons <- list(0)
        for (i in 1:length(sn1)) {
            sn1i <- sn1[i]
            if (sn1i == "") 
                next
            found <- which(sapply(operons, function(j) sn1i %in% 
                j))
            if (length(found) > 0) 
                operons[[found[1]]] <- c(operons[[found[1]]], 
                  sn2[i])
            else operons[[length(operons) + 1]] <- c(sn1i, sn2[i])
        }
        operons <- operons[-1]
        operons <- lapply(operons, unique)
        operons <- lapply(operons, function(i) i[i != ""])
        operons <- operons[sapply(operons, length) > 1]
        gns <- unlist(operons)
        search.names <- c(gns, as.character(genome.info$feature.names$id))
        if (exists("ratios")) 
            search.names <- c(attr(ratios, "rnames"), search.names)
        search.names <- unique(search.names)
        mc <- get.parallel(length(operons))
        all.ids <- unique(genome.info$feature.names$id)
        nms <- mc$apply(1:length(operons), function(i) {
            s <- get.synonyms(operons[[i]])
            s <- lapply(s, function(i) i[i %in% search.names])
            ids <- unlist(lapply(s, function(i) i[i %in% all.ids][1]))
            if (length(ids) <= 0) {
                warning(paste("No genome annotation for any genes in operon #", 
                  i, " -- don't know what to do!", call. = F))
                return("")
            }
            ids[is.na(ids)] <- names(ids)[is.na(ids)]
            vngs <- unlist(lapply(s, function(i) {
                out <- i[!i %in% all.ids]
                if (length(out) <= 0) 
                  out <- i[i %in% search.names]
                if (length(out) <= 0) 
                  out <- i[genome.info$feature.names$id %in% 
                    i & genome.info$feature.names$id == "primary"]
                if (length(out) <= 0) 
                  out <- i
                if (length(out) > 1 && any(out %in% attr(ratios, 
                  "rnames"))) 
                  out <- out[out %in% attr(ratios, "rnames")]
                out
            }))
            coos <- get.gene.coords(ids, op.shift = F)
            vngs <- vngs[ids %in% coos$names]
            if (is.null(coos) || nrow(coos) <= 0) {
                warning(paste("No genome annotation for any genes in operon #", 
                  i, " -- don't know what to do!", call. = F))
                return("")
            }
            if (mean(as.character(coos$strand) == "D") > 0.6) 
                head <- vngs[which.min(coos$start_pos)]
            else if (mean(as.character(coos$strand) == "R") > 
                0.6) 
                head <- vngs[which.max(coos$end_pos)]
            else {
                head <- ""
                warning(paste("About 50% of operon #", i, "are on opposite strands -- don't know what to do!", 
                  call. = F))
            }
            head
        })
        names(operons) <- unlist(nms)
        operons <- operons[names(operons) != ""]
        operons <- do.call(rbind, lapply(names(operons), function(h) data.frame(head = h, 
            gene = operons[[h]])))
        operons <- subset(operons, head != "")
    }
    if (!is.null(operons)) 
        attr(operons, "source") <- fetch.predicted.operons
    closeAllConnections()
    operons
}
get.parallel <-
function (X = k.clust, verbose = F, para.cores = get("parallel.cores")) 
{
    if (is.na(para.cores) || (is.logical(para.cores) && para.cores == 
        FALSE) || (is.numeric(para.cores) && para.cores <= 1)) {
        out <- list(mc = FALSE, par = para.cores, apply = lapply)
        if (verbose) 
            cat("NOT PARALLELIZING\n")
    }
    else {
        try(has.multi <- require(multicore, quietly = T), silent = T)
        if (!has.multi || (has.multi && multicore:::isChild())) {
            out <- list(mc = FALSE, par = para.cores, apply = lapply)
            if (verbose) 
                cat("NOT PARALLELIZING\n")
        }
        else {
            mc <- has.multi && !multicore:::isChild() && X > 
                1 && !is.na(para.cores) && (is.numeric(para.cores) && 
                para.cores > 1) || (is.logical(para.cores) && 
                para.cores == TRUE)
            par <- para.cores
            out.apply <- lapply
            if (mc) {
                if (is.logical(par) && par == TRUE) 
                  par <- multicore:::detectCores()
                par <- min(c(X, par, multicore:::detectCores()))
                if (verbose) 
                  cat("PARALLELIZING:", par, ": ")
                foreach.register.backend(par)
                if (verbose) 
                  cat(getDoParName(), getDoParWorkers(), "\n")
                out.apply <- function(list, FUN, ...) foreach(l = list) %dopar% 
                  {
                    FUN(l, ...)
                  }
            }
            else {
                par <- 1
                if (verbose) 
                  cat("NOT PARALLELIZING:", par, "\n")
            }
            out <- list(mc = mc, par = par, apply = out.apply)
        }
    }
    if (is.numeric(out$par) && !is.na(out$par)) 
        options(cores = out$par)
    else if (is.na(out$par) || (is.logical(out$par) && out$par == 
        TRUE)) 
        options(cores = NULL)
    else options(cores = 1)
    out
}
get.pv.ev.single <-
function (mast.out, rows) 
{
    pv.ev <- NULL
    if (length(grep("Error reading log-odds matrix file", mast.out)) <= 
        0 && class(mast.out) != "try-error" && length(mast.out) > 
        0) {
        pv.ev <- getMastPValuesAndEValues(mast.out, get.p.values = rows)
        attr(pv.ev, "mast.command.line") <- attr(mast.out, "mast.command.line")
        if (length(pv.ev) > 0 && nrow(pv.ev[[1]]) == 0 && nrow(pv.ev[[2]]) == 
            0) {
            pv.ev <- NULL
        }
        else {
            for (i in 1) {
                tmp <- as.matrix(pv.ev[[i]][, 2:ncol(pv.ev[[i]])])
                rownames(tmp) <- pv.ev[[i]][, 1]
                pv.ev[[i]] <- tmp
            }
        }
    }
    pv.ev
}
get.rows <-
function (k) 
clusterStack[[k]]$rows
get.row.scores <-
function (k, cols = get.cols(k), for.rows = "all", ratios = ratios[[1]], 
    ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (is.null(for.rows) || for.rows[1] == "all") 
        for.rows <- rownames(ratios)
    rows <- rows[rows %in% rownames(ratios)]
    cols <- cols[cols %in% colnames(ratios)]
    if (length(rows) < 1 || length(cols) < 1) 
        return(rep(NA, length(for.rows)))
    rats <- ratios[for.rows, cols, drop = F]
    rats.mn <- colMeans(rats[rows, , drop = F], na.rm = T)
    rats[, ] <- t(t(rats) - rats.mn)^2
    col.weights <- if (exists("get.col.weights")) 
        get.col.weights(rows, cols, ratios)
    else NA
    if (is.na(col.weights[1])) 
        rats <- rowMeans(rats, na.rm = T)
    else rats <- apply(rats, 1, weighted.mean, w = col.weights[cols], 
        na.rm = T)
    rats <- log(rats + 1e-99)
    return(rats)
}
get.sequences <-
function (k, seq.type = paste(c("upstream", "upstream.noncod", 
    "upstream.noncod.same.strand", "downstream", "gene")[1], 
    "meme"), verbose = F, filter = T, distance = motif.upstream.search[[seq.type]], 
    op.shift, ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else if (!is.null(genome.info$genome.seqs) && k %in% names(genome.info$genome.seqs)) 
        return(genome.info$genome.seqs[k])
    else rows <- k
    if (is.null(rows)) 
        return(NULL)
    start.stops <- NULL
    n.seq.type <- strsplit(seq.type, " ")[[1]][1]
    if (substr(n.seq.type, 1, 8) == "fstfile=") {
        if (!is.null(genome.info$all.upstream.seqs[[seq.type]]) && 
            length(genome.info$all.upstream.seqs[[seq.type]]) > 
                0) {
            seqs <- genome.info$all.upstream.seqs[[seq.type]]
        }
        else {
            seqs <- read.fasta(fname = gzfile(strsplit(n.seq.type, 
                "=")[[1]][2]))
        }
        seqs <- seqs[rows]
        names(seqs) <- toupper(rows)
    }
    else if (substr(n.seq.type, 1, 8) == "csvfile=") {
        if (!is.null(genome.info$all.upstream.seqs[[seq.type]]) && 
            length(genome.info$all.upstream.seqs[[seq.type]]) > 
                0) {
            seqs <- genome.info$all.upstream.seqs[[seq.type]]
        }
        else {
            tab <- read.csv(gzfile(strsplit(n.seq.type, "=")[[1]][2]), 
                head = F)
            seqs <- as.character(tab[, 2])
            names(seqs) <- toupper(as.character(tab[, 1]))
        }
        seqs <- seqs[rows]
        names(seqs) <- toupper(rows)
    }
    else {
        if (is.null(genome.info$feature.tab) || !"genome.seqs" %in% 
            names(genome.info) || is.null(genome.info$genome.seqs)) {
            if (!is.null(genome.info$all.upstream.seqs[[seq.type]])) {
                out.seqs <- genome.info$all.upstream.seqs[[seq.type]][rows]
                out.seqs <- out.seqs[!is.na(out.seqs) & out.seqs != 
                  ""]
                if (filter) 
                  out.seqs <- filter.sequences(out.seqs, NULL, 
                    seq.type, distance, verbose = verbose, ...)
                return(invisible(out.seqs))
            }
            else {
                stop("Motif searching is on but no ", seq.type, 
                  " sequences!")
            }
        }
        if (missing(op.shift)) {
            if (is.na(seq.type) || seq.type == "gene") 
                op.shift <- FALSE
            else op.shift <- operon.shift[seq.type]
            if (is.na(op.shift)) 
                op.shift <- operon.shift[1]
            if (n.seq.type %in% c("gene", "upstream.noncod", 
                "upstream.noncod.same.strand")) 
                op.shift <- FALSE
        }
        coos <- get.gene.coords(rows, op.shift = op.shift)
        if (is.null(coos) || nrow(coos) <= 0) 
            return(NULL)
        coos <- subset(coos, !is.na(start_pos) & !is.na(end_pos))
        if (is.null(coos) || nrow(coos) <= 0) 
            return(NULL)
        seqs <- character()
        if (n.seq.type %in% c("upstream.noncod", "upstream.noncod.same.strand")) {
            all.coos <- genome.info$feature.tab[, c("id", "name", 
                "contig", "strand", "start_pos", "end_pos")]
            all.coos <- subset(all.coos, grepl("^NP_", id, perl = T))
        }
        mc <- get.parallel(nrow(coos))
        tmp <- mc$apply(1:nrow(coos), function(i) {
            if (n.seq.type == "gene") {
                st.st <- coos[i, c("start_pos", "end_pos"), drop = F]
            }
            else if (n.seq.type == "upstream") {
                st.st <- if (coos$strand[i] == "D") 
                  c(coos$start_pos[i] - 1 - distance[2], coos$start_pos[i] - 
                    1 - distance[1])
                else c(coos$end_pos[i] + 1 + distance[1], coos$end_pos[i] + 
                  1 + distance[2])
            }
            else if (n.seq.type == "downstream") {
                st.st <- if (coos$strand[i] == "D") 
                  c(coos$end_pos[i] + 1 + distance[1], coos$end_pos[i] + 
                    1 + distance[2])
                else c(coos$start_pos[i] - 1 - distance[2], coos$start_pos[i] - 
                  1 - distance[1])
            }
            else if (n.seq.type %in% c("upstream.noncod", "upstream.noncod.same.strand")) {
                cc <- all.coos[as.character(all.coos$contig) == 
                  as.character(coos$contig[i]) & abs(all.coos$start_pos - 
                  coos$start_pos[i]) <= 1e+05, ]
                if (n.seq.type == "upstream.noncod.same.strand") 
                  cc <- all.coos[as.character(all.coos$strand) == 
                    as.character(coos$strand[i]), ]
                if (coos$strand[i] == "D") {
                  nearest <- max(cc$end_pos[cc$end_pos < coos$start_pos[i]])
                  st.st <- c(nearest, coos$start_pos[i] - distance[1] - 
                    1)
                }
                else if (coos$strand[i] == "R") {
                  nearest <- min(cc$start_pos[cc$start_pos > 
                    coos$end_pos[i]])
                  st.st <- c(coos$end_pos[i] + distance[1] + 
                    1, nearest)
                }
            }
            seq <- substr(genome.info$genome.seqs[[as.character(coos$contig[i])]], 
                st.st[1], st.st[2])
            if (coos$strand[i] == "R") 
                seq <- rev.comp(seq)
            if (seq.type != "gene" && nchar(seq) > abs(diff(distance))) {
                if (coos$strand[i] == "D") 
                  seq <- substr(seq, 1, abs(diff(distance)))
                else seq <- rev.comp(substr(rev.comp(seq), 1, 
                  abs(diff(distance))))
            }
            out <- list(seq = seq, name = as.character(coos$names[i]), 
                start.stops = data.frame(start = st.st[1], end = st.st[2], 
                  strand = as.character(coos$strand[i]), contig = as.character(coos$contig[i])))
            out
        })
        for (i in tmp) {
            seqs[i$name] <- i$seq
            start.stops <- rbind(start.stops, i$start.stops)
            rownames(start.stops)[nrow(start.stops)] <- i$name
        }
        rownames(start.stops) <- names(seqs) <- make.unique(rownames(start.stops))
        rows <- rows[rows %in% names(seqs)]
        start.stops <- start.stops[rows, , drop = F]
        seqs <- seqs[rows]
        names(seqs) <- rownames(start.stops) <- rows
    }
    if (any(is.na(seqs))) {
        warning("Warning: could not find '", n.seq.type, "' sequences for all input genes", 
            call. = F)
        if (!is.null(start.stops)) 
            start.stops <- start.stops[!is.na(seqs), ]
        seqs <- seqs[!is.na(seqs)]
    }
    if (filter) 
        seqs <- filter.sequences(seqs, start.stops, seq.type, 
            distance, verbose = verbose, ...)
    attr(seqs, "start.stops") <- start.stops
    invisible(seqs)
}
get.stats <-
function (mean.func = mean) 
{
    changed <- NA
    if (!exists("row.memb")) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        if (is.vector(row.memb)) 
            row.memb <- t(row.memb)
        rownames(row.memb) <- attr(ratios, "rnames")
        col.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "cnames") %in% get.cols(k))
        if (is.vector(col.memb)) 
            col.memb <- t(col.memb)
        rownames(col.memb) <- attr(ratios, "cnames")
    }
    if (!exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
    }
    cs <- as.list(clusterStack)
    resids <- sapply(cs, "[[", "resid")
    if (is.matrix(resids)) 
        resids <- apply(resids, 1, function(r) mean.func(r[r != 
            1], na.rm = T))
    else resids <- mean.func(resids[resids != 1], na.rm = T)
    p.clusts <- sapply(cs, "[[", "p.clust")
    if (is.list(p.clusts)) 
        p.clusts <- sapply(cs[!sapply(p.clusts, is.null)], "[[", 
            "p.clust")
    if (is.matrix(p.clusts)) 
        p.clusts <- apply(p.clusts, 1, mean.func, na.rm = T)
    else p.clusts <- mean.func(p.clusts, na.rm = T)
    out <- data.frame(iter = iter, changed = changed, row.scores = mean.func(row.scores[, 
        ][row.memb[, ]], na.rm = T), col.scores = mean.func(col.scores[, 
        ][col.memb[, ]], na.rm = T), mot.scores = if (!is.null(mot.scores)) 
        mean.func(mot.scores[, ][row.memb[, ]], na.rm = T)
    else NA, net.scores = if (!is.null(net.scores)) 
        mean.func(net.scores[, ][row.memb[, ]], na.rm = T)
    else NA, resid = weighted.mean(resids, row.weights, na.rm = T), 
        nrow = mean.func(sapply(cs, "[[", "nrows"), na.rm = T), 
        ncol = mean.func(sapply(cs, "[[", "ncols"), na.rm = T), 
        p.clust = if (!all(is.na(p.clusts))) 
            weighted.mean(p.clusts, mot.weights, na.rm = T)
        else NA)
    if (length(resids) > 1) 
        for (i in names(resids)) {
            out <- cbind(out, resids[i])
            names(out)[ncol(out)] <- paste("resid", i, sep = ".")
        }
    if (length(p.clusts) > 1) 
        for (i in names(p.clusts)) {
            out <- cbind(out, p.clusts[i])
            names(out)[ncol(out)] <- paste("p.clust", i, sep = ".")
        }
    if (length(networks) > 1) {
        for (i in names(net.weights)) {
            if (exists("cluster.net.scores") && i %in% colnames(cluster.net.scores)) 
                out <- cbind(out, weighted.mean(cluster.net.scores[, 
                  i], sapply(cs, "[[", "nrows"), na.rm = T))
            else out <- cbind(out, rep(NA, nrow(out)))
            names(out)[ncol(out)] <- paste("net", i, sep = ".")
        }
        if (exists("cluster.net.scores") && "net.scores" %in% 
            colnames(cluster.net.scores)) {
            out[, "net.scores"] <- weighted.mean(cluster.net.scores[, 
                "net.scores"], sapply(cs, "[[", "nrows"), na.rm = T)
        }
        else {
            out[, "net.scores"] <- mean.func(net.scores[, ][row.memb[, 
                ]], na.rm = T)
        }
    }
    out
}
get.STRING.links <-
function (org.id = genome.info$org.id$V1[1], detailed = T) 
{
    url <- sprintf("http://baliga.systemsbiology.net/cmonkey/data/STRING/%d_STRING.tsv.gz", 
        org.id)
    fname <- sprintf("data/%s/string_links.tab.gz", rsat.species)
    if ((!file.exists(fname) || file.info(fname)$size <= 10)) {
        err <- dlf(fname, url, paste("Fetching STRING protein links file", 
            url, "\n"))
    }
    if (file.exists(fname) && file.info(fname)$size > 10) {
        cat("Loading EMBL STRING interaction links from local file", 
            fname, "\n")
        string.links <- read.delim(gzfile(fname), sep = "\t", 
            head = F)
        colnames(string.links) <- c("protein1", "protein2", "combined_score")
    }
    else {
        string.links <- get.string.links.NEW(org.id)
    }
    closeAllConnections()
    invisible(string.links)
}
get.STRING.links.NEW <-
function (org.id = genome.info$org.id$V1[1], all.genes = attr(ratios, 
    "rnames"), score = "score", min.score = 2, string.url = "http://string-db.org/") 
{
    if (file.exists(paste("data/", rsat.species, "/string_links_FALSE_", 
        org.id, ".tab", sep = ""))) {
        string.links <- read.delim(paste("data/", rsat.species, 
            "/string_links_FALSE_", org.id, ".tab", sep = ""), 
            head = T, sep = " ")
        string.links$protein1 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein1)
        string.links$protein2 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein2)
        return(string.links)
    }
    file <- sprintf("data/%s/string_links_%s.tab", rsat.species, 
        org.id)
    proc.string.df <- function(file) {
        err <- try(tmp <- unique(read.delim(file, head = F, sep = "")))
        if ("try-catch" %in% class(err)) 
            return(NULL)
        if (!exists("tmp") || nrow(tmp) <= 0) 
            return(NULL)
        tmp2 <- strsplit(as.character(tmp$V15), "[:|]", perl = T)
        tmp2a <- sapply(tmp2, function(i) which(i == score))
        tmp2b <- sapply(1:length(tmp2), function(i) if (length(tmp2a[[i]]) == 
            0) 
            NA
        else as.numeric(tmp2[[i]][tmp2a[[i]] + 1]))
        string.links <- data.frame(protein1 = gsub(paste("string:", 
            org.id, ".", sep = ""), "", tmp$V1), protein2 = gsub(paste("string:", 
            org.id, ".", sep = ""), "", tmp$V2), combined_score = tmp2b)
        string.links <- unique(subset(string.links, !is.na(combined_score)))
        string.links
    }
    string.links <- NULL
    tried <- character()
    if (file.exists(file)) 
        string.links <- proc.string.df(file)
    if (file.exists(sprintf("%s.tried", file))) 
        tried <- readLines(sprintf("%s.tried", file))
    tmp2 <- all.genes %in% tried
    if (!file.exists(file) || (!is.null(string.links) && any(!tmp2))) {
        if (!is.null(string.links)) 
            all.genes <- all.genes[!tmp2]
        id.file <- tempfile()
        options(timeout = 300)
        for (i in seq(1, length(all.genes), by = 100)) {
            ids <- all.genes[i:min(i + 99, length(all.genes))]
            ids <- c(ids, paste(org.id, all.genes[i:min(i + 99, 
                length(all.genes))], sep = "."))
            cat(i, "of", length(all.genes), "\n")
            if (org.id == 3702) {
                url <- paste(string.url, "api/tsv/resolveList?caller_identity=cMonkey&identifiers=", 
                  URLencode(paste(ids, collapse = "\r"), reserved = T), 
                  sep = "")
                dlf(id.file, url, mode = "wb", quiet = T)
                ids <- unique(as.character(read.delim(id.file)$stringId))
                unlink(id.file)
            }
            url <- paste(string.url, "api/psi-mi-tab/interactionsList?required_score=", 
                min.score, "&caller_identity=cMonkey&network_graph=2&limit=99999&identifiers=", 
                URLencode(paste(ids, collapse = "\r"), reserved = T), 
                "&species=", org.id, sep = "")
            if (!file.exists(file)) 
                dlf(file, url, mode = "wb", msg = "Fetching STRING protein links (piecewise)... this may take a while...", 
                  quiet = F)
            else dlf(file, url, mode = "ab", quiet = F)
        }
        string.links <- proc.string.df(file)
        writeLines(unique(c(all.genes, tried)), sprintf("%s.tried", 
            file))
        options(timeout = 60)
    }
    invisible(string.links)
}
get.synonyms <-
function (gns, ft = genome.info$feature.names, ignore.case = T, 
    verbose = F, fast = F, force = F) 
{
    gns.input <- gns
    if (exists("no.genome.info") && no.genome.info) {
        out <- as.list(gns)
        names(out) <- gns
        return(out)
    }
    out <- list()
    if ((!force && exists("genome.info") && !is.null(genome.info$synonyms))) {
        gns.cached <- gns[gns %in% names(genome.info$synonyms)]
        out <- genome.info$synonyms[gns.cached]
        gns <- gns[!gns %in% names(genome.info$synonyms)]
        if (length(gns) <= 0 || (is.null(ft) && (!exists("translation.tab") || 
            is.null(translation.tab)))) 
            return(out)
    }
    tmp.out <- as.list(gns)
    names(tmp.out) <- gns
    if (is.null(ft) && (!exists("translation.tab") || is.null(translation.tab))) 
        return(c(out, tmp.out))
    gns.orig <- gns
    gns <- gsub("m$|_\\d$|\\-\\S$", "", gns, perl = T)
    gns <- gsub("([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\\\\\1", 
        gns, perl = T)
    gns <- gns[!is.na(gns) & gns != ""]
    ft <- ft[, c("id", "names")]
    if (exists("translation.tab") && !is.null(translation.tab)) 
        ft <- rbind(ft, data.frame(id = as.character(translation.tab$V1), 
            names = as.character(translation.tab$V2)))
    ft <- subset(ft, names != "")
    if (verbose) 
        ggggg <- gns[seq(1, length(gns), by = min(length(gns), 
            100))]
    mc <- get.parallel(length(gns), verbose = F)
    tmp <- mc$apply(gns, function(g) {
        if (verbose && g %in% ggggg) 
            cat(" ...", g)
        greg <- paste("^", g, "$", sep = "")
        tmp <- subset(ft, grepl(greg, id, perl = T, ignore = ignore.case) | 
            grepl(greg, names, perl = T, ignore = ignore.case))
        if (nrow(tmp) <= 0) 
            return(g)
        tmp2 <- unique(c(g, as.character(tmp$id), as.character(tmp$names)))
        if (!fast) {
            tmp2 <- subset(ft, id %in% tmp2 | names %in% tmp2)
            tmp2 <- unique(c(g, as.character(tmp2[, 1]), as.character(tmp2[, 
                2])))
        }
        tmp2 <- gsub("\\\\([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", 
            "\\1", tmp2, perl = T)
        tmp2
    })
    names(tmp) <- gns.orig
    if (verbose) 
        cat("\n")
    out <- c(tmp, out)
    out[gns.input[gns.input %in% names(out)]]
}
get.unpreprocessed.ratios <-
function (...) 
{
    return(ratios.raw)
}
get.updated.memberships <-
function (row.membership, col.membership, rr.scores, cc.scores) 
{
    rm <- t(apply(rr.scores, 1, order, decreasing = T)[1:n.clust.per.row, 
        , drop = F])
    rm <- t(apply(rm, 1, sort))
    if (n.clust.per.row == 1) 
        rm <- t(rm)
    if (ncol(rm) < ncol(row.membership)) 
        rm <- cbind(rm, matrix(0, nrow = nrow(rm), ncol = ncol(row.membership) - 
            ncol(rm)))
    for (i in 1:nrow(rm)) {
        if (all(rm[i, ] %in% row.membership[i, ])) 
            next
        mc <- max.changes["rows"]
        if (mc < 1 && mc > 0 && runif(1) > mc) 
            next
        for (ii in 1:mc) {
            if (sum(!rm[i, ] %in% row.membership[i, ]) >= mc) {
                if (any(row.membership[i, ] == 0)) {
                  col.change <- which(row.membership[i, ] == 
                    0)[1]
                }
                else {
                  ttmp <- tabulate(row.membership[i, ])
                  if (any(ttmp > 1)) {
                    col.change <- which(row.membership[i, ] %in% 
                      which(ttmp > 1))[1]
                  }
                  else {
                    delta <- rr.scores[i, rm[i, ]] - rr.scores[i, 
                      row.membership[i, ]]
                    if (any(row.membership[i, ] %in% rm[i, ])) 
                      delta[row.membership[i, ] %in% rm[i, ]] <- 0
                    if (all(is.na(delta) | delta <= 0)) 
                      next
                    col.change <- which.max(delta)
                  }
                }
                if (exists("maintain.seed") && !is.null(maintain.seed) && 
                  !is.null(maintain.seed$rows) && !is.null(maintain.seed$rows[[as.character(row.membership[i, 
                  col.change])]]) && rownames(row.membership)[i] %in% 
                  maintain.seed$rows[[as.character(row.membership[i, 
                    col.change])]]) 
                  next
                if (!rm[i, col.change] %in% row.membership[i, 
                  ]) 
                  row.membership[i, col.change] <- rm[i, col.change]
            }
        }
    }
    if (!is.null(cc.scores)) {
        cm <- t(apply(cc.scores, 1, order, decreasing = T)[1:n.clust.per.col, 
            , drop = F])
        if (ncol(cm) < ncol(col.membership)) 
            cm <- cbind(cm, matrix(0, nrow = nrow(cm), ncol = ncol(col.membership) - 
                ncol(cm)))
        for (i in 1:nrow(cm)) {
            mc <- max.changes["cols"]
            if (mc < 1 && mc > 0 && runif(1) > mc) 
                next
            for (ii in 1:mc) {
                if (sum(!cm[i, ] %in% col.membership[i, ]) >= 
                  mc) {
                  if (any(col.membership[i, ] == 0)) {
                    col.change <- which(col.membership[i, ] == 
                      0)[1]
                  }
                  else {
                    ttmp <- tabulate(col.membership[i, ])
                    if (any(ttmp > 1)) {
                      col.change <- which(col.membership[i, ] %in% 
                        which(ttmp > 1))[1]
                    }
                    else {
                      delta <- cc.scores[i, cm[i, ]] - cc.scores[i, 
                        col.membership[i, ]]
                      if (all(is.na(delta) | delta <= 0)) 
                        next
                      col.change <- which.max(delta)
                    }
                  }
                  if (exists("maintain.seed") && !is.null(maintain.seed) && 
                    !is.null(maintain.seed$cols) && !is.null(maintain.seed$cols[[as.character(col.membership[i, 
                    col.change])]]) && colnames(col.membership)[i] %in% 
                    maintain.seed$cols[[as.character(col.membership[i, 
                      col.change])]]) 
                    next
                  col.membership[i, col.change] <- cm[i, col.change]
                }
            }
        }
    }
    invisible(list(r = row.membership, c = col.membership))
}
id.duplicate.clusters <-
function (scores = r.scores, cor.cutoff = 0.9) 
{
    cors <- cor(scores[, ], use = "pairwise", method = "pearson")
    cors[lower.tri(cors, diag = T)] <- NA
    tmp <- which(cors >= cor.cutoff, arr = T)
    cbind(tmp, cors[tmp])
}
install.binaries <-
function (meme.version = "4.3.0", url = "ftp://ftp.ebi.edu.au/pub/software/MEME/4.3.0/meme_4.3.0.tar.gz", 
    make = "make -j 4", path = system.file(package = "cMonkey")) 
{
    cwd <- setwd(path)
    on.exit(setwd(cwd))
    if (!exists("progs")) 
        dir.create("progs")
    setwd("progs/")
    cMonkey:::dlf(sprintf("meme_%s.tar.gz", meme.version), url)
    system(sprintf("tar -xzf meme_%s.tar.gz", meme.version))
    unlink(sprintf("meme_%s.tar.gz", meme.version))
    setwd(sprintf("meme_%s", meme.version))
    dir.create("local")
    system(sprintf("./configure --prefix=%s/local/ --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial --enable-build-libxml2 --enable-build-libxslt --disable-shared --enable-static --with-gnu-ld", 
        getwd()))
    system(make)
    system("make install")
    setwd("..")
    system(sprintf("ln -s meme_%s/local/bin/meme", meme.version))
    system(sprintf("ln -s meme_%s/local/bin/mast", meme.version))
    system(sprintf("ln -s meme_%s/local/bin/dust", meme.version))
    setwd(cwd)
}
load.data <-
function (org.code = NULL) 
{
    try(dlf("data/NCBI/taxdump.tar.gz", sprintf("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")))
    system("cd data/NCBI; tar -xvzf taxdump.tar.gz names.dmp; gzip -v names.dmp")
    if (!is.null(org.code)) {
        lines <- system(sprintf("gunzip -c data/NCBI/names.dmp.gz | grep \"%s\"", 
            as.character(org.code)), intern = T)
    }
    lines <- do.call(rbind, lapply(strsplit(lines, "\t"), function(i) i[i != 
        "" & i != "|"]))[, 1:2]
    while (length(unique(lines[, 1])) > 1) {
    }
    taxon.id <- as.integer(unique(lines[, 1]))
    org.name <- unique(lines[, 2])
    try(dlf(sprintf("%s/genome.txt", data.dir), sprintf("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%s;export=genome", 
        taxon.id)))
    system(sprintf("gzip -v %s/genome.txt", data.dir))
    genome <- read.fasta(gzfile(sprintf("%s/genome.txt.gz", data.dir)))
    data.dir <- sprintf("data/%s/", gsub("[.,; ]", "_", org.name[1]))
    try(dlf(sprintf("%s/genome_info.tsv", data.dir), sprintf("http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab", 
        taxon.id)))
    system(sprintf("gzip -v %s/genome_info.tsv", data.dir))
    genome.info <- read.delim(gzfile(sprintf("%s/genome_info.tsv.gz", 
        data.dir)), head = T)
    invisible(list(taxon.id = taxon.id, org.name = org.name, 
        genome = genome, genome.info = genome.info))
}
load.genome.info.MicrobesOnline <-
function (id = taxon.id) 
{
    f <- sprintf("data/%s/microbesOnlineGenomeInfo_%d.tsv", rsat.species, 
        id)
    if (!file.exists(sprintf("%s.gz", f))) {
        try(dlf(f, sprintf("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab", 
            id)))
        system(sprintf("gzip -fv %s", f))
    }
    cat("Loading:", f, "\n")
    out <- read.delim(gzfile(sprintf("%s.gz", f)), row.names = 1, 
        sep = "\t", as.is = T, head = T)
    out
}
load.ratios <-
function (ratios) 
{
    if (is.null(ratios)) 
        return(NULL)
    if (is.character(ratios) && file.exists(ratios)) {
        cat("Loading ratios file", ratios, "\n")
        ratios <- read.delim(file = gzfile(ratios), sep = "\t", 
            as.is = T, header = T, check.names = F)
    }
    if (is.matrix(ratios) || is.data.frame(ratios)) {
        if (class(ratios[, 1]) == "character") {
            ratios <- ratios[!duplicated(ratios[, 1]), ]
            rownames(ratios) <- attr(ratios, "rnames") <- ratios[, 
                1]
            ratios <- ratios[, -1]
        }
        if (class(ratios[, 1]) == "character") 
            ratios <- ratios[, -1]
    }
    cat("Original ratios matrix is", paste(dim(ratios), collapse = "x"), 
        "\n")
    if (!is.matrix(ratios)) 
        ratios <- as.matrix(ratios)
    if (is.null(attr(ratios, "isPreProcessed")) || attr(ratios, 
        "isPreProcessed") == FALSE) {
        ratios <- preprocess.ratios(ratios)
        attr(ratios, "isPreProcessed") <- TRUE
    }
    cat("Processed ratios matrix is", paste(dim(ratios), collapse = "x"), 
        "\n")
    closeAllConnections()
    ratios
}
make.pv.ev.matrix <-
function (out.ms, make.ev = F) 
{
    mot.rows <- character()
    for (k in 1:k.clust) {
        if (is.null(out.ms[[k]]$pv.ev) || length(out.ms[[k]]$pv.ev) <= 
            0) 
            next
        mot.rows <- unique(c(mot.rows, rownames(out.ms[[k]]$pv.ev[[1]])))
    }
    mot.rows <- sort(mot.rows)
    out.pv <- out.ev <- NULL
    for (k in 1:k.clust) {
        m <- out.ms[[k]]
        if (is.null(m) || is.null(m$pv.ev) || length(m$pv.ev) <= 
            0) {
            out.pv <- cbind(out.pv, rep(NA, length(mot.rows)))
            if (make.ev) 
                out.ev <- cbind(out.ev, rep(NA, length(mot.rows)))
        }
        else {
            m.scores <- numeric(length = length(mot.rows))
            tmp <- m$pv.ev[[1]][, "p.value"]
            names(tmp) <- rownames(m$pv.ev[[1]])
            m.scores <- tmp[mot.rows]
            out.pv <- cbind(out.pv, m.scores)
            colnames(out.pv) <- NULL
            if (make.ev) {
                m.scores <- numeric(length = length(mot.rows))
                tmp <- m$pv.ev[[1]][, "e.value"]
                names(tmp) <- rownames(m$pv.ev[[1]])
                m.scores <- tmp[mot.rows]
                out.ev <- cbind(out.ev, m.scores)
                colnames(out.ev) <- NULL
            }
            out.ms[[k]]$pv.ev[[1]] <- NULL
        }
    }
    rownames(out.pv) <- mot.rows
    if (!is.null(out.pv)) 
        rownames(out.pv) <- mot.rows
    out.pv
}
matrix.reference <-
function (m, ...) 
{
    return(m)
}
meme.one.cluster <-
function (k, seq.type = names(mot.weights)[1], verbose = F, force = F, 
    keep.meme.out = F, keep.mast.out = F, ...) 
{
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    meme.out <- mast.out <- NULL
    cmd <- sprintf(meme.cmd[seq.type], n.motifs[[seq.type]][iter])
    pal.opt <- "non"
    if ("pal.opt" %in% names(list(...))) 
        pal.opt <- list(...)$pal.opt
    get.meme.consensus <- function(cmd, min.iter = 500, max.e.value = 0.1, 
        ...) {
        if (grepl("-cons $compute", cmd, fixed = T)) {
            if (iter > min.iter && !is.null(ms) && !is.null(ms$meme.out)) {
                e.val <- sapply(1:length(ms$meme.out), function(i) ms$meme.out[[i]]$e.value)
                if (min(e.val, na.rm = T) < max.e.value) {
                  best <- which.min(e.val)
                  consensus <- toupper(pssm.to.string(ms$meme.out[[best]]$pssm))
                  cmd <- gsub("$compute", consensus, cmd, fixed = T)
                }
            }
        }
        if (grepl("-cons $compute", cmd, fixed = T)) 
            cmd <- gsub("-cons $compute", "", cmd, fixed = T)
        cmd
    }
    ms <- NULL
    if (is.numeric(k)) {
        ms <- try(meme.scores[[seq.type]][[k[1]]])
        if ("try-error" %in% class(ms)) 
            ms <- NULL
    }
    if (!is.null(ms)) 
        cmd <- get.meme.consensus(cmd, ...)
    if (grepl("-cons $none", cmd, fixed = T)) 
        cmd <- gsub("-cons $none", "", cmd, fixed = T)
    if (grepl("-cons $compute", cmd, fixed = T)) 
        cmd <- gsub("-cons $compute", "", cmd, fixed = T)
    bg.list <- genome.info$bg.list[[seq.type]]
    bg.fname <- genome.info$bg.fname[seq.type]
    bgo <- bg.order[seq.type]
    if (TRUE && !force && !is.null(ms) && !is.null(ms$prev.run) && 
        length(ms$prev.run$rows) == length(rows) && all(ms$prev.run$rows %in% 
        rows) && all(rows %in% ms$prev.run$rows) && cmd == ms$prev.run$cmd && 
        (ms$prev.run$bg.order == bgo || sum(is.na(c(ms$prev.run$bg.order, 
            bgo) == 1))) && all(motif.upstream.scan[[seq.type]] == 
        ms$prev.run$m.u.scan)) {
        message("SKIPPING CLUSTER (UNCHANGED): ", k)
        ms$iter = iter
        ms$last.run = TRUE
        return(ms)
    }
    seqs <- get.sequences(rows, seq.type = seq.type, verbose = verbose, 
        ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k, last.run = FALSE))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k, last.run = FALSE))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    if (is.null(all.seqs)) 
        all.seqs <- get.sequences("all", seq.type = seq.type, 
            distance = motif.upstream.scan[[seq.type]], filter = F, 
            ...)
    if (!is.na(bgo) && grepl("-bfile $bgFname", cmd, fixed = T)) {
        if (!recalc.bg && !file.exists(bg.fname)) {
            genome.info$bg.fname[seq.type] <- bg.fname <- my.tempfile("meme.tmp", 
                suf = ".bg")
            capture.output(genome.info$bg.list[[seq.type]] <- mkBgFile(all.seqs, 
                order = bgo, bgfname = bg.fname, input.list = bg.list, 
                use.rev.comp = grepl("-revcomp", meme.cmd[seq.type])))
        }
        if (recalc.bg || is.null(bg.list)) {
            tmp.seqs <- all.seqs[!names(all.seqs) %in% rows]
            bg.fname <- my.tempfile("meme.tmp", suf = ".bg")
            capture.output(bg.list <- mkBgFile(tmp.seqs, order = bg.order[seq.type], 
                verbose = F, bgfname = bg.fname, use.rev.comp = grepl("-revcomp", 
                  meme.cmd[seq.type])))
            rm(tmp.seqs)
        }
    }
    if (is.na(bgo) || is.null(bg.list)) {
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
        bg.list <- NULL
    }
    psps <- NULL
    run.meme <- function(sgenes, seqs, cmd, seq.type, ...) {
        if (pal.opt == "non") {
            out <- runMeme(sgenes, seqs, cmd, seq.type = seq.type, 
                ...)
        }
        out
    }
    if (verbose) {
        meme.out <- run.meme(names(seqs), seqs, cmd, verbose = verbose, 
            bg.list = bg.list, bgfname = bg.fname, psps = psps, 
            seq.type = seq.type, ...)
    }
    else {
        capture.output(meme.out <- try(run.meme(names(seqs), 
            seqs, cmd, verbose = verbose, bg.list = bg.list, 
            bgfname = bg.fname, psps = psps, seq.type = seq.type, 
            ...)))
    }
    meme.out2 <- getMemeMotifInfo(meme.out)
    attr(meme.out2, "meme.command.line") <- attr(meme.out, "meme.command.line")
    attr(meme.out2, "is.pal") <- pal.opt == "pal"
    if (length(meme.out2) <= 0) 
        return(list(k = k, last.run = FALSE))
    if (is.na(bgo) || is.null(bg.list) || !file.exists(bg.fname)) 
        bg.fname <- NULL
    if (verbose) 
        mast.out <- try(runMast(meme.out, mast.cmd[seq.type], 
            names(all.seqs), all.seqs, verbose = verbose, seq.type = seq.type, 
            bg.list = bg.list, bgfname = bg.fname, ...))
    else capture.output(mast.out <- try(runMast(meme.out, mast.cmd[seq.type], 
        names(all.seqs), all.seqs, verbose = verbose, seq.type = seq.type, 
        bg.list = bg.list, bgfname = bg.fname, ...)))
    pv.ev <- get.pv.ev.single(mast.out, rows)
    if (recalc.bg && !is.null(bg.fname) && file.exists(bg.fname) && 
        !"unlink" %in% names(list(...))) 
        unlink(bg.fname)
    prev.run <- list(rows = rows, cmd = cmd, bg.order = bgo, 
        m.u.scan = motif.upstream.scan[[seq.type]])
    out <- list(k = k, last.run = FALSE, meme.out = meme.out2, 
        pv.ev = pv.ev, prev.run = prev.run)
    if (keep.meme.out) 
        out$meme.out.orig <- meme.out
    if (keep.mast.out) 
        out$mast.out.orig <- mast.out
    invisible(out)
}
mkBgFile <-
function (bgseqs = NULL, order = 0, bgfname = NULL, input.list = NULL, 
    use.rev.comp = T, verbose = T) 
{
    if (!is.null(input.list) && !is.null(bgfname)) {
        tmp <- unlist(input.list[2:length(input.list)])
        tmp2 <- sprintf("%.8f", tmp)
        names(tmp2) <- names(tmp)
        write.table(tmp2, row.names = names(tmp2), col.names = paste("#", 
            order, "th order Markov background model"), quote = F, 
            file = bgfname)
        return(input.list)
    }
    repl <- list(R = c("G", "A"), Y = c("T", "C"), K = c("G", 
        "T"), M = c("A", "C"), S = c("G", "C"), W = c("A", "T"), 
        N = c("G", "A", "T", "C"))
    bad.seqs <- grep("[^GATCX]", bgseqs, perl = T)
    if (length(bad.seqs) > 0) {
        if (verbose) 
            message(length(bad.seqs), " sequences with degenerate residues...fixing.")
        for (i in bad.seqs) {
            tmp <- strsplit(bgseqs[i], character(0))[[1]]
            inds <- grep("[^GATCX]", tmp, perl = T)
            for (ind in inds) tmp[ind] <- sample(repl[[tmp[ind]]], 
                1)
            bgseqs[i] <- paste(tmp, collapse = "")
        }
    }
    if (verbose) 
        cat("Calculating", order, "th order background Markov model from", 
            length(bgseqs), "sequences\n")
    if (use.rev.comp && verbose) 
        cat("Using reverse-complement too.\n")
    if (use.rev.comp) 
        bgseqs <- c(bgseqs, rev.comp(bgseqs))
    bgseqs <- bgseqs[!get.dup.seqs(bgseqs)]
    mc <- get.parallel(order + 1)
    apply.func <- lapply
    tmp <- mc$apply(0:order, function(ord, mc.cores) {
        out <- list()
        if (verbose) 
            cat("Calculating", ord, "th order part of background Markov model from", 
                length(bgseqs), "sequences\n")
        if (ord == 0) {
            all.substrings <- unlist(strsplit(bgseqs, character(0)), 
                use.names = F)
        }
        else {
            all.substrings <- sapply(1:(max(nchar(bgseqs)) - 
                ord), function(i) substr(bgseqs, i, i + ord))
            all.substrings <- as.vector(all.substrings)
        }
        all.substrings <- all.substrings[!is.na(all.substrings) & 
            all.substrings != "" & nchar(all.substrings) == ord + 
            1]
        counts <- table(as.factor(all.substrings))
        tmp <- names(counts)
        counts <- as.vector(counts)
        names(counts) <- tmp
        if (length(counts) < 4^(ord + 1)) {
            all.pairs <- unique(combn(rep(col.let, ord + 1), 
                ord + 1, FUN = paste, sep = "", collapse = ""))
            all.pairs <- all.pairs[!all.pairs %in% names(counts)]
            for (let in all.pairs) counts[let] <- 0.1
        }
        counts <- sort(counts)
        counts <- counts/length(all.substrings)
        counts <- counts[grep("N", names(counts), val = T, invert = T)]
        out <- as.list(counts)
        for (i in names(out)) {
            names(out[[i]]) <- NULL
            if (verbose && ord <= 3) 
                cat("FREQ:", i, "=", counts[i], "\n")
        }
        out
    }, mc.cores = min(order + 1, mc$par))
    out <- list()
    out$order <- order
    for (i in 1:length(tmp)) for (j in 1:length(tmp[[i]])) out[[names(tmp[[i]])[j]]] <- tmp[[i]][[j]]
    if (!is.null(bgfname) && !file.exists(bgfname)) {
        cat("Writing to file:", bgfname, "\n")
        tmp <- unlist(out)
        tmp <- tmp[2:length(tmp)]
        tmp2 <- sprintf("%.8f", tmp)
        names(tmp2) <- names(out)[2:length(out)]
        write.table(tmp2, row.names = names(tmp2), col.names = paste("#", 
            order, "th order Markov background model"), quote = F, 
            file = bgfname)
    }
    invisible(out)
}
mkTempMemeFiles <-
function (sgenes, seqs, fname = "meme.tmp.fst", bgseqs = NULL, 
    bgfname = NULL, filter.seqs = T, bg.list = NULL, force.overwrite = F, 
    seq.type = names(mot.weights)[1], seq.weights = NULL, psps = NULL, 
    ...) 
{
    if (!file.exists(fname) || file.info(fname)$size == 0 || 
        force.overwrite) {
        sgenes <- sgenes[!(is.na(seqs) | is.null(seqs) | seqs == 
            "")]
        seqs <- seqs[!(is.na(seqs) | is.null(seqs) | seqs == 
            "")]
        max.width <- as.integer(strsplit(meme.cmd[seq.type], 
            " ")[[1]][which(strsplit(meme.cmd[seq.type], " ")[[1]] == 
            "-maxw") + 1])
        if (filter.seqs) {
            sgenes <- sgenes[nchar(seqs) >= max.width]
            seqs <- seqs[nchar(seqs) >= max.width]
        }
        lengths <- sum(nchar(seqs)) + length(seqs) * 3
        if (!is.null(seq.weights)) {
            seq.weights <- seq.weights[sgenes]
            seq.weights[is.na(seq.weights)] <- 0
            cat(paste(">WEIGHTS", paste(seq.weights, collapse = " ")), 
                paste(">", sgenes, "\n", seqs, sep = ""), file = fname, 
                sep = "\n")
        }
        else {
            cat(paste(">", sgenes, "\n", seqs, sep = ""), file = fname, 
                sep = "\n")
        }
    }
    if (force.overwrite || (!is.null(bgfname) && (!file.exists(bgfname) || 
        file.info(bgfname)$size <= 0))) {
        if (!is.null(bg.list)) 
            mkBgFile(input.list = bg.list, order = bg.list$order, 
                bgfname = bgfname)
        else if (!is.null(bgseqs)) 
            mkBgFile(bgseqs, order = 0, bgfname = bgfname)
    }
    length(seqs)
}
motif.all.clusters <-
function (ks = 1:k.clust, seq.type = names(mot.weights)[1], verbose = T, 
    debug = F, ...) 
{
    out.ms <- meme.scores[[seq.type]]
    mc <- get.parallel(length(ks), verbose = T, para.cores = get("parallel.cores.motif"))
    if (any(grepl("foreach", deparse(mc$apply))) && getDoParName() == 
        "doMC") 
        mc$apply <- function(list, FUN, ...) foreach(l = list, 
            .options.multicore = list(preschedule = F, set.seed = T)) %dopar% 
            {
                FUN(l, ...)
            }
    if (!debug) {
        out.ms <- mc$apply(ks, FUN = function(k) try(motif.one.cluster(k, 
            seq.type = seq.type, verbose = F, ...)))
    }
    else {
        message("DEBUG MODE: NOT PARALLELIZING!\n")
        out.ms <- lapply(ks, FUN = function(k) motif.one.cluster(k, 
            seq.type = seq.type, verbose = T, ...))
    }
    out.ms[[k.clust + 1]] <- ""
    for (k in ks) {
        if (length(out.ms) < k || is.null(out.ms[[k]]) || class(out.ms[[k]]) == 
            "try-error" || out.ms[[k]]$k != k || (!is.null(out.ms[[k]]$iter) && 
            out.ms[[k]]$iter != iter)) {
            out <- try(motif.one.cluster(k, seq.type = seq.type, 
                verbose = T, ...))
        }
        else {
            out <- out.ms[[k]]
        }
        if (class(out) == "try-error") 
            out <- try(motif.one.cluster(k, seq.type = seq.type, 
                verbose = T, ...))
        if (class(out) == "try-error" || is.null(out) || out$k != 
            k) {
            message("ERROR on cluster ", k)
            out <- list()
        }
        else if (verbose) {
            cat(iter, k, length(get.rows(k)), seq.type, "\t")
        }
        if (verbose) {
            if (is.null(out) || is.null(out$meme.out)) 
                cat("Inf \n")
            else {
                ind <- 1
                if (!is.null(out$pv.ev)) {
                  if ("p.value" %in% colnames(out$pv.ev[[1]])) 
                    mn <- mean(log10(out$pv.ev[[ind]][rownames(out$pv.ev[[ind]]) %in% 
                      get.rows(k), "p.value"]), na.rm = T)
                  else mn <- mean(log10(out$pv.ev[[ind]]$pvals), 
                    na.rm = T)
                }
                else {
                  mn <- "Inf"
                }
                cat(k, if (attr(out$meme.out, "is.pal")) 
                  "pal"
                else "non", sapply(out$meme.out[1:min(3, length(out$meme.out))], 
                  "[[", "e.value"), mn, "\t", pssm.to.string(out$meme.out[[1]]$pssm), 
                  "\n")
            }
        }
        out$iter <- iter
        out$k <- k
        out.ms[[k]] <- out
    }
    out.ms$all.pv <- make.pv.ev.matrix(out.ms)
    if (FALSE) {
        for (k in 1:k.clust) {
            m <- out.ms[[k]]
            if (!is.null(m) && !is.null(m$pv.ev)) 
                out.ms[[k]]$pv.ev[[1]] <- NULL
        }
    }
    attr(out.ms, "seq.type") <- seq.type
    invisible(out.ms)
}
motif.one.cluster <-
function (k, seq.type = names(mot.weights)[1], verbose = F, force = F, 
    ...) 
{
    st <- strsplit(seq.type, " ")[[1]]
    out <- meme.scores[[seq.type]][[k]]
    if (st[2] == "meme") {
        out <- meme.one.cluster(k, seq.type = seq.type, verbose, 
            force = force, ...)
    }
    invisible(out)
}
my.tempfile <-
function (pattern = "file", tmpdir = tempdir(), suffix = "", 
    n.rnd.char = 20) 
{
    file.path(paste(tmpdir, "/", pattern, "_", paste(sample(c(LETTERS, 
        letters, 0:9, 0:9, 0:9, 0:9), n.rnd.char), collapse = ""), 
        suffix, sep = ""))
}
new.cluster.motifs <-
function (motifs = "ALL", p.cutoff = 1e-06, include.bad = F, 
    n.cutoff = 10, mcl.I = 1.2, mcl.pi = 2, distance.weight.cutoff = 0.999999, 
    mcl.cmd = paste("./progs/mcl-10-201/local/bin/mcl new.mcltemp -o %s --abc", 
        " -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1"), 
    meme.the.clusters = T, improve.the.clusters = F, plot.them = F, 
    in.clusts = NULL) 
{
    in.args <- c(mget(names(formals()), env = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    nc <- sapply(e$genome.info$genome.seqs, nchar)
    if (!is.null(in.clusts)) {
        clusts <- in.clusts
    }
    else {
        if (exists("fimo.out")) 
            p.scans <- subset(fimo.out, `p-value` <= p.cutoff)
        else if (exists("pssm.scans")) 
            p.scans <- subset(pssm.scans, pvals <= p.cutoff)
        else {
            load(sprintf("filehash/fimo_out_%s.RData", as.character(p.cutoff)))
            p.scans <- subset(fimo.out, `p-value` <= p.cutoff)
        }
        if (motifs == "ALL") {
            motifs <- unlist(lapply(1:nrow(motif.widths), function(i) {
                ii <- which(motif.widths[i, ] > 0)
                if (length(ii) <= 0) 
                  return(NULL)
                paste("MOT", i, ii, sep = "_")
            }))
        }
        cat(length(motifs), "motifs\n")
        if (file.exists(sprintf("filehash/new_motif_shadows_%s.RData", 
            as.character(p.cutoff)))) {
            cat("Loading motif distances pre-computed...\n")
            load(sprintf("filehash/new_motif_shadows_%s.RData", 
                as.character(p.cutoff)))
        }
        else {
            m <- mclapply(1:length(motifs), function(i) {
                m.tmp <- lapply(nc, integer)
                names(m.tmp) <- names(nc)
                bi <- as.integer(strsplit(motifs[i], "_")[[1]][2])
                mo <- as.integer(strsplit(motifs[i], "_")[[1]][3])
                scans <- p.scans[J(c(bi, bi), c(mo, -mo))]
                if (nrow(scans) <= 0) 
                  return(lapply(m.tmp, function(i) integer(0)))
                scans <- scans[!is.na(scans$Start), ]
                if (nrow(scans) <= 0) 
                  return(lapply(m.tmp, function(i) integer(0)))
                for (iii in 1:nrow(scans)) {
                  inds <- scans$Start[iii]:scans$Stop[iii]
                  chr <- levels(scans$Seq)[scans$Seq[iii]]
                  m.tmp[[chr]] <- c(m.tmp[[chr]], inds)
                }
                cat(i, length(motifs), motifs[i], nrow(scans), 
                  "\n")
                lapply(m.tmp, unique)
            }, mc.preschedule = F)
            rm(m.tmp, p.scans)
            gc()
            names(m) <- motifs
            save(m, file = sprintf("filehash/new_motif_shadows_%s.RData", 
                as.character(p.cutoff)))
        }
        cat(length(motifs), "motifs\n")
        if (!include.bad && exists("coding.fracs")) {
            frac.in.coding <- coding.fracs$all.fracs[motifs]
            motifs <- motifs[!is.na(frac.in.coding) & frac.in.coding < 
                coding.fracs$mean.fracs - 0.01]
            cat(length(motifs), "motifs\n")
            m <- m[motifs]
            rm(frac.in.coding)
        }
        nc.cumsum <- c(0, cumsum(nc))[1:length(nc)]
        names(nc.cumsum) <- names(nc)
        for (i in 1:length(m)) {
            x <- m[[i]]
            if (length(unlist(x)) <= 0) 
                next
            for (ii in 1:length(x)) x[[ii]] <- x[[ii]] + nc.cumsum[ii]
            m[[i]] <- unlist(x)
            names(m[[i]]) <- NULL
        }
        dirname <- sprintf("filehash/new_motif_shadows_%s_%d", 
            as.character(p.cutoff), length(motifs))
        dir.create(dirname)
        if (!file.exists(sprintf("%s.tsv.bz2", dirname))) {
            tmp <- mclapply(1:(length(m) - 1), function(i) {
                fname <- sprintf("%s/%08d.tsv.bz2", dirname, 
                  i)
                if (file.exists(fname)) 
                  return()
                print(i)
                x <- m[[i]]
                if (length(x) <= 0) 
                  return()
                out <- rep(1, length(m))
                for (j in (i + 1):length(m)) {
                  y <- m[[j]]
                  if (length(y) <= 0) 
                    next
                  tmp <- x %in% y
                  if (!any(tmp)) 
                    next
                  tmp <- (sum(!tmp) + sum(!(y %in% x)))/length(unique(c(x, 
                    y)))
                  out[j] <- tmp
                }
                tmp <- which(out < 1)
                if (length(tmp) <= 0) 
                  return()
                write.table(cbind(i, tmp, out[tmp]), quote = F, 
                  sep = "\t", row.names = F, col.names = F, file = bzfile(fname))
                return()
            }, mc.preschedule = F)
            system(sprintf("find ./%s/ -name '*.tsv.bz2' -print | sort | xargs bunzip2 -c | bzip2 -c >%s.tsv.bz2", 
                dirname, dirname))
        }
        cat("Going to run mcl now...\n")
        system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
            dirname, distance.weight.cutoff))
        system(sprintf("bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
            dirname, distance.weight.cutoff))
        outfile <- sprintf("new.mcltemp.I%s.pi%s", gsub(".", 
            "", sprintf("%.1f", mcl.I), fixed = T), gsub(".", 
            "", sprintf("%.1f", mcl.pi), fixed = T))
        mcl.cmd <- sprintf(mcl.cmd, outfile, mcl.I, mcl.pi)
        print(mcl.cmd)
        mcl.out <- system(mcl.cmd, intern = T, ignore.stderr = F)
        if (!improve.the.clusters) 
            unlink("new.mcltemp")
        print(mcl.out)
        system(sprintf("gzip -fv %s", outfile))
        clusts <- lapply(strsplit(readLines(gzfile(sprintf("%s.gz", 
            outfile))), "\t"), as.integer)
        cat("GOT", length(clusts), "motif clusters with", length(unlist(clusts)), 
            "motifs.\n")
        clusts <- clusts[sapply(clusts, length) >= 3]
        clusts <- lapply(clusts, function(i) motifs[i])
        cat("GOT", sum(sapply(clusts, length) >= 3), "motif clusters (length > 3) with", 
            length(unlist(clusts[sapply(clusts, length) >= 3])), 
            "motifs.\n")
        cat("GOT", sum(sapply(clusts, length) >= 10), "motif clusters (length > 10) with", 
            length(unlist(clusts[sapply(clusts, length) >= 10])), 
            "motifs.\n")
        mc.length <- max(which(sapply(clusts, length) >= n.cutoff))
        if (is.infinite(mc.length)) 
            mc.length <- max(which(sapply(clusts, length) >= 
                3))
        attr(clusts, "mcl.cmd") <- mcl.cmd
        attr(clusts, "mcl.out") <- mcl.out
        attr(clusts, "mc.length") <- mc.length
        attr(clusts, "in.args") <- in.args
        attr(clusts, "motifs") <- motifs
        save(clusts, file = sprintf("filehash/new_motif_shadows_%s_clusts.RData", 
            as.character(p.cutoff)))
    }
    if (improve.the.clusters) {
        out2 <- cluster.the.new.motif.clusters(clusts)
        out$motif.cluster.clusters <- out2$motif.cluster.clusters
    }
    if (meme.the.clusters) {
        full.genome <- paste(e$genome.info$genome.seqs, collapse = "", 
            sep = "")
        seq.type <- 1
        bg.list <- e$genome.info$bg.list[[seq.type]]
        bgo <- e$bg.order[seq.type]
        bg.fname <- e$my.tempfile("meme.tmp", suf = ".bg")
        tmp <- e$mkBgFile(e$genome.info$genome.seqs, order = bgo, 
            bgfname = bg.fname, input.list = bg.list, use.rev.comp = T)
        require(Matrix)
        if (!exists("attrs")) 
            attrs <- list()
        for (i in (length(attrs) + 1):mc.length) {
            print(i)
            if (length(unlist(m[clusts[[i]]])) <= 0) 
                next
            thresh.frac = 10
            tmp <- table(unlist(m[clusts[[i]]]))
            if (length(tmp) <= 1) 
                return(NULL)
            tmp2 <- tmp[tmp >= max(tmp)/thresh.frac]
            if (length(tmp2) <= 1) 
                return(NULL)
            m.tmp <- matrix(0, nrow = sum(nc), ncol = 1)
            m.tmp[as.integer(names(tmp2)), 1] <- as.integer(tmp2)
            starts <- which(m.tmp > 0 & c(0, m.tmp[-length(m.tmp)]) == 
                0)
            cat(i, length(starts), "\n")
            while (length(starts) > 300 || length(starts) < 5) {
                was.toobig <- length(starts) > 300
                if (length(starts) > 300) 
                  thresh.frac <- thresh.frac/1.1
                else if (length(starts) < 5) 
                  thresh.frac <- thresh.frac * 1.1
                tmp2 <- tmp[tmp >= max(tmp)/thresh.frac]
                m.tmp <- matrix(0, nrow = sum(nc), ncol = 1)
                m.tmp[as.integer(names(tmp2)), 1] <- as.integer(tmp2)
                starts <- which(m.tmp > 0 & c(0, m.tmp[-length(m.tmp)]) == 
                  0)
                if (was.toobig && length(starts) < 30) {
                  thresh.frac <- thresh.frac * 1.1
                  tmp2 <- tmp[tmp >= max(tmp)/thresh.frac]
                  m.tmp <- matrix(0, nrow = sum(nc), ncol = 1)
                  m.tmp[as.integer(names(tmp2)), 1] <- as.integer(tmp2)
                  starts <- which(m.tmp > 0 & c(0, m.tmp[-length(m.tmp)]) == 
                    0)
                  break
                  cat(i, thresh.frac, length(starts), "\n")
                }
                gc()
            }
            ends <- which(m.tmp > 0 & c(m.tmp[-1], 0) == 0)
            if (length(ends) != length(starts)) 
                next
            tmp <- which(ends - starts < 10)
            starts[tmp] <- starts[tmp] - 5
            ends[tmp] <- ends[tmp] + 5
            sseqs <- substring(full.genome, starts, ends)
            tmp <- table(nchar(sseqs))
            tmp <- as.integer(names(tmp)[c(min(which(tmp > 1)), 
                max(which(tmp > 3)))])
            if (any(is.na(tmp))) 
                tmp <- range(nchar(sseqs))
            if (tmp[2] > 40) 
                tmp[2] <- 40
            starts <- starts - 10
            ends <- ends + 10
            sseqs <- substring(full.genome, starts, ends)
            meme.out <- "Segmentation fault"
            cat(i, length(sseqs), "\n")
            while (any(grepl("Segmentation fault", meme.out))) {
                meme.out <- e$runMeme(as.character(1:length(sseqs)), 
                  sseqs, verbose = T, cmd = sprintf(paste("./progs/meme $fname -bfile %s -time 600 -dna -revcomp -maxsize", 
                    "9999999 -nmotifs 2 -evt 1e9 -minw %d -maxw %d -mod zoops -nostatus -text"), 
                    bg.fname, tmp[1], tmp[2]), filter = F)
                tmp[2] <- tmp[2] - 2
            }
            meme.out <- e$getMemeMotifInfo(meme.out)
            cat(i, length(sseqs), meme.out[[1]]$e.value, "\n")
            meme.out[[1]]$posns$abs.start <- starts[as.integer(as.character(meme.out[[1]]$posns$gene))] + 
                meme.out[[1]]$posns$start
            attrs[[i]] <- list(meme.out = meme.out, sseqs = sseqs)
        }
        for (i in 1:mc.length) if (!is.null(attrs[[i]]) && class(attrs[[i]]) != 
            "try-error") {
            attr(clusts[[i]], "meme.out") <- attrs[[i]]$meme.out
            attr(clusts[[i]], "sseqs") <- attrs[[i]]$sseqs
        }
        if (get.combined.pssms) {
            for (i in 1:length(clusts)) {
                tmp <- out$cluster.motifs("hclust", motifs = clusts[[i]], 
                  min.gene.overlap = 1, e.value.cutoff = Inf, 
                  p.value.cutoff = 1e-04, resid.cutoff = Inf, 
                  n.cutoff = 1, expand = F, include.bad = T, 
                  find.bad = NA, in.tt.out = NULL, k.cut = 1)
                attr(clusts[[i]], "tt.out") <- tmp$tt.out
                attr(clusts[[i]], "tt.out2") <- tmp$tt.out2[[1]]
            }
        }
    }
    if (plot.them && (meme.the.clusters || get.combined.pssms)) {
        pdf("Rplots.pdf")
        par(mfrow = c(4, 4))
        for (i in 1:mc.length) {
            print(i)
            if (!get.combined.pssms) {
                mo <- attr(clusts[[i]], "meme.out")[[1]]
                e$viewPssm(mo$pssm, main = paste(i, length(clusts[[i]]), 
                  length(attr(clusts[[i]], "sseqs")), mo$e.value))
            }
            else {
                pssm <- attr(attr(clusts[[i]], "tt.out2"), "combined.pssm")
                e$viewPssm(pssm, main = paste(i, length(clusts[[i]])))
            }
        }
        dev.off()
    }
    clusts
}
operon.list <-
function () 
{
    out <- list()
    ops <- as.matrix(genome.info$operons)
    for (i in 1:nrow(ops)) {
        if (!ops[i, 1] %in% names(out)) 
            out[[ops[i, 1]]] <- ops[i, 2]
        else out[[ops[i, 1]]] <- c(out[[ops[i, 1]]], ops[i, 2])
    }
    out
}
plotClust <-
function (k, cluster = NULL, w.motifs = T, all.conds = T, title = NULL, 
    o.genes = NULL, dont.plot = F, network = "all", short.names = organism == 
        "sce", seq.type = names(mot.weights), ...) 
{
    if (!dont.plot && names(dev.cur()) != "devSVG") 
        opar <- par(no.readonly = T)
    if (!is.null(cluster)) {
        if (!dont.plot) 
            plotCluster.motif(cluster, seqs = cluster$seqs, p.val.shade.cutoff = 1, 
                o.genes = o.genes, no.plotCluster = all.conds, 
                ...)
        return(invisible(cluster))
    }
    c <- get.clust(k, varNorm = F)
    rows <- get.rows(k)
    if (!is.null(o.genes)) 
        rows <- unique(c(rows, o.genes))
    if (length(rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    if (!w.motifs && !dont.plot) {
        if (all.conds) 
            plotCluster.all.conds(c, o.genes = o.genes, ...)
        else plotCluster(c, o.genes = o.genes, ...)
    }
    else {
        c$seq.type <- seq.type
        for (st in seq.type) {
            c[[st]] <- list()
            c[[st]]$motif.out <- meme.scores[[st]][[k]]
            tmp <- cluster.pclust(k, st)
            c[[st]]$e.val <- tmp$e.vals
            c[[st]]$p.clust <- tmp$p.clusts
            c[[st]]$motif.out$pssms <- lapply(c[[st]]$motif.out$meme.out, 
                "[[", "pssm")
            c[[st]]$motif.out$e.values <- c[[st]]$e.val
            if (!is.null(c[[st]]$motif.out$pv.ev)) {
                if ("gene" %in% colnames(c[[st]]$motif.out$pv.ev[[1]])) 
                  c[[st]]$motif.out$pv.ev[[2]] <- c[[st]]$motif.out$pv.ev[[1]]
                if (!is.null(meme.scores[[st]]$all.pv)) {
                  tmp <- cbind(p.value = meme.scores[[st]]$all.pv[, 
                    k], e.value = if ("all.ev" %in% names(meme.scores[[st]])) 
                    meme.scores[[st]]$all.ev[, k]
                  else NA)
                }
                else {
                  pv.ev <- meme.scores[[st]][[k]]$pv.ev[[1]]
                  if (ncol(pv.ev) <= 2) 
                    pv.ev <- meme.scores[[st]][[k]]$pv.ev[[2]]
                  tmp <- NULL
                  if (ncol(pv.ev) > 0) {
                    tmp <- as.matrix(pv.ev[, 2:ncol(pv.ev)])
                    rownames(tmp) <- pv.ev[, 1]
                    colnames(tmp) <- c("p.value", "posns", "mots")
                  }
                }
                c[[st]]$motif.out$pv.ev[[1]] <- tmp
                if (!is.null(tmp)) {
                  c[[st]]$motif.out$p.values <- log10(c[[st]]$motif.out$pv.ev[[1]][, 
                    "p.value"])
                  names(c[[st]]$motif.out$p.values) <- rownames(c[[st]]$motif.out$pv.ev[[1]])
                }
            }
        }
    }
    if (!is.na(mot.iters[1]) && !no.genome.info) {
        c$seqs <- get.sequences(rows, distance = motif.upstream.scan[[seq.type[1]]], 
            seq.type = seq.type[1], filter = T, uniq = F)
        tmp <- c$seqs[rows]
        if (!is.null(tmp)) 
            names(tmp) <- rows
        attr(tmp, "start.stops") <- attr(c$seqs, "start.stops")
        c$seqs <- tmp
        rm(tmp)
    }
    else c$seqs <- NULL
    if (!is.na(net.iters[1])) {
        if (network == "all") 
            network <- names(networks)
        for (i in network) {
            if (!i %in% names(networks)) 
                next
            tmp.net <- networks[[i]][networks[[i]]$protein1 %in% 
                rows & networks[[i]]$protein2 %in% rows, ]
            tmp.net <- cbind(tmp.net, net = rep(i, nrow(tmp.net)))
            c$network <- if (!is.null(c$network)) 
                rbind(c$network, tmp.net)
            else tmp.net
        }
    }
    c$gene.coords <- get.long.names(rows, short = short.names)
    if ("cog.code" %in% names(genome.info)) 
        c$cog.code <- genome.info$cog.code[rows]
    if (!is.null(title)) 
        c$name <- title
    if (!dont.plot) {
        plotCluster.motif(c, seqs = c$seqs, p.val.shade.cutoff = 1, 
            o.genes = o.genes, no.plotCluster = all.conds, ...)
        if (names(dev.cur()) != "devSVG" && !"layout" %in% names(list(...))) 
            par(opar)
    }
    invisible(c)
}
plotCluster <-
function (cluster, imag = F, cond.labels = F, o.genes = NULL, 
    col.func = if (imag) topo.colors else rainbow, rats.names = names(ratios), 
    main = NULL, range.r = NULL, no.par = F, sort = F, box.plot = F, 
    ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    if (is.null(main)) 
        main <- paste(sprintf("Cluster: %04d %s; resid: %s; r/c: %d/%d", 
            k, organism, paste(sprintf("%.2f", cluster$resid[rats.names]), 
                collapse = " "), length(cluster$rows), length(cluster$cols)))
    rats <- get.cluster.matrix(unique(c(cluster$rows, o.genes)), 
        cluster$cols, matrices = rats.names)
    cols.b <- colnames(rats)[colnames(rats) %in% cluster$cols]
    if (sort) {
        o1 <- order(apply(rats[cluster$rows, cols.b, drop = F], 
            2, mean, na.rm = T))
        cols.b <- cols.b[o1]
        rats <- rats[, cols.b, drop = F]
    }
    if (all(is.na(rats))) {
        plot(0, 0, typ = "n", min = main, ...)
        return()
    }
    if (is.vector(rats)) {
        rats <- t(rats)
        rownames(rats) <- cluster$rows
    }
    if (imag) {
        grey.image <- function(mat, n.gray = 32, x = 1:nrow(mat), 
            y = 1:ncol(mat), col = gray((0:n.gray)/n.gray), ...) image(x, 
            y, mat, col = col, ...)
        grey.image(t(rats), col = col.func(256))
        return()
    }
    if (is.null(range.r)) 
        range.r <- range(rats[rats != min(rats, na.rm = T) & 
            rats != max(rats, na.rm = T)], na.rm = T)
    if (cond.labels && cluster$ncols < 100) 
        range.r[1] <- range.r[1] * 1.5
    if (!no.par) 
        par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    plot(1:length(cols.b), ylim = range.r, xlab = NA, ylab = NA, 
        main = main, typ = "n", xaxs = "i", ...)
    if (length(rats.names) > 1) {
        ind <- 0.5
        rts <- NULL
        for (i in 1:length(rats.names)) {
            col <- sapply(col2rgb(i + 1)/255 + 0.9, function(cc) min(cc, 
                1))
            col <- rgb(col[1], col[2], col[3])
            rect(ind, range.r[1] + 0.05, ind + sum(colnames(rats) %in% 
                colnames(ratios[[rats.names[i]]])), range.r[2] - 
                0.05, col = col, dens = NA)
            ind <- ind + sum(colnames(rats) %in% colnames(ratios[[rats.names[i]]]))
            rts <- cbind(rts, rats[, cols.b[cols.b %in% colnames(ratios[[rats.names[i]]]), 
                drop = F]])
        }
        rats <- rts
        rm(rts)
        cols.b <- colnames(rats)
    }
    if (exists("col.rug")) {
        if (is.integer(col.rug)) 
            colmap <- col.func(max(col.rug))[col.rug[cols.b]]
        else colmap <- col.rug[cols.b]
    }
    else if (all(deparse(col.func) == deparse(rainbow))) {
        colmap <- col.func(length(cols.b))
    }
    else {
        colmap <- col.func(cols.b)
    }
    if (box.plot) {
        colMeans <- apply(rats[cluster$rows, , drop = F], 2, 
            mean, na.rm = T)
        colSd <- apply(rats[cluster$rows, , drop = F], 2, sd, 
            na.rm = T)
        matlines(1:length(cols.b), cbind(colMeans - 2 * colSd, 
            colMeans + 2 * colSd), lty = 1, col = "lightgrey")
        boxplot(as.data.frame(rats[cluster$rows, , drop = F]), 
            ylim = range.r, names = NA, main = main, col = colmap, 
            outline = FALSE, border = FALSE, add = T, xaxs = "i", 
            xaxt = "n", ...)
        if (sort) 
            lines(1:length(cols.b), colMeans, lty = 1, lwd = 1, 
                col = "red")
    }
    else {
        cmap <- col.func(cluster$nrows)
        matlines(1:length(cols.b), t(rats[cluster$rows, , drop = F]), 
            ylim = range.r, xlab = NA, ylab = NA, main = main, 
            col = cmap, lty = 1, ...)
        if (exists("col.rug")) 
            for (i in unique(col.rug)) rug(which(cols.b %in% 
                names(which(col.rug == i))), col = colmap[which(col.rug == 
                i)[1]])
    }
    if (cond.labels) {
        tmp.y <- rep(range.r[1] * 0.85, cluster$ncols)
        cols <- if (box.plot) 
            colmap
        else "black"
        text(1:cluster$ncols, tmp.y, cols.b, srt = 90, col = cols, 
            ...)
    }
    if (names(dev.cur()) == "devSVG") {
        par(family = "Arial")
        for (c in 1:length(cols.b)) {
            setSVGShapeToolTip(cols.b[c])
            rect(c, range.r[1], c + 1, range.r[2], col = NA, 
                border = NA)
        }
    }
    if (!is.null(o.genes)) {
        matlines(1:length(cols.b), t(rats[o.genes, , drop = F]), 
            lty = 1, lwd = 3, col = 2:6)
        legend("bottomright", legend = o.genes, lty = 1, lwd = 3, 
            col = 2:6, cex = 0.7, bty = "n")
    }
}
plotCluster.all.conds <-
function (cluster, imag = F, cond.labels = F, o.genes = NULL, 
    rats.names = names(ratios), range.r = NULL, sort = F, box.plot = F, 
    col.func = if (imag) topo.colors else rainbow, only.in.conds = F, 
    ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    main <- paste(sprintf("Cluster: %04d %s; resid: %s; r/c: %d/%d", 
        k, organism, sprintf("%.2f", weighted.mean(cluster$resid[rats.names], 
            row.weights[rats.names], na.rm = T)), length(cluster$rows), 
        length(cluster$cols)))
    rats <- get.cluster.matrix(unique(c(cluster$rows, o.genes)), 
        NULL, matrices = rats.names)
    cols.b <- c(colnames(rats)[colnames(rats) %in% cluster$cols], 
        colnames(rats)[!colnames(rats) %in% cluster$cols])
    if (only.in.conds) 
        cols.b <- colnames(rats)[colnames(rats) %in% cluster$cols]
    if (sort) {
        inClust <- colnames(rats)[colnames(rats) %in% cluster$cols]
        o1 <- order(apply(rats[cluster$rows, inClust, drop = F], 
            2, mean, na.rm = T))
        outClust <- colnames(rats)[!colnames(rats) %in% cluster$cols]
        o2 <- order(apply(rats[cluster$rows, outClust, drop = F], 
            2, mean, na.rm = T))
        cols.b <- c(inClust[o1], outClust[o2])
        if (only.in.conds) 
            cols.b <- inClust[o1]
    }
    len.b <- length(cols.b)
    rats <- rats[, cols.b, drop = F]
    par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    if (all(is.na(rats))) {
        plot(0, 0, typ = "n", main = main, ...)
        return()
    }
    if (is.vector(rats)) {
        rats <- t(rats)
        rownames(rats) <- cluster$rows
    }
    if (imag) {
        grey.image(t(rats), col = col.func(256))
        lines(rep(cluster$ncols + 0.5, 2), c(-999, 9999), col = 2, 
            lwd = 3, lty = 2)
        return()
    }
    if (is.null(range.r)) 
        range.r <- range(rats[rats != min(rats, na.rm = T) & 
            rats != max(rats, na.rm = T)], na.rm = T)
    if (cond.labels && len.b < 100) 
        range.r[1] <- range.r[1] * 1.5
    plot(1:len.b, xlim = c(0.95, len.b + 0.05), ylim = range.r, 
        xlab = NA, ylab = NA, main = main, typ = "n", xaxs = "i", 
        ...)
    if (length(ratios) > 1) {
        ind <- 0.5
        rts.in <- rts.out <- NULL
        for (in.out in 1:2) {
            if (in.out == 1) 
                cols <- cols.b[cols.b %in% cluster$cols]
            else if (in.out == 2) 
                cols <- cols.b[!cols.b %in% cluster$cols]
            for (i in 1:length(ratios)) {
                col <- col.func(length(ratios), s = 0.2)[i]
                rect(ind, range.r[1] + 0.01, ind + sum(cols %in% 
                  colnames(ratios[[i]])), range.r[2] - 0.01, 
                  col = col, dens = NA)
                ind <- ind + sum(cols %in% colnames(ratios[[i]]))
                if (in.out == 1) 
                  rts.in <- cbind(rts.in, rats[, cols[cols %in% 
                    colnames(ratios[[i]])], drop = F])
                else if (in.out == 2) 
                  rts.out <- cbind(rts.out, rats[, cols[cols %in% 
                    colnames(ratios[[i]])], drop = F])
            }
        }
        rats <- cbind(rts.in, rts.out)
        cols.b <- colnames(rats)
        len.b <- length(cols.b)
        rm(rts.in, rts.out)
    }
    if (exists("col.rug")) {
        if (is.integer(col.rug)) 
            colmap <- col.func(max(col.rug))[col.rug[cols.b]]
        else colmap <- col.rug[cols.b]
    }
    else if (length(rats.names) > 1) {
        colmap <- sapply(cols.b, function(col) which(sapply(ratios[rats.names], 
            function(i) col %in% colnames(i)))[1])
        if (is.list(colmap)) {
            colmap <- unlist(colmap)
            names(colmap) <- cols.b
        }
        colmap <- col.func(max(colmap))[colmap]
    }
    else if (all(deparse(col.func) == deparse(rainbow))) {
        colmap <- col.func(length(cols.b))
    }
    else {
        colmap <- col.func(cols.b)
    }
    if (box.plot) {
        colMeans <- apply(rats[cluster$rows, , drop = F], 2, 
            mean, na.rm = T)
        colSd <- apply(rats[cluster$rows, , drop = F], 2, sd, 
            na.rm = T)
        matlines(1:length(cols.b), cbind(colMeans - 2 * colSd, 
            colMeans + 2 * colSd), lty = 1, col = "lightgrey")
        boxplot(as.data.frame(rats[cluster$rows, , drop = F]), 
            ylim = range.r, names = NA, main = main, col = colmap, 
            outline = FALSE, border = FALSE, add = T, xaxs = "i", 
            xaxt = "n", ...)
        if (sort) 
            lines(1:length(cols.b), colMeans, lty = 1, lwd = 1, 
                col = "red")
    }
    else {
        cmap <- col.func(cluster$nrows)
        matlines(1:len.b, t(rats[cluster$rows, , drop = F]), 
            ylim = range.r, col = cmap, main = main, xlab = NA, 
            ylab = NA, lty = 1, ...)
        if (exists("col.rug")) 
            for (i in unique(col.rug)) rug(which(cols.b %in% 
                names(which(col.rug == i))), col = colmap[which(col.rug == 
                i)[1]])
    }
    cols.in <- colnames(rats)[colnames(rats) %in% cluster$cols]
    if (!only.in.conds) 
        lines(rep(length(cols.in) + 0.5, 2), range.r, col = 2, 
            lwd = 3, lty = 2)
    if (!is.null(o.genes)) {
        matlines(1:len.b, t(rats[o.genes, , drop = F]), lty = 1, 
            lwd = 3, col = 2:6)
        legend("bottomright", legend = o.genes, lty = 1, lwd = 3, 
            col = 2:6, cex = 0.7, bty = "n")
    }
    if (cond.labels) {
        tmp.y <- rep(range.r[1] * 0.85, len.b)
        cols <- if (box.plot) 
            colmap
        else "black"
        text(1:len.b, tmp.y, cols.b, srt = 90, col = cols, ...)
    }
    if (names(dev.cur()) == "devSVG") {
        par(family = "Arial")
        for (c in 1:length(cols.b)) {
            setSVGShapeToolTip(cols.b[c])
            rect(c, range.r[1], c + 1, range.r[2], col = NA, 
                border = NA)
        }
    }
}
plotCluster.motif <-
function (cluster, seqs = cluster$seqs, layout = NULL, colors = NULL, 
    motif.e.cutoff = Inf, no.plotCluster = T, addl.text = NULL, 
    ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    if (names(dev.cur()) == "devSVG") 
        par(family = "Arial")
    if (any(!cluster$rows %in% attr(ratios, "rnames"))) {
        cluster$rows <- cluster$rows[cluster$rows %in% attr(ratios, 
            "rnames")]
        cluster$nrows <- length(cluster$rows)
        warning(cluster$k, ": Some cluster rows are not in the ratios. Will plot without these rows.\n")
    }
    if (any(!cluster$cols %in% attr(ratios, "cnames"))) {
        cluster$cols <- cluster$cols[cluster$cols %in% attr(ratios, 
            "cnames")]
        cluster$ncols <- length(cluster$cols)
        warning(cluster$k, ": Some cluster cols are not in the ratios. Will plot without these cols.\n")
    }
    seq.types <- cluster$seq.type
    if (length(seq.types) == 1) {
        n.pssm.plot <- 3
    }
    else {
        n.pssm.plot <- 6
    }
    if (is.null(layout)) {
        if (length(seq.types) == 1) {
            layout <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 8, 2, 
                2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
                8, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 
                4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 
                6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7), ncol = 17, 
                byrow = T)
        }
        else {
            layout <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 
                2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
                11, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 
                1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 
                1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 
                3, 3, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 
                10, 10, 5, 5, 5, 5, 6, 6, 6, 6, 10, 10, 10, 10, 
                10, 10, 10, 10, 10, 7, 7, 7, 7, 9, 9, 9, 9, 10, 
                10, 10, 10, 10, 10, 10, 10, 10, 8, 8, 8, 8, 9, 
                9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10), 
                ncol = 17, byrow = T)
        }
        if (no.plotCluster) {
            layout[layout == 1 | layout == max(layout)] <- 2
            layout <- layout - 1
            layout[, 1][layout[, 1] == 1] <- max(layout) + 1
        }
    }
    layout(layout)
    k <- cluster$k
    if (!is.null(ratios)) {
        args <- list(...)
        args <- args[names(args) != "p.val.shade.cutoff"]
        args$cluster <- cluster
        do.call(plotCluster.all.conds, args)
        if (!no.plotCluster) 
            do.call(plotCluster, args)
    }
    rows <- cluster$rows
    if (is.null(colors) || !is.null(cluster$cog.code)) {
        tmp.lett <- 1:26
        names(tmp.lett) <- LETTERS
        if (!is.null(cluster$cog.code)) 
            coo <- cluster$cog.code[rows]
        else coo <- 1:length(rows)
        tmp <- unique(tmp.lett[coo])
        names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
        cols <- rainbow(length(tmp))
        names(cols) <- names(tmp)
        cols <- cols[names(tmp.lett[coo])]
        cols[is.na(names(cols))] <- "darkgrey"
        names(cols) <- rows
    }
    else {
        cols <- rainbow(length(rows))
        names(cols) <- rows
    }
    colors <- cols
    n.plotted <- 1
    for (seq.type in seq.types) {
        if (n.plotted > n.pssm.plot) 
            break
        if (is.null(cluster[[seq.type]]$e.val) || all(is.na(cluster[[seq.type]]$e.val)) || 
            is.null(cluster[[seq.type]]$motif.out) || is.null(cluster[[seq.type]]$motif.out$pssms)) 
            next
        pssm <- cluster[[seq.type]]$motif.out$pssms
        for (ppp in 1:min(floor(n.pssm.plot/length(seq.types)), 
            length(pssm))) {
            if (n.plotted > n.pssm.plot) 
                break
            if (cluster[[seq.type]]$motif.out$e.values[ppp] > 
                motif.e.cutoff) 
                next
            viewPssm(pssm[[ppp]], mot.ind = ppp, main.title = sprintf("%s PSSM #%d; e=%.3g", 
                seq.type, ppp, cluster[[seq.type]]$motif.out$e.values[ppp]), 
                cex.main = 0.9)
            n.plotted <- n.plotted + 1
        }
    }
    while (n.plotted <= n.pssm.plot) {
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
        n.plotted <- n.plotted + 1
    }
    suppressWarnings(cluster <- plotCluster.network(cluster, 
        ...))
    if (is.null(seqs)) {
        seqs <- rep("", length(cluster$rows))
        names(seqs) <- cluster$rows
    }
    if (!is.null(seq.type) && !is.null(seqs) && length(seqs) > 
        0) 
        plotClusterMotifPositions(cluster, seqs, colors = colors, 
            ...)
    else plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    try({
        par(mar = rep(0.5, 4), mgp = c(2, 1, 0) * 0.5)
        plot(c(0.5, 2.5), c(-1, 1), type = "n", tck = 0.01, cex.lab = 0.2, 
            cex.sub = 0.2, cex.axis = 0.2, axes = F)
        if (names(dev.cur()) == "devSVG") {
            par(family = "Arial")
            setSVGShapeToolTip(title = paste("Cluster:", sprintf("%04d", 
                k), organism, cmonkey.version), desc1 = sprintf("resid = %s; genes = %d; conds = %d", 
                paste(sprintf("%.2f", cluster$resid), collapse = " "), 
                length(cluster$rows), length(cluster$cols)))
            setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                paste(organism, ":", cluster$rows, sep = "", 
                  collapse = "+"), sep = ""))
            rect(0.5, -1, 3.25, +1, col = "lightgreen", border = NA)
        }
    })
    if (!is.null(cluster$name)) 
        text(0.65, 0, cluster$name, srt = 90, xpd = NA, cex = 1)
    else if (!is.null(addl.text)) 
        text(0.65, 0, addl.text, srt = 90, xpd = NA, cex = 1)
    text(1.5, 0, sprintf("%s iter=%d", date.run, iter), srt = 90, 
        xpd = NA, cex = 1)
    text(2.35, 0, paste("cMonkey Version", cmonkey.version, organism), 
        srt = 90, xpd = NA, cex = 1)
    invisible(cluster)
}
plotClusterMotifPositions <-
function (cluster, seqs = cluster$seqs, long.names = T, shade = T, 
    p.val.shade.cutoff = 999, colors = NULL, sort.by = "p.value", 
    o.genes = NULL, no.key = F, short.names = organism == "sce", 
    seq.type = cluster$seq.type[1], ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    rows <- cluster$rows
    if (!is.null(o.genes)) 
        rows <- unique(c(rows, o.genes))
    motif.out <- NULL
    if (!is.null(seq.type) && seq.type %in% names(cluster)) 
        motif.out <- cluster[[seq.type]]$motif.out
    is.dup.seq <- get.dup.seqs(cluster$seqs)
    p.clust <- cluster$p.clust
    e.clust <- cluster$e.val
    motif.info <- NULL
    if ((!all(is.na(p.clust)) || !all(is.na(e.clust))) && !is.null(motif.out) && 
        !is.null(motif.out$pv.ev) && length(motif.out$pv.ev) > 
        1) 
        motif.info <- subset(motif.out$pv.ev[[2]], gene %in% 
            rows)
    if (is.null(colors) || "cog.code" %in% names(genome.info)) {
        tmp.lett <- 1:26
        names(tmp.lett) <- LETTERS
        if (!is.null(cluster$cog.code)) 
            coo <- cluster$cog.code[rows]
        else coo <- 1:length(rows)
        tmp <- unique(tmp.lett[coo])
        names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
        cols <- rainbow(length(tmp))
        names(cols) <- names(tmp)
        cols <- cols[names(tmp.lett[coo])]
        cols[is.na(names(cols))] <- "darkgrey"
        names(cols) <- rows
    }
    else {
        cols <- rainbow(length(rows))
        names(cols) <- rows
    }
    cluster$colors <- colors <- cols
    no.motif <- FALSE
    p.values <- motif.widths <- pssm <- NULL
    if (!is.null(motif.out) && !is.null(motif.info) && nrow(motif.info) > 
        0 && (!all(is.na(p.clust)) || !all(is.na(e.clust)))) {
        p.values <- motif.out$p.values[rows]
        motif.widths <- sapply(motif.out$pssms, nrow, simplify = T)
        pssm <- motif.out$pssms
    }
    else {
        no.motif <- TRUE
        p.values <- numeric(length(rows))
        motif.widths <- 0
    }
    seqs <- seqs[rows]
    names(seqs) <- names(p.values) <- rows
    seq.lengths <- nchar(seqs)
    seq.lengths[seq.lengths == 2] <- NA
    if (any(seq.lengths[!is.na(seq.lengths)] > median(seq.lengths, 
        na.rm = T))) {
        seqs <- substr(seqs, 1, median(seq.lengths, na.rm = T))
        seq.lengths <- nchar(seqs)
    }
    maxlen <- max(seq.lengths, na.rm = T)
    if (!is.null(seq.type) && (maxlen == 0 || is.infinite(maxlen))) 
        maxlen <- diff(motif.upstream.search[[seq.type]])
    inds <- integer()
    if (no.motif && (sort.by == "p.value" || sort.by == TRUE)) 
        sort.by <- "gene.name"
    if (sort.by == "gene.name") 
        inds <- sort(rows, decreasing = T, index = T)$ix
    else if (sort.by == "p.value" || sort.by == TRUE) 
        inds <- order(p.values[rows], decreasing = T, na.last = F)
    else if (sort.by == "resid") 
        inds <- order(row.scores[rows, k], decreasing = T)
    if (length(inds) < length(rows)) 
        inds <- c((1:length(rows))[!1:length(rows) %in% inds], 
            inds)
    x.range <- c(-maxlen * 0.08, maxlen * 1.15)
    y.range <- c(0.5, length(rows) + 1)
    plot(x.range, y.range, type = "n", axes = F, xlab = "sequence position", 
        ylab = "")
    cexes <- 1
    if (!no.key) 
        axis(side = 1, pos = 0.6, tck = 0.01, mgp = c(0.1, 0.1, 
            0.1), labels = c(-1, seq(-100, -maxlen, -100)), at = seq(maxlen, 
            0, -100) + motif.upstream.scan[[seq.type]][1], ...)
    if (max(seq.lengths, na.rm = T) > 0) 
        sapply(maxlen - c(0, motif.upstream.search[[seq.type]]) + 
            motif.upstream.scan[[seq.type]][1], function(i) lines(rep(i, 
            2), c(-999, 999), col = "lightgray", lty = 2))
    colmap <- rainbow(length(rows))
    mots.used <- numeric()
    if (is.list(motif.widths)) {
        if (length(motif.widths) <= 0) 
            motif.widths <- 0
        else {
            for (i in 1:length(motif.widths)) if (is.null(motif.widths[[i]])) 
                motif.widths[[i]] <- 0
            motif.widths <- unlist(motif.widths)
        }
    }
    lwd <- 3
    if (length(rows) > 20) 
        lwd <- 1
    else if (length(rows) > 10) 
        lwd <- 2
    if (no.key) 
        lwd <- 1
    if (!no.motif) {
        tmp.mot.info <- subset(motif.info, gene %in% rows)
        tmp.mot.info <- subset(tmp.mot.info, posns <= diff(motif.upstream.scan[[seq.type]]))
        p.min <- quantile(log10(tmp.mot.info$pvals), 0.1, na.rm = T)
        if (is.na(p.min)) 
            p.min <- -5
        p.max <- quantile(log10(tmp.mot.info$pvals), na.rm = T, 
            0.9)
        if (is.na(p.max)) 
            p.max <- log10(p.val.shade.cutoff)
    }
    for (j in 1:length(rows)) {
        jj <- inds[j]
        cur.gene <- rows[jj]
        seq.len <- seq.lengths[jj]
        if (!is.null(rows)) {
            label <- rows[jj]
            if (!is.null(colors)) 
                rect(maxlen + 5, j - 0.18, maxlen * 1.195, j + 
                  0.18, col = colors[label], border = colors[label], 
                  lwd = 3)
        }
        if (!no.motif) {
            rects <- NULL
            mot.info <- subset(tmp.mot.info, gene == cur.gene)
            if (nrow(mot.info) > 0) {
                mots <- mot.info$mots
                starts <- mot.info$posns
                widths <- motif.widths[abs(mots)]
                for (i in 1:length(mots)) {
                  mot <- mots[i]
                  if (is.na(mot) || is.na(seq.len)) 
                    next
                  start <- starts[i]
                  if (start > seq.len) 
                    next
                  end <- start + widths[i]
                  mots.used <- unique(c(mots.used, abs(mot)))
                  col <- abs(mot) + 1
                  if (shade) {
                    if (!is.null(mot.info)) 
                      p.val <- mot.info$pvals[i]
                    else p.val <- 1e-05
                    if (is.na(p.val) || p.val > p.val.shade.cutoff || 
                      p.val > 1) 
                      next
                    else if (p.val <= 0) 
                      p.val <- 1e-05
                    p.val <- log10(p.val)
                    col <- col2rgb(palette()[col])/255
                    col[col > 0] <- 1
                    tmp <- if (p.val < 10) 
                      min(1, max(0, (p.val - p.min)/(p.max - 
                        p.min)))
                    else 0.99
                    if (names(dev.cur()) != "X11") {
                      alpha <- tmp
                    }
                    else {
                      col[col == 0] <- tmp
                      alpha <- 0
                    }
                    col[col < 0] <- 0
                    col[col > 1] <- 1
                    col <- rgb(col["red", 1], col["green", 1], 
                      col["blue", 1], 1 - alpha)
                  }
                  start.1 <- start + maxlen - seq.len
                  end.1 <- end + maxlen - seq.len
                  if (names(dev.cur()) == "devSVG") {
                    par(family = "Arial")
                    setSVGShapeToolTip(title = sprintf("Motif # %2d", 
                      abs(mot)), desc1 = paste(ifelse(mot < 0, 
                      "Rev.", "For."), "strand,", start - maxlen, 
                      "to", end - maxlen), desc2 = sprintf("p-value = %.2g", 
                      10^p.val))
                  }
                  if (!is.null(mot.info)) {
                    if (names(dev.cur()) == "devSVG") {
                      if (mot > 0) 
                        rect(start.1, j + 0.01, end.1, j + 0.3, 
                          col = col, border = col)
                      else if (mot < 0) 
                        rect(start.1, j - 0.3, end.1, j - 0.01, 
                          col = col, border = col)
                    }
                    else {
                      if (mot > 0) 
                        rects <- rbind(rects, c(start.1, j + 
                          0.01, end.1, j + 0.3, col = col, border = col))
                      else if (mot < 0) 
                        rects <- rbind(rects, c(start.1, j - 
                          0.3, end.1, j - 0.01, col = col, border = col))
                    }
                  }
                  else {
                    if (names(dev.cur()) == "devSVG") {
                      if (mot > 0) 
                        rect(start.1, j + 0.01, end.1, j + 0.3, 
                          border = col)
                      else if (mot < 0) 
                        rect(start.1, j - 0.3, end.1, j - 0.01, 
                          border = col)
                    }
                    else {
                      if (mot > 0) 
                        rects <- rbind(rects, c(start.1, j + 
                          0.01, end.1, j + 0.3, col = NA, border = col))
                      else if (mot < 0) 
                        rects <- rbind(rects, c(start.1, j - 
                          0.3, end.1, j - 0.01, col = NA, border = col))
                    }
                  }
                }
            }
            if (!is.null(rects)) 
                rect(rects[, 1], rects[, 2], rects[, 3], rects[, 
                  4], col = rects[, 5], border = rects[, 6])
        }
        slen <- seq.lengths[jj]
        if (all(seq.lengths[!is.na(seq.lengths)]) == 0) 
            slen <- 50
        lines(c(maxlen - slen, maxlen), c(j, j), lwd = lwd + 
            as.integer(rows[jj] %in% o.genes), col = colmap[jj])
        if (grepl("N", seqs[jj])) {
            locs <- which(strsplit(seqs[jj], "")[[1]] == "N")
            diff.locs <- c(diff(locs), 999)
            for (i in 1:sum(diff.locs > 1)) {
                if (i == 1) 
                  lines(c(locs[1], locs[diff.locs > 1][i]) + 
                    maxlen - slen, c(j, j), lwd = lwd + as.integer(rows[jj] %in% 
                    o.genes), col = "gray", lty = 2)
                else lines(c(locs[which(diff.locs > 1) + 1][i - 
                  1], locs[diff.locs > 1][i]) + maxlen - slen, 
                  c(j, j), lwd = lwd + as.integer(rows[jj] %in% 
                    o.genes), col = "gray", lty = 2)
            }
        }
        if (!is.null(rows)) {
            label <- rows[jj]
            col <- "black"
            if (exists("all.tfs") && label %in% all.tfs) 
                col <- "tomato3"
            if (names(dev.cur()) == "devSVG" || !long.names && 
                !no.key) {
                label <- substr(label, 1, 80)
                text(maxlen * 1.2, j, labels = label, adj = c(1, 
                  0.5), col = col, xpd = NA, ...)
            }
            if (long.names || names(dev.cur()) == "devSVG") {
                g.name <- toupper(label)
                if (!is.null(cluster$gene.coords)) 
                  g.name <- cluster$gene.coords[label]
                g.name[is.na(g.name)] <- label[is.na(g.name)]
                if (is.na(g.name)) 
                  g.name <- label
                if (names(dev.cur()) == "devSVG") {
                  par(family = "Arial")
                  setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                    organism, ":", label, sep = ""))
                  if (!is.na(g.name)) 
                    setSVGShapeToolTip(label, g.name)
                  else setSVGShapeToolTip(label)
                  rect(maxlen * 1.2, j - 0.18, maxlen, j + 0.18, 
                    col = NA, border = NA, xpd = NA)
                }
                else if (!is.na(g.name)) {
                  lab <- label
                  if (toupper(g.name) != toupper(label) && g.name != 
                    "") {
                    g.name <- gsub("^[:\\s+]+", "", gsub("\\s+$", 
                      "", g.name, perl = T), perl = T)
                    if (names(dev.cur()) == "X11") 
                      g.name <- strtrim(g.name, 40)
                    else gname <- strtrim(g.name, 60)
                    if (label != "") 
                      lab <- paste(g.name, ": ", label, sep = "")
                    else lab <- g.name
                  }
                  if (nchar(lab) > 60) 
                    lab <- substr(lab, nchar(lab) - 60, nchar(lab))
                  if (!no.key) 
                    text(maxlen * 1.2, j, labels = lab, adj = c(1, 
                      0.5), col = col, xpd = NA, ...)
                }
            }
            if ((!all(is.na(p.clust)) || !all(is.na(e.clust))) && 
                !no.key) 
                text(-maxlen * 0.07, j, labels = sprintf("%.2f", 
                  p.values[label]), xpd = NA, col = if (label %in% 
                  names(which(!is.dup.seq))) 
                  "black"
                else "blue", ...)
        }
    }
    if (!no.key && (!all(is.na(p.clust)) && !all(is.na(e.clust)))) {
        text(-maxlen * 0.15, length(rows) + 0.9, labels = sprintf("log10(P) %s", 
            seq.type), pos = 4, ...)
        mots.used <- sort(unique(mots.used))
        if (length(mots.used) > 1) {
            sapply(1:length(mots.used), function(j) text(maxlen * 
                0.24 + (j + 0) * maxlen * 0.03, length(rows) + 
                0.9, as.character(mots.used[j]), col = mots.used[j] + 
                1, xpd = NA, adj = c(0, 0.5), ...))
        }
        n.unique.seqs <- sum(!is.dup.seq)
        text(maxlen * 1.2, length(rows) + 0.9, sprintf("log10(P.clust)=%.2f; %d seqs; %d uniq", 
            p.clust[seq.type], length(seqs), n.unique.seqs), 
            xpd = NA, adj = c(1, 0.5), ...)
    }
}
plotCluster.network <-
function (cluster, network = "all", o.genes = NULL, colors = NULL, 
    cex = 0.7, no.legend = F, ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    require(igraph0)
    rows <- cluster$rows
    if (is.null(cluster$network)) {
        if (network == "all") 
            network <- names(networks)
        for (i in network) {
            if (!i %in% names(networks)) 
                next
            tmp.net <- networks[[i]][networks[[i]]$protein1 %in% 
                rows & networks[[i]]$protein2 %in% rows, ]
            tmp.net <- cbind(tmp.net, net = rep(i, nrow(tmp.net)))
            cluster$network <- if (!is.null(cluster$network)) 
                rbind(cluster$network, tmp.net)
            else tmp.net
        }
    }
    network <- cluster$network
    nrows <- rows
    if (!is.null(o.genes)) 
        nrows <- unique(c(nrows, o.genes))
    if (is.null(cluster$cog.code) && "cog.code" %in% names(genome.info)) 
        cluster$cog.code <- genome.info$cog.code[rows]
    if (is.null(cluster$colors)) {
        if (is.null(colors) || "cog.code" %in% names(genome.info)) {
            tmp.lett <- 1:26
            names(tmp.lett) <- LETTERS
            if (!is.null(cluster$cog.code)) 
                coo <- cluster$cog.code[rows]
            else coo <- 1:length(rows)
            tmp <- unique(tmp.lett[coo])
            names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
            cols <- rainbow(length(tmp))
            names(cols) <- names(tmp)
            cols <- cols[names(tmp.lett[coo])]
            cols[is.na(names(cols))] <- "darkgrey"
            names(cols) <- rows
            cluster$colors <- cols
        }
        else {
            cols <- rainbow(length(rows))
            names(cols) <- rows
            cluster$colors <- cols
        }
    }
    colors <- cluster$colors
    if (is.null(network) || nrow(network) <= 0) 
        network <- data.frame(protein1 = nrows, protein2 = nrows, 
            combined_score = jitter(rep(1/50, length(nrows))), 
            net = rep("none", length(nrows)))
    not.in <- nrows[!nrows %in% network$protein1 & !nrows %in% 
        network$protein2]
    for (i in not.in) network <- rbind(network, data.frame(protein1 = i, 
        protein2 = i, combined_score = 0, net = "none"))
    gr <- graph.edgelist(as.matrix(network[, 1:2]), directed = F)
    net.wts <- as.numeric(network$combined_score)
    names(net.wts) <- as.character(network$net)
    for (n in unique(names(net.wts))) {
        if (n == "none") 
            next
        net.wts[names(net.wts) == n] <- net.wts[names(net.wts) == 
            n]/max(net.wts[names(net.wts) == n], na.rm = T)
    }
    gr.layout <- layout.fruchterman.reingold(gr, niter = 3000, 
        weights = net.wts/5)
    gr.layout <- layout.norm(gr.layout, -1, 1, -1, 1)
    edge.colors <- character()
    curves <- rep(0, nrow(network))
    nets <- unique(as.character(network$net))
    if ("none" %in% nets) {
        nets <- unique(c("none", nets))
        inds <- c(1, 1:(length(nets) - 1))
        inds <- inds[inds != 0]
        inds <- inds[1:length(nets)]
        net.colors <- t(col2rgb(inds, T))/255
        net.colors[1, 4] <- 0
    }
    else {
        net.colors <- t(col2rgb(1:length(nets), T))/255
    }
    rownames(net.colors) <- nets
    for (i in 1:nrow(network)) {
        net <- as.character(network$net)[i]
        shade <- net.wts[i]
        nodes <- c(as.character(network$protein1)[i], as.character(network$protein2)[i])
        sub.net <- subset(network, protein1 %in% nodes & protein2 %in% 
            nodes)
        sub.nets <- unique(as.character(sub.net$net))
        curve.it <- max(0, nrow(sub.net) - 2)/2 * 0.33
        col <- net.colors[net, ]
        col2 <- col
        col2[col2 == 0] <- 1 - shade
        edge.colors[i] <- if (names(dev.cur()) != "X11") 
            rgb(col[1], col[2], col[3], shade)
        else rgb(col2[1], col2[2], col2[3])
        curves[i] <- curve.it * floor(which(sub.nets == net)/2) * 
            (if (which(sub.nets == net)%%2 == 0) 
                -1
            else 1)
    }
    if (all(curves == 0)) 
        curves <- FALSE
    if (!no.legend) {
        labels <- try(get.long.names(get.vertex.attribute(gr, 
            "name"), short = T), silent = T)
        if (class(labels) == "try-error") 
            labels <- get.vertex.attribute(gr, "name")
        labels[is.na(labels) | labels == ""] <- get.vertex.attribute(gr, 
            "name")[is.na(labels) | labels == ""]
    }
    else labels <- NA
    plot(gr, layout = gr.layout, margin = 0, rescale = F, edge.curved = curves, 
        vertex.color = colors[get.vertex.attribute(gr, "name")], 
        vertex.frame.color = colors[get.vertex.attribute(gr, 
            "name")], vertex.label.cex = cex, vertex.size = 7, 
        vertex.label = labels, vertex.label.family = if (names(dev.cur()) != 
            "pdf") 
            "Arial"
        else "sans", edge.color = edge.colors, edge.width = round(net.wts) + 
            1)
    if (!no.legend && length(nets[nets != "none"]) > 0) 
        legend("bottomright", legend = nets[nets != "none"], 
            col = 1:length(nets[nets != "none"]), lty = 1, lwd = 2, 
            bty = "n", cex = 0.5)
    if (names(dev.cur()) == "devSVG") {
        names <- cluster$gene.coords
        for (i in 1:nrow(gr.layout)) {
            gene <- get.vertex.attribute(gr, "name")[i]
            setSVGShapeToolTip(title = gene, desc1 = ifelse(is.na(names[gene]) || 
                is.null(names[gene]), "", names[gene]))
            setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                organism, ":", gene, sep = ""))
            points(gr.layout[i, 1], gr.layout[i, 2], col = "#FF000001", 
                cex = 10/3)
        }
    }
    cluster
}
plotScores <-
function (k, o.genes = NULL, b.genes = NULL, recompute = F) 
{
    opar <- par(no.readonly = T)
    rows <- get.rows(k)
    if (recompute || !exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores(k, force.row = T, force.col = T, 
                force.motif = T, force.net = T)
            rs <- tmp$r
            ms <- tmp$m
            ns <- tmp$n
            cs <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices(k)
            rs <- tmp$r[, 1]
            ms <- tmp$m[, 1]
            ns <- tmp$n[, 1]
            if (all(is.na(ms)) || all(ms == ms[1])) 
                rm(ms)
            if (all(is.na(ns)) || all(ns == ns[1])) 
                rm(ns)
        }
    }
    tmp.scale <- round(attr(ratios, "nrow")/length(rows)/4)
    layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3, byrow = T))
    if (!exists("rs")) 
        rs <- row.scores[, k, drop = T]
    rs[rs < -220] <- min(rs[rs > -220], na.rm = T)
    h <- try(hist(rs, breaks = 20, main = paste("Cluster", k), 
        xlab = "Ratios scores"), silent = T)
    if (class(h) != "try-error") {
        try(hist(rep(rs[rows], tmp.scale), breaks = h$breaks, 
            col = "red", border = "red", add = T), silent = T)
        try(hist(rs, breaks = h$breaks, add = T), silent = T)
    }
    if (exists("ms") || (!is.null(mot.scores) && !all(is.na(mot.scores[, 
        k])) && !no.genome.info)) {
        if (!exists("ms")) 
            ms <- mot.scores[, k, drop = T]
        ms[ms < -20] <- min(ms[ms > -20], na.rm = T)
        h <- try(hist(ms, breaks = 20, main = NULL, xlab = "Motif scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ms[rows], tmp.scale * 3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ms, breaks = h$breaks, add = T), silent = T)
        }
    }
    else {
        ms <- NULL
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
    }
    if (exists("ns") || (!is.null(net.scores) && !all(net.scores[, 
        k] == 0))) {
        if (!exists("ns")) {
            ns <- net.scores[, k, drop = T]
            ns[ns < -20] <- min(ns[ns > -20], na.rm = T)
            ns <- -log10(-ns)
        }
        ns[is.infinite(ns)] <- max(ns[!is.infinite(ns)], na.rm = T) + 
            0.1
        h <- try(hist(ns, breaks = 20, main = NULL, xlab = "-log10(-Network scores)"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ns[rows], tmp.scale/3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ns, breaks = h$breaks, add = T), silent = T)
        }
    }
    else {
        ns <= NULL
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
    }
    row.memb <- attr(ratios, "rnames") %in% rows
    if (!is.null(ms) && !all(is.na(ms)) && !all(ns == 0)) {
        plot(rs, ms, typ = "n", main = paste("Cluster", k), xlab = "Ratios scores", 
            ylab = "Mot scores")
        text(rs, ms, label = 1:length(rs), col = row.memb + 1, 
            cex = 0.5)
    }
    else if (!is.null(ns) && !all(ns == 0)) {
        plot(rs, ns, typ = "n", main = paste("Cluster", k), xlab = "Ratios scores", 
            ylab = "Net scores")
        text(rs, ns, label = 1:length(rs), col = row.memb + 1, 
            cex = 0.5)
    }
    else {
        plot(rs, jitter(rep(0, length(rs))), typ = "n", main = paste("Cluster", 
            k), xlab = "Ratios scores", ylab = "")
        text(rs, jitter(rep(0, length(rs))), label = 1:length(rs), 
            col = row.memb + 1, cex = 0.5)
    }
    if (!is.null(o.genes)) 
        text(rs[o.genes], ms[o.genes], label = which(attr(ratios, 
            "rnames") %in% o.genes), col = "green", cex = 0.5)
    if (!is.null(b.genes)) 
        text(rs[b.genes], ms[b.genes], label = which(attr(ratios, 
            "rnames") %in% b.genes), col = "blue", cex = 0.5)
    try({
        tmp <- get.combined.scores(quant = F)
        r.scores <- tmp$r
        c.scores <- tmp$c
        rr <- get.density.scores(ks = k, r.scores, c.scores, 
            plot = "rows")$r
        rr <- rr[, k, drop = T]
        h <- try(hist(log10(rr), breaks = 50, main = NULL, xlab = "Density (membership) scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(log10(rr[rows]), tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(log10(rr), breaks = h$breaks, add = T), 
                silent = T)
        }
    }, silent = T)
    par(opar)
}
plotStats <-
function (iter = stats$iter[nrow(stats)], plot.clust = NA, new.dev = F, 
    ...) 
{
    if (!exists("row.memb")) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        if (is.vector(row.memb)) 
            row.memb <- t(row.memb)
        rownames(row.memb) <- attr(ratios, "rnames")
        col.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "cnames") %in% get.cols(k))
        if (is.vector(col.memb)) 
            col.memb <- t(col.memb)
        rownames(col.memb) <- attr(ratios, "cnames")
    }
    if (!exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
            tmp <- get.combined.scores(quant = T)
            r.scores <- tmp$r
            c.scores <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
    }
    opar <- par(no.readonly = T)
    tmp.scale <- round(1/mean(row.memb, na.rm = T)/4)
    if (new.dev) {
        if (length(dev.list()) < 1) 
            dev.new()
        dev.set(2)
    }
    layout(matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 9, 
        11, 8, 10, 12), byrow = T, ncol = 3))
    par(mar = c(3, 3, 2, 0.1), mgp = c(3, 1, 0) * 0.5)
    stats <- stats[stats[, "iter"] <= iter, , drop = F]
    try(matplot(stats[, "iter"], stats[, grep("resid", colnames(stats), 
        val = T)], typ = "l", xlab = "iter", ylab = "Mean resid", 
        main = sprintf("Iter: %d", iter), lty = 1), silent = T)
    sapply(c(51, 101, 21), function(kwin) try(matlines(stats[, 
        "iter"], apply(stats[, grep("resid", colnames(stats), 
        val = T), drop = F], 2, function(i) runmed(i, k = min(length(i), 
        kwin))), lty = 2, lwd = 0.6), silent = T))
    if ((nn <- length(grep("resid", colnames(stats)))) > 1) 
        legend("bottomleft", legend = gsub("resid.", "", grep("resid", 
            colnames(stats), val = T)), lwd = 1, bty = "n", col = 1:nn, 
            lty = 1:nn, cex = 0.5)
    if (exists("row.scores") && !is.null(mot.scores) && !all(is.na(mot.scores[, 
        ]))) {
        rs <- row.scores[]
        rs[rs < -20] <- min(rs[rs > -20], na.rm = T)
        h <- try(hist(rs, breaks = 50, main = NULL, xlab = "Ratios scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(rs[row.memb == 1], tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(rs, breaks = h$breaks, add = T), silent = T)
        }
    }
    if (exists("mot.scores") && !is.null(mot.scores) && !all(is.na(mot.scores[, 
        ]))) {
        ms <- mot.scores[, ]
        ms[ms < -20] <- min(ms[ms > -20], na.rm = T)
        ms[ms >= 0] <- NA
        h <- try(hist(ms, breaks = 50, main = NULL, xlab = "Motif scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ms[row.memb == 1], tmp.scale * 3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ms, breaks = h$breaks, add = T), silent = T)
        }
    }
    tmp <- stats[, grep("p.clust", colnames(stats), val = T), 
        drop = F]
    if (!all(is.na(tmp))) {
        try(matplot(stats[, "iter"], tmp, typ = "l", xlab = "iter", 
            ylab = "Mean motif p-value", main = sprintf("Motif scaling: %.3f", 
                mot.scaling[max(1, iter - 1)]), lty = 1), silent = T)
        sapply(c(51, 101, 21), function(kwin) try(matlines(stats[!is.na(tmp), 
            "iter"], apply(tmp, 2, function(i) runmed(i[!is.na(i)], 
            k = min(sum(!is.na(i)), kwin))), lty = 2, lwd = 0.6), 
            silent = T))
        if ((nn <- length(grep("p.clust", colnames(stats)))) > 
            1) 
            legend("bottomleft", legend = gsub("p.clust.", "", 
                grep("p.clust", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5)
    }
    if (exists("net.scores") && !is.null(net.scores) && !all(net.scores[, 
        ] == 0)) {
        ns <- net.scores[, ]
        ns[ns < -20] <- min(ns[ns > -20], na.rm = T)
        ns[ns >= 0] <- NA
        ns[, ] <- -log10(-ns)
        tmp.scale <- ceiling(tmp.scale * mean(!is.na(ns), na.rm = T))
        h <- try(hist(ns, breaks = 50, main = NULL, xlab = "-log10(-Network scores)"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ns[row.memb == 1], tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ns, breaks = h$breaks, add = T), silent = T)
        }
        tmp <- stats[, grep("net.", colnames(stats), val = T, 
            fixed = T), drop = F]
        try(matplot(stats[, "iter"], tmp, typ = "l", xlab = "iter", 
            ylab = "Mean net-score", main = sprintf("Net scaling: %.3f", 
                net.scaling[max(1, iter - 1)]), lty = 1), silent = T)
        try(matlines(stats[, "iter"], apply(tmp, 2, function(i) runmed(i[!is.na(i)], 
            k = min(sum(!is.na(i)), 51))), lty = 2, lwd = 0.6), 
            silent = T)
        if ((nn <- length(grep("net.", colnames(stats)))) > 1) 
            try(legend("bottomleft", legend = gsub("net.", "", 
                grep("net.", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5), 
                silent = T)
    }
    clusterStack <- get.clusterStack(ks = 1:k.clust)
    resids <- sapply(as.list(clusterStack), "[[", "resid")
    try(hist(resids[resids <= 1.5], main = NULL, xlab = "Cluster Residuals", 
        xlim = if (all(resids > 0, na.rm = T)) 
            c(0, 1.5)
        else range(resids, na.rm = T), breaks = k.clust/4), silent = T)
    if (!exists("mot.scores") || is.null(mot.scores) || all(is.na(mot.scores[, 
        ]))) {
        pclusts <- sapply(as.list(clusterStack), "[[", "p.clust")
        try(hist(pclusts[pclusts <= 1], main = NULL, xlab = "Cluster Motif P-values", 
            xlim = if (all(pclusts > 0, na.rm = T)) 
                c(0, 1)
            else range(pclusts, na.rm = T), breaks = k.clust/4), 
            silent = T)
    }
    if (mot.scaling[iter] > 0) {
        plot.all.clusterMotifPositions <- function(ks = 1:k.clust, 
            mots = 1, e.cutoff = 1, p.cutoff = 0.05, seq.type = names(mot.weights)[1], 
            breaks = 100, ...) {
            if (seq.type == "ALL") 
                seq.type <- names(mot.weights)
            df <- NULL
            for (st in seq.type) {
                ms <- meme.scores[[st]][ks]
                ind <- 1
                if (!"posns" %in% colnames(ms[[1]]$pv.ev[[1]])) 
                  ind <- 2
                posns <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$posns)))
                pvals <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$pvals)))
                imots <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$mots)))
                genes <- as.vector(unlist(sapply(ms, function(i) as.character(i$pv.ev[[ind]]$gene))))
                slens <- nchar(genome.info$all.upstream.seqs[st][[st]][genes])
                clusts <- as.vector(unlist(sapply(ms, function(i) rep(i$k, 
                  if (is.null(i$pv.ev[[ind]])) 0 else nrow(i$pv.ev[[ind]])))))
                ms <- meme.scores[[st]]
                evals <- sapply(1:length(imots), function(i) ms[[clusts[i]]]$meme.out[[abs(imots[i])]]$e.value)
                df <- rbind(df, data.frame(clusts, posns, pvals, 
                  imots, evals, genes, slens, seq.type = rep(st, 
                    length(clusts))))
            }
            df2 <- subset(df, evals < e.cutoff & pvals < p.cutoff & 
                abs(imots) %in% mots)
            psns <- df2$posns - df2$slens - do.call(rbind, motif.upstream.scan[df2$seq.type])[, 
                1]
            h <- hist(psns, breaks = breaks, xlab = sprintf("%s %s", 
                paste(seq.type, collapse = " "), paste(mots, 
                  collapse = " ")), ...)
            dd <- density(psns, bw = 5)
            lines(dd$x, dd$y * max(h$counts)/max(dd$y) * 0.9, 
                col = "red")
            invisible(df)
        }
        try(plot.all.clusterMotifPositions(main = "Positions of motif #1", 
            ...), silent = T)
        if (any(sapply(meme.scores$upstream[1:k.clust], function(i) length(i$meme.out)) == 
            2)) 
            try(plot.all.clusterMotifPositions(mots = 2, main = "Positions of motif #2", 
                ...), silent = T)
    }
    n.rows <- sapply(1:k.clust, function(k) length(get.rows(k)))
    try(hist(n.rows, main = NULL, xlab = "Cluster Nrows", breaks = k.clust/4, 
        xlim = c(-5, max(n.rows, na.rm = T))), silent = T)
    n.cols <- sapply(1:k.clust, function(k) length(get.cols(k)))
    try(hist(n.cols, main = NULL, xlab = "Cluster Ncols", breaks = k.clust/4, 
        xlim = c(-5, attr(ratios, "ncol"))), silent = T)
    nr <- table(unlist(sapply(clusterStack, "[[", "rows")))
    if (length(nr) < attr(ratios, "nrow")) 
        nr <- c(nr, rep(0, attr(ratios, "nrow") - length(nr) + 
            1))
    try(hist(nr, breaks = seq(-0.5, max(nr, na.rm = T) + 0.5, 
        by = 1), main = NULL, xlab = "NClust per gene"), silent = T)
    if (!is.na(plot.clust)) {
        if (new.dev) {
            if (length(dev.list()) < 2) 
                dev.new()
            dev.set(3)
        }
        try(plotClust(plot.clust, w.motifs = T, cex = 0.7), silent = T)
        if (new.dev) {
            if (length(dev.list()) < 3) 
                dev.new()
            dev.set(4)
        }
        try(plotScores(plot.clust), silent = T)
    }
    par(opar)
}
preprocess.ratios <-
function (ratios, filter = T, normalize = T, col.groups = NULL, 
    frac.cutoff = 0.98) 
{
    if (is.null(col.groups)) 
        col.groups <- rep(1, ncol(ratios))
    if (is.null(names(col.groups))) 
        names(col.groups) <- colnames(ratios)
    if (filter) {
        cat("Filtering out nochange rows/cols from ratios matrix...\n")
        tmp1 <- apply(ratios, 1, function(i) mean(is.na(i) | 
            abs(i) <= 0.17)) < frac.cutoff
        tmp2 <- apply(ratios, 2, function(i) mean(is.na(i) | 
            abs(i) <= 0.1)) < frac.cutoff
        ratios <- ratios[tmp1, , drop = F]
        ratios <- ratios[, tmp2, drop = F]
        cat("Filtered ratios matrix is", paste(dim(ratios), collapse = "x"), 
            "\n")
        col.groups <- col.groups[tmp2]
    }
    if (normalize) {
        for (cg in unique(col.groups)) {
            cols <- names(which(col.groups == cg))
            cat("Converting ratios matrix", cg, "to z-scores...\n")
            ratios[, cols] <- t(scale(t(ratios[, cols, drop = F]), 
                center = apply(ratios[, cols, drop = F], 1, median, 
                  na.rm = T), scale = apply(ratios[, cols, drop = F], 
                  1, sd, na.rm = T)))
        }
    }
    ratios
}
pssm.to.string <-
function (pssm, cutoff.1 = 0.7, cutoff.2 = 0.4) 
{
    maxes <- max.col(pssm)
    letters <- col.let[maxes]
    values <- pssm[cbind(1:nrow(pssm), maxes)]
    letters[letters == "A" & values < cutoff.1] <- "a"
    letters[letters == "C" & values < cutoff.1] <- "c"
    letters[letters == "G" & values < cutoff.1] <- "g"
    letters[letters == "T" & values < cutoff.1] <- "t"
    letters[values < cutoff.2] <- "n"
    return(paste(letters, collapse = ""))
}
quantile.normalize.scores <-
function (scores, weights = NULL, keep.nas = F) 
{
    if (!is.list(scores) || sum(!sapply(scores, is.null)) <= 
        1) 
        return(scores)
    scores <- scores[!sapply(scores, is.null)]
    d <- dim(scores[[1]])
    dn <- dimnames(scores[[1]])
    tmp2 <- sapply(scores, function(i) sort(i[, ], na.last = T))
    if (is.null(weights)) 
        tmp2.mn <- rowMeans(tmp2, na.rm = T)
    else {
        for (i in 1:ncol(tmp2)) tmp2[, i] <- tmp2[, i] * weights[i]
        tmp2.mn <- rowMeans(tmp2, na.rm = T)/sum(weights, na.rm = T)
    }
    rm(tmp2)
    out <- list()
    for (n in names(scores)) {
        tmp <- rank(scores[[n]][, ], ties = "min", na = "keep")
        z <- matrix(tmp2.mn[tmp], nrow = d[1], ncol = d[2])
        dimnames(z) <- dn
        if (keep.nas) 
            z[is.na(scores[[n]][, ])] <- NA
        out[[n]] <- z
    }
    out
}
read.fasta <-
function (fname, lines = NULL) 
{
    if (is.null(lines)) 
        lines <- readLines(fname)
    lines <- lines[lines != ""]
    starts <- grep("^>", lines, perl = T)
    if (length(starts) > 1) 
        stops <- c(starts[2:length(starts)], length(lines) + 
            1)
    else stops <- length(lines) + 1
    seqs <- sapply(1:length(starts), function(i) paste(lines[(starts[i] + 
        1):(stops[i] - 1)], collapse = "", sep = ""))
    names(seqs) <- gsub("^>", "", lines[starts], perl = T)
    seqs
}
remove.low.complexity <-
function (seqs, length = 8, entropy.cutoff = 0.6, repl = "N", 
    use.dust = T, seq.type = names(mot.weights)[1]) 
{
    write.fasta <- function(seqs, fname) writeLines(paste(paste(">", 
        names(seqs), sep = ""), seqs, sep = "\n"), con = fname)
    if (use.dust) {
        seqs <- seqs[!is.null(seqs) & !is.na(seqs)]
        max.width <- as.integer(strsplit(meme.cmd[seq.type], 
            " ")[[1]][which(strsplit(meme.cmd[seq.type], " ")[[1]] == 
            "-maxw") + 1])
        seqs <- seqs[nchar(seqs) >= max.width]
        if (length(seqs) > 0) {
            fname <- my.tempfile("dust", suf = ".fst")
            write.fasta(seqs, fname)
            cmd <- gsub("$fname", fname, dust.cmd, fixed = T)
            fst <- system.time.limit(paste(cmd, "2>/dev/null"), 
                tlimit = 60)
            unlink(fname)
            if (length(fst) <= 1) 
                cat("WARNING: you probably don't have 'dust' installed.\n")
            else seqs <- read.fasta(NULL, fst)
            return(seqs)
        }
    }
}
rev.comp <-
function (seqs) 
{
    sapply(seqs, function(seq) paste(rev(strsplit(toupper(chartr("ATCG", 
        "tagc", seq)), "")[[1]]), collapse = ""))
}
row.col.membership.from.clusterStack <-
function (cs) 
{
    row.memb <- col.memb <- NULL
    for (k in 1:length(cs)) {
        row.memb <- cbind(row.memb, rep(0, attr(ratios, "nrow")))
        if (ncol(row.memb) == 1) 
            rownames(row.memb) <- attr(ratios, "rnames")
        rows <- cs[[k]]$rows
        rows <- rows[!is.na(rows)]
        row.memb[rows, k] <- k
        col.memb <- cbind(col.memb, rep(0, attr(ratios, "ncol")))
        if (ncol(col.memb) == 1) 
            rownames(col.memb) <- attr(ratios, "cnames")
        cols <- cs[[k]]$cols
        cols <- cols[!is.na(cols)]
        col.memb[cols, k] <- k
    }
    row.memb <- t(apply(row.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    row.memb <- row.memb[, apply(row.memb, 2, sum) != 0, drop = F]
    colnames(row.memb) <- NULL
    col.memb <- t(apply(col.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    col.memb <- col.memb[, apply(col.memb, 2, sum) != 0, drop = F]
    colnames(col.memb) <- NULL
    if (ncol(row.memb) < n.clust.per.row) 
        row.memb <- cbind(row.memb, rep(0, nrow(row.memb)))
    if (ncol(col.memb) < n.clust.per.col) 
        col.memb <- cbind(col.memb, rep(0, nrow(col.memb)))
    list(r = row.memb, c = col.memb)
}
runMast <-
function (memeOut, mast.cmd, genes, seqs, bgseqs = NULL, bg.list = NULL, 
    bgfname = NULL, unlink = T, verbose = F, ...) 
{
    fname <- my.tempfile("mast.tmp", suf = ".fst")
    if (is.null(bgfname)) {
        bgfname <- my.tempfile("mast.tmp", suf = ".bg")
        on.exit(unlink(bgfname))
    }
    memeOutFname <- my.tempfile("meme.tmp", suf = ".out")
    cat(memeOut, sep = "\n", file = memeOutFname)
    tmp <- mkTempMemeFiles(genes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname, seq.weights = NULL, 
        psps = NULL, ...)
    if (tmp <= 0) 
        return(NULL)
    cmd <- mast.cmd
    if (is.null(bgfname) || !file.exists(bgfname)) 
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
    else cmd <- gsub("$bgFname", bgfname, cmd, fixed = T)
    cmd <- gsub("$memeOutFname", memeOutFname, cmd, fixed = T)
    cmd <- gsub("$fname", fname, cmd, fixed = T)
    if (verbose) 
        cat(cmd, "\n")
    output <- system.time.limit(cmd)
    attr(output, "mast.command.line") <- cmd
    if (unlink) 
        unlink(c(memeOutFname, fname))
    output
}
runMeme <-
function (sgenes, seqs, cmd = meme.cmd[names(mot.weights)[1]], 
    bgseqs = NULL, bgfname = NULL, bg.list = NULL, nmotif = 1, 
    unlink = T, verbose = T, seq.weights = NULL, psps = NULL, 
    ...) 
{
    fname <- my.tempfile("meme.tmp", suf = ".fst")
    if (is.null(bgfname)) {
        bgfname <- my.tempfile("meme.tmp", suf = ".bg")
        on.exit(unlink(bgfname))
    }
    tmp <- mkTempMemeFiles(sgenes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname, seq.weights = seq.weights, 
        psps = psps, ...)
    if (tmp <= 0) 
        return(NULL)
    if (is.null(bgfname) || !file.exists(bgfname)) 
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
    else cmd <- gsub("$bgFname", bgfname, cmd, fixed = T)
    if (is.null(psps)) 
        cmd <- gsub("-psp $pspFname", "", cmd, fixed = T)
    else cmd <- gsub("$pspFname", sprintf("%s.psp", fname), cmd, 
        fixed = T)
    cmd <- gsub("$fname", fname, cmd, fixed = T)
    if (verbose) 
        cat(cmd, "\n")
    output <- system.time.limit(cmd)
    attr(output, "meme.command.line") <- cmd
    if (unlink) 
        unlink(c(fname, sprintf("%s.psp", fname)))
    return(output)
}
save.cmonkey.env <-
function (env = NULL, file = NULL, verbose = T) 
{
    if (is.null(env)) {
        for (i in ls(env = globalenv())) if (is.environment(get(i, 
            globalenv())) && "cmonkey" %in% class(get(i, globalenv()))) 
            save.cmonkey.env(get(i, globalenv()), file, verbose)
        return(invisible())
    }
    if (is.null(file)) 
        file <- paste(env$cmonkey.filename, ".RData", sep = "")
    if (verbose) 
        message("Saving environment to ", file)
    save(env, file = file)
    invisible(env)
}
seed.clusters <-
function (k.clust, seed.method = "rnd", col.method = "rnd") 
{
    if (seed.method == "custom" && exists("seed.clusters.custom")) 
        return(seed.clusters.custom(k.clust, col.method))
    if (substr(seed.method, 1, 3) == "net" && length(networks) <= 
        0) {
        cat("Seed method is", seed.method, ", but no networks -- using 'kmeans' instead.\n")
        seed.method <- "kmeans"
    }
    if (seed.method == "rnd") {
        rm <- t(sapply(1:attr(ratios, "nrow"), function(i) sample(1:k.clust, 
            n.clust.per.row[1], replace = n.clust.per.row[1] > 
                attr(ratios, "nrow"))))
    }
    else if (substr(seed.method, 1, 5) == "list=") {
        rm <- matrix(0, nrow = attr(ratios, "nrow"), ncol = n.clust.per.row[1])
        rownames(rm) <- attr(ratios, "rnames")
        fname <- strsplit(seed.method, "=")[[1]][2]
        if (exists(fname)) 
            lists <- get(fname)
        else if (file.exists(fname)) 
            lists <- strsplit(readLines(fname), split = "[ \t,]", 
                perl = T)
        for (k in 1:min(c(k.clust, length(lists)))) {
            probes <- unlist(lapply(get.synonyms(lists[[k]]), 
                function(i) i %in% rownames(rm)))
            rm[probes[rm[probes, 1] == 0], 1] <- k
            rm[probes[rm[probes, 1] != 0], 2] <- k
        }
        if (length(lists) < k.clust) {
            for (k in (length(lists) + 1):k.clust) {
                rnames <- attr(ratios, "rnames")[!attr(ratios, 
                  "rnames") %in% unlist(lists)]
                rows <- sample(rnames, 5)
                rm[rows[rm[rows, 1] == 0], 1] <- k
                rm[rows[rm[rows, 1] != 0 & rm[rows, 2] == 0], 
                  2] <- k
            }
        }
    }
    else if (substr(seed.method, 1, 4) == "rnd=") {
        n.samp <- as.integer(strsplit(seed.method, "=")[[1]][2])
        rm <- matrix(0, nrow = attr(ratios, "nrow"), ncol = n.clust.per.row[1])
        rownames(rm) <- attr(ratios, "rnames")
        for (i in 1:n.clust.per.row) {
            sampled <- rep(FALSE, attr(ratios, "nrow"))
            names(sampled) <- attr(ratios, "rnames")
            for (k in 1:k.clust) {
                g <- sample(attr(ratios, "rnames")[!sampled], 
                  n.samp)
                rm[g, 1] <- k
                sampled[g] <- TRUE
            }
        }
    }
    else if (seed.method == "kmeans") {
        if (!exists("ratios")) 
            stop("kmeans seed method but no ratios")
        tmp.rat <- get.cluster.matrix()
        tmp.rat[is.na(tmp.rat)] <- 0
        km <- kmeans
        tmp.rat[is.na(tmp.rat)] <- rnorm(sum(is.na(tmp.rat)))
        rm <- km(tmp.rat, centers = k.clust, iter.max = 20, nstart = 2)$cluster
        names(rm) <- attr(ratios, "rnames")
        if (n.clust.per.row[1] > 1) 
            rm <- cbind(rm, matrix(rep(0, attr(ratios, "nrow") * 
                (n.clust.per.row[1] - 1)), ncol = n.clust.per.row[1] - 
                1))
    }
    if (is.vector(rm)) 
        rm <- t(rm)
    if (nrow(rm) == 1) 
        rm <- t(rm)
    if (col.method == "rnd") {
        cm <- t(sapply(1:attr(ratios, "ncol"), function(i) sample(1:k.clust, 
            n.clust.per.col[1], replace = n.clust.per.col[1] > 
                k.clust)))
    }
    else if (col.method == "best") {
        if (!exists("ratios")) 
            stop("best col seed method but no ratios")
        all.rats <- get.cluster.matrix()
        attr(all.rats, "all.colVars") <- apply(all.rats, 2, var, 
            use = "pair", na.rm = T)
        col.scores <- -sapply(1:k.clust, function(k) if (sum(rm == 
            k, na.rm = T) <= 0) 
            rep(NA, attr(ratios, "ncol"))
        else get.col.scores(k = rownames(which(rm == k, arr = T)), 
            ratios = all.rats, method = "orig"))
        cm <- t(apply(col.scores, 1, function(i) order(i, decreasing = T)[1:n.clust.per.col[1]]))
    }
    rownames(rm) <- attr(ratios, "rnames")
    rownames(cm) <- attr(ratios, "cnames")
    list(rm = rm, cm = cm)
}
set.param <-
function (name, val, env = cmonkey.params, override = F, quiet = F) 
{
    if (!exists(name, envir = env) || override) {
        if (!quiet) 
            try({
                cat(name, "-> ")
                str(val, digits.d = 9, no.list = T)
            })
        assign(name, val, envir = env)
    }
    else {
        val <- get(name, envir = env)
        if (!quiet) 
            try({
                cat(name, "= ")
                str(val, digits.d = 9, no.list = T)
            })
        assign(name, val, envir = env)
    }
    assign(name, val, envir = parent.frame())
}
system.time.limit <-
function (cmd, tlimit = 600) 
{
    out <- readLines(pipe(cmd, open = "rt"))
    out
}
update.all.clusters <-
function (env, dont.update = F, ...) 
{
    tmp <- env$row.col.membership.from.clusterStack(env$clusterStack)
    row.membership <- tmp$r
    col.membership <- tmp$c
    tmp <- get.all.scores(...)
    env$row.scores <- tmp$r
    env$mot.scores <- tmp$m
    env$net.scores <- tmp$n
    env$col.scores <- tmp$c
    env$meme.scores <- tmp$ms
    if (!is.null(tmp$cns)) 
        env$cluster.net.scores <- tmp$cns
    tmp <- get.combined.scores(quant = F)
    r.scores <- tmp$r
    c.scores <- tmp$c
    if (length(tmp$scalings) > 0) {
        env$row.scaling[iter] <- tmp$scalings["row"]
        env$mot.scaling[iter] <- tmp$scalings["mot"]
        env$net.scaling[iter] <- tmp$scalings["net"]
    }
    rm(tmp)
    row.memb <- sapply(1:k.clust, function(k) attr(ratios, "rnames") %in% 
        get.rows(k))
    if (is.vector(row.memb)) 
        row.memb <- t(row.memb)
    rownames(row.memb) <- attr(ratios, "rnames")
    col.memb <- sapply(1:k.clust, function(k) attr(ratios, "cnames") %in% 
        get.cols(k))
    if (is.vector(col.memb)) 
        col.memb <- t(col.memb)
    rownames(col.memb) <- attr(ratios, "cnames")
    if (row.scaling[iter] > 0 && fuzzy.index[iter] > 1e-05) {
        r.scores[, ] <- r.scores[, ] + rnorm(length(r.scores[, 
            ]), sd = sd(r.scores[, ][row.memb[, ] == 1], na.rm = T) * 
            fuzzy.index[iter])
        if (!is.null(c.scores)) 
            c.scores[, ] <- c.scores[, ] + rnorm(length(c.scores[, 
                ]), sd = sd(c.scores[, ][col.memb[, ] == 1], 
                na.rm = T) * fuzzy.index[iter])
    }
    tmp <- get.density.scores(ks = 1:k.clust, r.scores, col.scores)
    rr.scores <- tmp$r
    cc.scores <- tmp$c
    rm(tmp)
    size.compensation.func.rows <- function(n) exp(-n/(attr(ratios, 
        "nrow") * n.clust.per.row/k.clust))
    size.compensation.func.cols <- function(n) exp(-n/(attr(ratios, 
        "ncol") * n.clust.per.col/k.clust))
    for (k in 1:k.clust) {
        tmp <- sum(row.memb[, k])
        rr.scores[, k] <- rr.scores[, k] * size.compensation.func.rows(max(tmp, 
            cluster.rows.allowed[1]))
        if (!is.null(cc.scores)) {
            tmp <- sum(col.memb[, k])
            cc.scores[, k] <- cc.scores[, k] * size.compensation.func.cols(max(tmp, 
                attr(ratios, "ncol")/10))
        }
    }
    tmp <- get.updated.memberships(row.membership, col.membership, 
        rr.scores, cc.scores)
    row.membership <- tmp$r
    col.membership <- tmp$c
    rm(tmp)
    if (env$post.adjust == TRUE && env$iter == env$n.iter) {
        env$adjust.all.clusters(env, expand.only = F)
        gc()
    }
    if (!dont.update) {
        env$clusterStack <- lapply(1:env$k.clust, function(k) list(rows = rownames(which(row.membership == 
            k, arr = T)), cols = rownames(which(col.membership == 
            k, arr = T))))
        env$clusterStack <- env$get.clusterStack(ks = 1:k.clust)
    }
    env
}
update.cmonkey.env <-
function (object, ...) 
{
    if (file.exists("cmonkey-funcs.R")) {
        tmp.e <- new.env()
        sys.source("cmonkey.R", envir = tmp.e)
    }
    else {
        tmp.e <- environment(cMonkey:::cmonkey)
    }
    for (i in ls(tmp.e)) {
        if (i %in% c("DATE", "VERSION")) 
            next
        f <- try(get(i, envir = tmp.e))
        f2 <- try(get(paste("super", i, sep = "."), envir = object), 
            silent = T)
        if (class(f) == "function") {
            environment(f) <- object
            if (class(f2) != "function") 
                assign(i, f)
            else assign(paste("super", i, sep = "."), f)
        }
    }
    rm(f, f2, tmp.e, i)
    for (i in ls()) {
        if (i %in% c("i", "object")) 
            next
        f <- get(i)
        if (is.function(f)) 
            assign(i, f, object)
    }
}
viewPssm <-
function (pssm, e.val = NA, mot.ind = NA, use.char = T, main.title = NA, 
    no.par = F, scale.e = NA, boxes = F, new = T, xoff = 0, yoff = 0, 
    no.axis.labels = F, min.height.drawn = 1e-05, no.render = F, 
    ...) 
{
    if (is.null(pssm)) 
        return()
    getEntropy <- function(pssm) {
        pssm[pssm == 0] <- 1e-05
        entropy <- apply(pssm, 1, function(i) -sum(i * log2(i)))
        return(entropy)
    }
    char.coords = list(T = list(x = c(0.45, 0.55, 0.55, 1, 1, 
        0, 0, 0.45), y = c(0, 0, 0.9, 0.9, 1, 1, 0.9, 0.9), color = 2), 
        A = list(x = c(0, 0.1, 0.28, 0.72, 0.68, 0.32, 0.5, 0.9, 
            1, 0.55, 0.45, 0), y = c(0, 0, 0.4, 0.4, 0.5, 0.5, 
            0.9, 0, 0, 1, 1, 0), color = 3), C = list(x = c(1, 
            1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 
            0.85, 1, 1, 0.9, 0.9, 0.8, 0.55, 0.45, 0.2, 0.1, 
            0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9), y = c(0.6, 
            0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 
            0.35, 0.4, 0.4, 0.35, 0.2, 0.1, 0.1, 0.2, 0.42, 0.58, 
            0.8, 0.9, 0.9, 0.8, 0.65, 0.6), color = 4), G = list(x = c(1, 
            1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 
            0.85, 1, 1, 0.7, 0.7, 0.9, 0.8, 0.55, 0.45, 0.2, 
            0.1, 0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9), y = c(0.6, 
            0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 
            0.35, 0.5, 0.5, 0.4, 0.4, 0.2, 0.1, 0.1, 0.2, 0.42, 
            0.58, 0.8, 0.9, 0.9, 0.8, 0.65, 0.6), color = "orange"))
    draw.char <- function(char = col.let, rect = c(0, 0, 1, 1), 
        border = NULL, ...) {
        if (rect[4] <= min.height.drawn) 
            return()
        x <- char.coords[[char]]$x * rect[3] + rect[1]
        y <- char.coords[[char]]$y * rect[4] + rect[2]
        color <- char.coords[[char]]$color
        if (is.null(border)) 
            border <- color
        if (!boxes) {
            polygon(x + xoff, y + yoff, col = color, border = border, 
                density = NA, ...)
        }
        else {
            rect[1] <- rect[1] + xoff
            rect[2] <- rect[2] + yoff
            rect(rect[1], rect[2], rect[1] + rect[3], rect[2] + 
                rect[4], col = color, border = border, density = NA, 
                ...)
        }
    }
    win.size <- nrow(pssm)
    if (!no.par) 
        par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 0.75)
    if (any(pssm <= 0)) 
        pssm <- pssm + 1e-10
    if (any(pssm > 1)) 
        pssm <- t(apply(pssm, 1, function(i) i/(sum(i) + 1e-10)))
    entr <- getEntropy(pssm)
    if (is.na(scale.e)) 
        scale.e <- (2 - entr)/2
    x.range <- c(0.5, win.size + 0.5)
    y.range <- c(0, max(scale.e))
    if (new) 
        plot(x.range, y.range, type = "n", tck = 0.01, cex.lab = 0.2, 
            cex.sub = 0.2, cex.axis = 0.2, axes = F)
    if (!is.na(main.title[1])) {
        if (!is.na(mot.ind)) 
            title(main.title, col.main = mot.ind + 1, xpd = NA, 
                ...)
        else title(main.title, xpd = NA, ...)
    }
    else if (!is.na(mot.ind) || !is.na(e.val)) {
        if (!is.na(e.val)) 
            tmp.tit <- sprintf("PSSM #%d; E=%.3g", mot.ind, e.val)
        else tmp.tit <- tmp.tit <- sprintf("PSSM #%d", mot.ind)
        title(tmp.tit, col.main = mot.ind + 1, xpd = NA, ...)
    }
    if (no.render) 
        return(invisible(y.range))
    pssm.sc <- scale.e * pssm
    for (j in 1:win.size) {
        inds <- sort(pssm.sc[j, ], index = T)$ix
        for (i in 1:4) {
            ind <- inds[i]
            if (i == 1) {
                if (!use.char) {
                  rect((j - 0.5), 0, (j + 0.5), pssm.sc[j, ind], 
                    col = colMap[ind])
                  if (pssm[j, ind] > 0.05) 
                    text(j, 0 + pssm.sc[j, ind]/2, colLet[ind])
                }
                else {
                  draw.char(col.let[ind], c((j - 0.4), 0, 0.9, 
                    pssm.sc[j, ind] - 0.001), ...)
                }
                prev.h <- pssm.sc[j, ind]
            }
            else {
                if (!use.char) {
                  rect((j - 0.5), prev.h, (j + 0.5), (pssm.sc[j, 
                    ind] + prev.h), col = colMap[ind])
                  if (pssm.sc[j, ind] > 0.05) {
                    if (i == 2) 
                      text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                        colLet[ind], col = 8)
                    else text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                      colLet[ind])
                  }
                }
                else {
                  draw.char(col.let[ind], c((j - 0.4), prev.h, 
                    0.9, pssm.sc[j, ind] - 0.001), ...)
                }
                prev.h <- prev.h + pssm.sc[j, ind]
            }
        }
    }
    if (!no.axis.labels) {
        if (win.size < 10) 
            text(1:win.size, rep(-0.01, win.size), as.character(1:win.size), 
                cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else if (win.size < 20) 
            text(seq(1, win.size, 2), rep(-0.01, win.size), as.character(seq(1, 
                win.size, 2)), cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else if (win.size < 50) 
            text(seq(1, win.size, 5), rep(-0.01, win.size), as.character(seq(1, 
                win.size, 5)), cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else text(seq(1, win.size, 25), rep(-0.01, win.size), 
            as.character(seq(1, win.size, 25)), cex = 0.7, adj = c(0.5, 
                1), xpd = NA)
    }
    invisible(y.range)
}
write.project <-
function (ks = sapply(as.list(clusterStack), "[[", "k"), para.cores = 1, 
    out.dir = NULL, gaggle = T, seq.type = names(mot.weights)[1], 
    gzip = T, output = c("svg", "pdf", "png", "html", "main", 
        "rdata"), ...) 
{
    if (is.null(out.dir)) {
        out.dir <- cmonkey.filename
        if (iter != n.iter) 
            out.dir <- sprintf("%s_%04d", out.dir, iter)
    }
    cat("Outputing to", out.dir, "\n")
    if (!file.exists(out.dir)) 
        dir.create(out.dir, recursive = T, showWarnings = F)
    clusterStack <- clusterStack[ks]
    mc <- get.parallel(length(ks), para.cores = para.cores)
    if (!file.exists(paste(out.dir, "/svgs", sep = ""))) 
        dir.create(paste(out.dir, "/svgs", sep = ""), showWarnings = F)
    if (!file.exists(paste(out.dir, "/pdfs", sep = ""))) 
        dir.create(paste(out.dir, "/pdfs", sep = ""), showWarnings = F)
    if (!file.exists(paste(out.dir, "/htmls", sep = ""))) 
        dir.create(paste(out.dir, "/htmls", sep = ""), showWarnings = F)
    if ("svg" %in% output) {
        require(RSVGTipsDevice)
        if (!file.exists(sprintf("%s/svgs/stats.svg", out.dir)) && 
            !file.exists(sprintf("%s/svgs/stats.svgz", out.dir))) {
            cat("STATS...\n")
            devSVGTips(sprintf("%s/svgs/stats.svg", out.dir), 
                toolTipMode = 2, title = "Biclustering statistics", 
                xmlHeader = T)
            par(family = "Arial")
            plotStats(new.dev = F)
            dev.off()
        }
        require(igraph0)
        cat("SVGS: ")
        for (qqq in 1:3) {
            for (k in ks) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (file.exists(sprintf("%s/svgs/cluster%04d.svg", 
                  out.dir, k)) || file.exists(sprintf("%s/svgs/cluster%04d.svgz", 
                  out.dir, k))) 
                  next
                tmp.cl <- try(plotClust(k, w.motifs = T, seq.type = seq.type, 
                  dont.plot = T, ...))
                devSVGTips(sprintf("%s/svgs/cluster%04d.svg", 
                  out.dir, k), toolTipMode = 2, title = sprintf("Bicluster %04d", 
                  k), xmlHeader = T)
                try(plotClust(cluster = tmp.cl, w.motifs = T, 
                  seq.type = seq.type, ...))
                dev.off()
            }
        }
        cat("\n")
    }
    if ("pdf" %in% output) {
        require(igraph0)
        cat("PDFS: ")
        lapply(ks, function(k) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/pdfs/cluster%04d.pdf", 
                out.dir, k))) 
                next
            pdf(sprintf("%s/pdfs/cluster%04d.pdf", out.dir, k))
            plotClust(k, w.motifs = T, seq.type = seq.type, ...)
            dev.off()
        })
        cat("\n")
    }
    if (gaggle && "html" %in% output) {
        require(hwriter)
        cat("HTMLS: ")
        lapply(ks, function(k, ...) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d.html", 
                out.dir, k))) 
                return()
            rows <- sort(get.rows(k))
            if (length(rows) <= 0) 
                return()
            short.names <- try(get.long.names(rows, short = T))
            if (class(short.names) == "try-error") 
                short.names <- rep("", length(rows))
            short.names <- cbind(rows, short.names)
            rownames(short.names) <- colnames(short.names) <- NULL
            long.names <- try(get.long.names(rows, short = F))
            if (class(long.names) == "try-error") 
                long.names <- rep("", length(rows))
            long.names <- cbind(rows, long.names)
            rownames(long.names) <- colnames(long.names) <- NULL
            refseq.names <- unique(unlist(get.synonyms(rows)))
            refseq.names <- grep("^NP_", refseq.names, val = T)
            upstream.seqs <- try(get.sequences(k, filter = F, 
                uniq = F), silent = T)
            if (class(upstream.seqs) == "try-error" || is.null(upstream.seqs) || 
                length(upstream.seqs) == 0) {
                upstream.seqs <- rep("", length(rows))
                names(upstream.seqs) <- rows
            }
            upstream.seqs <- cbind(names(upstream.seqs), upstream.seqs)
            rownames(upstream.seqs) <- colnames(upstream.seqs) <- NULL
            htmltext <- paste(c("<html><head><title>Bicluster %K (%FILE)</title>", 
                "<style type=\"text/css\">", "  .hidden {", "     display: none;", 
                "   }", "  .gaggle-data {", "     color: green;", 
                "     font-size: xx-small;", "   }", "   p {", 
                "     color: red;", "     font-size: x-small;", 
                "   }", "</style>", "<script type=\"text/javascript\">", 
                "   function toggleVisible(id){", "      if (document.getElementById){", 
                "         obj = document.getElementById(id);", 
                "         if (obj) {", "            if (obj.style.display == 'none'){", 
                "               obj.style.display = 'block';", 
                "            } else {", "               obj.style.display = 'none';", 
                "            }", "         }", "      }", "   }", 
                "</script>", "</head>", "<table><tr><td>", "<iframe src=\"../svgs/cluster%K03d%K.svg\" width=\"600\" height=\"520\" frameborder=\"0\"></iframe>", 
                "</td><td>", "<p><a href=\"#bicluster%K03d%K\" onclick=\"toggleVisible('bicluster%K03d%K'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows and columns.</p>", 
                "<div id=\"bicluster%K03d%K\" style=\"display:none;\" class=\"gaggle-data bicluster\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                  length(rows), length(get.cols(k))), "   <div class=\"gaggle-cluster\">", 
                "      <ol class=\"gaggle-rowNames\">", paste("<li>", 
                  sort(rows), "</li>", sep = "", collapse = ""), 
                "      </ol>", "   <ol class=\"gaggle-columnNames\">", 
                paste("<li>", sort(get.cols(k)), "</li>", sep = "", 
                  collapse = ""), "      </ol>", "   </div>", 
                "</div>", "<p><a href=\"#bicluster%K03d%K_genes\" onclick=\"toggleVisible('bicluster%K03d%K_genes'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (genes).</p>", 
                "<div id=\"bicluster%K03d%K_genes\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K genes</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(rows)), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(rows), "</li>", 
                  sep = "", collapse = ""), "      </ol>", "   </div>", 
                "</div>", "<p><a href=\"#bicluster%K03d%K_short_names\" onclick=\"toggleVisible('bicluster%K03d%K_short_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (short gene names).</p>", 
                "<div id=\"bicluster%K03d%K_short_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K short names</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(short.names)), "   <span class=\"gaggle-namelist-tag hidden\">short_name</span>", 
                hwrite(short.names, table.class = "toc", col.class = list(NA, 
                  "short_name"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color: green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_long_names\" onclick=\"toggleVisible('bicluster%K03d%K_long_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (long gene names).</p>", 
                "<div id=\"bicluster%K03d%K_long_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K long names</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(long.names)), "   <span class=\"gaggle-namelist-tag hidden\">long_name</span>", 
                hwrite(long.names, table.class = "toc", col.class = list(NA, 
                  "long_name"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_refseq_names\" onclick=\"toggleVisible('bicluster%K03d%K_refseq_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (NCBI RefSeq gene IDs).</p>", 
                "<div id=\"bicluster%K03d%K_refseq_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K NCBI RefSeq IDs</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(refseq.names)), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(refseq.names), 
                  "</li>", sep = "", collapse = ""), "      </ol>", 
                "   </div>", "</div>", "<p><a href=\"#bicluster%K03d%K_upstream_seqs\" onclick=\"toggleVisible('bicluster%K03d%K_upstream_seqs'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K gene upstream sequences.</p>", 
                "<div id=\"bicluster%K03d%K_upstream_seqs\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K upstream sequences</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(upstream.seqs)), "   <span class=\"gaggle-namelist-tag hidden\">upstream</span>", 
                hwrite(upstream.seqs, table.class = "toc", col.class = list(NA, 
                  "upstream"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_arrays\" onclick=\"toggleVisible('bicluster%K03d%K_arrays'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K columns (arrays; conditions).</p>", 
                "<div id=\"bicluster%K03d%K_arrays\" style=\"display:none;\" class=\"gaggle-data arrays\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K arrays</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(get.cols(k))), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(get.cols(k)), 
                  "</li>", sep = "", collapse = ""), "      </ol>", 
                "   </div>", "</div>", "<p><a href=\"#bicluster%K03d%K_ratios\" onclick=\"toggleVisible('bicluster%K03d%K_ratios'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K ratios.</p>", "<div id=\"bicluster%K03d%K_ratios\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K ratios</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                  length(rows), length(get.cols(k))), "   <div class=\"gaggle-matrix-tsv\">", 
                "        RATIOS", "   </div>", "</div>", if (!is.null(seq.type) && 
                  !is.null(meme.scores[[seq.type]][[k]]$meme.out) && 
                  !is.null(meme.scores[[seq.type]][[k]]$meme.out[[1]])) paste("<p><a href=\"#bicluster%K03d%K_pssm1\" onclick=\"toggleVisible('bicluster%K03d%K_pssm1'); return false;\">[+]</a>", 
                  "Show/hide bicluster #%K motif PSSM #1.</p>", 
                  "<div id=\"bicluster%K03d%K_pssm1\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                  "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #1</span>", 
                  "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                  sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                    nrow(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm), 
                    ncol(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm)), 
                  "   <div class=\"gaggle-matrix-tsv\">", "           MOTIF1", 
                  "   </div>", "</div>") else "", if (!is.null(seq.type) && 
                  length(meme.scores[[seq.type]][[k]]$meme.out) >= 
                    2 && !is.null(meme.scores[[seq.type]][[k]]$meme.out[[2]])) paste("<p><a href=\"#bicluster%K03d%K_pssm2\" onclick=\"toggleVisible('bicluster%K03d%K_pssm2'); return false;\">[+]</a>", 
                  "Show/hide bicluster #%K motif PSSM #2.</p>", 
                  "<div id=\"bicluster%K03d%K_pssm2\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                  "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #2</span>", 
                  "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                  sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                    nrow(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm), 
                    ncol(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm)), 
                  "   <div class=\"gaggle-matrix-tsv\">", "           MOTIF2", 
                  "   </div>", "</div>") else "", "</td></table>", 
                if ("pdf" %in% output) sprintf("<a href=\"../pdfs/cluster%04d.pdf\">View PDF version</a>", 
                  k) else "", "</html>"), collapse = "\n")
            rm(short.names, long.names, refseq.names, upstream.seqs)
            htmltext <- gsub("%K03d%K", sprintf("%04d", k), htmltext)
            htmltext <- gsub("%K", k, htmltext)
            htmltext <- gsub("%FILE", cmonkey.filename, htmltext)
            htmltext <- gsub("%SPECIES", gsub("_", " ", rsat.species), 
                htmltext)
            tmp <- as.data.frame(get.cluster.matrix(rows, get.cols(k)))
            tmp <- cbind(GENES = rownames(tmp), tmp)
            tf <- tempfile()
            write.table(tmp, file = tf, sep = "\t", quote = F, 
                row.names = F)
            rm(tmp)
            htmltext <- sub("RATIOS", paste(readLines(tf), collapse = "\n"), 
                htmltext)
            unlink(tf)
            if (!is.null(seq.type) && !is.null(meme.scores[[seq.type]][[k]]$meme.out)) {
                if (!is.null(meme.scores[[seq.type]][[k]]$meme.out[[1]])) {
                  tmp <- as.data.frame(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm)
                  if (!is.null(tmp) && nrow(tmp) > 0) {
                    tmp <- cbind(1:nrow(tmp), tmp)
                    colnames(tmp) <- c("POSITION", "A", "C", 
                      "G", "T")
                    write.table(tmp, file = tf, sep = "\t", quote = F, 
                      row.names = F)
                    htmltext <- sub("MOTIF1", paste(readLines(tf), 
                      collapse = "\n"), htmltext)
                    unlink(tf)
                  }
                  rm(tmp)
                }
                if (length(meme.scores[[seq.type]][[k]]$meme.out) >= 
                  2 && !is.null(meme.scores[[seq.type]][[k]]$meme.out[[2]])) {
                  tmp <- as.data.frame(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm)
                  if (!is.null(tmp) && nrow(tmp) > 0) {
                    tmp <- cbind(1:nrow(tmp), tmp)
                    colnames(tmp) <- c("POSITION", "A", "C", 
                      "G", "T")
                    write.table(tmp, file = tf, sep = "\t", quote = F, 
                      row.names = F)
                    htmltext <- sub("MOTIF2", paste(readLines(tf), 
                      collapse = "\n"), htmltext)
                    unlink(tf)
                  }
                  rm(tmp)
                }
            }
            rm(tf)
            cat(htmltext, file = sprintf("%s/htmls/cluster%04d.html", 
                out.dir, k), sep = "\n")
            rm(htmltext)
        })
        cat("\n")
    }
    if ("png" %in% output) {
        mc <- get.parallel(length(ks), para = 1)
        cat("PROFILES: ")
        for (k in ks) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d_profile.png", 
                out.dir, k))) 
                return()
            try({
                c <- get.clust(k)
                png(sprintf("%s/htmls/cluster%04d_profile.png", 
                  out.dir, k), width = 128, height = 64, antialias = "subpixel")
                par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 
                  0.75)
                try(plotCluster(c, main = "", no.par = T, ...))
                dev.off()
            }, silent = T)
        }
        cat("\n")
        cat("NETWORKS: ")
        require(igraph0)
        for (k in ks) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d_network.png", 
                out.dir, k))) 
                return()
            try({
                png(sprintf("%s/htmls/cluster%04d_network.png", 
                  out.dir, k), width = 64, height = 64, antialias = "subpixel")
                par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 
                  0.75)
                c <- get.clust(k)
                try(plotCluster.network(c, cex = 0.3, no.legend = T, 
                  ...))
                dev.off()
            }, silent = T)
        }
        cat("\n")
        if (!is.null(seq.type)) {
            cat("MOTIFS: ")
            for (k in ks) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                e.vals <- lapply(meme.scores[[seq.type]][[k]]$meme.out, 
                  "[[", "e.value")
                pssms <- lapply(meme.scores[[seq.type]][[k]]$meme.out, 
                  "[[", "pssm")
                if (length(pssms) < 2) {
                  for (i in (length(pssms) + 1):2) {
                    pssms[[i]] <- matrix(0.25, nrow = 6, ncol = 4)
                    e.vals[[i]] <- Inf
                  }
                }
                for (pp in 1:length(pssms)) {
                  if (file.exists(sprintf("%s/htmls/cluster%04d_pssm%d.png", 
                    out.dir, k, pp))) 
                    return()
                  try({
                    png(sprintf("%s/htmls/cluster%04d_pssm%d.png", 
                      out.dir, k, pp), width = 128, height = 64, 
                      antialias = "subpixel")
                    if (is.matrix(pssms[[pp]])) 
                      try(viewPssm(pssms[[pp]], e.val = NA, mot.ind = pp, 
                        main.title = sprintf("e=%.3g", e.vals[[pp]]), 
                        cex.main = 0.7), silent = T)
                    dev.off()
                  }, silent = T)
                }
            }
            cat("\n")
            cat("MOTIF POSITIONS: ")
            for (k in ks) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (file.exists(sprintf("%s/htmls/cluster%04d_mot_posns.png", 
                  out.dir, k))) 
                  return()
                try({
                  png(sprintf("%s/htmls/cluster%04d_mot_posns.png", 
                    out.dir, k), width = 128, height = 12 + 6 * 
                    length(get.rows(k)), antialias = "subpixel")
                  par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 
                    0) * 0.75)
                  c <- plotClust(k, dont.plot = T, ...)
                  try(plotClusterMotifPositions(c, cex = 0.4, 
                    no.key = T, ...))
                  dev.off()
                }, silent = F)
            }
            cat("\n")
        }
    }
    if ("main" %in% output) {
        mc <- get.parallel(length(ks))
        cat("WRITING MAIN HTML TABLE...")
        require(hwriter)
        dlf(paste(out.dir, "hwriter.css", sep = "/"), "http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css")
        dlf(paste(out.dir, "sorttable.js", sep = "/"), "http://www.kryogenix.org/code/browser/sorttable/sorttable.js")
        cat("...")
        cluster.summ <- cluster.summary(e.cutoff = NA, nrow.cutoff = 2)
        write.table(cluster.summ, file = paste(out.dir, "/cluster.summary.tsv", 
            sep = ""), quote = F, sep = "\t")
        cat("...")
        html <- openPage(paste(out.dir, "/index.html", sep = ""), 
            link.javascript = "sorttable.js", title = paste("cMonkey bicluster summary for run", 
                cmonkey.filename), link.css = "hwriter.css")
        hwrite(paste("<h2>cMonkey bicluster summary for run", 
            cmonkey.filename, "</h2>"), html)
        hwrite("<ul><li>Download a tab-delimited version of this table", 
            html, link = "cluster.summary.tsv", style = "font-size:75%")
        hwrite("<li>Download a list of each bicluster's gene members", 
            html, link = "cluster.members.genes.txt", style = "font-size:75%")
        hwrite("<li>Download a list of each bicluster's array/condition members", 
            html, link = "cluster.members.arrays.txt", style = "font-size:75%")
        cat("...")
        hwrite("<li>Plots of summary statistics of biclustering run", 
            html, link = "svgs/stats.svg", style = "font-size:75%")
        hwrite("<li>Saved cMonkey R session file", html, link = "cm_session.RData", 
            style = "font-size:75%")
        hwrite("<li>Summary of cMonkey input parameters</ul>", 
            html, link = "cm.params.txt", style = "font-size:75%")
        hwrite("<br><center><b>Bicluster summary</b></center><br>", 
            html)
        hwrite("<br><center><b>Sort the table by a given column by clicking on the column's header.<br>Click on bicluster link in first column for more info.</b></center><br>", 
            html, style = "font-size:60%")
        cat("...")
        himg0 <- hwriteImage(sprintf("htmls/cluster%04d_profile.png", 
            as.integer(rownames(cluster.summ))), table = F)
        himg0 <- hwrite(paste(himg0, sprintf("Residual = %.3f", 
            cluster.summ$resid), sep = "<br>"), center = TRUE, 
            table = F)
        himg0a <- hwriteImage(sprintf("htmls/cluster%04d_network.png", 
            as.integer(rownames(cluster.summ))), table = F)
        if (!is.null(seq.type)) {
            e.val.1 <- lapply(meme.scores[[seq.type]][as.integer(rownames(cluster.summ))], 
                function(i) {
                  mo <- i$meme.out
                  if (length(mo) >= 1) 
                    mo[[1]]$e.value
                  else Inf
                })
            for (i in 1:length(e.val.1)) if (is.null(e.val.1[[i]])) 
                e.val.1[[i]] <- NA
            himg1 <- hwriteImage(sprintf("htmls/cluster%04d_pssm1.png", 
                as.integer(rownames(cluster.summ))), table = F, 
                title = sprintf("E-val = %.3g", unlist(e.val.1)))
            himg1 <- hwrite(paste(himg1, as.character(cluster.summ$consensus1), 
                sep = "<br>"), center = TRUE, table = F)
            if (!is.null(seq.type)) 
                e.val.2 <- lapply(meme.scores[[seq.type]][as.integer(rownames(cluster.summ))], 
                  function(i) {
                    mo <- i$meme.out
                    if (length(mo) > 1) 
                      mo[[2]]$e.value
                    else Inf
                  })
            else e.val.2 <- as.list(rep(NA, k.clust))
            for (i in 1:length(e.val.2)) if (is.null(e.val.2[[i]])) 
                e.val.2[[i]] <- NA
            himg2 <- hwriteImage(sprintf("htmls/cluster%04d_pssm2.png", 
                as.integer(rownames(cluster.summ))), table = F, 
                title = sprintf("E-val = %.3g", unlist(e.val.2)))
            himg2 <- hwrite(paste(himg2, as.character(cluster.summ$consensus2), 
                sep = "<br>"), center = TRUE, table = F)
            himg2a <- hwriteImage(sprintf("htmls/cluster%04d_mot_posns.png", 
                as.integer(rownames(cluster.summ))), table = F)
            e.val.1[is.na(e.val.1)] <- 9e+09
            e.val.2[is.na(e.val.2)] <- 9e+09
        }
        else {
            e.val.1 <- e.val.2 <- as.list(rep(NA, k.clust))
            himg1 <- himg2 <- himg2a <- NULL
        }
        cluster.summ$score <- sprintf("%.3f", cluster.summ$score)
        rn <- rownames(cluster.summ)
        cat("...")
        cluster.summ.orig <- cluster.summ
        cluster.summ <- cbind(bicluster = cluster.summ$k, n.genes = cluster.summ$nrow, 
            n.arrays = sapply(as.integer(rownames(cluster.summ)), 
                function(i) length(get.cols(i))), score = cluster.summ$score, 
            residual = sprintf("%.3f", cluster.summ$resid))
        if ("score.norm" %in% colnames(cluster.summ.orig)) 
            cluster.summ <- cbind(cluster.summ, score.norm = sprintf("%.3f", 
                cluster.summ.orig$score.norm))
        rownames(cluster.summ) <- rn
        rows <- list()
        for (k in as.integer(rn)) rows[[k]] <- sort(get.rows(k))
        himg3 <- hwrite(sapply(as.integer(rn), function(k) paste(rows[[k]], 
            collapse = " ")), table = F)
        cat("...\n")
        if (!no.genome.info) {
            himg4 <- hwrite(unlist(mc$apply(as.integer(rn), function(k) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (length(rows[[k]]) <= 0) 
                  return()
                tmp <- get.long.names(rows[[k]], short = T)
                tmp <- unique(tmp[!tmp %in% rows[[k]] & tmp != 
                  ""])
                paste(tmp, collapse = " ")
            })), table = F)
            cat("\n")
            himg5 <- hwrite(unlist(mc$apply(as.integer(rn), function(k) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (length(rows[[k]]) <= 0) 
                  return()
                tmp <- get.long.names(rows[[k]], short = F)
                tmp <- unique(tmp[!tmp %in% rows[[k]] & tmp != 
                  ""])
                paste(tmp, collapse = " | ")
            })), table = F)
            cat("\n")
        }
        else {
            himg4 <- himg5 <- NULL
        }
        nas <- rep(NA, nrow(cluster.summ))
        hwrite(cbind(cluster.summ[, 1:min(ncol(cluster.summ), 
            6)], profile = himg0, network = himg0a, motif1 = himg1, 
            motif2 = himg2, motif.posns = himg2a, probe.names = himg3, 
            short.names = himg4, long.names = himg5), html, row.names = F, 
            table.style = "text-align:center;font-size:70%;font-family:Arial", 
            table.class = "sortable", row.style = list("font-weight:bold;text-align:center;font-size:70"), 
            col.style = list(probe.names = "font-size:70%", orf.names = "font-size:50%", 
                short.names = "font-size:50%", long.names = "font-size:50%", 
                motif1 = "font-size:50%", motif2 = "font-size:50%"), 
            col.sorttable_customkey = list(residual = sprintf("%.3f", 
                cluster.summ.orig$residual), score.norm = if ("score.norm" %in% 
                colnames(cluster.summ.orig)) sprintf("%.3f", 
                cluster.summ.orig$score.norm) else NULL, profile = sprintf("%.3f", 
                cluster.summ.orig$resid), motif1 = sprintf("%.30f", 
                unlist(e.val.1)), e.val1 = sprintf("%.30f", unlist(e.val.1)), 
                motif2 = sprintf("%.30f", unlist(e.val.2)), e.val2 = sprintf("%.30f", 
                  unlist(e.val.2))), col.class = list(network = c("sorttable_nosort", 
                nas), motif.posns = c("sorttable_nosort", nas)), 
            col.link = list(sprintf("htmls/cluster%04d.html", 
                as.integer(rownames(cluster.summ)))))
        closePage(html, splash = F)
        for (i in sapply(1:k.clust, function(k) c(k, sort(get.rows(k))))) cat(i, 
            "\n", file = paste(out.dir, "/cluster.members.genes.txt", 
                sep = ""), append = T)
        for (i in sapply(1:k.clust, function(k) c(k, sort(get.cols(k))))) cat(i, 
            "\n", file = paste(out.dir, "/cluster.members.arrays.txt", 
                sep = ""), append = T)
        tmp <- capture.output(for (name in ls(cmonkey.params)) {
            cat(name, "= ")
            str(get(name, envir = cmonkey.params), no.list = T)
        })
        cat(tmp, file = paste(out.dir, "/cm.params.txt", sep = ""), 
            sep = "\n", collapse = "\n")
    }
    if (gzip) {
        rpl <- function(find, replace, file, ...) {
            f <- readLines(file)
            f <- gsub(find, replace, f, ...)
            writeLines(f, con = file)
        }
        system(sprintf("gzip -v %s/svgs/*.svg", out.dir))
        for (f in list.files(paste(out.dir, "/svgs", sep = ""), 
            full = T)) if (grepl(".svg.gz", f, fixed = T)) 
            system(sprintf("mv -v %s %s", f, sub(".svg.gz", ".svgz", 
                f, fixed = T)))
        lapply(c(list.files(sprintf("%s/htmls", out.dir), pattern = glob2rx("*.html"), 
            full = T), list.files(out.dir, pattern = glob2rx("*.html"), 
            full = T)), function(f) {
            cat(f, "\n")
            rpl(".svg\"", ".svgz\"", f, fixed = T)
        })
    }
    if ("rdata" %in% output) 
        save.cmonkey.env(file = paste(out.dir, "/cm_session.RData", 
            sep = ""))
    out.dir
}
