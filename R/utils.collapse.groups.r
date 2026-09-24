# Groups of populations joined by fixed difference counts <= tpop, following
# chains (connected components). Returns a list of sorted name vectors,
# one per group, ordered by their first name. NA distances do not join.
utils.collapse.groups <- function(mat, tpop) {
    pops <- rownames(mat)
    grp <- seq_along(pops)
    repeat {
        changed <- FALSE
        for (i in seq_along(pops)) {
            for (j in seq_along(pops)) {
                if (!is.na(mat[i, j]) && mat[i, j] <= tpop &&
                    grp[i] != grp[j]) {
                    lo <- min(grp[i], grp[j])
                    grp[grp == max(grp[i], grp[j])] <- lo
                    changed <- TRUE
                }
            }
        }
        if (!changed) break
    }
    groups <- lapply(split(pops, grp), sort)
    groups <- unname(groups[order(vapply(groups, `[`, "", 1))])
    groups
}
