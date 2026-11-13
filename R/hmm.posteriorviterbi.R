#########################
#Viterbi most probable track
##########################

###########################################################
#This code finds the most probable sequence of grid cells through the posteriors of HMMoce
###########################################################
#' Collaborative development by Drs. Julie Nielsen, Martin Arostegui, Camrin Braun, and Ben Galuardi
#' Cite: 
#' Nielsen JK, Bryan DR, Rand KM, Arostegui MC, Braun CD, Galuardi B, McDermott SF (2023) Geolocation of a demersal fish (Pacific cod) in a high-latitude island chain (Aleutian Islands, Alaska). Animal Biotelemetry 11: 29.
#'
#' Parameters
#' @param pars HMM parameters output from opt.params()
#' @param s HMM backward smoothed output of hmm.smoother()
#' @param g original grid on which HMM was run, assigned from L.res$g
#' @param iniloc dataframe of tag and popup locations/dates (2 rows, 5 columns)
#' @param dateVec timestep sequence used in HMM
#' @param lims.lat vector of c(min,max) desired latitude for reduced spatial extent in which to run algorithm
#' @param lims.lon vector of c(min,max) desired longitude for reduced spatial extent in which to run algorithm
#' @param percentile uncertainty percentile of the daily posterior to include in admissible paths (default is 99.9% of daily posterior)
#' @param approach for two-state outputs, use a priori most likely behavior state sequence with separate posteriors and kernels ('separate'), or joint posteriors and kernels weighted by relative behavior state probability ('joint')
#' 
hmm.posteriorviterbi <- function(pars, s, g, iniloc, dateVec, lims.lat, lims.lon, percentile = 99.9, approach = "separate") {
  if (!approach %in% c("joint", "separate")) return(NULL)

  thresh <- 1 - percentile / 100
  d_s <- dim(s)
  n_behav <- d_s[1L]
  T_len <- d_s[2L]
  n_col <- d_s[3L]
  n_row <- d_s[4L]
  n_pix <- n_col * n_row

  if (approach == "joint") {
    if (n_behav > 1L) {
      state_sums <- rowSums(s, dims = 2L)
      p.resid <- state_sums[2L, ]
      p.migrant <- 1 - p.resid
    } else {
      p.resid <- rep(NA, length.out = length(dateVec))
    }

    sp <- aperm(s, c(3L, 4L, 2L, 1L))
    if (n_behav == 2L) {
      M1 <- matrix(sp[, , , 1L], nrow = n_pix, ncol = T_len)
      M2 <- matrix(sp[, , , 2L], nrow = n_pix, ncol = T_len)
      p1 <- state_sums[1L, ]
      p2 <- state_sums[2L, ]
      M <- matrix(0, nrow = n_pix, ncol = T_len)
      for (i in seq_len(T_len)) {
        c1 <- M1[, i]
        c2 <- M2[, i]
        n1 <- (c1 / max(c1, na.rm = TRUE)) * p1[i]
        n2 <- (c2 / max(c2, na.rm = TRUE)) * p2[i]
        n1[is.na(n1)] <- 0
        n2[is.na(n2)] <- 0
        M[, i] <- n1 + n2
      }
      M <- .pv_threshold_slices(M, thresh)
      sp <- array(M, dim = c(n_col, n_row, T_len, 1L))
    } else {
      M <- .pv_threshold_slices(matrix(sp, nrow = n_pix, ncol = T_len), thresh)
      sp <- array(M, dim = c(n_col, n_row, T_len, 1L))
    }
  } else {
    sp <- aperm(s, c(3L, 4L, 2L, 1L))
    M <- .pv_threshold_slices(matrix(sp, nrow = n_pix, ncol = T_len * n_behav), thresh)
    sp <- array(M, dim = c(n_col, n_row, T_len, n_behav))
  }

  keep_lat <- which(g$lat[, 1L] >= lims.lat[1L] & g$lat[, 1L] <= lims.lat[2L])
  keep_lon <- which(g$lon[1L, ] >= lims.lon[1L] & g$lon[1L, ] <= lims.lon[2L])

  L.vit <- aperm(sp[keep_lon, rev(keep_lat), , , drop = FALSE], c(2L, 1L, 3L, 4L))

  longs.vit <- g$lon[1L, keep_lon]
  lats.vit <- g$lat[keep_lat, 1L]
  rows.vit <- dim(L.vit)[1L]

  ilo.vit <- which.min(abs(iniloc[1L, "lon"] - longs.vit))
  ila.vit <- which.min(abs(iniloc[1L, "lat"] - rev(lats.vit)))

  mult <- if (approach == "joint") 6 else 4
  sigmas <- if (length(pars) == 4L) pars[1:2] else pars[1L]
  sizes <- rep(ceiling(sigmas[1L] * mult), 2L)
  if (sizes[1L] %% 2L == 0L) sizes[1L] <- sizes[1L] + 1L
  if (sizes[2L] %% 2L == 0L) sizes[2L] <- sizes[2L] + 1L

  K1 <- gausskern.pg(sizes[1L], sigmas[1L], muadv = 0)
  if (length(pars) == 4L) {
    K2 <- gausskern.pg(sizes[2L], sigmas[2L], muadv = 0)
    K <- list(K1, K2)
  } else {
    K <- list(K1)
  }

  if (approach == "joint") {
    sv_or_weights <- if (is.na(p.resid[1L])) rep(1L, T_len) else list(p.migrant = p.migrant, p.resid = p.resid)
    mpt <- .pv_viterbi_core(L.vit, if (is.na(p.resid[1L])) list(K[[1L]]) else K, sv_or_weights, sigmas[1L], ila.vit, ilo.vit)
  } else {
    if (n_behav == 2L) {
      st_sums <- matrix(colSums(M), nrow = T_len, ncol = 2L)
      sv <- (-(st_sums[, 1L] > st_sums[, 2L]) + 1L) + 1L
      p.resid <- st_sums[, 2L]
    } else {
      sv <- rep(1L, T_len)
      p.resid <- rep(NA, length.out = length(dateVec))
    }
    mpt <- .pv_viterbi_core(L.vit, if (n_behav == 2L) K else list(K[[1L]]), sv, sigmas[1L], ila.vit, ilo.vit)
  }

  data.frame(
    "date" = dateVec,
    "lat" = lats.vit[(rows.vit + 1L) - mpt[, 1L]],
    "lon" = longs.vit[mpt[, 2L]],
    "p.resid" = p.resid
  )
}

.pv_threshold_slices <- function(M, thresh) {
  for (j in seq_len(ncol(M))) {
    col_j <- M[, j]
    rg <- range(col_j, na.rm = TRUE)
    dr <- rg[2L] - rg[1L]
    below <- if (dr != 0) {
      ((col_j - rg[1L]) / dr) < thresh
    } else if (!is.na(col_j[1L]) && col_j[1L] > 1) {
      (col_j * 0 + 1) < thresh
    } else if (!is.na(col_j[1L]) && col_j[1L] < 0) {
      (col_j * 0) < thresh
    } else {
      col_j < thresh
    }
    idx <- which(below)
    if (length(idx) > 0L) {
      col_j[idx] <- 0
      M[, j] <- col_j
    }
  }
  M
}

.pv_prep_kern <- function(kern, allow_sep = TRUE) {
  kdim <- nrow(kern)
  h <- kdim %/% 2L
  if (allow_sep && all(kern > 0)) {
    hc <- 0.5 * log(kern[h + 1L, h + 1L])
    list(
      sep = TRUE, kdim = kdim, h = h, min_pos = min(kern),
      lr = log(kern[, h + 1L]) - hc,
      lc = log(kern[h + 1L, ]) - hc,
      kern = kern
    )
  } else {
    pos <- kern > 0
    lk <- matrix(-Inf, kdim, kdim)
    lk[pos] <- log(kern[pos])
    list(
      sep = FALSE, kdim = kdim, h = h,
      min_pos = if (any(pos)) min(kern[pos]) else 0,
      log_kern = lk, pos = pos, kern = kern
    )
  }
}

.pv_viterbi_core <- function(L.vit, K_list, sv_or_weights, sigmas1, ila.vit, ilo.vit) {
  rows.vit <- dim(L.vit)[1L]
  cols.vit <- dim(L.vit)[2L]
  timesteps.vit <- dim(L.vit)[3L]
  behavs.vit <- dim(L.vit)[4L]

  bp_tgt <- vector("list", timesteps.vit)
  bp_src <- vector("list", timesteps.vit)

  start_lin <- as.integer(ila.vit + (ilo.vit - 1L) * rows.vit)
  bp_tgt[[1L]] <- start_lin
  bp_src[[1L]] <- start_lin

  prev_prob <- matrix(-Inf, nrow = rows.vit, ncol = cols.vit)
  prev_prob[ila.vit, ilo.vit] <- 0.0
  upd_lin <- start_lin
  r_min <- ila.vit; r_max <- ila.vit
  c_min <- ilo.vit; c_max <- ilo.vit

  is_joint_2state <- is.list(sv_or_weights)
  K_prep <- if (!is_joint_2state) lapply(K_list, .pv_prep_kern, allow_sep = TRUE) else NULL

  if (timesteps.vit >= 2L) {
    for (i in 2:timesteps.vit) {
      if (is_joint_2state) {
        kern <- K_list[[1L]] * sv_or_weights$p.migrant[i] + K_list[[2L]] * sv_or_weights$p.resid[i]
        kp <- .pv_prep_kern(kern, allow_sep = FALSE)
        L_day <- L.vit[, , i, 1L]
      } else {
        st <- if (length(K_prep) == 1L) 1L else sv_or_weights[i]
        kp <- K_prep[[st]]
        L_day <- if (behavs.vit == 1L || st == 1L) L.vit[, , i, 1L] else L.vit[, , i, 2L]
      }

      h <- kp$h
      kdim <- kp$kdim
      cur_prob <- matrix(-Inf, nrow = rows.vit, ncol = cols.vit)
      any_updated <- FALSE

      if (length(upd_lin) > 0L && kp$min_pos > 0) {
        tr0 <- max(1L, r_min - h)
        tr1 <- min(rows.vit, r_max + h)
        tc0 <- max(1L, c_min - h)
        tc1 <- min(cols.vit, c_max + h)

        L_sub <- L_day[tr0:tr1, tc0:tc1, drop = FALSE]
        pos_L <- !is.na(L_sub) & (L_sub > 0)

        if (any(pos_L)) {
          ry <- which(rowSums(pos_L) > 0L)
          cx <- which(colSums(pos_L) > 0L)
          ry0 <- ry[1L]; ry1 <- ry[length(ry)]
          cx0 <- cx[1L]; cx1 <- cx[length(cx)]

          if (ry0 > 1L || ry1 < nrow(L_sub) || cx0 > 1L || cx1 < ncol(L_sub)) {
            tr0 <- tr0 + ry0 - 1L
            tr1 <- tr0 + (ry1 - ry0)
            tc0 <- tc0 + cx0 - 1L
            tc1 <- tc0 + (cx1 - cx0)
            L_sub <- L_sub[ry0:ry1, cx0:cx1, drop = FALSE]
            pos_L <- pos_L[ry0:ry1, cx0:cx1, drop = FALSE]
          }

          nr <- tr1 - tr0 + 1L
          nc <- tc1 - tc0 + 1L
          stride <- nr + 2L * h

          P_pad <- matrix(-Inf, nrow = stride, ncol = nc + 2L * h)
          sr0 <- max(1L, tr0 - h)
          sr1 <- min(rows.vit, tr1 + h)
          sc0 <- max(1L, tc0 - h)
          sc1 <- min(cols.vit, tc1 + h)
          pr0 <- sr0 - (tr0 - h) + 1L
          pr1 <- pr0 + (sr1 - sr0)
          pc0 <- sc0 - (tc0 - h) + 1L
          pc1 <- pc0 + (sc1 - sc0)
          P_pad[pr0:pr1, pc0:pc1] <- prev_prob[sr0:sr1, sc0:sc1]

          no_underflow <- (min(L_sub[pos_L]) * kp$min_pos) > 0

          if (kp$sep && no_underflow) {
            nc_pad <- nc + 2L * h
            tmp_prob <- matrix(-Inf, nrow = nr, ncol = nc_pad)
            tmp_ky <- matrix(0L, nrow = nr, ncol = nc_pad)
            lr <- kp$lr
            for (ky in kdim:1L) {
              dy <- ky - (h + 1L)
              cand <- P_pad[(h + 1L - dy):(h + nr - dy), , drop = FALSE] + lr[ky]
              b <- which(cand > tmp_prob)
              if (length(b) > 0L) {
                tmp_prob[b] <- cand[b]
                tmp_ky[b] <- ky
              }
            }

            best_trans <- matrix(-Inf, nrow = nr, ncol = nc)
            best_ky <- matrix(0L, nrow = nr, ncol = nc)
            best_kx <- matrix(0L, nrow = nr, ncol = nc)
            lc <- kp$lc
            for (kx in kdim:1L) {
              dx <- kx - (h + 1L)
              c_seq <- (h + 1L - dx):(h + nc - dx)
              cand <- tmp_prob[, c_seq, drop = FALSE] + lc[kx]
              b <- which(cand > best_trans)
              if (length(b) > 0L) {
                best_trans[b] <- cand[b]
                best_kx[b] <- kx
                best_ky[b] <- tmp_ky[, c_seq, drop = FALSE][b]
              }
            }

            win <- which(pos_L & (best_kx > 0L))
            if (length(win) > 0L) {
              any_updated <- TRUE
              loc_r <- ((win - 1L) %% nr) + 1L
              loc_c <- ((win - 1L) %/% nr) + 1L
              ty <- tr0 + loc_r - 1L
              tx <- tc0 + loc_c - 1L
              sy <- ty - (best_ky[win] - (h + 1L))
              sx <- tx - (best_kx[win] - (h + 1L))
              tgt_lin <- as.integer(ty + (tx - 1L) * rows.vit)
              bp_tgt[[i]] <- tgt_lin
              bp_src[[i]] <- as.integer(sy + (sx - 1L) * rows.vit)
              cur_prob[tgt_lin] <- best_trans[win] + log(L_sub[win])
            }
          } else if (no_underflow) {
            P_vec <- as.vector(cbind(-Inf, P_pad, -Inf))
            total_len <- stride * nc
            best_vec <- rep.int(-Inf, total_len)
            best_off_vec <- integer(total_len)
            lk <- kp$log_kern
            pos_k <- kp$pos
            for (kx in kdim:1L) {
              dx <- kx - (h + 1L)
              base_off <- (h + 1L - dx) * stride
              for (ky in kdim:1L) {
                if (!pos_k[ky, kx]) next
                off <- base_off - (ky - (h + 1L))
                cand <- P_vec[(off + 1L):(off + total_len)] + lk[ky, kx]
                b <- which(cand > best_vec)
                if (length(b) > 0L) {
                  best_vec[b] <- cand[b]
                  best_off_vec[b] <- ky + (kx - 1L) * kdim
                }
              }
            }
            inner_rows <- (h + 1L):(h + nr)
            best_trans <- matrix(best_vec, nrow = stride, ncol = nc)[inner_rows, , drop = FALSE]
            best_off <- matrix(best_off_vec, nrow = stride, ncol = nc)[inner_rows, , drop = FALSE]
            win <- which(pos_L & (best_off > 0L))
            if (length(win) > 0L) {
              any_updated <- TRUE
              off_w <- best_off[win]
              ky_w <- ((off_w - 1L) %% kdim) + 1L
              kx_w <- ((off_w - 1L) %/% kdim) + 1L
              loc_r <- ((win - 1L) %% nr) + 1L
              loc_c <- ((win - 1L) %/% nr) + 1L
              ty <- tr0 + loc_r - 1L
              tx <- tc0 + loc_c - 1L
              sy <- ty - (ky_w - (h + 1L))
              sx <- tx - (kx_w - (h + 1L))
              tgt_lin <- as.integer(ty + (tx - 1L) * rows.vit)
              bp_tgt[[i]] <- tgt_lin
              bp_src[[i]] <- as.integer(sy + (sx - 1L) * rows.vit)
              cur_prob[tgt_lin] <- best_trans[win] + log(L_sub[win])
            }
          } else {
            best_prob <- matrix(-Inf, nrow = nr, ncol = nc)
            best_off <- matrix(0L, nrow = nr, ncol = nc)
            kern <- kp$kern
            pos_k <- kp$pos
            L_clean <- ifelse(pos_L, L_sub, 0)
            for (kx in kdim:1L) {
              dx <- kx - (h + 1L)
              P_col <- P_pad[, (h + 1L - dx):(h + nc - dx), drop = FALSE]
              for (ky in kdim:1L) {
                if (!pos_k[ky, kx]) next
                dy <- ky - (h + 1L)
                prod_sub <- L_clean * kern[ky, kx]
                pm <- prod_sub > 0
                if (!any(pm)) next
                cand <- log(prod_sub) + P_col[(h + 1L - dy):(h + nr - dy), , drop = FALSE]
                b <- which(pm & (cand > best_prob))
                if (length(b) > 0L) {
                  best_prob[b] <- cand[b]
                  best_off[b] <- ky + (kx - 1L) * kdim
                }
              }
            }
            win <- which(best_off > 0L)
            if (length(win) > 0L) {
              any_updated <- TRUE
              off_w <- best_off[win]
              ky_w <- ((off_w - 1L) %% kdim) + 1L
              kx_w <- ((off_w - 1L) %/% kdim) + 1L
              loc_r <- ((win - 1L) %% nr) + 1L
              loc_c <- ((win - 1L) %/% nr) + 1L
              ty <- tr0 + loc_r - 1L
              tx <- tc0 + loc_c - 1L
              sy <- ty - (ky_w - (h + 1L))
              sx <- tx - (kx_w - (h + 1L))
              tgt_lin <- as.integer(ty + (tx - 1L) * rows.vit)
              bp_tgt[[i]] <- tgt_lin
              bp_src[[i]] <- as.integer(sy + (sx - 1L) * rows.vit)
              cur_prob[tgt_lin] <- best_prob[win]
            }
          }
        }
      }

      # Fallback if bounded kernel support did not reach any positive cell in L_day
      if (!any_updated) {
        if (length(upd_lin) > 0L) {
          best_prev_lin <- upd_lin[which.max(prev_prob[upd_lin])]
          bp_r <- ((best_prev_lin - 1L) %% rows.vit) + 1L
          bp_c <- ((best_prev_lin - 1L) %/% rows.vit) + 1L
          pos_l <- which(!is.na(L_day) & L_day > 0, arr.ind = TRUE)
          if (nrow(pos_l) == 0L) pos_l <- matrix(c(bp_r, bp_c), nrow = 1L)
          d2 <- (pos_l[, 1L] - bp_r)^2 + (pos_l[, 2L] - bp_c)^2
          sig_eff <- max(sigmas1, 1)
          lp <- log(pmax(L_day[pos_l], 1e-300)) - 0.5 * (d2 / (sig_eff^2)) + prev_prob[bp_r, bp_c]
          lin_l <- as.integer(pos_l[, 1L] + (pos_l[, 2L] - 1L) * rows.vit)
          bp_tgt[[i]] <- lin_l
          bp_src[[i]] <- rep.int(as.integer(best_prev_lin), length(lin_l))
          cur_prob[lin_l] <- lp
          any_updated <- TRUE
        }
      }

      prev_prob <- cur_prob
      if (any_updated) {
        upd_lin <- bp_tgt[[i]]
        upd_r <- ((upd_lin - 1L) %% rows.vit) + 1L
        upd_c <- ((upd_lin - 1L) %/% rows.vit) + 1L
        r_min <- min(upd_r); r_max <- max(upd_r)
        c_min <- min(upd_c); c_max <- max(upd_c)
      } else {
        upd_lin <- integer(0L)
        r_min <- NA_integer_; r_max <- NA_integer_
        c_min <- NA_integer_; c_max <- NA_integer_
      }
    }
  }

  # Reconstruct MPT from grid cell at last time step with highest cumulative probability
  mpt <- matrix(NA_integer_, nrow = timesteps.vit, ncol = 2L)
  end_lin <- upd_lin[which.max(prev_prob[upd_lin])]
  mpt[timesteps.vit, 1L] <- ((end_lin - 1L) %% rows.vit) + 1L
  mpt[timesteps.vit, 2L] <- ((end_lin - 1L) %/% rows.vit) + 1L

  prev.loc <- bp_src[[timesteps.vit]][match(end_lin, bp_tgt[[timesteps.vit]])]
  if (timesteps.vit >= 2L) {
    for (i in (timesteps.vit - 1L):1L) {
      mpt[i, 1L] <- ((prev.loc - 1L) %% rows.vit) + 1L
      mpt[i, 2L] <- ((prev.loc - 1L) %/% rows.vit) + 1L
      prev.loc <- bp_src[[i]][match(prev.loc, bp_tgt[[i]])]
    }
  }
  mpt
}
