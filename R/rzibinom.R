#'
#'
#' @author VGAM package

qzibinom<-function (p, size, prob, pstr0 = 0)
{
  LLL <- max(length(p), length(size), length(prob), length(pstr0))
  p <- rep_len(p, LLL)
  size <- rep_len(size, LLL)
  prob <- rep_len(prob, LLL)
  pstr0 <- rep_len(pstr0, LLL)
  ans <- p
  ans[p <= pstr0] <- 0
  ans[p > pstr0] <- qbinom((p[p > pstr0] - pstr0[p > pstr0])/(1 -
                                                                pstr0[p > pstr0]), size[p > pstr0], prob[p > pstr0])
  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0/(1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposbinom((p[pindex] - Pobs0)/(1 - Pobs0),
                             size = size[pindex], prob = prob[pindex])
  }
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN
  ans
}


rzibinom<-function (n, size, prob, pstr0 = 0)
{
  qzibinom(runif(n), size, prob, pstr0 = pstr0)
}


qposbinom<-function (p, size, prob){
  ans <- qbinom(pbinom(0, size, prob, lower.tail = FALSE) *
                  p + dbinom(0, size, prob), size, prob)
  ans[p == 1] <- size[p == 1]
  ans[p == 0] <- 1
  ans[prob == 0] <- NaN
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans
}
