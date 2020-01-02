## direct methods for sparse linear systems
## episode 03 - chapter 2

# multiplication of sparse matrix by dense vector
# Y = Ax, where A is a sparse matrix and Y & x is are dense vectors

A <- matrix(round(runif(100, 0, 10)), 10, 10)
A <- matrix(c(2, 0, 0, -1, 0, -2, 0, 3, 8), 3, 3)
x <- c(1, 2, 3)


# standard matrix multiplication

`%*s%` <- function(x, y){
  if(is.vector(x)) x <- matrix(x, ncol = 1) # convert x to column vector
  if(is.vector(y)) y <- matrix(y, ncol = 1) # convert y to column vector
  
  xi <- nrow(x)
  xj <- ncol(x)
  yi <- nrow(y)
  yj <- ncol(y)
  
  out <- matrix(0, nrow = xi, ncol = yj)
  if(xj != yi) stop(gettext("non-conformable arguments"))
  
  for(j in seq_len(yj)){
    for(i in seq_len(xi)){
      out[i, j] <- sum(x[i, ] * y[, j])
    }
  }
  out
}
A%*%x
A%*s%x

## coerce to sparse matrix
setClass("CCsparse", slots = c(i = "numeric", p = "numeric", Dim = "numeric", x = "numeric"))
setClass("TRIPsparse", slots = c(i = "numeric", j = "numeric", Dim = "numeric", x = "numeric"))


## cumulative sum function -- to be effective for column pointers, must be
## executed on concatenation of row counts with 1 (e.g., c(1, row_counts))
sp_cumsum <- function(x){
  nz <- 0
  for(i in 1:length(x)){
    x[i] <- nz + x[i]
    nz <- x[i]
  }
  x
}



## function for producing sparse matrices
## initializes all matrices in triplet form, then converts to column-compressed if type = "cc"
mksparse <- function(x, type = "cc", print = FALSE){
  if(type != "cc" & type != "trip") stop(gettext("type must be column-compressed 'cc' or triplet 'trip'"))
  
  Dim  <- c(nrow(x), ncol(x))
  xx <- as.vector(x)
  xj <- as.vector(col(x))[xx != 0]
  xi <- as.vector(row(x))[xx != 0]
  xx <- xx[xx != 0]
  
  nz   <- length(xx)

  ## preparing column-compressed form -- the default object
  if(type == "cc"){
    p    <- vector(mode = "numeric", length = Dim[2])
    for(i in 1:nz){
      p[xj[i]] <- p[xj[i]] + 1
    }
    p <- sp_cumsum(c(1, p))
    
    out <- list(i = xi, p = p, Dim = Dim, x = xx)
  } else {
    out <- list(i = xi, j = xj, Dim = Dim, x = xx)
  }
  
  if(print){
    x.mat <- matrix(nrow = Dim[1], ncol = Dim[2])
    x.mat[cbind(xi, xj)] <- xx
    x.mat[is.na(x.mat)] <- '.'
    print(noquote(x.mat))
  }
  
  out
}
mksparse(x=A)


mkTRIPsparse <- function(x){
  if(class(x) == "matrix" | class(x) == "data.frame"){
    nz   <- length(x[x != 0])
    Dim  <- c(nrow(x), ncol(x))
    xi   <- vector(mode = "numeric", length = nz)
    xx   <- vector(mode = "numeric", length = nz)
    ind  <- 1L
    
    xj <- vector(mode = "numeric", length = nz)
    for(j in seq_len(Dim[2])){
      for(k in seq_len(Dim[1])){
        if(x[k, j] != 0){
          xx[ind] <- x[k, j]
          xi[ind] <- k
          xj[ind] <- j
          ind <- ind + 1
        }
      }
    }
    out <- list(i = xi, j = xj, Dim = Dim, x = xx)
  } else {
    nz  <- length(x$x)
    Dim <- x$Dim
    p   <- x$p
    xi  <- x$i
    xj  <- vector(mode = "numeric", length = nz)
    ind <- 1L
    
    for(j in 1:Dim[2]){
      if(p[j+1] > p[j]){
        xj[p[j]:(p[j+1]-1)] <- j
      }
    }
    out <- list(i = xi, j = xj, Dim = Dim, x = x$x)
  }
  out
}

## how to adapt to sparse system?
# block matrix multiplication - you can take matrix A and matrix B, slice them into conformable 
# chunks--rows and columns of A and rows and columns of B--such that multiplying the submatrices
# produces a partitioned matrix C
A <- matrix(round(runif(100, 0, 10)), 10, 10)
A[cbind(sample(A, size = 60), sample(A, 60))] <- 0
B <- matrix(round(runif(30, 0, 10)), 10, 3)
B[cbind(sample(B, size = 15), c(2, 3, 1, 2, 1))] <- 0

# sparse matrix times a dense vector 'gaxpy'
`%*s%` <- function(x, y, plot = F){
  if(is.vector(x)) x <- matrix(x, ncol = 1) # convert x to column vector
  if(is.vector(y)) y <- matrix(y, ncol = 1) # convert y to column vector
  
  if(class(x) != "sparsematrix") x <- mksparse(x)
  n.row <- x$Dim[1]
  n.col <- x$Dim[2]
  xi    <- x$i
  xp    <- x$p
  x     <- x$x
  yi    <- nrow(y)
  yj    <- ncol(y)
  
  out <- matrix(0, nrow = xi, ncol = yj)
  if(n.col != yi) stop(gettext("non-conformable arguments"))
  
  y.new <- vector(mode = "numeric", length = n.col)
  for(j in 1:n.col){
    for(p in xp[j]:(xp[j+1]-1)){
      y.new[xi[p]] <- y.new[xi[p]] + x[p]*y[j]
    }
  }
  y.new
}
A%*%c(1,2,3,4,5,6,7,8,9,0)
A%*s%c(1,2,3,4,5,6,7,8,9,0)



## sparse transpose -- lecture 06 on youtube, section 2.5 in book
sp_t <- function(x, print = FALSE){
  
  ## operations on original matrix
  Ax  <- x$x # extract x values 
  Ai  <- x$i # extract row indices
  Ap  <- x$p # extract column pointers
  nz  <- length(Ax) # number of nonzeros
  Dim <- x$Dim # extract row and column dimensions
  
  ## initializing new matrix
  xx  <- vector(mode = "numeric", length = nz) # x vector
  xi  <- vector(mode = "numeric", length = nz) # row index vector
  xp  <- vector(mode = "numeric", length = Dim[1] + 1) # column pointer vector

  ## initialize workspace vector
  w <- vector(mode = "numeric", length = Dim[1])
  for(p in 1:nz){             # loop over all entries of row index vector of old matrix 
    w[Ai[p]] <- w[Ai[p]] + 1  # add 1 each time a row is observed to produce row counts
  }
  
  w <- xp <- sp_cumsum(c(1, w)) # take cumulative sum of row counts to produce column pointers
                                # stored twice--one to operate over and one to return in output
  
  for(j in 1:Dim[2]){               # loop over columns
    if(Ap[j] < Ap[j+1]){            # skip zero columns
      for(p in Ap[j]:(Ap[j+1]-1)){  # loop over rows within columns
        xi[w[Ai[p]]] <- j           # place column index 'j' in xi where row index is present
        xx[w[Ai[p]]] <- Ax[p]       # same as above, but for x values
        w[Ai[p]] <- w[Ai[p]] + 1    # postincrement indicator
      }
    }
  }

  if(print){
    x.mat <- matrix(nrow = Dim[2], ncol = Dim[1])
    x.mat[cbind(xi, sort(Ai))] <- xx
    x.mat[is.na(x.mat)] <- '.'
    print(noquote(x.mat))
  }
  out <- list(i = xi, p = xp, Dim = c(Dim[2], Dim[1]), x = xx)
  out
}
sp_t(E2)
sp_t(A2)


## flexible function for preserving certain values in matrix
## accepts as input a sparse matrix and a function that can be evaluated to a logical
sp_keep <- function(x, what = NULL){
  if(!is.function(what)) stop(gettext("'what' must be a function that evaluates to a logical"))
  
  An <- x$Dim[2]
  Ap <- x$p
  Ai <- x$i
  Ax <- x$x
  nz <- 1
  
  for(j in 1:An){
    p <- Ap[j]
    Ap[j] <- nz
    for(k in p:(Ap[j+1]-1)){
      if(what(Ax = Ax[k], Ai = Ai[k], j = j)){
        Ax[nz] <- Ax[k]
        Ai[nz] <- Ai[k]
        nz <- nz + 1
      }
    }
  }
  Ax <- Ax[1:(nz-1)]
  Ai <- Ai[1:(nz-1)]
  Ap[An+1] <- nz
  out <- list(i = Ai, p = Ap, Dim = x$Dim, x = Ax)
  out
}
sp_keep(A2, what = function(Ax, Ai, j){Ax > 5})



## helper function that operates on values in a blank workspace--hence scatter
## values must be drawn back into desired matrix
sp_scatter <- function(A, j, beta, w1, w2, mark, ci, nz){
  Ap <- A$p
  Ai <- A$i
  Ax <- A$x

  for(p in Ap[j]:(Ap[j+1]-1)){
    i <- Ai[p]
    if(w1[i] < mark){
      w1[i] <- mark
      ci[nz] <- i
      w2[i] <- beta * Ax[p]
      nz <- nz + 1
    } else {
      w2[i] <- w2[i] + beta * Ax[p]
    }
  }
  out <- list(w1 = w1, w2 = w2, ci = ci, nz = nz)
}



## matrix multiplication function
## operates over columns of y and values within rows of x
sp_multiply <- function(x, y, sort = FALSE){
  nz <- 1               # initialize number of nonzero entries (assume 1)
  
  m <- x$Dim[1]         # initialize nrow of new matrix 'c' according to nrow(x)
  xp <- x$p             # extract column pointers for 'x'
  xi <- x$i             # extract row indices for 'x'
  xx <- x$x             # extract numeric values from 'x'
  xnz <- length(xx)     # get number of nonzero entries in 'x'
  
  n <- y$Dim[2]         # initialize ncol of new matrix 'c' according to ncol(y)
  yp <- y$p             # same extraction methods as above but for 'y'
  yi <- y$i
  yx <- y$x
  ynz <- length(yx)
  
  w1 <- vector(mode = "numeric", length = m)   # initialize two workspace vectors
  w2 <- vector(mode = "numeric", length = m)
  
  ci <- vector(mode = "numeric")               # initialize row inidicator vector for 'c'
  cp <- vector(mode = "numeric", length = n+1) # initialize column pointer vector for 'c'
  cx <- vector(mode = "numeric")               # initialize numeric value vector for 'c'
  
  for(j in 1:n){                               # outer loop operates over columns of 'c'
    cp[j] <- nz                                # first pointer gets initial value for nonzeros
    for(p in yp[j]:(yp[j+1]-1)){               # second loop operates over rows present in column j of matrix 'y'
      # out <- sp_scatter(x, yi[p], yx[p], w1, w2, j, ci, nz)
      for(k in xp[yi[p]]:(xp[yi[p]+1]-1)){     # third loop operates over rows present in column of 'x' corresponding to rows present in column j of 'y'
        i <- xi[k]                             # row entry i recorded as row present in x
        if(w1[i] < j){                         # tests whether this row has been visited before
          w1[i] <- j                           # if not, update value
          ci[nz] <- i                          #    update row indicator
          w2[i] <- yx[p] * xx[k]               #    update numeric value for c
          nz <- nz + 1                         #    increment nonzero inidcator
        } else {                               # if it has been visited
          w2[i] <- w2[i] + yx[p] * xx[k]       #    add new value to value already present
        }
      }
    }
    for(p in cp[j]:(nz-1)){                    # align x values with their row indicators
      cx[p] <- w2[ci[p]]
    }
  }
  cp[n+1] <- nz                                # assign final column pointer and prepare output
  out <- list(i = ci, p = cp, Dim = c(m, n), x = cx)
  if(sort) out <- sp_sort(out)
  out
}
sp_multiply(A2, B2)


sp_crossprod <- function(x, y = NULL, sort = FALSE){
  if(is.null(y)) y <- x
  x <- sp_t(x)
  
  if(x$Dim[2] != y$Dim[1]) stop(gettext("non-conformable arguments"))
  
  sp_multiply(x, y, sort = sort)
}


sp_tcrossprod <- function(x, y = NULL, sort = FALSE){
  if(is.null(y)) y <- x
  y <- sp_t(y)
  
  if(x$Dim[2] != y$Dim[1]) stop(gettext("non-conformable arguments"))
  
  sp_multiply(x, y, sort = sort)
}


## matrix addition
## related to matrix multiplication, except it does not iterate over rows present in column j of y
sp_add <- function(x, y, alpha = 1, beta = 1, sort = FALSE){
  nz <- 1
  
  m <- x$Dim[1]
  xp <- x$p
  xi <- x$i
  xx <- x$x
  xnz <- length(xx)
  
  n <- y$Dim[2]
  yp <- y$p
  yi <- y$i
  yx <- y$x
  ynz <- length(yx)
  
  w1 <- vector(mode = "numeric", length = m)
  w2 <- vector(mode = "numeric", length = m)
  
  ci <- vector(mode = "numeric")
  cp <- vector(mode = "numeric", length = n+1)
  cx <- vector(mode = "numeric")
  
  for(j in 1:n){
    cp[j] <- nz
    
    for(k in xp[j]:(xp[j+1]-1)){
      i <- xi[k]
      if(w1[i] < j){
        w1[i] <- j
        ci[nz] <- i
        w2[i] <- alpha * xx[k]
        nz <- nz + 1
      } else {
        w2[i] <- w2[i] + alpha * xx[k]
      }
    }
    
    for(k in yp[j]:(yp[j+1]-1)){
      i <- yi[k]
      if(w1[i] < j){
        w1[i] <- j
        ci[nz] <- i
        w2[i] <- beta * yx[k]
        nz <- nz + 1
      } else {
        w2[i] <- w2[i] + beta * yx[k]
      }
    }
    cp[n+1] <- nz
    for(p in cp[j]:(nz-1)){
      cx[p] <- w2[ci[p]]
    }
  }
  
  out <- list(i = ci, p = cp, Dim = c(m, n), x = cx)
  if(sort) out <- sp_sort(out)
  out
}
sp_add(A2, A2)



sp_pvec <- function(x, p){
  n <- length(x)
  if(length(p) != n) stop(gettext("'p' and 'n' must be equal in length"))
  if(any(p > n)) stop(gettext("max(p) must be <= length(x)"))
  b <- vector(mode = "numeric", length = n)
  
  for(k in 1:n) b[k] <- x[p[k]]
  b
}
sp_pvec(x = c(5,6,7), p = c(2,3,1))


sp_ipvec <- function(x, p){
  n <- length(x)
  if(length(p) != n) stop(gettext("'p' and 'n' must be equal in length"))
  if(any(p > n)) stop(gettext("max(p) must be <= length(x)"))
  b <- vector(mode = "numeric", length = n)
  
  for(k in 1:n) b[p[k]] <- x[k]
  b
}
sp_ipvec(x = c(5,6,7), p = c(2,3,1))



sp_pinv <- function(pvec){
  n <- length(pvec)
  pinv <- vector(mode = "numeric", length = n)
  for(i in 1:n) pinv[pvec[i]] <- i
  pinv
}



sp_sort <- function(x){
  Ai <- x$i; Ap <- x$p; Ax <- x$x
  nz <- length(Ax)
  xi <- xx <- vector(mode = "numeric", length = nz)
  w <- vector(mode = "numeric", length = x$Dim[1])
  
  for(p in 1:nz){
    w[Ai[p]] <- w[Ai[p]] + 1
  }
  w <- xp <- sp_cumsum(c(1, w))
  
  for(j in 1:x$Dim[2]){
    for(p in Ap[j]:(Ap[j+1]-1)){
      xi[w[Ai[p]]] <- j
      xx[w[Ai[p]]] <- Ax[p]
      w[Ai[p]] <- w[Ai[p]] + 1
    }
  }
  
  for(j in 1:x$Dim[1]){
    for(p in xp[j]:(xp[j+1]-1)){
      Ai[Ap[xi[p]]] <- j
      Ax[Ap[xi[p]]] <- xx[p]
      Ap[xi[p]] <- Ap[xi[p]] + 1
    }
  }
  
 out <- list(i = Ai, p = c(1, Ap[1:x$Dim[2]]), Dim = x$Dim, x = Ax)
 out
}
sp_sort(AB)
sp_sort(ztz)


sp_permute <- function(x, iiperm = NULL, jperm = NULL, sort = F){
  if(is.null(iiperm)) iiperm <- seq_len(x$Dim[1])
  if(length(iiperm) != x$Dim[1]) stop(gettext("row permutation vector must equal nrow(x)"))
  
  if(is.null(jperm)) jperm <- seq_len(x$Dim[2])
  if(length(jperm) != x$Dim[2]) stop(gettext("column permutation vector must equal ncol(x)"))
  
  nz <- 1
  Dim <- x$Dim
  Ai <- x$i
  Ap <- x$p
  Ax <- x$x
  
  xp <- vector(mode = "numeric", length = Dim[2])
  xi <- vector(mode = "numeric", length = length(Ai))
  xx <- vector(mode = "numeric", length = length(Ax))
  
  for(j in 1:Dim[2]){
    xp[j] <- nz
    k <- jperm[j]
    for(i in Ap[k]:(Ap[k+1]-1)){
      xx[nz] <- Ax[i]
      xi[nz] <- iiperm[Ai[i]]
      nz <- nz + 1
    }
  }
  xp[Dim[2]+1] <- nz
  out <- list(i = xi, p = xp, Dim = Dim, x = xx)
  if(sort) out <- sp_sort(out)
  out
}
sp_permute(B2, iiperm = sp_pinv(c(3, 1, 2, 4, 5, 6, 7, 8, 9, 10)), sort = T)



## checks symmetry of SORTED matrix 'x'
## if not sorted, specify sorted = FALSE for preprocessing
sp_isSymmetric <- function(x, sorted = TRUE){
  if(!sorted) x <- sp_sort(x)
  Dim <- x$Dim
  if(Dim[1] != Dim[2]) return(FALSE)

  Ai <- x$i; Ap <- w <- x$p; Ax <- x$x
 
  for(j in 1:Dim[2]){
    for(p in Ap[j]:(Ap[j+1]-1)){
      if(Ax[p] != Ax[w[Ai[p]]]){
        return(FALSE)
      }
      w[Ai[p]] <- w[Ai[p]] + 1
    }
  }
  TRUE
}
sp_isSymmetric(ztz, sorted = F)



sp_crossprod <- function(x, y = NULL, sort = FALSE){
  if(is.null(y)) y <- x
  nz <- 1               # initialize number of nonzero entries (assume 1)
  
  m <- x$Dim[1]         # initialize nrow of new matrix 'c' according to nrow(x)
  xp <- x$p             # extract column pointers for 'x'
  xi <- x$i             # extract row indices for 'x'
  xx <- x$x             # extract numeric values from 'x'
  xnz <- length(xx)     # get number of nonzero entries in 'x'
  
  n <- y$Dim[2]         # initialize ncol of new matrix 'c' according to ncol(y)
  yp <- y$p             # same extraction methods as above but for 'y'
  yi <- y$i
  yx <- y$x
  ynz <- length(yx)
  
  w1 <- vector(mode = "numeric", length = m)   # initialize two workspace vectors
  w2 <- vector(mode = "numeric", length = m)
  
  ci <- vector(mode = "numeric")               # initialize row inidicator vector for 'c'
  cp <- vector(mode = "numeric", length = n+1) # initialize column pointer vector for 'c'
  cx <- vector(mode = "numeric")               # initialize numeric value vector for 'c'
  
  # for(j in 1:n){                               # outer loop operates over columns of 'c'
  #   cp[j] <- nz                                # first pointer gets initial value for nonzeros
    # for(p in yp[j]:(yp[j+1]-1)){               # second loop operates over rows present in column j of matrix 'y'
      # for(k in xp[yi[p]]:(xp[yi[p]+1]-1)){     # third loop operates over rows present in column of 'x' corresponding to rows present in column j of 'y'
        # i <- xi[k]                             # row entry i recorded as row present in x
        # if(w1[i] < j){                         # tests whether this row has been visited before
        #   w1[i] <- j                           # if not, update value
        #   ci[nz] <- i                          #    update row indicator
        #   w2[i] <- yx[p] * xx[k]               #    update numeric value for c
        #   nz <- nz + 1                         #    increment nonzero inidcator
        # } else {                               # if it has been visited
        #   w2[i] <- w2[i] + yx[p] * xx[k]       #    add new value to value already present
        # }
      # }
    # }
    # for(p in cp[j]:(nz-1)){                    # align x values with their row indicators
    #   cx[p] <- w2[ci[p]]
    # }
  # }
  
  for(j in 1:x$Dim[2]){
    for(p in xp[j]:(xp[j+1]-1)){
      w1[xi[p]] <- w1[xi[p]] + xx[p] * yx[xi[p]]
      print(w1)
    }
    # w[xi[p]] <- w[xi[p]]+1
  }
  # cp[n+1] <- nz                                # assign final column pointer and prepare output
  # out <- list(i = ci, p = cp, Dim = c(m, n), x = cx)
  # if(sort) out <- sp_sort(out)
  # out
  # list(ci, cx)
}
sp_crossprod(Z2)
sp_crossprod(sp_t(Z2), Z2)
sp_crossprod(B2, A2)

# last completed lecture 10



