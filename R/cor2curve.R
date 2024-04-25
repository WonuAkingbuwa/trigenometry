#' Interpolate the (non)linear shape of a relation between two traits based on genetic correlations
#'
#' Based on genetic correlations between trait y (can be binary or continuous trait) and a number of GWASs of trait x, where x is binned  into segments of the distributions of x, and apira of bins are used as case.control in a GWAS
#'
#' @param rg genetic correlation genetic correlation between outcome y and a pair of bins of x
#' @param b1 value assigned to bin b1 used in the GWAS of "x" usually a count, or median/mean phenotype for participants in the bin, used to compute dx
#' @param b2 value assigned to bin b2 used in the GWAS of "x" , usually a count, or median/mean phenotype for participants in the bin, used to compute dx
#' @param method either "polynomial' or "cubic-spline" specifies the method used for interpolation
#'
#' @return plot of the interpolated curves (optional)
#' @return plot of re-sampled version of  interpolated curves to indicate statistical uncertainty (optional)
#' @return a list with any data pruning undertaken, parameters, optimizer status and the method used for interpolation
cor2curve <- function(rg,b1,b2,se,method = "polynomial",q1=NULL,q2=NULL,plot=TRUE,boot=FALSE,x.trait="",y.trait=""){
  # Arguments:
  # rg: genetic correlation genetic correlation between this GWAS of a pair of bins of x, and the outcome y.
  # b1:value assigned to bin b1 used in the GWAS of "x" usually a count, or median/mean phenotype for participants in the bin, used to compute dx
  # b2: value assigned to bin b2 used in the GWAS of "x" , usually a count, or median/mean phenotype for participants in the bin, used to compute dx
  # method, either "polynomial' or "cubic-spline"

  # q1 and q2 are cut points for the cubic spline, default to .4 and .6 of sorted c(b1,b2)
  # plot, do we want a plot of the curve?

  # boot, do you need bootstrapped alternate curves?


  # dy (test argument!)
  dx <- abs(b1 - b2)
  int <- mean(c(max(b1),min(b1)))

  # check and omit NA's
  dat <- cbind(rg,b1,b2,se)
  dat2 <- na.omit(dat)

  rg <- dat2[,1]
  b1 <- dat2[,2]
  b2 <- dat2[,3]
  se <- dat2[,4]

  mess <- "No missing data"
  if( nrow(dat2) < nrow(dat)){
    print(paste0("Had to omit ",nrow(dat) - nrow(dat2)," rows with missing values in rg, b1, b2, or se" ))
    mess <- paste0("Had to omit ",nrow(dat) - nrow(dat2)," rows with missing values in rg, b1, b2, or se" )

  }


  # warn if rg our of bounds
  mess2 <- "All genetic correlations between -1 and 1"
  if(sum(abs(rg)  > 1) > 0){
    print(paste0("Had to omit ",sum(abs(rg)  > 1)," rows with rg our of bounds (-1 to 1)" ))
    mess2 <- paste0("Had to omit ",sum(abs(rg)  > 1)," rows with rg our of bounds (-1 to 1)" )
  }
  # remove out of bounds
  b1 <- b1[abs(rg) <= 1]
  b2 <- b2[abs(rg) <= 1]
  se <- se[abs(rg) <= 1]
  rg <- rg[abs(rg) <= 1]

  #if(sum(dx==0) > 0){stop("Some distances between bin values b1 and b2 are 0, which shouldn't be posible")}

  # the correlation, divided by distance, as an angle
  angle <- 90 - acos(rg) * 180/pi

  # distance on y given know distance on x and the angle
  dy <- tan(angle/ (180/pi))



  if(method == "polynomial"){

    solve_fn <-  function(par){
      a <- par[1]
      b <- par[2]
      c <- par[3]
      d <- par[4]
      e <- par[5]


      # real distance, minus implied distance
      ak <- (dy - ((a*b1^5 + b*b1^4 + c*b1^3 +d*b1^2 + e*b1^1 )  - (a*b2^5 + b*b2^4 + c*b2^3 +d*b2^2 + e*b2^1 ) )) * (1/(se^2))

      # minimize that distance squared
      sum(ak^2)/ sum(1/(se^2))

    }

    # this is where we really minimize...
    fit <- optim(par = c(0,0,0,0,0),fn = solve_fn,method = "BFGS")

    out <- fit # update


    if(plot==TRUE){
      min <- min(c(b1,b2))
      max <- max(c(b1,b2))

      beta <- unique(b1)

      ys <- fit$par[1]*beta ^5 + fit$par[2]*beta ^4 + fit$par[3]*beta ^3 + fit$par[4]*beta ^2 + fit$par[5]*beta

      # build an intercept
      int <- mean(ys)

      # automatically pick sensible y limits
      ylim <-  c(min(ys) - int,max(ys) - int)

      # Margins aroudn the y-limits
      ylim[1] <- (ylim[1] - .2*(abs(ylim[1] - ylim[2])))
      ylim[2] <- (ylim[2] + .2*(abs(ylim[1] - ylim[2])))

      curve(fit$par[1]*x^5 + fit$par[2]*x^4 + fit$par[3]*x^3 + fit$par[4]*x^2 + fit$par[5]*x  - int,from = min, to = max,ylim=ylim,ylab=y.trait,xlab=x.trait)

    }


    sep <- rep(NA,5)

    if(boot==T){

      sev <- matrix(NA,200,5)

      for(i in 1:200){

        rgi <- rnorm(length(rg),rg,se)

        rgi[abs(rgi) > .9999] <- .9999

        # the correlation, divided by distance, as an angle
        angle <- 90 - acos(rgi) * 180/pi

        # distance on y given know distance on x and the angle
        dy <- tan(angle/ (180/pi))



        if(method == "polynomial"){
          solve_fn <-  function(par){
            a <- par[1]
            b <- par[2]
            c <- par[3]
            d <- par[4]
            e <- par[5]


            # real distance, minus implied distance
            ak <- (dy - ((a*b1^5 + b*b1^4 + c*b1^3 +d*b1^2 + e*b1^1 )  - (a*b2^5 + b*b2^4 + c*b2^3 +d*b2^2 + e*b2^1 ) )) * (1/(se^2))

            # minimize that distance squared
            sum(ak^2)

          }

          # this is where we really minimize...
          sefit <- optim(par = c(0,0,0,0,0),fn = solve_fn,method = "BFGS")
        }
        sev[i,] <- sefit$par
        if(plot==TRUE){

          seys <- sefit$par[1]*beta^5 + sefit$par[2]*beta^4 + sefit$par[3]*beta^3 + sefit$par[4]*beta^2 + sefit$par[5]*beta

          seint <- mean(seys)

          curve(sefit$par[1]*x^5 + sefit$par[2]*x^4 + sefit$par[3]*x^3 + sefit$par[4]*x^2 + sefit$par[5]*x - seint,add=T,col="grey",lty="dashed")

        }

      }

      sep <- apply(sev,2,sd)

    }

    if(plot==T){
      curve(fit$par[1]*x^5 + fit$par[2]*x^4 + fit$par[3]*x^3 + fit$par[4]*x^2 + fit$par[5]*x - int,from = min, to = max,ylim=ylim,col="red",lwd=2,add=T)
      abline(v=beta,lty="dashed",col="lightblue",lwd=.5)
    }

    out <-  list(method=method,
                 messages = c(mess,mess2),
                 observations=length(rg),
                 coefficients = cbind.data.frame(term = c("5th","4th","3rd","2nd","1st"),
                                                 parameter = fit$par,
                                                 se = sep,
                                                 z = fit$par/sep),
                 convergence= fit$convergence,
                 value =fit$value)
    return(out)
  }

  if(is.null(q1)){
    lq1 <- round(.4*length(c(b1,b2)))
    q1 <- sort(c(b1,b2))[lq1]

    lq2 <- round(.6*length(c(b1,b2)))
    q2 <- sort(c(b1,b2))[lq2]


  }

  if(method == "cubic-spline"){

    spline_fn <- function(par){

      ## NOTATION GUIDANCE (SORRY SORRY FUTURE WONU/MICHEL)
      # a_*, b_* and c_* where * i 1-3 are parameters of the cubic spline functions
      # b1 and b2 are the bin mid values
      # q1 and q2 are the spline cut points (default to avlues within b1/b2)


      a_1 <- par[1]
      b_1 <- par[2]
      c_1 <- par[3]
      d_1 <- par[4]

      a_2 <- par[5]
      b_2 <- par[6]
      c_2 <- par[7]
      d_2 <- par[8]

      a_3 <- par[9]
      b_3 <- par[10]
      c_3 <- par[11]
      d_3 <- par[12]




      constr1 <- (a_1*q1^3 + b_1*q1^2 + c_1*q1^1 + d_1) - (a_2*q1^3 + b_2*q1^2 + c_2*q1^1 + d_2)
      constr1d <- (3*a_1*q1^2 + 2*b_1*q1 + c_1) - (3*a_2*q1^2 + 2*b_2*q1 + c_2)

      constr2 <- (a_2*q2^3 + b_2*q2^2 + c_2*q2^1 + d_2) - (a_3*q2^3 + b_3*q2^2 + c_3*q2^1 + d_3)
      constr2d <- (3*a_2*q2^2 + 2*b_2*q2 + c_2) - (3*a_3*q2^2 + 2*b_3*q2 + c_3)

      intr <- matrix(NA,nrow=length(b1),ncol=2)

      intr[b2  <= q1,1] <- (a_1*b2[b2 <=q1]^3 + b_1*b2[b2 <=q1]^2 + c_1*b2[b2 <=q1]^1 + d_1)
      intr[b2 > q1 & b2  <= q2 ,1] <- (a_2*b2[b2  > q1 & b2  <= q2]^3 + b_2*b2[b2  > q1 & b2  <= q2]^2 + c_2*b2[b2  > q1 & b2  <= q2 ]^1 + d_2)
      intr[ b2  > q2 ,1] <- (a_3*b2[b2  > q2]^3 + b_3*b2[b2  > q2]^2 + c_3*b2[b2  > q2]^1 + d_3)

      intr[b1  <=q1,2] <- (a_1*b1[b1  <=q1]^3 + b_1*b1[b1  <=q1]^2 + c_1*b1[b1  <=q1]^1 + d_1)
      intr[b1  > q1 & b1  <= q2 ,2] <- (a_2*b1[b1  > q1 & b1  <= q2]^3 + b_2*b1[b1  > q1 & b1  <= q2]^2 + c_2*b1[b1  > q1 & b1  <= q2]^1 + d_2)
      intr[ b1  > q2 ,2] <- (a_3*b1[b1  > q2 ]^3 + b_3*b1[b1  > q2]^2 + c_3*b1[b1  > q2]^1 + d_3)

      ak <- mean((dy - (intr[,2] - intr[,1]))^2 * 1/((se^2)))
      ak + 1000*(constr1^2 + constr1d^2 + constr2^2 + constr2d^2)
    }


    # this is where we really minimize...
    fit <- optim(par = rep(0,12),fn = spline_fn,method = "BFGS")


    if(plot==TRUE){
      min <- min(c(b1,b2))
      max <- max(c(b1,b2))



      beta <- unique(b1)

      ys <- ifelse(beta < q1,fit$par[1]*beta^3 + fit$par[2] * beta^2 + fit$par[3] * beta + fit$par[4],
                   ifelse(beta >= q2,fit$par[9]*beta^3 + fit$par[10] * beta^2 + fit$par[11] * beta + fit$par[12],fit$par[5]*beta^3 + fit$par[6] * beta^2 + fit$par[7] * beta + fit$par[8]))

      int <- mean(ys)



      ylim <- c(min(ys) - int,max(ys) - int)

      ylim[1] <- ylim[1] - .2*(abs(ylim[1] - ylim[2]))
      ylim[2] <- ylim[2] + .2*(abs(ylim[1] - ylim[2]))


      plot(NULL,ylim=ylim,xlim=c(min,max),ylab=y.trait,xlab=x.trait)
      curve(fit$par[1]*x^3 + fit$par[2] * x^2 + fit$par[3] * x + fit$par[4]  - int ,from =min,to=q1,add=T)

      curve(fit$par[5]*x^3 + fit$par[6] * x^2 + fit$par[7] * x + fit$par[8] - int ,from = q1,to=q2,add=T)

      curve(fit$par[9]*x^3 + fit$par[10] * x^2 + fit$par[11] * x + fit$par[12] - int ,from = q2,to=max,add=T)

    }

    sep <- rep(NA,12)

    if(boot==T){

      sev <- matrix(NA,100,12)

      for(i in 1:100){

        rgi <- rnorm(length(rg),rg,se)

        rgi[abs(rgi) > .9999] <- .9999

        # the correlation, divided by distance, as an angle
        angle <- 90 - acos(rgi) * 180/pi

        # distance on y given know distance on x and the angle
        dy <- tan(angle/ (180/pi))


        # this is where we really minimize...
        sefit <- optim(par = rep(0,12),fn = spline_fn,method = "BFGS")

        sev[i,] <- sefit$par


        if(plot==TRUE){

          seys <- ifelse(beta < q1,sefit$par[1]*beta^3 + sefit$par[2] * beta^2 + sefit$par[3] * beta + sefit$par[4],
                         ifelse(beta >= q2,sefit$par[9]*beta^3 + sefit$par[10] * beta^2 + sefit$par[11] * beta + sefit$par[12],sefit$par[5]*beta^3 + sefit$par[6] * beta^2 + sefit$par[7] * beta + sefit$par[8]))

          seint <- mean(seys)


          curve(sefit$par[1]*x^3 + sefit$par[2] * x^2 + sefit$par[3] * x + sefit$par[4]  - seint,from =min,to=q1,add=T,col="grey",lty="dashed")

          curve(sefit$par[5]*x^3 + sefit$par[6] * x^2 + sefit$par[7] * x + sefit$par[8] - seint,from = q1,to=q2,add=T,col="grey",lty="dashed")

          curve(sefit$par[9]*x^3 + sefit$par[10] * x^2 + sefit$par[11] * x  + sefit$par[12] - seint,from = q2,to=max,add=T,col="grey",lty="dashed")

        }
      }



      sep <- apply(sev,2,sd)

    }


    if(plot==TRUE){

      curve(fit$par[1]*x^3 + fit$par[2] * x^2 + fit$par[3] * x + fit$par[4] - int,from =min,to=q1,add=T,col="red",lwd=2)

      curve(fit$par[5]*x^3 + fit$par[6] * x^2 + fit$par[7] * x + fit$par[8]  - int,from = q1,to=q2,add=T,col="red",lwd=2)

      curve(fit$par[9]*x^3 + fit$par[10] * x^2 + fit$par[11] * x + fit$par[12] - int,from = q2,to=max,add=T,col="red",lwd=2)
      abline(v=beta,lty="dashed",col="lightblue",lwd=.5)
    }



    out <-  list(method=method,
                 messages = c(mess,mess2),
                 observations=length(rg),
                 cubic1 = cbind.data.frame(term = c("3rd","2nd","1st","int"),
                                           parameters = fit$par[1:4],
                                           se = sep[1:4]),
                 cubic2 = cbind.data.frame(term = c("3rd","2nd","1st","int"),
                                           parameters = fit$par[5:8],
                                           se = sep[5:8]),
                 cubic3 = cbind.data.frame(term = c("3rd","2nd","1st","int"),
                                           parameters = fit$par[9:12],
                                           se = sep[9:12]),

                 cuts = c(q1,q2),
                 convergence= fit$convergence,
                 value =fit$value,
                 int=int)


    return(out)

  }



}
