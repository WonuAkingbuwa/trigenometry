#' Interpolate the (non)linear shape of a relation between two traits based on genetic correlations
#'
#' Based on genetic correlations between trait y (can be binary or continuous trait) and a number of GWASs of trait x, where x is binned  into segments of the distributions of x, and a pair of bins are used as case.control in a GWAS
#'
#' @param rg genetic correlation genetic correlation between outcome y and a pair of bins of x
#' @param b1 value assigned to bin b1 used in the GWAS of "x" usually a count, or median/mean phenotype for participants in the bin, used to compute dx
#' @param b2 value assigned to bin b2 used in the GWAS of "x" , usually a count, or median/mean phenotype for participants in the bin, used to compute dx
#' @param method either "polynomial' or "cubic-spline" specifies the method used for interpolation
#'
#' @param q1 value (in terms of the scale of b1/b2) at which to place the first of two knows for the method 'cubic-spline', ignored if method is 'polynomial', automatically initiated in left at the default value NULL
#' @param q1 value (in terms of the scale of b1/b2) at which to place the second of two knows for the method 'cubic-spline', ignored if method is 'polynomial', automatically initiated in left at the default value NULL
#' @param boot logical used to indicate whether the genetic correlations should be re-sampled 200 times (based on point estimate and s.e.) to estimate parameter uncertainty (standard errors) and introduce curves to a plot to visualize uncertainty  (optional)
#'
#'
#' @return a list with any data pruning undertaken, parameters, optimizer status and the method used for interpolation
#' @return (optionally) a curve reflection the estimated non-linear genetic relationship.
cor2curve <- function(rg,b1,b2,se,method = "polynomial",q1=NULL,q2=NULL,boot=FALSE,se.filter=.30,x.trait="",y.trait=""){
  # Arguments:
  # rg: genetic correlation genetic correlation between this GWAS of a pair of bins of x, and the outcome y.
  # b1:value assigned to bin b1 used in the GWAS of "x" usually a count, or median/mean phenotype for participants in the bin, used to compute dx
  # b2: value assigned to bin b2 used in the GWAS of "x" , usually a count, or median/mean phenotype for participants in the bin, used to compute dx
  # method, either "polynomial' or "cubic-spline"

  # q1 and q2 are cut points for the cubic spline, default to .4 and .6 of sorted c(b1,b2)
  # plot, do we want a plot of the curve?

  # boot, do you need bootstrapped alternate curves?
  sev <- NULL

  # dy (test argument!)
  dx <- abs(b1 - b2)
  int <- mean(c(max(b1),min(b1)))

  # check and omit NA's
  dat <- cbind(rg,b1,b2,se)
  dat2 <- na.omit(dat)

  # se filter
  dat3 <-  dat2[dat2[,4] < se.filter,]  

  rg <- dat3[,1]
  b1 <- dat3[,2]
  b2 <- dat3[,3]
  se <- dat3[,4]

  mess <- "No missing data"
  if( nrow(dat2) < nrow(dat)){
    print(paste0("Had to omit ",nrow(dat) - nrow(dat2)," rows with missing values in rg, b1, b2, or se" ))
    mess <- paste0("Had to omit ",nrow(dat) - nrow(dat2)," rows with missing values in rg, b1, b2, or se" )

  }
 mess <- "no se's pout of bounds"
  if( nrow(dat2) < nrow(dat)){
    print(paste0("Had to omit ",nrow(dat2) - nrow(dat3)," rows with se's out of (user specified) bounds" ))
    mess <- paste0("Had to omit ",nrow(dat2) - nrow(dat3)," rows with se's out of (user specified) bound" )

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


  ### Filter for the bootstrap outliers:
  IQR.outliers <- function(x) {
  if(any(is.na(x)))
    stop("x is missing values")
  if(!is.numeric(x))
    stop("x is not numeric")
  Q3<-quantile(x,0.75)
  Q1<-quantile(x,0.25)
  IQR<-(Q3-Q1)
  left<- (Q1-(4*IQR))
  right<- (Q3+(4*IQR))
  out <- x[x >left  & x < right]
  if(length(out) < length(x)){
    print(paste('Had to omit',length(x) - length(out),' bootstrap itterations where the optimizer value was far outside the norm!'))
  }
   c(left,right)
  }




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


    sep <- rep(NA,5)

    if(boot==T){

      sev <- matrix(NA,200,5)
      value <- matrix(NA,200,1)
      for(i in 1:200){

        rgi <- rnorm(length(rg),rg,se)

        rgi[rgi > .9999] <-  .9999
        rgi[rgi < -.9999] <-  -.9999
        
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
        value[i] <- sefit$value
        
      }
      lr <-  IQR.outliers(value)
      sev <-  sev[value > lr[1] & value < lr[2],]
      sep <- apply(sev,2,sd)

    }



    out <-  list(method=method,
                 messages = c(mess,mess2),
                 observations=length(rg),
                 coefficients = cbind.data.frame(term = c("5th","4th","3th","2nd","1st"),
                                                 parameter = fit$par,
                                                 se = sep,
                                                 z = fit$par/sep),
                 convergence= fit$convergence,
                 value =fit$value,
                 boot=boot,
                 b1 = b1,
                 b2 = b2,
                 par=fit$par,
                 bpar=sev,
                 x.trait=x.trait,
                 y.trait=y.trait)

    class(out) <- c("list", "cor2curve")
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
    fit <- optim(par = rep(0,12),fn = spline_fn,method = "BFGS",control=list(maxit=500))




    sep <- rep(NA,12)
    
    if(boot==T){

      sev <- matrix(NA,200,12)
      value <- matrix(NA,200,1)
      for(i in 1:200){

        rgi <- rnorm(length(rg),rg,se)

        rgi[rgi > .9999] <-  .9999
        rgi[rgi < -.9999] <-  -.9999

        # the correlation, divided by distance, as an angle
        angle <- 90 - acos(rgi) * 180/pi

        # distance on y given know distance on x and the angle
        dy <- tan(angle/ (180/pi))


        # this is where we really minimize...
        sefit <- optim(par = rep(0,12),fn = spline_fn,method = "BFGS",control=list(maxit=500))
        sev[i,] <- sefit$par
        value[i] <- sefit$value
        

      }


      lr <-  IQR.outliers(value)
      sev <-  sev[value > lr[1] & value < lr[2],]
      sep <- apply(sev,2,sd)

    }





    out <-  list(method=method,
                 messages = c(mess,mess2),
                 observations=length(rg),
                 cubic1 = cbind.data.frame(term = c("3th","2nd","1st","int"),
                                           parameters = fit$par[1:4],
                                           se = sep[1:4]),
                 cubic2 = cbind.data.frame(term = c("3th","2nd","1st","int"),
                                           parameters = fit$par[5:8],
                                           se = sep[5:8]),
                 cubic3 = cbind.data.frame(term = c("3th","2nd","1st","int"),
                                           parameters = fit$par[9:12],
                                           se = sep[9:12]),

                 cuts = c(q1,q2),
                 convergence= fit$convergence,
                 boot=boot,
                 value =fit$value,
                 int=int,
                 par=fit$par,
                 b1=b1,
                 b2=b2,
                 boot=boot,
                 x.trait=x.trait,
                 y.trait=y.trait,
                 bpar=sev)

    class(out) <- c("list", "cor2curve")
    return(out)

  }
}



#' Plot a curve based on cor2curve
#'
#' This function plot a curve based on the input parameters.
#'
#' @param x an object of class 'cor2curve' generated with the cor2curve() function.
#' @param xlim numeric vector, optional. The x-axis limits for the plot, defaults to sensible values.
#' @param ylim numeric vector, optional. The y-axis limits for the plot, defaults to sensible values.
#' @param main character string, optional. The main title for the plot.
#' @param sub character string, optional. The subtitle for the plot.
#' @param xlab character string, optional. The x-axis label for the plot, as a default its inherited from cor2curve object.
#' @param ylab character string, optional. The y-axis label for the plot, as a default its inherited from cor2curve object.
#' @param boot logical, optional. Whether to include bootstrapped lines in the plot, as a default its inherited from cor2curve object.
#'
#' @return A plot of the curve.
#'
plot.cor2curve <- function(x,xlim=NULL,ylim=NULL,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,col.curve="black", col.boot="grey",boot=NULL,add=FALSE, ... ){

  if(is.null(xlab)){xlab <- x$x.trait}
  if(is.null(ylab)){ylab <- x$y.trait}
  if(is.null(boot)){boot <- x$boot}


  if(x$method=="polynomial"){

    par <-x$par
    bpar <- x$bpar
    min <- min(c(x$b1,x$b2)) 
    max <- max(c(x$b1,x$b2)) 

    if(is.null(xlim)){ xlim <- c(min,max)}

    beta <- unique(x$b1)
    # "draw rthe values on y given x:
    ys <- par[1]*beta^5 + par[2]*beta^4 + par[3]*beta^3 + par[4]*beta^2 + par[5]*beta
    # build an intercept
    int <- mean(ys)

    # pick smart ylims...
    if(is.null(ylim)){
      ylim <- c(min(ys) - int,max(ys) - int)

      ylim[1] <- ylim[1] - .4*(abs(ylim[1] - ylim[2]))
      ylim[2] <- ylim[2] + .4*(abs(ylim[1] - ylim[2]))
    }
    if(add==FALSE){
    plot(NULL,xlim=xlim,ylim=ylim,main=main,sub=sub,ylab=ylab,xlab=xlab, ... )
    }
    
    #if bootstrapped lines are to be drawn:
    if(x$boot){
      for(i in 1:nrow(bpar)){
        ys <- bpar[i,1]*beta ^5 + bpar[i,2]*beta ^4 + bpar[i,3]*beta ^3 + bpar[i,4]*beta ^2 + bpar[i,5]*beta

        # build an intercept
        int <- mean(ys)

        curve(bpar[i,1]*x^5 + bpar[i,2]*x^4 + bpar[i,3]*x^3 + bpar[i,4]*x^2 + bpar[i,5]*x - int,from = min, to = max,add=T,col=col.boot,lty="dashed",...)
      }
    }

    # reset ys and int for main plot:
    ys <- par[1]*beta^5 + par[2]*beta^4 + par[3]*beta^3 + par[4]*beta^2 + par[5]*beta

    # build an intercept
    int <- mean(ys)

    curve(par[1]*x^5 + par[2]*x^4 + par[3]*x^3 + par[4]*x^2 + par[5]*x  - int,from = xlim[1], to = xlim[2],add=T, col = col.curve)


  }



  if(x$method =="cubic-spline"){

    bpar <- x$bpar
    par <- x$par
    q1 <- x$cuts[1]
    q2 <- x$cuts[2]

    min <- min(c(x$b1,x$b2)) 
    max <- max(c(x$b1,x$b2)) 

    if(is.null(xlim)){ xlim <- c(min,max)}

    beta <- unique(x$b1)

    ys <- ifelse(beta < q1,x$par[1]*beta^3 + x$par[2] * beta^2 + x$par[3] * beta + x$par[4],
                 ifelse(beta >= q2,x$par[9]*beta^3 + x$par[10] * beta^2 + x$par[11] * beta + x$par[12],x$par[5]*beta^3 + x$par[6] * beta^2 + x$par[7] * beta + x$par[8]))

    int <- mean(ys)


    if(is.null(ylim)){
      ylim <- c(min(ys) - int,max(ys) - int)

      ylim[1] <- ylim[1] - .4*(abs(ylim[1] - ylim[2]))
      ylim[2] <- ylim[2] + .4*(abs(ylim[1] - ylim[2]))
    }

    if(add==FALSE){
    plot(NULL,xlim=xlim,ylim=ylim,main=main,sub=sub,ylab=ylab,xlab=xlab, ... )
    }
  
    ### Draw bootstrapped lines
    if(x$boot){
      for(i in 1:nrow(bpar)){
        seys <- ifelse(beta < q1,bpar[i,1]*beta^3 + bpar[i,2] * beta^2 + bpar[i,3] * beta + bpar[i,4],
                       ifelse(beta >= q2,bpar[i,9]*beta^3 + bpar[i,10] * beta^2 + bpar[i,11] * beta + bpar[i,12],x$bpar[i,5]*beta^3 + bpar[i,6] * beta^2 + bpar[i,7] * beta + bpar[i,8]))

        seint <- mean(seys)


        curve(bpar[i,1]*x^3 + bpar[i,2] * x^2 + bpar[i,3] * x + bpar[i,4]  - seint,from =min,to=q1,add=T,col=col.boot,lty="dashed")

        curve(bpar[i,5]*x^3 + bpar[i,6] * x^2 + bpar[i,7] * x + bpar[i,8] - seint,from = q1,to=q2,add=T,col=col.boot,lty="dashed")

        curve(bpar[i,9]*x^3 + bpar[i,10] * x^2 + bpar[i,11] * x  + bpar[i,12] - seint,from = q2,to=max,add=T,col=col.boot,lty="dashed")

      }
    }

    # Plot the curve...

    curve(par[1]*x^3 + par[2] * x^2 + par[3] * x + par[4]  - int ,from =min,to=q1,add=T,col=col.curve,...)

    curve(par[5]*x^3 + par[6] * x^2 + par[7] * x + par[8] - int ,from = q1,to=q2,add=T,col=col.curve,...)

    curve(par[9]*x^3 + par[10] * x^2 + par[11] * x + par[12] - int ,from = q2,to=max,add=T,col=col.curve,...)

  }


}


#' Print a cor2curve object
#'
#' @param x A cor2curve object
#' @param ... Additional arguments (not used)
#' @param quote A logical indicating whether to quote the output (default: FALSE)
#'
#' @description
#' This function prints a summary of a cor2curve object, which is a custom object
#' containing the results of a cor2curve analysis. The output includes the method
#' used, messages, number of observations, coefficients or cubic spline parameters,
#' cuts, convergence information, value, and intercept.
#'
#' @details
#' The function prints different sets of information depending on the method used
#' to fit the cor2curve. If the method is "polynomial", it prints the coefficients
#' of the polynomial. If the method is "cubic-spline", it prints the parameters
#' of the cubic spline.
#'
#' @return
#' A list containing the summary information
#' @examples
#' # Create a cor2curve object (not shown)
#' print.cor2curve(x)
print.cor2curve = function(x){

  if(x$method == "polynomial"){

    print(list(method=x$method,
               messages = x$message,
               observations=x$observaions,
               coefficients = x$coefficients,
               cuts = x$cuts,
               convergence= x$convergence,
               value =x$value,
               int=x$int))
  }

  if(x$method == "cubic-spline"){

    print(list(method=x$method,
               messages = x$message,
               observations=x$observaions,
               cubic1 = x$cubic1,
               cubic2 = x$cubic2,
               cubic3 = x$cubic3,
               cuts = x$cuts,
               convergence= x$convergence,
               value =x$value,
               int=x$int))
  }

}



