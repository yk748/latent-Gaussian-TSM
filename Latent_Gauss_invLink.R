#################################################################################
# Main function for computing inverse link functions
#################################################################################
Latent_Gauss_invLink <- function(d,p,Cov_X,Link){
  
  # ------------------------------------------------- #
  # Compute inverse link functions from the estimated link functions
  Cov_Z <- array(NA,dim(Cov_X))
  
  # Grid setting:
  grid_u <- c(seq(-1,-0.72,0.01),seq(-0.7,0,0.2),seq(0.1,0.7,0.2),seq(0.72,1,0.01))
  
  for (i in 1:d){
    for (j in 1:d){
      Cov_Z[i,j,] <- interpolation(Link[,i,j],grid_u,Cov_X[i,j,])
    }
  }
  
  return(Cov_Z)
}

#################################################################################
# Main function for interpolation: Find the interval where the values are belonging to
#################################################################################
interpolation <- function(coef,u,v){
  
  # ------------------------------------------------- #
  # Result from natural spline
  Spline_result <- nat_spline(coef,u);
  c <- Spline_result[[1]];
  d <- Spline_result[[2]];
  h <- Spline_result[[3]];
  
  # ------------------------------------------------- #
  # Compute the knots
  knot <- vector("numeric",length=length(u))
  for (i in 1:length(u)){
    pow <- 1:length(coef)
    knot[i] <- coef[1:length(coef)]%*%(u[i]^pow)[1:length(coef)]
  }
  n <- length(knot);
  
  # ------------------------------------------------- #
  # main loop
  L_inv_v <- vector("numeric",length(v));
  for (i in 1:length(v)){
    if ( v[i] < knot[1] |  v[i] > knot[n]){# Cutoff L(u) < -1 or L(u) > 1:
      
      if ( v[i] < knot[1] ){# Cutoff L(u) < -1
        L_inv_v[i] <- -1;  
        next;
      }
      else {# Cutoff L(u) > 1
        L_inv_v[i] <- 1; 
        next;
      }
  
    }else{
      idx <- 1;
      while(v[i] > knot[idx]){
        idx <- idx+1;
      }
      
      idx <- idx-1;
      First_term <- d[(idx+1)]*(v[i]-knot[idx])^3/(6*h[(idx+1)])
      Second_term <- d[idx]*(knot[(idx+1)]-v[i])^3/(6*h[(idx+1)])
      Third_term <- c[1,idx]*(v[i]-knot[idx])
      Fourth_term <- c[2,idx]*(knot[(idx+1)]-v[i])
      
      L_inv_v[i] <- First_term + Second_term + Third_term + Fourth_term
      next;
    }
  } 
  return(L_inv_v)
}

#################################################################################
# Function for cubic spline
#################################################################################
nat_spline <- function(coef,u){
  
  # ------------------------------------------------- #
  # compute the knots from the given polynomials
  v <- vector("numeric",length=length(u))
  for (i in 1:length(u)){
    pow <- 1:length(coef)
    v[i] <- coef[1:length(coef)]%*%(u[i]^pow)[1:length(coef)]
  }
  
  # ------------------------------------------------- #
  # set up
  L_inv_v <- u;
  n <- length(v);
  
  h <- vector("numeric",length=n);
  LHS <- matrix(0,nrow=n,ncol=n);
  RHS <- vector("numeric",length=n);
  
  h[1] <- 0; h[n] <- v[n]-v[(n-1)];
  LHS[1,1] <- 1; LHS[n,n] <- 1;
  RHS[1] <- 0; RHS[n] <- 0;
  
  # ------------------------------------------------- #
  # main loop
  for (i in 2:(n-1)){
    h[i] <- v[i] - v[(i-1)];
    LHS[i,(i-1)] <- (v[i] - v[(i-1)])/6;
    LHS[i,i] <- ( (v[i] - v[(i-1)]) + (v[(i+1)] - v[i]) )/3;
    LHS[i,(i+1)] <- (v[(i+1)] - v[i])/6;
    RHS[i] <- ( L_inv_v[(i+1)]-L_inv_v[i] ) / (v[(i+1)] - v[i]) - ( L_inv_v[i]-L_inv_v[(i-1)] ) / (v[i] - v[(i-1)]);
  }
  
  # ------------------------------------------------- #
  # compute the piecewise ploynomial
  c <- matrix(0,nrow=2,ncol=(n-1));
  d <- solve(LHS)%*%RHS;
  
  for (i in 1:(n-1)){
    c[1,i] <- L_inv_v[(i+1)]/h[(i+1)]-d[(i+1)]*h[(i+1)]/6;
    c[2,i] <- L_inv_v[i]/h[(i+1)]-d[i]*h[(i+1)]/6;
  }
  
  return(list(c,d,h));
}
