#!/usr/bin/env Rscript
library(partitions)
#################################################################################
# The packages above are downloadable from:
# http://cran.r#project.org/web/packages/partitions/
# http://cran.r#project.org/web/packages/RColorBrewer/
#
# This code is to generate results used for TABLE.txt
# Escpecially for data where coverage, error rate, edit distances
# that are not specified in that file.
#
# To run it use this command:
# $ ./theor_est_general.r
#
# Copyright 2013 (Michiaki Hamada and Edward Wijaya)
#################################################################################

# Modify these 3 variables to suit your need
cv  <- c(100);
edis <- c(0,1,2)
error <- c(1,2);


# When the L is large, it has no effect to the 
# the estimation of theoretical mean
the_L <- 117080 # ref genome length
the_k <- 0:4


###################################################
# Subroutine 
################################################### 
compute_rk <- function(err,k,l) {
  rk <- ((1-err)^(l-k)) *((err/3)^k)
  return(rk); 
}

compute_rk_prime <- function(err,k,d,l) {
   # Equation 12
   # Example compute_rk_prime(2/100,2,4,36)
   # k has to be less than d
   rkp <- 0;
   for (kp in (k-d):k) {
      tmp_rk <- choose(k,kp)* (((1-err)^(l-k)) * ((1 - (err/3)))^(k-kp)) * (err/3)^kp
      rkp <- rkp + tmp_rk;
   }
  return(rkp); 
}

pk_dge0_ege0_d_less_k <- function(n,err,l,L,c,d,krange) {
    # Equation 11
    # k has to be less than d
    rkprime <- compute_rk_prime(err,krange,d,l);
    all_pk_dge0_ege0_d_less_k  <- 0;
    max_m <- find_max_m(n,c,l,L);

    for (m in 0:max_m) {
        tmp <- p_d0_e0(n+m,c,l,L) * choose(n+m,m) * (rkprime)^n * ((1-rkprime)^m)
        all_pk_dge0_ege0_d_less_k <- all_pk_dge0_ege0_d_less_k+tmp;
    }
    
    return(all_pk_dge0_ege0_d_less_k);
}

p_dge0_e0 <- function(n,c,l,L,d) {
    # Equation 3
    N <- floor(c * L / l)
    prob <- pbinom(n,size=N,prob=((1+2 * floor(d/2))/L)) - pbinom(n-1,size=N,prob=((1+2 * floor(d/2))/L));

    return (prob);
}

p_dge0_e0_norm <- function(n,c,l,L,d) {
    # Only used when N large enough e.g. (200x)
    m <- (c/l) * ((2 * floor(d/2)) + 1) 
    sgsq <- (c/l) * ((2 * floor(d/2)) + 1) 
    prob <- dnorm(n,mean=m,sd=sqrt(sgsq));
    return (prob);
}

p_d0_ege0 <- function(n,err,l,L,c,d,krange) {
   # Equation 6
   all_pn_p_d0_ege0 <- c();

   for(ni in 1:length(n)){
    pn_p_d0_ege0 <- 0;
    sum_of_all_k <- 0;
    for (ki in 1:length(krange)){
        pk <- compute_pk(l,err,krange[ki]);
        d_min_k <- d-krange[ki];
        p_k_pk_d0_ege0 <- 0;

        #cat("K ",krange[ki]," d-k= ",d_min_k, " ni_n: ",n[ni],", pk ",pk,"\n",sep="");

        if (d_min_k >= 0){
             p_k_pk_d0_ege0 <- pk * pk_d0_ege0(n[ni],err,l,L,c,d,krange[ki])
        }
        else if (d_min_k < 0   )  {
             p_k_pk_d0_ege0 <- pk * pk_d0_ege0(n[ni],err,l,L,c,d,krange[ki])
        }

         sum_of_all_k <- sum_of_all_k + p_k_pk_d0_ege0;
     }

      pn_p_d0_ege0 <- sum_of_all_k;
      all_pn_p_d0_ege0 <-c(all_pn_p_d0_ege0,pn_p_d0_ege0);
   }
   return(all_pn_p_d0_ege0)
}

###
find_max_m <- function(n,c,l,L) {
    # only take 'n' as single value NOT vector
    # create a function
    # given value of n, find max of 'm' where eq8 start giving 0
    # m = number of reads that contain errorneous nucleotide

    max_m <- 0;

    for (m in 1:100) {
        pb <- p_d0_e0(n+m,c,l,L);
        #cat(m,",",n,",",pb,"\n",sep="");

        # Must NOT set to 0
        if (pb <= 1e-2 ) {
         max_m <- m;
         break;
        }
    }

   return(max_m);
}
p_d0_e0 <- function(n,c,l,L) {
    # n = vector of degree
    # L = reference length
    # N = number of reads
    # c = coverage
    # l = length of read

    # N = cL/l

      N <- floor(c * L / l)

    # compute p(n|d=0,e=0) = NCn (1/L)^n ((L-1/N))^(N-n)
    # return probability value given degree 'n'
     prob <- choose(N,n) * (1/L)^n * ((L - 1 )/ L)^(N-n);

    # Return 0 when 'n' is sufficiently large
    # eg. n=66 for c=200
    return (prob);
}


compute_pk <- function(l,err,k) {
    # compute the probability that a read of length l contain k 
    # errorneous nucleotides
    pk <- choose(l,k) * (1-err)^(l-k) * (err^k);
    return(pk)
}


find_max_k <- function(l,err) {

    max_k <- 0
    for (ki in 0:l) {
        pk <- compute_pk(l,err,ki);
        #cat(ki,"\t",pk,"\n",sep="");
        if (pk <= 1e-2) {
           max_k <- ki;
           break;
        }
    }
  
   return(max_k);
}


In_k <- function(n,d){
    # expand the value of n based on edit distance
    # to determine all the possible neighbors of 'n'

    nof_range <- (2*floor(d/2)) +1; 
    mat <- t(as.matrix(compositions(n, nof_range)))
    #print(mat)
    
    return(mat)
}

q_nk <- function(n,err,l,L,c,k_up){
  # takes 'n' as vector
  # compute 
  # q(n,k) = sum_{m>=0} [ p(n+m|d=0,e=0) * (n+m)Cm * pdk{n} * (1-pdk)^m]


  # m = number of reads that contain errorneous nucleotide
  # k = errorneous nucleotides

   pdk <- sum(compute_pk(l,err,0:k_up));

   prob_all <- c();

   for (n_i in 1:length(n)) {
       max_m <- find_max_m(n_i,c,l,L);

       sum_of_all_prob <- 0
       for (maxm_i in 0:max_m) {
            p_p_d0_e0_ni <- p_d0_e0(n[n_i]+maxm_i,c,l,L);
            #cat("p_eq8_ni ",p_eq8_ni,n[n_i],maxm_i,"\n", sep =",");

            prod_all <- p_p_d0_e0_ni * choose((n[n_i]+maxm_i),maxm_i) * (pdk)^(n[n_i]) * (1-pdk)^maxm_i
            #cat("Eq16. PROD ALL ",prod_all,pdk,p_p_d0_e0_ni,n[n_i],maxm_i,"\n",sep="\t");
            sum_of_all_prob <- sum_of_all_prob + prod_all;
       }

      prob_all <- c(prob_all,sum_of_all_prob);

   }

  return (prob_all);
}

pk_d0_ege0 <- function(n,err,l,L,c,d,krange) {
   # Equation 8
   # n here is not in RANGE format
   # but individual value.
    rk <- compute_rk(err,krange,l);
    all_pk_d0_ege0  <- 0;
    max_m <- find_max_m(n,c,l,L);

    for (m in 0:max_m) {
        tmp <- p_d0_e0(n+m,c,l,L) * choose(n+m,m) * (rk)^n * ((1-rk)^m)
        all_pk_d0_ege0 <- all_pk_d0_ege0 + tmp;
    }


 return(all_pk_d0_ege0)

}


pk_dge0_ege0 <- function(n,err,l,L,c,d,krange) {
   # Equation 10
   # n here is not in RANGE format
   # but individual value.

   #compute p^(k)(n|d>=0,e>=0) = sum_{n_i\inIn}[\prod_i q(n_i,d-2|i]-k)]

    all_prod_sum <- 0

    for (ki in 1:length(krange)){
        i_range<-c(-floor((d-krange[ki])/2):floor((d-krange[ki])/2));
        d_min_k <- d-krange[ki];

        if (d_min_k >= 0){
            # To set range for loop
            nof_col <- (2*floor(d_min_k/2)) +1;           #range for n_i 
            nof_row<-choose(n+(2*floor(d_min_k/2)), n)    #range for set member of (In)
            In_mat <- In_k(n,d_min_k);

            #print(In_mat)

                for (ri in 1:nof_row){
                    all_in_prod <- 1;

                    #Inner loop for product
                    for (ci in 1:length(i_range)) {
                        ni<-In_mat[ri,ci]
                            for (ki in 1:length(krange)){
                                k_fact2 <- d-(2*abs(i_range[ci]))-krange[ki];
                                prod_in <- q_nk(ni,err,l,L,c,k_fact2);

                                #cat("PROD IN ",prod_in,"\n",sep="");

                                all_in_prod <- all_in_prod*prod_in
                            }
                    }

                    all_prod_sum <- all_prod_sum + all_in_prod;
                }

        }
    }

 
 #cat(all_prod_sum,"\n",sep="");
 return(all_prod_sum)

}


p_dge0_ege0 <- function(n,err,l,L,c,d,krange) {
   # Equation 9
   all_pn_p_dge0_ege0 <- c();

   for(ni in 1:length(n)){
    pn_p_dge0_ege0 <- 0;
    sum_of_all_k <- 0;
    for (ki in 1:length(krange)){
        pk <- compute_pk(l,err,krange[ki]);
        d_min_k <- d-krange[ki];
        p_k_pk_dge0_ege0 <- 0;

        #cat("K ",krange[ki]," d-k= ",d_min_k, " ni_n: ",n[ni],", pk ",pk,"\n",sep="");

        if (d_min_k >= 0){
             p_k_pk_dge0_ege0 <- pk * pk_dge0_ege0(n[ni],err,l,L,c,d,krange[ki])
        }
        else if (d_min_k < 0   )  {
             p_k_pk_dge0_ege0 <- pk * pk_dge0_ege0_d_less_k(n[ni],err,l,L,c,d,krange[ki])
        }

         sum_of_all_k <- sum_of_all_k + p_k_pk_dge0_ege0;
     }

      #cat("SUM OF ALL K ",sum_of_all_k,"\n",sep="");
      pn_p_dge0_ege0 <- sum_of_all_k;
      all_pn_p_dge0_ege0 <-c(all_pn_p_dge0_ege0,pn_p_dge0_ege0);
   }
   return(all_pn_p_dge0_ege0)
}


###################################################
# End of subroutines 
################################################### 

cat("#Cv\tdis\terror\tthe.avg\tth.sd\temp.avg\temp.std")
for (i in 1:length(cv)) {

    the_N<- floor(cv[i]*the_L/36);

    for (j in 1:length(edis)) {
        for (k in 1:length(error)) {

              if ((cv[i] > 50 || cv[i]==10) && edis[j] == 6) {
                next;
              }


              title <- paste(cv[i],"x, e=",error[k],"%, d=",edis[j],sep="");
                
               X_max <- 50

              ######################################################
              # compute theoretical  densities
              ######################################################

              if (edis[j] >= 0 & error[k] ==0 & cv[i] >= 200) {
                 # Case 1 when N is large
                 # Use normal distributions, N must be large
                 new_theor_dens <- p_dge0_e0_norm(0:X_max,cv[i],36,the_L,as.numeric(edis[j]))
              }
              else if (edis[j] >= 0 & error[k] ==0 & cv[i] < 200) {
                 # Case 1 when N is small
                 new_theor_dens <- p_dge0_e0(0:X_max,cv[i],36,the_L,as.numeric(edis[j]))
              }
              else if (edis[j] == 0 & error[k] >=0) {
                 # Case 2 
                new_theor_dens <- p_d0_ege0(0:X_max,error[k]/100,36,the_L,cv[i],as.numeric(edis[j]),the_k);
              }
              else {
                # Case 3 most general
                new_theor_dens <- p_dge0_ege0(0:X_max,error[k]/100,36,the_L,cv[i],as.numeric(edis[j]),the_k);
              }

              ts <- sum(new_theor_dens);
              cat("\n")


              ###########################################################
              # Compute theoretical means and SD based on densities above
              ###########################################################

              new_theor_dat <- floor(new_theor_dens * the_N);
              new_theor_datrep <- rep(0:X_max,new_theor_dat);

              theor_avg <- round(mean(new_theor_datrep),digit=2)
              theor_sd   <- round(sd(new_theor_datrep),digit=2);

              cat(cv[i],edis[j],error[k],theor_avg, theor_sd, theor_avg, theor_sd,"\t",sep="\t");
 
        }	
    }  


}


