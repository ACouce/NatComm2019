# .. This routine simulates a single serial culture experiment. It is intended to be 
# .. called by a script 'single_timeseries_X.R'

threshold=1

while ((sum(x)/(sum(x)+sum(y)) <= threshold) && (t<time) && (sum(x)>0)) { 
  t=t+1
  
  # .. this loop represents one day in a flask 
  while ( sum(x)+sum(y) < Nmax ) {

    # .. deterministic growth
    g=g+1
    x[]=rx[]*x
    y[]=ry[]*y

	for (i in 1:n_del) {

	    # .. stochastic allocation of deleterious mutations
	    mdx<- rpois(1,lambda=x[1,i]*mu_del)
        if (is.na(mdx)) { # .. Normal approximation for values > .Machine$integer.max
            mdx<-round(rnorm(1,mean=x[1,i] * mu_del, sd=sqrt(x[1,i] * mu_del)))
        }
        
	    mdy<- rpois(1,lambda=y[1,i]*mu_del*m)
        if (is.na(mdy)) { # .. Normal approximation for values > .Machine$integer.max
            mdy<-round(rnorm(1,mean=y[1, i] * mu_del * m, sd=sqrt(y[1, i] * mu_del * m)))
        }

	    # .. stochastic allocation of lethal mutations
        mlx<- rpois(1,lambda=x[1,i]*mu_let)
        if (is.na(mlx)) { # .. Normal approximation for values > .Machine$integer.max
            mlx<-round(rnorm(1,mean=x[1, i] * mu_let, sd=sqrt(x[1, i] * mu_let)))
        }
	    mly<- rpois(1,lambda=y[1,i] * mu_let * k * m)
        if (is.na(mly)) { # .. Normal approximation for values > .Machine$integer.max
             mly<-round(rnorm(1,mean=y[1, i] * mu_let * k * m,sd=sqrt(y[1, i] * mu_let * k * m)))
        }
	    #mly<- rpois(1,lambda=y[1,i] * mu_let * m) # .. activate for spectrum-independent lethal mutation 
        # if (is.na(mly)) { # .. Normal approximation for values > .Machine$integer.max
        #     mly<-round(rnorm(1,mean=y[1, i] * mu_let * m,sd=sqrt(y[1, i] * mu_let * m)))
        #}

	    # .. the next if is to handle the boundaries of the genotypic matrix
	    a=1
	    if (i==n_del) { a=0 }

	    x[1,i]=x[1,i]-mdx-mlx
	    x[1,i+a]=x[1,i+a]+mdx
	
	    y[1,i]=y[1,i]-mdy-mly
	    y[1,i+a]=y[1,i+a]+mdy
	
	    x[x<0]=0
	    y[y<0]=0
    
	}

  }
    ratiox=x/(sum(x)+sum(y))
    ratioy=y/(sum(x)+sum(y))

	x_t[t,1]<-sum(x)/sum(y)
    g_t[t,1]<-g

    Xfinal<-sum(x)
    Yfinal<-sum(y)

    if (sum(x)/(sum(x)+sum(y)) >= threshold) {  #if the threshold is hit, record final population sizes
	    Dpt[1,rep]<-Xfinal
	    Dqt[1,rep]<-Yfinal
	    break
	}

   # .. the stochastic bottleneck
	for (i in 1:n_del) {

	  if ((ratiox[1,i]>0.1)) {
	    x[1,i]=rbinom (1, Btlnck, ratiox[1,i]) 
	  }
	  else if ((ratiox[1,i]>0)) { # .. Poisson approximation for small p

	    x[1,i]=rpois (1, lambda=Btlnck*ratiox[1,i])
	  } 
	  else {
            x[1,i]=0
	  }

	  if ((ratioy[1,i]>0.1)) {
	    y[1,i]=rbinom (1, Btlnck, ratioy[1,i])
	  }
	  else if ((ratioy[1,i]>0)) { # .. Poisson approximation for small p

	    y[1,i]=rpois (1, lambda=Btlnck*ratioy[1,i])
	  } 
	  else {
            y[1,i]=0
	  }

	}

    if (sum(x)/(sum(x)+sum(y)) >= threshold) {  #in case the threshold is hit as a product of the bottleneck
	    Dpt[1,rep]<-Xfinal
	    Dqt[1,rep]<-Yfinal

	    break
	}

}

if (sum(x)==0) { #in case the antimutator allele is lost by drift
	g<-NA
}
