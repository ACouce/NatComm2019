# .. This program simulates the serial passage of an asexual population, initially dominated by mutator cells ('y'), in an environment to which it is fully adapted. At the start of each simulation, a single antimutator allele, restoring the mutation rate to wild-type levels, is introduced. Deleterious mutations (with coefficient 's_del)' ocurr in wild-type cells ('x') at a rate 'mu_del'. Mutators not only exhibit a deleterious mutation rate increased by a factor of 'm', but also change the effect of these deleterious mutations by a factor of 'kappa'. The simulation ends when the population is dominated by wild-type cells. Reference: Couce, A. & Tenaillon, O. (2019). Mutation bias and GC content shape antimutator invasions'. TBD.

  rm(list=ls())

# .. start the stopwatch
  clock <- Sys.time()

  # .. PARAMETERS
  replicas=25
  
  # .. population sizes
  Nmax=1e9
  Btlnck=1e7

  x_ini=1
  time=5000
  
  # .. mutation rates 
  m=300
  mu_ben=0
  mu_del=2e-4
  mu_let=1e-5

  # .. fitness effects   
  s_del=-0.008
  kappa=c(0.25,1, 4)

  # .. VARIABLES
  # .. genotypic matrices
  n_del=50
  range=3
  rx=matrix(nrow = 1, ncol = n_del)
  ry=matrix(nrow = 1, ncol = n_del)

  # .. data matrices

  eff=matrix(nrow = replicas, ncol = range)
  prom_eff=matrix(0,nrow = range, ncol = range)
  sd_eff=matrix(0,nrow = range, ncol = range)
  Dpt=matrix(0,1,replicas)
  Dqt=matrix(0,1,replicas)

  # .. ALGORITHM

  x11(width = 8.5, height = 6)
  par(mfrow=c(1,2))

  plot(1:800,1:800, col='white', log='y', ylim=c(5e-8,1e7), type='l', xlab="generations", ylab="ratio antimutator to mutator")
  print(c('m','k','s_del', 's_eff'))

  counter_k=0

  for (k in kappa){ 
      counter_k=counter_k+1
      counter_c=0
  
      for (c in s_del) {
            counter_c=counter_c+1

            # .. VARIABLES
            x=matrix(nrow = 1, ncol = n_del)
            y=matrix(nrow = 1, ncol = n_del)

            ratiox=matrix(nrow = 1, ncol = n_del)
            x_t=matrix(nrow = time,  ncol = 1)

            ratioy=matrix(nrow = 1, ncol = n_del)
            y_t=matrix(nrow = time,  ncol = 1) 

            data_eff=matrix(0,1,replicas)
            g_t=matrix(nrow = time,  ncol = 1)

            # .. fitness matrices

            for (i in 1:n_del) { 
                rx[1,i]=2+(i-1)*c
                ry[1,i]=2+(i-1)*(c*k)
         			 }

            # .. the actual experiments
            rep=0
            while (rep<replicas) {
                rep<-rep+1
                t=0
                g=0
                x[]=0
                x[1,1]=x_ini # .. a single antimutator allele enters the population
                y[]=0
                y[1,1]=Btlnck-x_ini
                ratiox=x/(sum(x)+sum(y))
                ratioy=y/(sum(x)+sum(y))

                x_t[]=0

                source("serial_culture_antimut.R")
                
                # .. restart if the antimutator allele is lost by drift
                if (is.na(g)) {
                    rep<-rep-1
                    next
                }
                data_eff[1,rep]<-g
                lines(g_t[1:800], x_t[1:800], col=counter_k, type='o', lwd=0.5)

            }

        # .. calculation effective selection coefficient as s=log(Rt/R0)/g, R=p/q
        p0<-x_ini
        q0<-Btlnck-x_ini
        eff[,counter_k]<-t(log((Dpt/Dqt)/(p0/q0))/data_eff) 
        print(c(m,k,c, mean(eff[,counter_k], na.rm='T')))


    }
    eff[!is.finite(eff)] <- NA # to handle Inf and NaN
    prom_eff[counter_k,]<-apply(eff,2,mean, na.rm='T')
    sd_eff[counter_k,]<-apply(eff,2,sd, na.rm='T')


}
    boxplot(eff, names=c(0.25,1,4), ylim=c(0,0.1), xlab="spectrum effect (k)", ylab="s_eff")

# .. stop the stopwatch
  print(Sys.time()-clock)


