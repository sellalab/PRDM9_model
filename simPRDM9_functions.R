### These functions are designed to implement a model for the co-evolution of PRDM9 and its binding sites
### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence

### calculate.Pf - solves for the free concentration of PRDM9.
### Exact solutions are used if 1 or 2 heats are considered.
### When more than 2 heats are considered, the solution is derived numerically.
calculate.pf = function(Si, Ki, Pt=5000, Pf0 = 0, acc = -5:11) {

  # Exact solution for Pf when 1 heat is considered
  if (length(Si) == 1) {
    Pf = (Pt - Ki - 4*Si + ((Ki + 4*Si - Pt)^2 + (4*Ki*Pt))^0.5)/2
  }

  # Exact solution for Pf when 2 heats are considered
  if (length(Si) == 2) {
    b = sum(Ki) + 4*sum(Si) - Pt
    c = Ki[1]*Ki[2] + 4*Si[1]*Ki[2] + 4*Si[2]*Ki[1] - Pt*sum(Ki)
    d = -Pt*Ki[1]*Ki[2]

    Q = (3*c - b^2)/9
    R = (9*b*c - 27*d - 2*(b^3))/54

    Q = as.complex(Q)
    R = as.complex(R)

    S = (R + ((Q^3) + R^2)^(1/2))^(1/3)
    t = (R - ((Q^3) + R^2)^(1/2))^(1/3)

    # i = complex(real = 0, imaginary = 1)
    # x1 = S + t - b/(3*a)
    # x2 = -(S+t)/2 - b/(3*a) + (S-t)*i*(3^0.5)/2
    # x3 = -(S+t)/2 - b/(3*a) - (S-t)*i*(3^0.5)/2
    Pf = Re(S + t - b/3)
  }

  # Numerical solution for Pf when more than 2 heats are considered
  if (length(Si) > 2) {
    Pf = Pf0
    for (i in 1:length(acc)) {
      accuracy = acc[i]
      while (Pt - sum(Si*4*Pf/(Pf+Ki)) > Pf) {
        Pf = Pf+1*10^(-accuracy)
      }
      Pf = Pf-1*10^(-accuracy)
    }
    dif = Pf - (Pt - sum(Si*4*Pf/(Pf+Ki)))
    Pf = Pf+1*10^(-acc[length(acc)])
    dif2 = Pf - (Pt - sum(Si*4*Pf/(Pf+Ki)))
    if (abs(dif) == min(abs(c(dif,dif2)))) {
      Pf = Pf-1*10^(-acc[length(acc)])
    }
  }
  Pf
}

### calculate.W.hom - calculates the fitness of a homozygote
calculate.W.hom = function(Si, Ki, Pt=5000, D=300, r=1/40) {
  Pf = calculate.pf(Si=Si,Ki=Ki,Pt=Pt)
  Hi = Pf/(Pf+Ki)
  c = D/(Pt-Pf)
  ai = 1-(1-c*(Hi^2)*(2-Hi)*(2-c*Hi))^2
  W = 1-prod((1-ai)^(Si*r))
  W
}

### calculate.W.het - calculates the fitness of a heterozygote
calculate.W.het = function(SiA, SiB, Ki, Pt=5000, D=300, r=1/40) {
  PfA = calculate.pf(Si=SiA,Ki=Ki,Pt=Pt/2) # half PRDM9 expression
  PfB = calculate.pf(Si=SiB,Ki=Ki,Pt=Pt/2) # half PRDM9 expression
  HiA = PfA/(PfA+Ki)
  HiB = PfB/(PfB+Ki)
  c = D/(Pt-PfA-PfB) # need to consider Pf from both alleles
  aiA = 1-(1-c*(HiA^2)*(2-HiA)*(2-c*HiA))^2
  aiB = 1-(1-c*(HiB^2)*(2-HiB)*(2-c*HiB))^2
  W = 1-prod((1-aiA)^(SiA*r))*prod((1-aiB)^(SiB*r)) # twice as many 'attemPts' in terms of sites
  W
}

### calc.hom.details - calculates for a homozygote, fitness (W) and the number
### of sites lost to gene conversion for both classes of sites (dS1 and dS2)
### Note: number of sites lost needs to be scaled by frequency of homozygotes
### i.e., the number of sites lost is given under an assumption that the considered allele is fixed
calc.hom.details = function(Si, Ki, Pt = 5000, pf = 0, D = 300, r = 1/40, B = 0.7, N=1e3, u = 1.25e-7) {
  if (pf == 0) {
    pf = calculate.pf(Si = Si, Ki = Ki, Pt = Pt)
  }
  H = pf/(pf+Ki)
  c = D/(Pt-pf)
  a = 1-(1-(c*(H^2))*(2-H)*(2-c*H))^2
  W = 1-prod((1-a)^(r*Si))
  dS = -4*N*u*c*H*B*Si
  c(W,dS)
}

### calculate.pf.matrix - alculates the free concentration of PRDM9 given a particular binding distribution for both homozygotes and heterozygotes
calculate.pf.matrix = function(Si, Ki, Pt = 5000) {

  # homozygotes
  b = sum(Ki) + 4*apply(Si,1,sum) - Pt
  c = Ki[1]*Ki[2] + 4*Si[,1]*Ki[2] + 4*Si[,2]*Ki[1] - Pt*sum(Ki)
  d = -Pt*Ki[1]*Ki[2]

  Q = (3*c - b^2)/(9)
  R = (9*b*c - 27*d - 2*(b^3))/(54)

  Q = as.complex(Q)
  R = as.complex(R)

  S = (R + ((Q^3) + R^2)^(1/2))^(1/3)
  t = (R - ((Q^3) + R^2)^(1/2))^(1/3)

  pf.hom = as.numeric(Re(S + t - b/(3)))

  # heterozygotes
  Pt = Pt/2

  b = sum(Ki) + 4*apply(Si,1,sum) - Pt
  c = Ki[1]*Ki[2] + 4*Si[,1]*Ki[2] + 4*Si[,2]*Ki[1] - Pt*sum(Ki)
  d = -Pt*Ki[1]*Ki[2]

  Q = (3*c - b^2)/(9)
  R = (9*b*c - 27*d - 2*(b^3))/(54)

  Q = as.complex(Q)
  R = as.complex(R)

  S = (R + ((Q^3) + R^2)^(1/2))^(1/3)
  t = (R - ((Q^3) + R^2)^(1/2))^(1/3)

  pf.het = as.numeric(Re(S + t - b/(3)))

  rbind(pf.hom, pf.het)
}


if (FALSE) {
  N = 1e3
  u = 1.25e-7
  v = 1e-5
  Pt = 5000
  D = 300
  B = 0.7
  r = 1/40
  k1 = 5
  S1 = 1:5000
  k2 = NA
  S2 = 200000
  prop.S2.bound = 0.99
  num.gen = 5500000
  starting.alleles = 3
  record.per.gen = TRUE
  record.per.allele = TRUE
  time.run = TRUE
  sim.name = "1e3_5"
}

### run.2heat.model.sim - runs a simulation of the two heat models and records various details
run.2heat.model.sim = function(N = 1e6, u = 1.25e-7, v = 1e-5,
                               Pt = 5000, D = 300, B = 0.7, r = 1/40,
                               k1 = 50, S1 = 1:5000,
                               k2 = NA, S2 = 200000, prop.S2.bound = 0.99,
                               num.gen = 1e3,
                               starting.alleles = 3,
                               record.per.gen = TRUE, record.per.allele = TRUE, time.run = TRUE,
                               sim.name = "1e6_50") {

  ############################################
  ### setting things up for the simulation ###
  ############################################

  if (time.run) {
    start.time = Sys.time()
  }

  ### If k2 is not provided, calculate under an assumPtion that prop.S2.bound
  ### describes the proportion bound when S1 = 0
  if (is.na(k2)) {
    Pb = Pt*prop.S2.bound
    Pf = Pt-Pb
    H2 = Pb/(4*S2)
    k2 = (Pf/H2) - Pf
  }

  ### Establishing table of segregating PRDM9 alleles and their details
  prdm9.table = data.frame(AN = 1:starting.alleles,
                           AF = as.numeric(table(sample(1:starting.alleles,2*N,TRUE))/(2*N)),
                           initial.S1 = sample(S1,starting.alleles),
                           S1 = 0,
                           S2 = S2,
                           initial.gen = 0,
                           gen = 0,
                           t.time = 0)
  prdm9.table$S1 = prdm9.table$initial.S1
  max.AN = starting.alleles

  ### initialize result tables
  if (record.per.gen) {
    per.gen.results = NULL
  }
  if (record.per.allele) {
    per.allele.results = NULL
  }

  #######################################
  ### actually running the simulation ###
  #######################################
  gen = 1
  pb = txtProgressBar(0, num.gen, style = 3)
  for (gen in 1:num.gen) {
    setTxtProgressBar(pb,gen)

    ### generate new PRDM9 alleles (to appear by mutation) in batches of 1000 generations worth
    if ((gen %% 1000) == 1) {
      num.mut = rbinom(1000,2*N,v)

      unused.alleles = data.frame(AN = (max.AN+1):(max.AN+sum(num.mut)),
                                  AF = 1/(2*N),
                                  initial.S1 = sample(S1, sum(num.mut), replace = TRUE),
                                  S1 = 0,
                                  S2 = S2,
                                  initial.gen = unlist(sapply(1:length(num.mut), FUN = function(x) {
                                    rep(x,num.mut[x])
                                  }))+(gen-1),
                                  gen = 0,
                                  t.time = 1/(2*N))
      unused.alleles$S1 = unused.alleles$initial.S1
      unused.alleles$gen = unused.alleles$initial.gen
      max.AN = max(unused.alleles$AN)
    }

    ### mutate new PRDM9 alleles
    t.gen = (gen %% 1000)
    if (t.gen == 0) {t.gen = 1000}
    if (num.mut[t.gen] > 0) {
      for (i in 1:num.mut[t.gen]) {
        if (dim(prdm9.table)[1] == 1) {
          mut.allele = prdm9.table$AN
        } else {
          mut.allele = sample(prdm9.table$AN, 1, prob = prdm9.table$AF)
        }
        prdm9.table$AF[which(prdm9.table$AN == mut.allele)] = prdm9.table$AF[which(prdm9.table$AN == mut.allele)] - 1/(2*N)
      }
      prdm9.table = rbind(prdm9.table, unused.alleles[which(unused.alleles$initial.gen == gen),])
      unused.alleles = unused.alleles[which(unused.alleles$initial.gen > gen),]
    }

    prdm9.table$t.time = prdm9.table$t.time + prdm9.table$AF

    ### remove alleles at frequency of zero
    prdm9.table$gen = gen
    to.remove = which(prdm9.table$AF == 0)
    if (length(to.remove) > 0) {
      ### record details for each allele removed
      if (record.per.allele) {
        write.table(prdm9.table[to.remove,],paste("per.allele.results",sim.name,"txt",sep="."),
                    append = TRUE,col.names = FALSE,row.names = FALSE)
      }
      prdm9.table = prdm9.table[which(prdm9.table$AF != 0),]
    }

    ### If there is only 1 segregating PRDM9 allele, only need to calculate values in homoyzgotes
    if (dim(prdm9.table)[1] == 1) {

      temp = calc.hom.details(c(prdm9.table$S1[1], prdm9.table$S2[1]), c(k1,k2), Pt = Pt, D=D, r=r, B=B, N=N, u=u)
      prdm9.table$S1 = prdm9.table$S1 + temp[2]
      prdm9.table$S2 = prdm9.table$S2 + temp[3]

      # record per generation results
      if (record.per.gen) {
        # mean S1, mean fitness, diversity at PRDM9
        write.table(t(c(prdm9.table$S1[1], temp[1], 0)),
                    paste("per.gen.results",sim.name,"txt",sep="."),
                    append = TRUE,col.names = FALSE,row.names = FALSE)
      }

    } else {
    ### If there are multiple segregating PRDM9 alleles, need to calculate values for all possible genotypes
      # this involves calculating matrices for fitness / number of sites lost in each

      # this specifies which positions of such matrices correspond to homozygotes
      diag = (1:dim(prdm9.table)[1])*dim(prdm9.table)[1] - ((dim(prdm9.table)[1]:1)-1)

      # generate a genotype frequency matrix
      freq.table = matrix(prdm9.table$AF, dim(prdm9.table)[1], dim(prdm9.table)[1], byrow = TRUE)*prdm9.table$AF

      # calculate pf values for all PRDM9 alleles
      # (first row is pf in homozygotes, second row is pf in heterozygotes)
      pf = calculate.pf.matrix(prdm9.table[,c(4,5),],c(k1,k2), Pt = Pt)

      # calculate values of H1 and H2 (in both homs and hets, as above)
      H1 = pf/(pf+k1)
      H2 = pf/(pf+k2)

      # generate matrices for values of g and alpha - 0.1953246 [4.683175]
      # 0.02245263 [0.4920017]
      c.table = D/(Pt - (matrix(pf[2,], dim(prdm9.table)[1], dim(prdm9.table)[1], byrow = TRUE)+pf[2,]))
      c.table[diag] = D/(Pt-pf[1,])

      # 0.03465947 [0.7724192]
      H1.table = matrix(H1[2,], dim(prdm9.table)[1], dim(prdm9.table)[1])
      H2.table = matrix(H2[2,], dim(prdm9.table)[1], dim(prdm9.table)[1])
      H1.table[diag] = H1[1,]
      H2.table[diag] = H2[1,]

      # 0.0876864 [1.722615]
      g1.table = c.table*H1.table*B
      g2.table = c.table*H2.table*B

      # 0.1428395 [2.742673]
      a1.table = 1 - (1 - (c.table*(H1.table^2))*(2-H1.table)*(2-(c.table*H1.table)))^2
      a2.table = 1 - (1 - (c.table*(H2.table^2))*(2-H2.table)*(2-(c.table*H2.table)))^2

      # generate fitness matrix - 0.1525953 [3.6722]
      W.table = ((1-a1.table)^(r*prdm9.table$S1))*((1-a2.table)^(r*prdm9.table$S2))
      W.table = 1 - (W.table*t(W.table))
      W.table[diag] = 1 - ((1-a1.table[diag])^(r*prdm9.table$S1))*((1-a2.table[diag])^(r*prdm9.table$S2))

      #record per generation results - 0.02938221 [0.9043908]
      if (record.per.gen) {
        pop.W = sum(freq.table*W.table)
        write.table(t(c(sum(prdm9.table$AF*prdm9.table$S1), pop.W, 1-sum(prdm9.table$AF^2))),
                    paste("per.gen.results",sim.name,"txt",sep="."),
                    append = TRUE,col.names = FALSE,row.names = FALSE)
      }

      # calculate parental/expected PRDM9 allele frequencies - 0.09690132 [1.935955]
      freq.table = (freq.table*W.table)/pop.W
      exp.freqs = apply(freq.table,1,sum)

      # calculate new PRDM9 allele frequencies - 0.004789788 [0.08204269]
      prdm9.table$AF = as.numeric(rmultinom(1:dim(prdm9.table)[1], 2*N, exp.freqs)/(2*N))

      # calculate number of sites lost per PRDM9 allele - 0.2087681 [4.643561]
      freq.table = 2*freq.table
      freq.table[diag] = freq.table[diag]/2
      prdm9.table$S1 = prdm9.table$S1 - apply(4*N*u*g1.table*freq.table*prdm9.table$S1,1,sum)
      prdm9.table$S2 = prdm9.table$S2 - apply(4*N*u*g2.table*freq.table*prdm9.table$S2,1,sum)

    }
  }
  close(pb)

  if (time.run) {
    end.time = Sys.time()
    time.elapsed = as.numeric(end.time-start.time)
    print(paste("Simulation took",time.elapsed))
  }


}






