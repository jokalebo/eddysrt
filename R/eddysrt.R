#' eddysrt
#'
#' @name eddysrt.R#1
#' @title: Surface Renewal Parameters
#'
#' @description:  Calculating Surface Renewal Parameters for a 30 minute - 20 Hz Dataset
#'
#' Read in the file
#'
#' This function will take in the 30 minute -20 Hz dataset and calculate tau, amplitude and r
#'
#' @param z A value of 20 is utilized (derived from 20 Hz per second of data measurement)
#' @param df The dataframe consisting of  30 minute - 20 Hz file (30 * 20 * 60)=36000 values
#' @return a numeric matrix vector with 3 values ( rvalue, tau and amplitude)
#' @details The matrix vector:rvalue (seconds) - maximum ratio of 3rd moment to r, tau (total ramp duration) and amplitude of ramp
#' @seealso{NISTunits}
#' @export
#'
#'
eddysrt <-function(z, df){
  dat <- data.frame(col1=numeric(z), col2=numeric(z), col3=numeric(z),stringsAsFactors = FALSE)
##require(NISTunits)
df$H2O<-1000*df$H2O/18.01528

tempTot3<-0
p=20
#### Third Moment ####
for (p in 1:p)
{

  tempTot2<-0
  tempTot3<-0
  tempTot5<-0

  dat$col1[p]<-p
  dat$col2[p]<-p/20

  tem2<-0
  tem3<-0
  tem5<-0
  totnum<-dim(df)[1]
  for (i in (p+1):totnum)
  {
    tem2<-(df$H2O[i]-df$H2O[i-p])^2
    tempTot2<-tem2+tempTot2
    tem3<-(df$H2O[i]-df$H2O[i-p])^3
    tempTot3<-tem3+tempTot3
    tem5<-(df$H2O[i]-df$H2O[i-p])^5
    tempTot5<-tem5+tempTot5
  }

  dat$col3[p]<-(tempTot2/(totnum-dat$col1[p]))
  dat$col4[p]<-(tempTot3/(totnum-dat$col1[p]))
  dat$col5[p]<-(tempTot5/(totnum-dat$col1[p]))

}

dat$moment2<-abs(dat$col3)
dat$moment3<-abs(dat$col4)
dat$moment5<-abs(dat$col5)

dat$p1<-10*dat$col3 - (dat$col5/dat$col4)
dat$p<-10*dat$col3 - (dat$col5/dat$col4)
dat$diff<-dat$p - dat$p1
dat$q <- 10*dat$col4
dat$D<-(((dat$q)^2)/4)+(((dat$p)^3)/27)

Math.cbrt <- function(k)
{
  return(sign(k) * abs(k)^(1/3))

}

for (i in 1:20) {
  if (dat$D[i]>0)
  {
    one1<-(dat$D[i])^0.5
    one2<--0.5*dat$q[i]
    dat$x1[i]<-Math.cbrt(one2+one1)
    dat$x2[i]<-Math.cbrt(one2-one1)
    dat$a[i]<-dat$x1[i]+dat$x2[i]
    rm(one1)
    rm(one2)
  }
  else {
    if (dat$D[i]<0)
    {
      one1<-(dat$D[i])^0.5
      one2<--0.5*dat$q[i]
      dat$x1[i]<-Math.cbrt(one1+one2)
      dat$x2[i]<-Math.cbrt(one1-one2)
      rm(one1)
      rm(one2)

      part1<-((-1)*dat$p[i]/3)^3
      b<-NISTradianTOdeg(acos((-0.5*dat$q[i])/((part1)^0.5)))

      a1=2*((-dat$p[i]/3)^0.5)*(cos(NISTdegTOradian(b/3)))
      a2=2*((-dat$p[i]/3)^0.5)*(cos(NISTdegTOradian(120+b/3)))
      a3=2*((-dat$p[i]/3)^0.5)*(cos(NISTdegTOradian(240+b/3)))
      c<-max(a1, a2, a3)
      dat$a[i]<-c
      rm(c, a1, a2, a3, part1)

    }
  }
}

dat$answer<--(1*dat$col2*(dat$a)^3)/dat$col4

dat$ratio <- dat$a/dat$answer

##require(polynom)

for (i in 1:20){

  p<-polynomial(coef = c(dat$q[i],dat$p[i], 0, 1))
  p

  pz <- solve(p)
  pz
  length(pz)
  for (j in 1:length(pz))
  {
    assign(paste("sol", j, sep = ""), pz[j]) }
  dat$sol1[i]<-sol1
  dat$sol2[i]<-sol2
  dat$sol3[i]<-sol3
  rm(pz)
  rm(sol1, sol2, sol3)
}

for(i in 1:20)
{
  z1<-dat$sol1[i]
  z2<-dat$sol2[i]
  z3<-dat$sol3[i]

  dat$realsol1[i]<-Re(z1)
  dat$Imsol1[i]<-Im(z1)
  dat$realsol1[i]
  dat$Imsol1[i]


  dat$realsol2[i]<-Re(z2)
  dat$Imsol2[i]<-Im(z2)
  dat$realsol2[i]
  dat$Imsol2[i]

  dat$realsol3[i]<-Re(z3)
  dat$Imsol3[i]<-Im(z3)
  dat$realsol3[i]
  dat$Imsol3[i]

  max(dat$realsol1[i],dat$realsol2[i],dat$realsol3[i])
  dat$valA[i]<-max(dat$realsol1[i],dat$realsol2[i],dat$realsol3[i])
  dat$DandS[i]<-((-1*dat$col2[i])*(dat$valA[i])^3)/dat$col4[i]
  dat$newRatio[i]<-dat$valA[i]/dat$DandS[i]
  dat$calc[i]<-abs(dat$col4[i])/dat$col2[i]
  }


best<-which.max(dat$calc)
rvalue<-c(dat$col2[best])
tau<-c(dat$DandS[best])
amplitude<-c(dat$valA[best])
ak<-data.frame(rvalue, tau, amplitude)
return(ak)

}

