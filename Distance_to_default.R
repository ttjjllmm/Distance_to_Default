# Project exercises in Quantitative Finance course

# Merton's distance-to-defualt model. Analysing 4 stocks.

# Defining starting and ending dates
mystart<-'2023-01-31'
myend<-'2024-02-01'
data <- read.csv("ymky2h1xbelngour.csv") # Tickers of 4 stocks and corresponding accounting values, data from WRDS data library

stockdata <- data[, c(8,10,11,12,13)] # Making a new data frame with only tic and the accounting values
colnames(stockdata)[1] <- "ticker" # Renaming tic-column to ticker
View(stockdata)

myticker <- stockdata[, 1] # defining the values in column 1 (tickers) as a list called myticker
myticker

nstocks <- as.numeric(length(myticker)) #Needed for the coming loops. Counts the number of companies we want to analyze and thus the number of iterations in the loops



# Looping to get prices each ticker
# NOTE: The individual loops can be attached to the same for-loop, but for clarity, they are split into different parts

prices <- NULL
for( i in 1:nstocks) {
  tempdata <- getSymbols(Symbols = myticker[i], from = mystart, to =myend, auto.assign = FALSE)
  prices <- cbind(prices, tempdata) #cbind() adds them as columns
}
View(prices) # Contains price data for the tickers in "mytickers" list

# Looping to get sigmas for each ticker based on adjusted closing price (columns 6, 12, 18, 24...)
sigmas <- NULL
for( i in 1:nstocks) {
  tempsigma <- sd(diff(log(prices[ ,6*i])), na.rm=TRUE)*sqrt(252) # Each 6th column has the adjusted close price for the stocks
  sigmas <- c(sigmas, tempsigma)
}
round(sigmas, 2) # Stock price volatilities rounded to 2 decimals in a list



# Based on accounting data (annual) for the previous year, in 1000000s of USD
# Looping
shortterm <- NULL
for (i in 1:nstocks) {
  tempshort <- (stockdata[i,2] + stockdata[i, 4])
  shortterm <- c(shortterm, tempshort) #c() makes them a list
}
shortterm # Accounts payable + short/current long-term debt	(from balance sheet) --> ap + dd1


longterm <- NULL
for (i in 1:nstocks) {
  templong <- (stockdata[i, 5])
  longterm <- c(longterm, templong)
}
longterm # Long-term debt from the balance sheet, millions of USD


# Looping to get a list with debts for each firm
D <- NULL
for (i in 1:nstocks) {
  dtemp <- shortterm[i] + 0.5 * longterm[i]
  D <- c(D, dtemp)
}
D # Simplified way to attain debt in 1 year

shares <- stockdata[ ,3]	  # the 3rd column contains 'shares outstanding' in 1000000s 


# P0 (= Initial price) for each stock (this time grouping them)
# E0 (= Initial Equity value)
P0 <- NULL
E0 <- NULL
for (i in 1:nstocks) {
  tempp0 <- as.numeric(tail(prices[ ,6*i],1))
  P0 <- c(P0, tempp0)
  
  tempe0 <- round(P0[i] * shares[i], 2)
  E0 <- c(E0, tempe0)
}
P0
E0


T<-1
myrate<-getSymbols(Symbols='DGS3MO',src='FRED',auto.assign = FALSE) # extracting the risk-free rate
myrate <- myrate[time(myrate) =='2024-01-31']

rfree<-as.numeric(myrate*0.01) # Unit-less
rfree


# A0 = Initial Asset Value
A0ini <- NULL
#Looping
for (i in 1:nstocks) {
  
  tempA0ini<-E0[i] + D[i] 
  A0ini <- c(A0ini, tempA0ini)
}
A0ini

# Defining the Merton formula
Merton<-function(par,E0,sigmaE,r,T,D,lower,upper)
{
  #function to be minimized
  F2plusG2<-function(par,E0,sigmaE,r,T,D)
  {
    A0<-par[1]
    sigmaA<-par[2]
    
    d1<-(log(A0/D)+(r+sigmaA^2/2)*T)/(sigmaA*sqrt(T))
    d2<-d1-sigmaA*sqrt(T)
    
    return((E0-A0*pnorm(d1)+exp(-r*T)*D*pnorm(d2))^2
           +(sigmaE*E0-pnorm(d1)*sigmaA*A0)^2)
  }
  # par contains initial values
  result<-optim(par=par, fn=F2plusG2, gr = NULL,
                E0=E0, sigmaE=sigmaE, r=r, T=T, D=D,
                method= "L-BFGS-B", lower=lower, upper=upper)
  return(result)
}


sigmaAini<-0.1 # Volatility of Assets
myresult <- NULL
for (i in 1:nstocks) {
  
  tempmyresult <- Merton(par = c(A0ini[i], sigmaAini), E0 = E0[i], sigmaE= sigmas[i], r=rfree, T=T, D=D[i], lower=c(A0ini[i]*0.001, 0.001), upper=c(A0ini[i]*10,4))
  myresult <- cbind(myresult, tempmyresult) 
}

myresult # Merton-list for each stock




# Computing A0 and sigmaA for each stock from the list of lists
A0 <- NULL
for (i in 1:nstocks) {
  
  tempA0<-round(x=as.numeric(myresult[[1,i]][1]),digits=2) # In the [[]] defines as usual, first value is the row and second is the column. In the second [], it defines the element in the row
  A0 <- c(A0, tempA0)
}
A0  

sigmaA <- NULL
for (i in 1:nstocks) {
  
  tempsigmaA<-round(x=as.numeric(myresult[[1,i]][2]),digits=2) #Based on the above info, this extracts the second element from each column's first row
  sigmaA <- c(sigmaA, tempsigmaA)
}
sigmaA  


# Calculating d2 for each stock
d1 <- NULL
d2 <- NULL
prob <- NULL
for (i in 1:nstocks) {
  
  tempd1<-(log(A0[i]/D[i])+(rfree+sigmaA[i]^2/2)*T)/(sigmaA[i]*sqrt(T))# As in the BS formula
  d1 <- c(d1, tempd1)
  
  tempd2<-round(d1[i]-sigmaA[i]*sqrt(T),2) # Distance-to-default according to the BS formula.
  d2 <- c(d2, tempd2)
  
  tempprob<-pnorm(-d2[i])
  prob <- c(prob, tempprob)
  
}

d2 # Distance-to-default
round(prob,20) # Probability of default. Values are really close to 0. Thus, rounded to 20 decimals.



