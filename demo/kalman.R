## Comparing particle filtering with the Kalman filter
# Generate data from a local level linear model
N = 1000; W = 1^2; V = 1; m0 = 0; C0 = 1
true.x = rep(NA,N); true.x[1] = rnorm(1,m0,sqrt(C0))
for (i in 2:N) true.x[i] = rnorm(1,true.x[i-1],sqrt(W)) # Evolve x
y = rnorm(N,true.x,sqrt(V))                             # Noisy data

# Run a particle filter
J = 1e2
x  = matrix(NA,N,J); x[1,]  = rnorm(J,m0,C0) # Sample from the prior for x
ws = matrix(NA,N,J); ws[1,] = renormalize.weights(dnorm(y[1],x[1,],sqrt(V),log=TRUE), log=T)

# Run a Kalman filter
m = rep(NA,N); m[1] = m0 # Kalman filter expectation
M = rep(NA,N); M[1] = C0 # Kalman filter variance

for (i in 2:N) {
  # Particle filter
  component   = resample(ws[i-1,],J,"stratified","ess",0.8*J)
  x[i,]       = rnorm(J,x[i-1,component$indices],sqrt(W))
  log.weights = log(component$weights)+dnorm(y[i],x[i,],sqrt(V),log=TRUE)
  ws[i,]      = renormalize.weights(log.weights,log=TRUE)

  # Kalman filter
  K    = (M[i-1]+W)/(M[i-1]+W+V) # Adaptive coefficient
  m[i] = K*y[i]+(1-K)*m[i-1]
  M[i] = K*M[i-1]
}

pf.m = apply(x*ws,1,sum)
plot(m,type='l',ylim=range(pf.m,m),xlab='t',ylab='x')
lines(pf.m,col='red')
legend("bottomleft",inset=0.01,c("Kalman filter mean","Particle filter mean"),
       col=c("black","red"),lty=rep(1,2),bg="white")

