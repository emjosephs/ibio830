```{r LinearRegressionMCMC,echo=FALSE,cache=TRUE}
ln.post.prob <- function(y,x,alpha,beta,n){
	return(
		ln.L(y,x,alpha,beta,n) + 
		dnorm(alpha,mean=0,sd=1e6,log=TRUE) + 
		dnorm(beta,mean=0,sd=1e6,log=TRUE)
	)
}

n.iter <- 1e6
samp.freq <- 2e3
alpha.sample <- rep(NA,1+n.iter/samp.freq)
beta.sample <- rep(NA,1+n.iter/samp.freq)
post.prob.sample <- rep(NA,1+n.iter/samp.freq)
alpha <- rnorm(1,0,1)
beta <- rnorm(1,0,1)
post.prob <- ln.post.prob(y=occ.lfs,x=plant.height,alpha=alpha,beta=beta,n=n.leaves)
alpha.sample[1] <- alpha
beta.sample[1] <- beta
post.prob.sample[1] <- post.prob
for(i in 2:n.iter){
	alpha.prime <- alpha + rnorm(1,0,0.1)
	beta.prime <- beta + rnorm(1,0,0.1)
	post.prob.prime <- ln.post.prob(y=occ.lfs,x=plant.height,alpha=alpha.prime,beta=beta.prime,n=n.leaves)
	R <- exp(post.prob.prime - post.prob)
	if(R > runif(1)){
		alpha <- alpha.prime
		beta <- beta.prime
		post.prob <- post.prob.prime
	}
	if(i %% samp.freq == 0){
		alpha.sample[1+i/samp.freq] <- alpha
		beta.sample[1+i/samp.freq] <- beta
		post.prob.sample[1+i/samp.freq] <- post.prob
	}
}
```


## Binomial regression

```{r,echo=FALSE,fig.width=7.5,fig.height=5.5}
	par(mfrow=c(1,2),mar=c(4,2,1,1))
		plot(alpha.sample,lwd=2,col=1,type='l',xlab="",ylab="") ; abline(h=ph.alpha,col=2,lwd=2)
			text(x=250,y=-0.75,labels="alpha",cex=1.5)
		plot(beta.sample,lwd=2,col=1,type='l',xlab="",ylab="") ; abline(h=ph.beta,col=2,lwd=2)
			text(x=250,y=1.05,labels="beta",cex=1.5)
		mtext(side=1,adj=-4.25,padj=3,text="sampled mcmc iterations",cex=1.5)
```


## Binomial regression {#slideID} 

<style> 
  #slideID > p { 
    margin-top: -50px; 
  } 
</style>

<div class="centered">
```{r,echo=FALSE,fig.width=7,fig.height=5.5}

plot(plant.height,occ.lfs,pch=20,col=adjustcolor(1,0.3),ylim=c(0,n.leaves),xlab="plant height",ylab="occupied leaves")
```
</div>

## Binomial regression
<div class="centered">
```{r bin.reg.samples,echo=FALSE,fig.width=7,fig.height=5.5,cache=TRUE}
tmp.seq <- seq(0,10,length.out=100)
plot(plant.height,occ.lfs,pch=20,col=adjustcolor(1,0.3),ylim=c(0,n.leaves),xlab="plant height",ylab="occupied leaves")
	invisible(
		sapply(sample(1:501,50),
			function(i){
				lines(tmp.seq,n.leaves*inv.logit.link(alpha.sample[i]+beta.sample[i]*tmp.seq),lwd=3,col=adjustcolor(2,0.2))
			})
	)
	lines(tmp.seq,n.leaves*inv.logit.link(ph.alpha+ph.beta*tmp.seq),lwd=3,col="blue")
	legend(x="bottomright",col=c(4,2),lwd=3,legend=c("simulated","estimated"))
```
</div>
