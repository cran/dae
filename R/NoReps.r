"power.exp" <- function(rm=5, df.num=1, df.denom=10, delta=1, sigma=1, 
                        alpha=0.05, print=FALSE)
{
	fcrit <- qf(1-alpha, df.num, df.denom)
	lambda <- rm / 2 * (delta / sigma)^2
	powr <- 1 - pf(fcrit, df.num, df.denom, lambda)
	if (print == TRUE)
		{ print.data.frame(data.frame(rm, df.num, df.denom, alpha, delta, sigma, lambda, powr))}
	powr
}

"power.diff" <- function(r=5, multiple=1, df.num=1, df.denom=expression((df.num+1)*(r-1)), 
                         delta=1, sigma=1, alpha=0.05, power=0.8, print=FALSE)
{
	rm <- r * multiple
	powrdiff <- abs(power-power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), 
	                                delta=delta, sigma=sigma, alpha=alpha, print=print))
}

"no.reps" <- function(multiple=1, df.num=1, df.denom=expression((df.num+1)*(r-1)), 
                        delta=1, sigma=1, alpha=0.05, power=0.8, tol = 0.025, 
                        print=FALSE)
{
	rmini <- 2
	rmaxi <- 50
	local <- TRUE
#
# loop while minimum is local
	while (local)
	{
		local <- FALSE
#
# get a minimum
		minimum <- optimize(power.diff, interval=c(rmini,rmaxi), tol=0.1, print=print,  
					 multiple=multiple, df.num=df.num, df.denom=df.denom, 
					 delta=delta, sigma=sigma, alpha=alpha, power=power)
#
# check for local minimum				
		if (minimum$objective > tol)
		{
			local <- TRUE
			r <- minimum$minimum
			rm <- r * multiple
			if (power < power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), 
										delta=delta, sigma=sigma))
			{
				rmaxi <- floor(minimum$minimum)
				rmini <- rmaxi / 2
			}
			else
			{
				rmini <- ceiling(minimum$minimum)
				rmaxi <- 2 * rmaxi
			}
		}
	}
#
# compute integer number of pure reps and the corresponding power
	r <- ceiling(minimum$minimum)
	rm <- r * multiple
	power <- power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), delta=delta, sigma=sigma)
	list(nreps=r, power=power)
}
