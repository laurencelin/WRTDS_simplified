
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R')
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_misc.r')


WRTDS = function(obsData, predData, replicationN=50, Yhalfwin=10, Shalfwin=0.5, Qhalfwin=2){ 
	
	prediction <- sapply(seq_len(dim(predData)[1]), function(i){
		if(is.na(predData$logQ[i]) ){
			return <- rep(NA, replicationN)
		}else{
			Yhalfwin_ = Yhalfwin
			Ydist = abs(obsData$ydecimal - predData$ydecimal[i])
			Yweight = (Ydist<=Yhalfwin)*(1-(Ydist/Yhalfwin)^3)^3
				while(sum(Yweight>0)<100){
					Yhalfwin_ = Yhalfwin_*1.1
					Yweight = (Ydist<= Yhalfwin_)*(1-(Ydist/Yhalfwin_)^3)^3
				}#while
			
			Shalfwin_ = Shalfwin
			Sdist = vectorMin( ceiling(Ydist)-Ydist, Ydist-floor(Ydist))
			Sweight = (Sdist <= Shalfwin)*(1-(Sdist/Shalfwin)^3)^3
				while(sum(Yweight>0 & Sweight>0)<100){
					Shalfwin_ = Shalfwin_*1.1
					Sweight = (Sdist <= Shalfwin_)*(1-(Sdist/Shalfwin_)^3)^3
				}#while
				
			Qhalfwin_ = Qhalfwin 
			Qdist = abs(obsData$logQ - predData$logQ[i])
			Qweight = (Qdist <= Qhalfwin_)*(1-(Qdist/Qhalfwin_)^3)^3
				while(sum(Yweight>0 & Sweight>0 & Qweight>0)<100){
					Qhalfwin_ = Qhalfwin_*1.1
					Qweight = (Qdist <= Qhalfwin_)*(1-(Qdist/Qhalfwin_)^3)^3
				}# while
				
			Tweight = Yweight * Sweight * Qweight * (!is.na(obsData$logC)) * (!is.na(obsData$logQ)); # sum(Tweight>0)
			cond = Tweight>0; #sum(cond); sum(Yweight>0); sum(Sweight>0); sum(Qweight>0)
			
			dailyReplication <- sapply(seq_len(replicationN), function(gg){
					checkresults = rep(NA,4)
					while( sum(is.na(checkresults))>0 ){
						tryCatch({
							booststrap = rep(1,sum(cond))
							while( max(table(booststrap))>5 ){
								booststrap = sample(obsData$index[cond],sum(cond),replace=T)
							}#while
							result = lm(logC~ydecimal+logQ+sin2pit+cos2pit, data=obsData[booststrap,], weights= Tweight[booststrap])
							# log(obsData$no3) = beta0 + beta1*obsData$ydecimal + beta2*obsData$logQ + beta3*obsData$sin2pit + beta4*obsData$cos2pit
							# new features coming: visual beta2 A:: {y:Q X:time}; countour beta2 by colors
							# B:: monthly boxplot of beta2; C::percentile Q boxplot of beta2
						
							alpha = sum(exp(result$residuals)*Tweight[booststrap])/sum(Tweight[booststrap])
							predConc = alpha * exp( result$coefficients[1] + sum(result$coefficients[2:5]*predData[i,c('ydecimal','logQ','sin2pit','cos2pit')]) )
							if(is.infinite(predConc) | predConc > 100*max(obsData$C[booststrap]) |  predConc < 0.01 *min(obsData$C[booststrap]) ){
								print(paste(
									alpha,min(result$residuals),max(result$residuals),
									result$coefficients[1],
									result$coefficients[2],
									result$coefficients[3],
									result$coefficients[4],
									result$coefficients[5]));
								checkresults <- c(NA, NA, NA, NA);
							}else{
							checkresults <- c(
								predConc,
								result$coefficients[3], # coefficent for the log(Q)
								sum(result$coefficients[4:5]*predData[i,c('sin2pit','cos2pit')]),
								summary(result)$r.squared);
							}
						},warning = function(w){
						    checkresults <- c(NA, NA, NA, NA)		
						},error=function(e){
						    checkresults <- c(NA, NA, NA, NA)
						})# tryCatch
					}#while	
					return <- checkresults	
				})#sapply: every day and replication
			# dailyReplication is a matrix: row is return and col is replication
			#row[1] # WRTDS prediction
			#row[2] # coefficent beta2
			#row[3] # r2
			
			return <- c(dailyReplication[1,], dailyReplication[2,], dailyReplication[3,], dailyReplication[4,])
		}#ifelse
	})# sapply
	# prediction is a matrix: col is daily; row is [1:replication] [1:replication] [1:replication]
	
	commonDates = intersectDate(list(obsData$date, predData$date))
	obs.dtsMatch = match(commonDates, obsData$date)
	pred.dtsMatch = match(commonDates, predData$date)
	nse = sapply(1:replicationN,function(ii){ NSE(obsData$C[obs.dtsMatch], prediction[ii,pred.dtsMatch]) })
	return <- list(
		conc = prediction[1:replicationN,],
		beta2 = prediction[(1:replicationN)+replicationN,],
		seasonal = prediction[(1:replicationN)+2*replicationN,],
		r2 = prediction[(1:replicationN)+3*replicationN,],
		NSE = nse,
		NSEmean = mean(nse),
		NSEsd = sd(nse))
}# function

 		
	




