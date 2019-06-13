
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
					booststrap = sample(obsData$index[cond],sum(cond),replace=T)
					result = lm(logC~ydecimal+logQ+sin2pit+cos2pit, data=obsData[booststrap,], weights= Tweight[booststrap])
					# log(obsData$no3) = beta0 + beta1*obsData$ydecimal + beta2*obsData$logQ + beta3*obsData$sin2pit + beta4*obsData$cos2pit
					# new features coming: visual beta2 A:: {y:Q X:time}; countour beta2 by colors
					# B:: monthly boxplot of beta2; C::percentile Q boxplot of beta2
				
					alpha = sum(exp(result$residuals)*Tweight[booststrap])/sum(Tweight[booststrap])
					return <- c(
						alpha * exp( result$coefficients[1] + sum(result$coefficients[2:5]*predData[i,c('ydecimal','logQ','sin2pit','cos2pit')]) ),
						result$coefficients[3], # coefficent for the log(Q)
						sum(result$coefficients[4:5]*predData[i,c('sin2pit','cos2pit')]),
						summary(result)$r.squared)
				})#sapply: every day and replication
			# dailyReplication is a matrix: row is return and col is replication
			#row[1] # WRTDS prediction
			#row[2] # coefficent beta2
			#row[3] # r2
			
			return <- c(dailyReplication[1,], dailyReplication[2,], dailyReplication[3,], dailyReplication[4,])
		}#ifelse
	})# sapply
	# prediction is a matrix: col is daily; row is [1:replication] [1:replication] [1:replication]
	return <- list(
		conc = prediction[1:replicationN,],
		beta2 = prediction[(1:replicationN)+replicationN,],
		seasonal = prediction[(1:replicationN)+2*replicationN,],
		r2 = prediction[(1:replicationN)+3*replicationN,])
}# function

 		
	




