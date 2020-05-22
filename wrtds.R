
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R')
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_misc.r')


WRTDS = function(obsData, predData, replicationN=50, Yhalfwin=10, Shalfwin=0.5, Qhalfwin=2){ 
	
    # obsData and predData are data frames with column names
    # A. must-have columns are:
    # 1) date   :  R Date object
    # 2) Q      :  discharge in m3/s
    # 3) Conc   :  concentration in mg/L = g/m3 (obsData only)
    # B. obsData should not have empty rows or missing-data rows.
    # C. users are responsible to pass quality data for the calculation below
    
    print("forming data ...")
    obsData$index = seq_len(dim(obsData)[1]);
    obsData$ydecimal = dailyTimeSeries(obsData$date)$ydecimal
    obsData$sin2pit = sin(2*pi* obsData$ydecimal)
    obsData$cos2pit = cos(2*pi* obsData$ydecimal)
    obsData$logQ = log(obsData$Q)
    obsData$logConc = log(obsData$Conc)
    
    predData$index = seq_len(dim(predData)[1]);
    predData$ydecimal = dailyTimeSeries(predData$date)$ydecimal
    predData$sin2pit = sin(2*pi* predData$ydecimal)
    predData$cos2pit = cos(2*pi* predData$ydecimal)
    predData$logQ = log(predData$Q)
    
    pred2ObsCOND = match(obsData$date, predData$date)
    predData$ObsQ = rep(NA,dim(predData)[1])
    predData$logObsQ = rep(NA,dim(predData)[1])
    predData$ObsQ[pred2ObsCOND] = obsData$Q
    predData$logObsQ[pred2ObsCOND] = obsData$logQ
    
    thresholdNUM = 100
    
    
    ## scanning model results
    print("begin to analysis ...")
    pb <- txtProgressBar(min = 0, max = dim(predData)[1], style = 3)
	prediction <- lapply(seq_len(dim(predData)[1]), function(i){ #
        
        # i is every day
        # every day we select a Yhalfwin_, Shalfwin_, and Qhalfwin_
        Yhalfwin_ = Yhalfwin
        Ydist = abs(obsData$ydecimal - predData$ydecimal[i])
        Yweight = (Ydist<=Yhalfwin)*(1-(Ydist/Yhalfwin)^3)^3
            while(sum(Yweight>0)<thresholdNUM){
                Yhalfwin_ = Yhalfwin_*1.1
                Yweight = (Ydist<= Yhalfwin_)*(1-(Ydist/Yhalfwin_)^3)^3
            }#while
        
        Shalfwin_ = Shalfwin
        Sdist = vectorMin( ceiling(Ydist)-Ydist, Ydist-floor(Ydist))
        Sweight = (Sdist <= Shalfwin)*(1-(Sdist/Shalfwin)^3)^3
            while(sum(Yweight>0 & Sweight>0)<thresholdNUM){
                Shalfwin_ = Shalfwin_*1.1
                Sweight = (Sdist <= Shalfwin_)*(1-(Sdist/Shalfwin_)^3)^3
            }#while
            
        Qhalfwin_ = Qhalfwin
        Qdist = abs(obsData$logQ - predData$logQ[i])
        Qweight = (Qdist <= Qhalfwin_)*(1-(Qdist/Qhalfwin_)^3)^3
            while(sum(Yweight>0 & Sweight>0 & Qweight>0)<thresholdNUM){
                Qhalfwin_ = Qhalfwin_*1.1
                Qweight = (Qdist <= Qhalfwin_)*(1-(Qdist/Qhalfwin_)^3)^3
            }# while
            
        Tweight = Yweight * Sweight * Qweight * (!is.na(obsData$logC)) * (!is.na(obsData$logQ)); # sum(Tweight>0)
        cond = Tweight>0; #sum(cond); sum(Yweight>0); sum(Sweight>0); sum(Qweight>0)
        
        dailyReplication <- do.call(rbind,lapply(seq_len(replicationN), function(gg){
                checkresults = NULL
                while( is.null(checkresults) ){
                    tryCatch({
                        #booststrap = rep(1,sum(cond))
                        #while( max(table(booststrap))>4 ){
                        #    # no data can be repeated more than three times
                        #    booststrap = sample(obsData$index[cond],sum(cond),replace=T)
                        #}#while
                        booststrap = sample(obsData$index[cond],size = 80)
                        result = lm(logConc~ydecimal+logQ+sin2pit+cos2pit, data=obsData[booststrap,], weights= Tweight[booststrap])
                        # log(obsData$no3) = beta0 + beta1*obsData$ydecimal + beta2*obsData$logQ + (beta3*obsData$sin2pit + beta4*obsData$cos2pit)
                        # new features coming: visual beta2 A:: {y:Q X:time}; countour beta2 by colors
                        # B:: monthly boxplot of beta2; C::percentile Q boxplot of beta2
                    
                        alpha = sum(exp(result$residuals)*Tweight[booststrap])/sum(Tweight[booststrap])
                        
                        predConc = alpha * exp( result$coefficients[1] + sum(result$coefficients[2:5]*predData[i,c('ydecimal','logQ','sin2pit','cos2pit')]) )# mg/L
                        predFlux = predConc * predData[i,'Q']*1000 # g/s -> mg/s
                        
                        predObsConc = alpha * exp( result$coefficients[1] + sum(result$coefficients[2:5]*predData[i,c('ydecimal','logObsQ','sin2pit','cos2pit')]) )# mg/L
                        predObsFlux = predObsConc * predData[i,'ObsQ']*1000 # g/s -> mg/s
                        
                        if(is.infinite(predConc) | predConc > 100*max(obsData$C[booststrap]) |  predConc < 0.01 *min(obsData$C[booststrap]) ){
                            # print(paste(
                                # alpha ,min(result$residuals),max(result$residuals),
                                # result$coefficients[1],
                                # result$coefficients[2],
                                # result$coefficients[3],
                                # result$coefficients[4],
                                # result$coefficients[5]));
                            checkresults <- NULL;
                        }else{
                            checkresults <- c(
                                conc = as.numeric(predConc),
                                flux = as.numeric(predFlux),
                                ycoef = as.numeric(result$coefficients[2]), # coefficent for year
                                qcoef = as.numeric(result$coefficients[3]), # coefficent for the log(Q)
                                season = as.numeric(sum(result$coefficients[4:5]*predData[i,c('sin2pit','cos2pit')])),
                                r2 = as.numeric(summary(result)$r.squared),
                                predObsC = as.numeric(predObsConc),
                                predObsF = as.numeric(predObsFlux));
                        }# end of if
                    },warning = function(w){
                        checkresults <- NULL
                    },error=function(e){
                        checkresults <- NULL
                    })# tryCatch
                }#while
                return <- checkresults
            }))#sapply: every day and replication
        # dailyReplication is a matrix: row is return and col is replication
        #row[1] # WRTDS conc prediction
        #row[2] # WRTDS flux predictions
        #row[3] # coefficent beta2
        #row[4] # seasonal signal
        #row[5] # r2
        
        setTxtProgressBar(pb, i)
        return <- list(
                    predictedConc = dailyReplication[,'conc'], # conc
                    predictedFlux = dailyReplication[,'flux'],
                    Ycoef = dailyReplication[,'ycoef'], # beta1
                    Qcoef = dailyReplication[,'qcoef'], # beta2
                    SeasonalSignal = dailyReplication[,'season'], # seasonal
                    CQr2 = dailyReplication[,'r2'], # r2
                    predictedObsConc = dailyReplication[,'predObsC'],
                    predictedObsFlux = dailyReplication[,'predObsF'],
                    actualYhalfwin = Yhalfwin_,
                    actualShalfwin = Shalfwin_,
                    actualQhalfwin = Qhalfwin_)
	})# lapply - daily loop -- prediction
	# prediction is a matrix: col is daily; row is [1:replication] [1:replication] [1:replication] 1 1 1
	# prediction is a list of daily "return list"
    
    # fittness: if "prediction" is based on "predData", which is not the same as "obsData", then the fittness is meaning less
    extracted_predictedObsConc = do.call(rbind, lapply(predData$index[pred2ObsCOND],function(jj){
        return <- prediction[[jj]]$predictedConc # vector of replication values
    }))
    extracted_predictedObsFlux = do.call(rbind, lapply(predData$index[pred2ObsCOND],function(jj){
        return <- prediction[[jj]]$predictedObsFlux # vector of replication values
    }))
	concNSE = sapply(1:replicationN,function(ii){ NSE(obsData$C, extracted_predictedObsConc[,ii]) })
    fluxNSE = sapply(1:replicationN,function(ii){ NSE(obsData$C*obsData$Q*1000, extracted_predictedObsFlux[,ii]) })
	
	# prediction is a list of daily list
	# report mean daily conc and upper and lower bound
	return <- list(
		ConcPrediction = do.call(rbind,lapply(seq_along(prediction),function(ii){
			return <- c(
				mConc = mean(prediction[[ii]]$predictedConc), 
				lConc = quantile(prediction[[ii]]$predictedConc,0.025),
				uConc = quantile(prediction[[ii]]$predictedConc,1-0.025) )})),
				
		FluxPrediction = do.call(rbind,lapply(seq_along(prediction),function(ii){
			return <- c(
				mFlux = mean(prediction[[ii]]$predictedFlux),
				lFlux = quantile(prediction[[ii]]$predictedFlux,0.025),
				uFLux = quantile(prediction[[ii]]$predictedFlux,1-0.025) )})),
		
		modelParam = do.call(rbind,lapply(seq_along(prediction),function(ii){
			return <- c(
				Ycoef = mean(prediction[[ii]]$Ycoef),
				Qcoef = mean(prediction[[ii]]$Qcoef),
				Season = mean(prediction[[ii]]$SeasonalSignal),
				Yhalfwin = mean(prediction[[ii]]$actualYhalfwin), 
				Shalfwin = mean(prediction[[ii]]$actualShalfwin), 
				Qhalfwin = mean(prediction[[ii]]$actualQhalfwin),
				dailyCQr2 = mean(prediction[[ii]]$CQr2) )})),
		
		modelFittiness = c(concNSE, fluxNSE)			
	)# end of return list 
	
}# function

 		
	




