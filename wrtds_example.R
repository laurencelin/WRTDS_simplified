
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R')
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_misc.r')
source('https://raw.githubusercontent.com/laurencelin/WRTDS_simplified/master/wrtds.R')


 	
##-----------------------------------------------------------------------------------##
##	for examples (below), using BES water chemistry data and USGS flow data
##-----------------------------------------------------------------------------------##
	
	site = c('BARN','DRKR','GFCP','GFGB','GFGL','GFuGR','GFVN','GRGF','MAWI','MCDN','POBR','RGHT')
	
	## read in water chemistry table; caution of missing and unexpected data entries.
	## make sure it has columns : year, month, and day.  you can make these easily on excel
	w = read.csv(paste('waterChem_',site[5],'.csv',sep=''))
	w$date = as.Date(paste(w[,'day'], w[,'month'], w[,'year'],sep="-"),format="%d-%m-%Y");
	w$goodEntryCond = sapply(w$no3,can.be.numeric) ## <<---- use [NO3] here;
	w$targetChem = rep(NA,dim(w)[1])
		w$targetChem[w$goodEntryCond] = as.numeric(w$no3[w$goodEntryCond])
		w$goodEntryCond = !is.na(w$targetChem) & w$targetChem >0 # has to be larger than zero for log transformation
		
	## read in USGS flow data (cfs); caution of missing and unexpected data entries.
	## make sure it has columns : year, month, and day.  you can make these easily on excel
	usgs = read.csv(paste('waterChem_',site[5],'_q.csv',sep=''), na.strings = "NA", stringsAsFactors=F)
	usgs$date = as.Date(paste(usgs[,'day'], usgs[,'month'], usgs[,'year'],sep="-"),format="%d-%m-%Y");
	usgs$goodEntryCond = sapply(usgs$cfs,can.be.numeric)
	usgs$m3ps = rep(NA,dim(usgs)[1])
		usgs$m3ps[usgs$goodEntryCond] = as.numeric(usgs$cfs[usgs$goodEntryCond])*0.0283168
		usgs$goodEntryCond = !is.na(usgs$m3ps) & usgs$m3ps>0
		
	## match flow dates to the water chemistry dates
	commonDates = intersectDate(list(w$date[w$goodEntryCond], usgs$date[usgs$goodEntryCond]))
	w.dtsMatch = match(commonDates, w$date)
	usgs.dtsMatch = match(commonDates, usgs$date)
	
	## combine chemistry and flow data together to form "obsData" to build the WRTDS model
	obsData = data.frame(date=w$date[w.dtsMatch])
	obsData$C = w$targetChem[w.dtsMatch]
	obsData$Q = usgs$m3ps[usgs.dtsMatch]
	obsData$logC = log(w$targetChem[w.dtsMatch])
	obsData$logQ = log(usgs$m3ps[usgs.dtsMatch])
	obsData$index = seq_len(length(w.dtsMatch))
	obsData$ydecimal = dailyTimeSeries(obsData $date)$ydecimal
	obsData$sin2pit = sin(2*pi* obsData$ydecimal)
	obsData$cos2pit = cos(2*pi* obsData$ydecimal)
	
	## "predData" is the flow data as inputs to the built-model to generate predictions (concentrations)
	## make sure it has columns : year, month, and day.  you can make these easily on excel
	## "predData" could be the same site that used to build the model, i.e., predData = usgs(above)
	## "predData" could be a different site with similar LULC char. and flow char.
	targetSite = read.csv('usgs01589290.csv') # SLB tributary
	targetSite$goodEntryCond = sapply(targetSite$siteLs,can.be.numeric)
	targetSite$m3ps = rep(NA,dim(targetSite)[1])
		targetSite$m3ps[targetSite$goodEntryCond] = as.numeric(targetSite$siteLs[targetSite$goodEntryCond])*0.001 # L -> m3
		targetSite$goodEntryCond = !is.na(targetSite$m3ps) & targetSite$m3ps >0
		
	predData = targetSite[targetSite$goodEntryCond,]
	predData$logQ = log(predData$m3ps) 
	predData$date = as.Date(paste(predData[,'day'], predData[,'month'], predData[,'year'],sep="-"),format="%d-%m-%Y"); 
	predData$ydecimal = dailyTimeSeries(predData$date)$ydecimal
	predData$sin2pit = sin(2*pi* predData$ydecimal)
	predData$cos2pit = cos(2*pi* predData$ydecimal)

	## run the model: output is a matrix of replication	(row) X data length
	output = WRTDS(obsData, predData)	
	
		## write output to file
		hold = data.frame(year=targetSite $year)
		hold$month = targetSite$month
		hold$day = targetSite$day
		for(jj in 1:50){
			hold[,paste('no3_',jj,sep='')] = output[jj,]
		}#i
		write.csv(hold,'usgs01589290_WRTDS_GFGL.csv',row.names=F)
		
	## plot predicted concentration over time	
	plot(predData$date, output[1,], type='n')
	for( jj in 1:50){
		lines(predData$date, output[jj,], col='gray')
	}#jj
	lines(predData$date, colMeans(output))
	
	
	
	




