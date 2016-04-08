require(plyr)

#
# Functions
#

# Read Sample CSV and do some fixing on the Data Frame
ReadBioanalyzerCSV = function(x) {
	csv <- read.csv(x, skip = 17, stringsAsFactors = F)
	csv <- csv[seq(1,nrow(csv)-2,1),]
	csv <- as.data.frame(apply(csv, 2, as.numeric))
	csv <-csv[!is.na(csv$Time), ]
	return(csv)
}

# Fix that nasty comma in the quote value
FixAndSplit = function(x) { 
		
		strsplit(
			gsub("\"(\\d+),(\\d+.*)\"","\\1\\2", x, perl=TRUE),
		",") 
}

# Read Results file, use a Regex to extract the right table, 
# split by lines and return a table
ReadLadderScale = function(f) {
	ladder.list = strsplit(
		gsub(".*Sample Name,Ladder.*Peak Table\r\n.*Time corrected area\r\n(.*)\r\n \r\nOverall.*","\\1",
			readChar(f, file.info(f)$size)
		), "\r\n"
		)[[1]]
	# Convert into a Data Frame after fixing the ugly comma-in-quotes issue
	ladder.table = ldply(lapply(llply(ladder.list, FixAndSplit), ldply))
	# Add colnames (not parsed due to encoding issues)
	colnames(ladder.table) = c("Size","Concentration",
							"Molarity","Observations",
							"Area","Aligned_Migration_Time",
							"Peak_Height","Peak_Width","PC_of_Total",
							"Time_corrected_area")
	
	ladder.table[,c(1,2,3,5,6,7,8,9,10)] = apply(ladder.table[,c(1,2,3,5,6,7,8,9,10)], 2, as.numeric)
	
	return(ladder.table)
}

# Plot original and hypothetical scale
plotScale = function(l.p, l.s, fit.scale=300) {
	par(pty="s", mfrow=c(1,2))

	### Plot time against size. 
	plot(l.p, y=l.s, xlab="Ladder Positions", ylab="Ladder Sizes")

	### The nonlinear relationship has to be found.
	### In this case logistic function gives a good fit
	plot(l.p, y=plogis(l.s, scale=fit.scale), xlab="Ladder Positions", ylab="Logistic function")
	cr = cor(plogis(l.s, scale=fit.scale), l.p)
	legend("topleft", paste("Corr = ", cr, sep=""), bty="n")

}

# Fit logistic function to the data
fitScale = function(ladder.positions, ladder.size, fit.scale=300) {		
	### Fit a model to transform between sizes and time
	fit <- glm(plogis(ladder.size, scale=fit.scale) ~ ladder.positions) 
	return(fit)
}

# Adjust a sample using the fitting
adjustSample = function(sp, fit, r, fit.scale=300){
	sample.subset   <- sp[sp$Time >= r[1] & sp$Time <= r[2], ]
	sample.df <- data.frame("ladder.positions" = sample.subset$Time)
	fragment.sizes <- round(as.numeric(qlogis(predict(fit, sample.df), scale=fit.scale)))
	return ( data.frame(Positions=fragment.sizes, Value=sample.subset$Value) )
}

# Wrapper using previous functions to read and parse bioanalyzer data
ReadBioanalyzerData = function(d) {
	# Read all files of the directory. 
	# Must contain: 1 *Results.csv
	# Must contain: at least 1 *Sample\d.csv	
	input.files = list.files(d, pattern="*.csv", full.names=T)
	
	# Get sample files and names
	sample.files = input.files[grepl("Sample", input.files)]
	sample.tags   = gsub(".*_(.*).csv","\\1", sample.files)

	# Parse the ladder data from the *Results.csv
	ladder.table = ReadLadderScale( input.files[grepl("Results", input.files)] )
	
	# Load Sample files into a named list
	sample.df = setNames(llply(sample.files, ReadBioanalyzerCSV), nm=sample.tags)

	return(list(Ladder=ladder.table, Samples=sample.df))
}

