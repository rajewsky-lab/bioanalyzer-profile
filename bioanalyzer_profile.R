require(plyr)

#
# Functions
#
ReadBioanalyzerCSV = function(x) {
	# Read CSV and do some fixing on the Data Frame
	csv <- read.csv(x, skip = 17, stringsAsFactors = F)
	csv <- csv[seq(1,nrow(csv)-2,1),]
	csv <- as.data.frame(apply(csv, 2, as.numeric))
	csv <-csv[!is.na(csv$Time), ]
	return(csv)
}
FixAndSplit = function(x) { 
				# Fix that nasty comma in the quote value
				strsplit(
					gsub("\"(\\d+),(\\d+.*)\"","\\1\\2", x, perl=TRUE),
			    ",") 
				}

ReadLadderScale = function(f) {
	# Read file, use a Regex to find the right table, split by lines
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

CheckScale = function(l.p, l.s) {
	
# 	x11(width=20, height=20, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')
	par(pty="s", mfrow=c(1,2))

	### Plot time against size. 
	plot(l.p, y=l.s, xlab="Ladder Positions", ylab="Ladder Sizes")

	### The nonlinear relationship has to be found.
	### In this case logistic function gives a good fit
	plot(l.p, y=plogis(l.s, scale=300), xlab="Ladder Positions", ylab="Logistic function")
	cr = cor(plogis(l.s, scale=300), l.p)
	legend("topleft", paste("Corr = ", cr, sep=""), bty="n")
	
	### Fit a model to transform between sizes and time
	fit <- glm(plogis(l.s, scale=300) ~ l.p) 
	return(fit)
}

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

