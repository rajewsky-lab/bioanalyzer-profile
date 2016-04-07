source("~/bioanalyzer-profile/bioanalyzer_profile.R")

# Get Input files from a folder.
# it must contain at least one *Results.csv and 1 *Sample*.csv files
data = ReadBioanalyzerData("/data/win/rajewsky/lab_organization/sequencing_info/bioanalyzer_files/temp_nikos/")

# Analyze fitting of the scale
plotScale(data$Ladder$Aligned_Migration_Time, data$Ladder$Size)
# If it fits, it sitz
fit = fitScale(data$Ladder$Aligned_Migration_Time, data$Ladder$Size)


# Use the fit to adjust one sample
sample = data$Samples$Sample1

# Adjust the sample (Range is the subset of the distribution we want to adjust
# For instance = fragment sizes under the ~99% of the distribution
sample.adjusted = adjustSample(sample, fit, r=c(64.5,96.5))

### Plot original bioanalyzer distribution vs adjusted (subset only)
x11(width=17.5, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')
par(pty="s", mfrow=c(1,2))
plot(sample$Time, sample$Value, type='l', col='red', 
	 xlab="Aligned Migration Time", ylab="Value")

### Plot them to see how they distribute
plot(sample.adjusted$Position, sample.adjusted$Value, type='l', col='black', 
	 xlab="Fragment Size", ylab="Value")

