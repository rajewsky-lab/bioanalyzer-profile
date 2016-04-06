### Read the raw data for the sample
bio.ds12 <- read.csv("/data/win/rajewsky/lab_organization/sequencing_info/bioanalyzer_files/temp_nikos/2100 expert_High Sensitivity DNA Assay_DE72901937_2016-03-10_17-21-25_Sample1.csv", skip = 17, stringsAsFactors = F)
bio.ds12$Time <- as.numeric(bio.ds12$Time)
bio.ds12$Value <- as.numeric(bio.ds12$Value)
bio.ds12 <- bio.ds12[!is.na(bio.ds12$Time), ]

### Plot to happily reproduce the bioanalyzer profile created by the machine
plot(bio.ds12$Time, bio.ds12$Value, type='l', col='red')

### Introduce ladder positions and sizes seen in the bioanalyzer pdf
ladder.positions <- c(43, 45.35, 50.95, 55.75, 60.40, 69.55, 77.6, 
                      83.35, 88.05, 91.25, 95.35, 101.60)
ladder.size <- c(35, 50, 100, 150, 200, 300, 400, 500, 600, 700, 1000, 2000)

### Plot time against size. 
plot(ladder.positions, y=ladder.size)

### The nonlinear relationship has to be found.
### In this case logistic function gives a good fit
plot(ladder.positions, y=plogis(ladder.size, scale=300))
cor(plogis(ladder.size, scale=300), ladder.positions)

### Fit a model to transform between sizes and time
fit <- glm(plogis(ladder.size, scale=300) ~ ladder.positions) 

### For instance "predict" the fragment sizes under the ~99% of the distribution
bio.ds12.subset <- bio.ds12[bio.ds12$Time >= 64.5 & bio.ds12$Time <= 96.5, 'Time']
newdata <- data.frame("ladder.positions"=bio.ds12.subset)
fragment.sizes <- round(as.numeric(qlogis(predict(fit, newdata), scale=300)))

### Plot them to see how they distribute
plot(fragment.sizes, bio.ds12[bio.ds12$Time %in% bio.ds12.subset, 'Value']  )









