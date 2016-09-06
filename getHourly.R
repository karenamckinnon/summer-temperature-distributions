# Note: this code was modified from one found online, but I can no longer find the 
# original source. Please let me know if you are the original source, so I can credit
# your work! Thanks! 
 
rm(list = ls())

yearStart = 1980
yearEnd = 2015
 
isdDir <- '/n/regal/huybers_lab/mckinnon/isd'
rawDir <- '/n/regal/huybers_lab/mckinnon/isd/raw'
csvDir <- '/n/regal/huybers_lab/mckinnon/isd/csv'


file <- "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv"
destFile = paste(isdDir, "/isd-history.csv", sep = "")
if (!file.exists(destFile)) {
  download.file(file, destFile)
}

namesList <- c("TOTAL.CHAR", "USAFID", "WBAN", "YR", "M",
	"D", "HR", "MIN", "FLAG", "LAT", "LONG", "CODE", "ELEV", "CALL",
	"QCTYPE", "WIND.DIR", "WIND.DIR.FLAG", "WIND.OBS.TYPE", "WIND.SPD",
	"WIND.SPD.FLAG", "CEIL.HGT", "CEIL.HGT.FLAG", "CEIL.HGT.METHOD",
	"CAVOK.CODE", "VISIBILITY.DIST", "VISIBILITY.DIST.FLAG", "VISIBILITY.VARIABLE",
	"VISIBILITY.VARIABLE.FLAG", "TEMP", "TEMP.FLAG", "DEWPOINT", "DEWPOINT.FLAG",
	"SLP", "SLP.FLAG")

varIndexSave <- c(2:11, 13, 29:34)
 
st <- read.csv(destFile)
names(st)[c(3, 9)] <- c("NAME", "ELEV")

st0 <- st

st <- st[!is.na(st$LAT), ]
st <- st[st$LAT > 30, ]

# st <- st[st$CTRY == "US", ] # only use US
st$BEGIN <- as.numeric(substr(st$BEGIN, 1, 4))
st$END <- as.numeric(substr(st$END, 1, 4))
st.full.time <- st[st$BEGIN <= yearStart & st$END >= yearEnd & !is.na(st$BEGIN), ]

# make a list of stations
if (RUNNUMBER == 1) {
	write.csv(st.full.time, file = paste(isdDir, "/isdList.csv", sep = ""), row.names = FALSE)
}

nstations <- dim(st.full.time)[1]
# calculate intervals
intervals <- ceiling(seq(1, nstations, length = TOTALCORES + 1))

for (s in intervals[RUNNUMBER]:(intervals[RUNNUMBER + 1] - 1)) { # 1:dim(st.full.time)[1]
  print(paste("Running station ", s, " of ", (intervals[RUNNUMBER + 1] - 1), sep = ""))
  for (y in yearStart:yearEnd) {
	# print(paste("Year ", y, sep = ""))
    outputs <- as.data.frame(matrix(NA, 1, 2))
    names(outputs) <- c("FILE", "STATUS")

	# name of file (USAF-WBAN)
	outputs[1] <- paste(sprintf("%06d", st.full.time[s, 1]), "-", sprintf("%05d", st.full.time[s, 2]), "-", y, sep = "")

	csvSaveName <- paste(csvDir, "/", outputs[1], ".csv", sep = "")

	if (!file.exists(csvSaveName)) {
		print(paste("Calculating ", outputs[1], sep = ""))
		# make wget command
		wget <- paste("wget -P ", rawDir, "/ ",  "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/", y, "/", outputs[1], ".gz", sep = "")

		# declare where file should go
		destFile <- paste(rawDir, "/", outputs[1], ".gz", sep = "")

		# name of unzipped file
		unzipFile <- paste(rawDir, "/", outputs[1], sep = "")

		if (!( file.exists(destFile) | file.exists(unzipFile)) ) {
		  outputs[2] <- try(system(wget, intern = FALSE, ignore.stderr = TRUE))
		} else {
		  outputs[2] <- 0 # file already there
		} 

		has.file <- outputs[2] == 0 & !is.na(outputs[2])

		if (sum(has.file) > 0) {
			print("Has data")
			files <- outputs[1]

			if (!file.exists(unzipFile)) {
		  		system(paste("gunzip ", destFile, sep = ""), intern = FALSE, ignore.stderr = TRUE)
		  		# file.remove(destFile)
			}

			column.widths <- c(4, 6, 5, 4, 2, 2, 2, 2, 1, 6,
							   7, 5, 5, 5, 4, 3, 1, 1, 4, 1, 5, 1, 1, 1, 6,
							   1, 1, 1, 5, 1, 5, 1, 5, 1)
			stations <- as.data.frame(matrix(NA, length(files), 5))
			names(stations) <- c("USAFID", "WBAN", "LAT", "LONG", "ELEV")
			for (i in 1:length(files)) {
			  data <- read.fwf(paste(rawDir, "/", files[i], sep = ""), column.widths)
			  data <- data[ , varIndexSave]
			  names(data) = namesList[varIndexSave]

			  data$LAT <- data$LAT/1000
			  data$LONG <- data$LONG/1000
			  # data$WIND.SPD <- data$WIND.SPD/10
			  data$TEMP <- data$TEMP/10
			  data$DEWPOINT <- data$DEWPOINT/10
			  data$SLP <- data$SLP/10

			  write.csv(data, file = csvSaveName, row.names = FALSE)
			}

		} else {
			print("no data")
		}
	} else {
		print(paste("Loading ", outputs[1], sep = ""))
	}

  }
}

