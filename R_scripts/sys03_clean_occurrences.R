# SYS-03  Cleaning specimen occurrence records: coordinates, duplicates, and names. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript sys03_clean_occurrences.R [records.csv reference_names.txt xmin xmax ymin ymax]
# records.csv columns: species, lon, lat, date, collector.  reference_names.txt: one accepted species name per line.
# xmin..ymax is the bounding box of your study region (longitude and latitude, decimal degrees).
# With no arguments the script builds a messy SYNTHETIC table (bounding box 116 to 127 E, 5 to 20 N) so it runs offline.
set.seed(404)

make_messy <- function() {
  ref <- c("Rattus tanezumi", "Rattus exulans", "Sus philippensis", "Cynopterus brachyotis", "Macaca fascicularis", "Tarsius syrichta", "Crocidura beatus", "Suncus murinus")
  n <- 300; d <- data.frame(species = sample(ref, n, TRUE), lon = runif(n, 118, 126), lat = runif(n, 7, 18), date = sample(seq(as.Date("2000-01-01"), as.Date("2020-12-31"), by = "day"), n, TRUE),
                            collector = sample(c("Cruz", "Reyes", "Santos", "Garcia"), n, TRUE), stringsAsFactors = FALSE)
  d$date <- as.character(d$date)
  d[sample(n, 8), c("lon", "lat")] <- 0                                            # zero coordinates
  sw <- sample(n, 6); tmp <- d$lon[sw]; d$lon[sw] <- d$lat[sw]; d$lat[sw] <- tmp    # swapped lon and lat
  d[sample(n, 10), "lon"] <- runif(10, 60, 80)                                     # far outside the region
  d[sample(n, 5), c("lon", "lat")] <- NA
  dup <- d[sample(n, 20), ]; d <- rbind(d, dup)                                     # exact duplicates
  typo <- sample(nrow(d), 25); d$species[typo] <- sub("Rattus", "Ratus", d$species[typo]); d$species[sample(nrow(d), 6)] <- "Sus philipensis"   # misspelled names
  list(data = d, ref = ref) }

clean_records <- function(d, ref, box) {
  log <- data.frame(step = "start", removed = 0, remaining = nrow(d))
  add <- function(step, keep) { log <<- rbind(log, data.frame(step = step, removed = sum(!keep), remaining = sum(keep))); d <<- d[keep, ] }
  add("missing coordinates", !is.na(d$lon) & !is.na(d$lat))
  add("zero coordinates (0, 0)", !(d$lon == 0 & d$lat == 0))
  sw <- d$lon >= box[3] & d$lon <= box[4] & d$lat >= box[1] & d$lat <= box[2]   # would fit if swapped
  d[!(d$lon >= box[1] & d$lon <= box[2] & d$lat >= box[3] & d$lat <= box[4]) & sw, c("lon", "lat")] <- d[!(d$lon >= box[1] & d$lon <= box[2] & d$lat >= box[3] & d$lat <= box[4]) & sw, c("lat", "lon")]
  log <- rbind(log, data.frame(step = "swapped lon/lat corrected (not removed)", removed = 0, remaining = nrow(d)))
  add("outside study region", d$lon >= box[1] & d$lon <= box[2] & d$lat >= box[3] & d$lat <= box[4])
  add("exact duplicates", !duplicated(d[, c("species", "lon", "lat", "date", "collector")]))
  # Name matching: exact, otherwise the closest reference name if one is within 10% edit distance.
  fix <- function(nm) { if (nm %in% ref) return(nm); dd <- adist(nm, ref)[1, ] / pmax(nchar(nm), nchar(ref)); if (min(dd) <= 0.10) ref[which.min(dd)] else NA }
  d$species_clean <- sapply(d$species, fix); n_fixed <- sum(d$species != d$species_clean, na.rm = TRUE)
  log <- rbind(log, data.frame(step = sprintf("names corrected by approximate matching (%d records changed)", n_fixed), removed = 0, remaining = nrow(d)))
  add("name not matched to the reference list", !is.na(d$species_clean))
  list(clean = d, log = log) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 6) { d <- read.csv(args[1], stringsAsFactors = FALSE); ref <- readLines(args[2]); box <- as.numeric(args[3:6]) } else { m <- make_messy(); d <- m$data; ref <- m$ref; box <- c(116, 127, 5, 20); cat("Practice data: synthetic table with planted errors\n") }
  raw_S <- length(unique(d$species)); r <- clean_records(d, ref, box)
  cat("\nRecords removed at each step:\n"); print(r$log, row.names = FALSE)
  cat(sprintf("\nSpecies count before cleaning (as written): %d; after cleaning: %d\n", raw_S, length(unique(r$clean$species_clean))))
  rng <- function(x) c(lon = diff(range(x$lon)), lat = diff(range(x$lat))); cat("Extent (degrees) before cleaning, non-missing coordinates:\n"); print(round(rng(d[complete.cases(d[, c("lon", "lat")]), ]), 1))
  cat("Extent after cleaning:\n"); print(round(rng(r$clean), 1))
  cat(sprintf("Share of the original records kept: %.1f%%\n", 100 * nrow(r$clean) / nrow(d)))
}
