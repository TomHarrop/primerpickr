#!/usr/bin/env Rscript

# first header column names
header_cn <- c("Labware",
               "Src.Barcode",
               "Src.List Name",
               "Dest.Barcode",
               "Dest.List Name",
               NA,
               NA,
               NA)

# second header column names
header_main <- c("Barcode ID",
                 "Labware",
                 "Source",
                 "Labware",
                 "Destination",
                 "Volume",
                 "Tool",
                 "Name")

# The plate names are irrelevant so leave them blank (NA)
# `cbind` this header to the table after position matching
full_header <- transpose(data.table(header_cn,
                                    c(1, NA, NA, NA, NA, NA, NA, NA),
                                    c(2, NA, NA, NA, NA, NA, NA, NA),
                                    c(3, NA, NA, NA, NA, NA, NA, NA),
                                    c(4, NA, NA, NA, NA, NA, NA, NA),
                                    c(5, NA, NA, NA, NA, NA, NA, NA),
                                    c(6, NA, NA, NA, NA, NA, NA, NA),
                                    header_main))
