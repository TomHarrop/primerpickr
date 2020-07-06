#' MakePrimerPickingFile
#'
#' Make an epMotion-compatible csv for pipetting combinatorial primers into
#' library prep plates
#'
#' @param library_layout A \code{data.table} with columns \code{library_id},
#'   \code{i7_primer}, \code{"i5_primer"} and EITHER (\code{lib_row} AND
#'   \code{lib_col}) OR \code{lib_well}. If there are multiple plates of library
#'   preps, specify them in the \code{lib_plate} column. A separate file will be
#'   written for each plate.
#' @param primer_set Which set of barcodes to use. Only supports NEBNext
#'   Multiplex Oligos for Illumina so far (\code{neb}).
#' @param output_path Folder to write the output files
#' @param output_file Prefix for the output files. Output files will be named
#'   \code{output_file}.csv, or if there are multiple plates,
#'   \code{output_file}_{plate_name}.csv.
#'
#' @return Writes files to the \code{output_path}
#'
#' @export

MakePrimerPickingFile <- function(
    library_layout,
    primer_set = "neb",
    primer_volume = 5/3,
    output_path = ".",
    output_file) {

    # don't modify the original
    my_library_layout <- copy(library_layout)

    # which primer uses which layout
    primer_layouts <- list(
        "neb" = neb_multiplex_oligos_for_illumina)

    # check if we have a file for the primerset
    if (!primer_set %in% names(primer_layouts)) {
        stop(glue::glue("Primer {primer_set} not in primer_layouts"))
    } else {
        my_primer_layout <- primer_layouts[[primer_set]]
    }

    # check if we have lib_well and primer columns
    required_colnames <- c("library_id", "i7_primer", "i5_primer")
    for (colname in required_colnames) {
        if (!colname %in% names(my_library_layout)) {
            stop(glue::glue("library_layout missing column {colname}"))
        }
    }

    # set up lib well if missing
    position_colnames <- c("lib_row", "lib_col")
    if (!"lib_well" %in% names(my_library_layout)) {
        for (colname in position_colnames) {
            if (!colname %in% names(my_library_layout)) {
                stop(glue::glue(
                    "library_layout missing columns lib_well and {colname}"))
            }
        }
        message(glue::glue(
            "Detected {position_colnames[[1]]}, generating lib_well"))
        my_library_layout[, lib_well := paste0(lib_row, lib_col)]
    }

    # check if we have lib_plate column and if so how many plates
    if ("lib_plate" %in% names(my_library_layout)) {
        n_plates <- my_library_layout[, length(unique(lib_plate))]
    } else {
        n_plates <- 1
        my_library_layout[, lib_plate := 1]
    }

    message(glue::glue("Laying out {n_plates} library plate(s)"))

    # check for duplicate rows
    # this would result in two primer combos in the same well
    if (any(duplicated(my_library_layout,
                       by = c("lib_col", "lib_row", "lib_plate")))) {
        stop(glue::glue("Duplicated rows in library_layout"))
    }

    # melt the library_layout by id
    melt_vars <- c("library_id", "lib_well", "lib_plate")
    my_libs_long <- melt(my_library_layout,
                         id.vars = melt_vars,
                         measure.vars = c("i7_primer", "i5_primer"),
                         variable.name = "primer_colname",
                         value.name = "primer")

    # get the primer positions
    my_libs_with_primers <- merge(my_libs_long,
                                  my_primer_layout[, .(primer, primer_well = well)],
                                  by = c("primer"))

    # make sure we got every lib
    if (my_libs_with_primers[, any(is.na(primer_well))]) {
        print(my_libs_with_primers[is.na(primer_well)])
        stop("Unable to match primers")
    }

    # re-organise by primer position
    my_libs_with_primers[, sort_order := factor(
        primer_well,
        levels = gtools::mixedsort(unique(primer_well)))]

    # coerce lib_plate to factor if necessary
    if (my_libs_with_primers[, any(is.na(as.integer(lib_plate)))]) {
        my_libs_with_primers[, lib_plate := factor(
            lib_plate,
            levels = gtools::mixedsort(unique(lib_plate)))]
    }
    my_libs_with_primers[, lib_plate_int := as.integer(lib_plate)]

    # final sort is on destination position
    my_libs_with_primers[, final_sort := factor(
        lib_well,
        levels = gtools::mixedsort(unique(lib_well)))]
    setorder(my_libs_with_primers, sort_order, lib_plate_int, final_sort)

    # check for output directory
    if (!dir.exists(output_path)) {dir.create(output_path, recursive = TRUE)}

    # generate the layout ONCE PER PLATE FIXME
    all_plates <- my_libs_with_primers[, unique(lib_plate_int)]
    for (plate_num in all_plates) {
        my_name <- my_libs_with_primers[lib_plate_int == plate_num][1, lib_plate]
        my_filename <- glue::glue("{output_path}/{output_file}_{my_name}.csv")
        my_layout <- my_libs_with_primers[lib_plate_int == plate_num][
            ,
            .(`Barcode ID` = NA,
              Labware = 1,
              Source = primer_well,
              Labware = 1,
              Destination = lib_well,
              Volume = primer_volume,
              Tool = "TS_10",
              Name = library_id)]
        full_layout <- rbind(full_header, my_layout, use.names = FALSE)
        message(glue::glue(
            "Writing {my_filename}"
        ))
        fwrite(full_layout, my_filename, col.names = FALSE, eol = "\r\n")
    }

}
