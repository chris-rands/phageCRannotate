library(Gviz)
library(stringr)

#Input directory format:
#$ head col_files_nonfiltered/id_#0N.NC_001416.1.N.colors 
#191|736|+|pink|terminase
#711|2636|+|pink|terminase
#2836|4437|+|maroon1|portal
#4418|5737|+|blue|head 
#6135|7160|+|blue|capsid
#7612|7965|+|purple|head-tail
#7977|8555|+|light blue| tail
#8552|8947|+|light blue| tail
#8955|9695|+|light blue| tail
#9711|10133|+|light blue| tail

# Get taxonomy and linked colour data
#source('R/gviz_colourData.R')

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
dir_glob <- args#[1]  # "col_files/*.colors"
#out_file <- args[-1]  # "gviz_Negatvicute_phages"

#get_genus_col <- function(in_file) {
#    return('grey')
#    genus <- rev(strsplit(in_file,"[.]")[[1]])[[2]]
#    order <- genus_lst[[genus]]
#    genus_col <- genus_colors_lst[[order]]
#    return(genus_col)
#}

get_annot_track <- function(in_file, strand, cex, names)
{
    # Get params from file name
    #name <- strsplit(basename(strsplit(in_file, "[.]")[[1]][1]), "[#]")[[1]][2]
    #genus_col <- get_genus_col(in_file)

    # Read data
    all_data <- read.table(in_file, sep = "|", blank.lines.skip = TRUE, col.names = c("Start", "End", "Strand", "Colour", "Annot"))
    
    # Get params
    if (strand) {
        strand <- all_data$Strand
    } else {
        strand <- "*"
    }
    print(all_data)
    # Create track
    aTrack <- AnnotationTrack(start = all_data$Start, end=all_data$End, strand=strand, stacking="dense", name="phage",
                              fill=as.vector(all_data$Colour), id=all_data$Annot, showFeatureId=names,
                              rotation.item=90, fontcolor.item="black", cex=cex)#), background.title=genus_col)
    return(aTrack)
}

run_gviz_plot <- function(dir_glob, out_file)
{
    gTrack <- GenomeAxisTrack()

    # Filter empty files out
    info <- file.info(dir_glob)
    files <- rownames(info[info$size != 0, ])

    # Sort files
    #genus_cols <- unlist(lapply(files, function(x) get_genus_col(x)))
    data_frame <- data.frame(files=files)#, genus_cols=genus_cols)
    #data_frame <- data_frame[order(genus_cols, decreasing=TRUE),]
    files <- as.vector(data_frame$files)

    # Set parameters based on number of genome maps
    len <- length(files)
    if (len <= 8) {
        strand <- TRUE
        width <- 14
        height <- 4
        names <- TRUE
        cex <- 0.6
    } else if (len <= 30) {
        strand <- TRUE
        width <- 14
        height <- 18
        names <- FALSE
        cex <- NULL
    } else if (len <= 50){
        strand <- TRUE
        width <- 14
        height <- 22
        names <- FALSE
        cex <- NULL
    } else if (len <= 100) {
        strand <- FALSE
        width <- 14
        height <- 18
        names <- FALSE
        cex <- NULL
    } else {
        strand <- FALSE
        width <- 14
        height <- 54
        names <- FALSE
        cex <- NULL
    }

    # Apply track functions
    tracks <- lapply(files, function(a, b, c, d) get_annot_track(a, b, c, d),
                     b=strand, c=cex, d=names)

    # Sort tracks
       
 
    # Combine with genome track
    tracks <- c(tracks, gTrack)

    pdf(file=paste(out_file, "pdf", sep="."), width=width, height=height)
    plotTracks(tracks, legend=TRUE)
    dev.off()
}

run_gviz_plot(Sys.glob(file.path(dir_glob)), "GENOME_MAP")

