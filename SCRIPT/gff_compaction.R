# Check for optparse library to load arguments from command line ----------
if(suppressMessages(!require("optparse"))) {
  stop("R Package optparse was not found. Exiting.")
}

# Load parser -------------------------------------------------------------
library("optparse")
#Serie de instrucciones para crear el parseador de instrucciones y argumentos.
opt_list = list(make_option(opt_str = c("-f", "--file"), 
                            type = "character", 
                            default = NULL, 
                            help = "GFF file as obtained by blast6_to_gff.py", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-s", "--source"),
                            type = "character", 
                            default = NULL, 
                            help = "Method used to obtain the GFF file", 
                            metavar = "[STRING]"))

opt_parser = OptionParser(option_list = opt_list)
opt = parse_args(opt_parser)

# Check for arguments -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("A file must be submitted", call.=FALSE)
}

if (is.null(opt$source)) {
  print_help(opt_parser)
  stop("A source must be submitted", call.=FALSE)
}

# Check for packages ------------------------------------------------------
cat("-Loading packages...\n")

if(suppressMessages(!require("GenomicRanges"))) {
  stop("R Package GenomicRanges was not found. Exiting.")
} else {
  cat("-GenomicRanges loaded\n")
}

if(nzchar(system.file(package = "genomation"))) {
  cat("-genomation found\n")
} else {
  stop("genomation was not found. Exiting\n")
}

if(nzchar(system.file(package = "rtracklayer"))) {
  cat("-rtracklayer found\n")
} else {
  stop("rtracklayer was not found. Exiting\n")
}

# Program -----------------------------------------------------------------

cat("Proccesing data...\n")

filename <- options[[1]]
source <- options[[2]]
type <- "gene"

outname <- strsplit(filename, split = "\\.")[[1]]
outname <- paste0(append(outname, "reduced", after = length(outname) - 1), collapse = ".")

#Se lee el gff. Se utiliza la paquetería genomation
gr <- genomation::gffToGRanges(filename)

######Se los nombres únicos de scaffolds.
seqnames_unique <- unique(as.vector(seqnames(gr)))

list_temp_gr <- list()
c = 1
for (SCAFFOLD in seqnames_unique) {
  gr_scaffold <- gr[which(seqnames(gr) == SCAFFOLD)] #Se filtra por scaffold
  query_unique <- unique(mcols(gr_scaffold)[["Query"]]) #Se obtienen los Query únicos presentes en el scaffold de la iteracion
  for (QUERY in query_unique) {
    #Se genera un gff temporal donde se deposita los registros de cada query, sorteados por cadena y por rango
    gr_filtered <- sort(gr_scaffold[mcols(gr_scaffold)[["Query"]] == QUERY]) #Se obtiene los registros verdaderos para la columna Query y se ordenan por rango y strand
    gr_filtered <- reduce(gr_filtered)
    mcols(gr_filtered)$"source" <- source
    mcols(gr_filtered)$"type" <- type
    mcols(gr_filtered)$"score" <- NA
    mcols(gr_filtered)$"phase" <- NA
    mcols(gr_filtered)$"Query" <- QUERY
    list_temp_gr[[c]] <- gr_filtered
    c = c + 1
  }
}

gr_final <- unlist(as(list_temp_gr, "GRangesList"))

rtracklayer::export.gff(gr_final, outname, version = "3")


