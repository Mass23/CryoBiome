#!/usr/bin/Rscript

# https://github.com/jlw-ecoevo/gRodon

############################## LOG
sink(file=file(snakemake@log$out, open="wt"), type=c("output", "message"))

############################## LIBS
suppressMessages(library(devtools))
if(!"gRodon" %in% installed.packages()[,"Package"]){
    devtools::install_github("jlw-ecoevo/gRodon", dependencies=FALSE)
}
write(sprintf("Done: %s", Sys.time()), file=snakemake@output$done, append=TRUE)
