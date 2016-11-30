# INLCUDE ALL GENES
arg <- as.character(unique(read.table("/Volumes/Valeria/CRUK/AR_downstream_analysis/ARGs/list_args.txt",header=TRUE))[,1])
length(arg)

arg[which(arg=="ASB7i")] <- "ASB7"


# Which args passed the filters of the pipeline ----
args_passed <- intersect(arg,VCapProbes_genesTF_passed$VCap.Gene_symbol)
length(args_passed)/length(arg)


# Which args that passed the filters of the pipeline are returned by RCade ----

args_1500 <- intersect(args_passed,symbol_1500_ordered$SymbolReannotated)
length(args_1500)/length(arg)
args_5000 <- intersect(args_passed,symbol_5000_ordered$SymbolReannotated)
length(args_5000)/length(arg)
args_10000 <- intersect(args_passed,symbol_10000_ordered$SymbolReannotated)
length(args_10000)/length(arg)
args_25000 <- intersect(args_passed,symbol_25000_ordered$SymbolReannotated)
length(args_25000)/length(arg)

length(args_1500)/length(symbol_1500_ordered$SymbolReannotated)
length(args_5000)/length(symbol_5000_ordered$SymbolReannotated)
length(args_10000)/length(symbol_10000_ordered$SymbolReannotated)
length(args_25000)/length(symbol_25000_ordered$SymbolReannotated)
