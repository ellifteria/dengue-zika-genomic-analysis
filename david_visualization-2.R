library(ggplot2)

dupli_df <- subset(read.csv("gene-ontology/dupl_genes_fach.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & !grepl("PTM", Category, fixed = TRUE))

dupli_df <- dupli_df[order(dupli_df$PValue),]

for(n in 1:nrow(dupli_df)){
    dupli_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",dupli_df[n, "Term"])))
    dupli_df[n, "Category"] = switch(
        dupli_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        dupli_df[n, "Category"]
    )
}

ggplot(dupli_df[1:min(5, nrow(dupli_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Overlapping upregulated RNAs: Most\nsignificant gene functions\n(ordered by PValue)"
    ) +
    scale_fill_manual(values=c("#003f5c",
                                "#bc5090",
                                "#ffa600"))

# ~~~~~~~~ BP ONLY ~~~~~~~~

dupli_df <- subset(read.csv("gene-ontology/dupl_genes_fach.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & (grepl("BP", Category, fixed = TRUE ) | grepl("BIOLOGICAL_PROCESS", Category, fixed = TRUE)))

dupli_df <- dupli_df[order(dupli_df$PValue),]

for(n in 1:nrow(dupli_df)){
    dupli_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",dupli_df[n, "Term"])))
    dupli_df[n, "Category"] = switch(
        dupli_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        dupli_df[n, "Category"]
    )
}


ggplot(dupli_df[1:min(5, nrow(dupli_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment)) +
    geom_bar(stat = "identity", fill="#003f5c") +
    coord_flip() +
    labs(
        x = "Term",
        y = "Fold enrichment",
        title = "Overlapping upregulated RNAs: Most\nsignificant gene biological processes\n(ordered by PValue)"
    )
