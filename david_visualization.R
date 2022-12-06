library(ggplot2)

DENV_downreg_df <- subset(read.csv("gene-ontology/DENV_infected_downreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & !grepl("PTM", Category, fixed = TRUE))

DENV_downreg_df <- DENV_downreg_df[order(DENV_downreg_df$PValue),]

for(n in 1:nrow(DENV_downreg_df)){
    DENV_downreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",DENV_downreg_df[n, "Term"])))
    DENV_downreg_df[n, "Category"] = switch(
        DENV_downreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        DENV_downreg_df[n, "Category"]
    )
}

ggplot(DENV_downreg_df[1:min(5, nrow(DENV_downreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Dengue downregulated RNAs: Most\nsignificant gene functions\n(ordered by PValue)"
    ) +
    scale_fill_manual(values=c("#003f5c",
                                "#bc5090",
                                "#ffa600"))

DENV_upreg_df <- subset(read.csv("gene-ontology/DENV_infected_upreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & !grepl("PTM", Category, fixed = TRUE))

DENV_upreg_df <- DENV_upreg_df[order(DENV_upreg_df$PValue),]

for(n in 1:nrow(DENV_upreg_df)){
    DENV_upreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",DENV_upreg_df[n, "Term"])))
    DENV_upreg_df[n, "Category"] = switch(
        DENV_upreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        DENV_upreg_df[n, "Category"]
    )
}

ggplot(DENV_upreg_df[1:min(5, nrow(DENV_upreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Dengue upregulated RNAs: Most\nsignificant gene functions\n(ordered by PValue)"
    ) +
    scale_fill_manual(values=c("#003f5c",
                                "#bc5090",
                                "#ffa600"))

ZIKV_downreg_df <- subset(read.csv("gene-ontology/ZIKV_infected_downreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & !grepl("PTM", Category, fixed = TRUE))

ZIKV_downreg_df <- ZIKV_downreg_df[order(ZIKV_downreg_df$PValue),]

for(n in 1:nrow(ZIKV_downreg_df)){
    ZIKV_downreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",ZIKV_downreg_df[n, "Term"])))
    ZIKV_downreg_df[n, "Category"] = switch(
        ZIKV_downreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        ZIKV_downreg_df[n, "Category"]
    )
}

ggplot(ZIKV_downreg_df[1:min(5, nrow(ZIKV_downreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Zika downregulated RNAs: Most\nsignificant gene functions\n(ordered by PValue)"
    ) +
    scale_fill_manual(values=c("#ffa600",
                                "#bc5090"))

ZIKV_upreg_df <- subset(read.csv("gene-ontology/ZIKV_infected_upreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & !grepl("PTM", Category, fixed = TRUE))

ZIKV_upreg_df <- ZIKV_upreg_df[order(ZIKV_upreg_df$PValue),]

for(n in 1:nrow(ZIKV_upreg_df)){
    ZIKV_upreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",ZIKV_upreg_df[n, "Term"])))
    ZIKV_upreg_df[n, "Category"] = switch(
        ZIKV_upreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        ZIKV_upreg_df[n, "Category"]
    )
}

ggplot(ZIKV_upreg_df[1:min(5, nrow(ZIKV_upreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Zika upregulated RNAs: Most\nsignificant gene functions\n(ordered by PValue)"
    ) +
    scale_fill_manual(values=c("#003f5c",
                                "#ffa600",
                                "#bc5090"))

# ~~~~~~~~ BP ONLY ~~~~~~~~

DENV_downreg_df <- subset(read.csv("gene-ontology/DENV_infected_downreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & (grepl("BP", Category, fixed = TRUE ) | grepl("BIOLOGICAL_PROCESS", Category, fixed = TRUE)))

DENV_upreDENV_downreg_dfg_df <- DENV_downreg_df[order(DENV_downreg_df$PValue),]

for(n in 1:nrow(DENV_downreg_df)){
    DENV_downreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",DENV_downreg_df[n, "Term"])))
    DENV_downreg_df[n, "Category"] = switch(
        DENV_downreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        DENV_downreg_df[n, "Category"]
    )
}


ggplot(DENV_downreg_df[1:min(5, nrow(DENV_downreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment)) +
    geom_bar(stat = "identity", fill="#003f5c") +
    coord_flip() +
    labs(
        x = "Term",
        y = "Fold enrichment",
        title = "Dengue downregulated RNAs: Most\nsignificant gene biological processes\n(ordered by PValue)"
    )

DENV_upreg_df <- subset(read.csv("gene-ontology/DENV_infected_upreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & (grepl("BP", Category, fixed = TRUE ) | grepl("BIOLOGICAL_PROCESS", Category, fixed = TRUE)))

DENV_upreg_df <- DENV_upreg_df[order(DENV_upreg_df$PValue),]

for(n in 1:nrow(DENV_upreg_df)){
    DENV_upreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",DENV_upreg_df[n, "Term"])))
    DENV_upreg_df[n, "Category"] = switch(
        DENV_upreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        DENV_upreg_df[n, "Category"]
    )
}

ggplot(DENV_upreg_df[1:min(5, nrow(DENV_upreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment)) +
    geom_bar(stat = "identity", fill="#003f5c") +
    coord_flip() +
    labs(
        x = "Term",
        y = "Fold enrichment",
        title = "Dengue upregulated RNAs: Most\nsignificant gene biological processes\n(ordered by PValue)"
    )

ZIKV_downreg_df <- subset(read.csv("gene-ontology/ZIKV_infected_downreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & (grepl("BP", Category, fixed = TRUE ) | grepl("BIOLOGICAL_PROCESS", Category, fixed = TRUE)))

ZIKV_downreg_df <- ZIKV_downreg_df[order(ZIKV_downreg_df$PValue),]

for(n in 1:nrow(ZIKV_downreg_df)){
    ZIKV_downreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",ZIKV_downreg_df[n, "Term"])))
    ZIKV_downreg_df[n, "Category"] = switch(
        ZIKV_downreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        ZIKV_downreg_df[n, "Category"]
    )
}

ggplot(ZIKV_downreg_df[1:min(5, nrow(ZIKV_downreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment)) +
    geom_bar(stat = "identity", fill="#003f5c") +
    coord_flip() +
    labs(
        x = "Term",
        y = "Fold enrichment",
        title = "Zika downregulated RNAs: Most\nsignificant gene biological processes\n(ordered by PValue)"
    )

ZIKV_upreg_df <- subset(read.csv("gene-ontology/ZIKV_infected_upreg_gene_names_FACh.txt", sep="\t"), (grepl("GOTERM", Category, fixed = TRUE ) | grepl("UP_KW", Category, fixed = TRUE)) & (grepl("BP", Category, fixed = TRUE ) | grepl("BIOLOGICAL_PROCESS", Category, fixed = TRUE)))

ZIKV_upreg_df <- ZIKV_upreg_df[order(ZIKV_upreg_df$PValue),]

for(n in 1:nrow(ZIKV_upreg_df)){
    ZIKV_upreg_df[n, "Term"] = gsub("(.{10,}?)\\s", "\\1\n", sub(".*:", "", sub(".*~","",ZIKV_upreg_df[n, "Term"])))
    ZIKV_upreg_df[n, "Category"] = switch(
        ZIKV_upreg_df[n, "Category"],
        "GOTERM_BP_DIRECT" = "Biological process",
        "GOTERM_MF_DIRECT" = "Molecular function",
        "KEGG_PATHWAY" = "KEGG pathway",
        "GOTERM_CC_DIRECT" = "Cellular component",
        "UP_SEQ_FEATURE" = "UniProt sequence feature",
        "UP_KW_BIOLOGICAL_PROCESS" = "Biological process",
        "UP_KW_CELLULAR_COMPONENT" = "Cellular component",
        "UP_KW_MOLECULAR_FUNCTION" = "Molecular function",
        ZIKV_upreg_df[n, "Category"]
    )
}

ggplot(ZIKV_upreg_df[1:min(5, nrow(ZIKV_upreg_df)),], aes(x = reorder(Term, -1*PValue), y = Fold.Enrichment)) +
    geom_bar(stat = "identity", fill="#003f5c") +
    coord_flip() +
    labs(
        x = "Term",
        y = "Fold enrichment",
        title = "Zika upregulated RNAs: Most\nsignificant gene biological processes\n(ordered by PValue)"
    )


