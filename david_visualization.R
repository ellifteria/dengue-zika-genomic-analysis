library(ggplot2)

d_u_FAC = data.frame(
    term = c(
        "response to endoplasmic\nreticulum stress",
        "aminoacyl-tRNA\nsynthetase",
        "ATP binding",
        "ubiquitin protein\nligase binding",
        "protein processing in\nendoplasmic reticulum"
    ),
    enrichment = c(
        39.40408163265306,
        52.14269788182832,
        3.4894128525196484,
        8.80905695611578,
        11.017543859649123
    ),
    Category = c(
        "Biological process",
        "Molecular function",
        "Molecular function",
        "Molecular function",
        "KEGG pathway"
    ),
    order = c(
        1,
        2,
        3,
        4,
        5
    )
)

ggplot(d_u_FAC, aes(x = reorder(term, -1* order), y = enrichment, fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        legend = "Category",
        x = "Term",
        y = "Fold enrichment",
        title = "Dengue upregulated RNAs:\nMost enriched functional annotation clusters"
    )
