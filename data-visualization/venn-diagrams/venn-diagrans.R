library(ggVennDiagram)
library(ggplot2)
library(ggvenn)

data <- c(
    dd = read.delim("./gene-ontology/DENV_infected_downreg_gene_names.txt", header=FALSE),
    du = read.delim("./gene-ontology/DENV_infected_upreg_gene_names.txt", header=FALSE),
    zd = read.delim("./gene-ontology/ZIKV_infected_downreg_gene_names.txt", header=FALSE),
    zu = read.delim("./gene-ontology/ZIKV_infected_upreg_gene_names.txt", header=FALSE)
)

names(data) <- c(
    "Dengue downregulated",
    "Dengue upregulated",
    "Zika downregulated",
    "Zika upregulated")
    

ggvenn(
    data,
    show_percentage = FALSE,
    columns = c(
        "Dengue downregulated",
        "Zika downregulated"),
        fill_color = c("#bc5090", "#ffa600"))
ggvenn(
    data,
    show_percentage = FALSE,
    columns = c(
        "Dengue upregulated",
        "Zika upregulated"),
        fill_color = c("#bc5090", "#ffa600"))

ggvenn(
    data,
    show_percentage = FALSE,
    columns = c(
        "Dengue downregulated",
        "Zika upregulated"),
        fill_color = c("#bc5090", "#ffa600"))

ggvenn(
    data,
    show_percentage = FALSE,
    columns = c(
        "Dengue upregulated",
        "Zika downregulated"),
        fill_color = c("#bc5090", "#ffa600"))

if(file.exists("./Rplots.pdf")){
	file.rename(from="./Rplots.pdf", to="./data-visualization/venn-diagrams/Rplots.pdf")
}