using DelimitedFiles, CSV, DataFrames

gene_df = DataFrame(
    file = [
        "./gene-ontology/DENV_infected_downreg_gene_names.txt",
        "./gene-ontology/DENV_infected_upreg_gene_names.txt",
        "./gene-ontology/ZIKV_infected_downreg_gene_names.txt",
        "./gene-ontology/ZIKV_infected_upreg_gene_names.txt"
    ],
    infection = ["DENV", "DENV", "ZIKV", "ZIKV"],
    regulated = ["down", "up", "down", "up"],
    num_genes = [0, 0, 0, 0],
    genes = [Matrix{Any}(undef, 0, 0), Matrix{Any}(undef, 0, 0), Matrix{Any}(undef, 0, 0), Matrix{Any}(undef, 0, 0)]
)

for row in eachrow(gene_df)
    row.genes = readdlm(row.file)
    row.num_genes = length(row.genes)
end

println(gene_df[:, ["infection", "regulated", "num_genes"]])

for row1 in eachrow(gene_df)
    for row2 in eachrow(gene_df)
        if row1 != row2
            overlap = row1.num_genes + row2.num_genes - length(union(row1.genes, row2.genes))
            println("Overlap $(row1.infection) + $(row1.regulated) & $(row2.infection) + $(row2.regulated): $overlap")
        end
    end
end