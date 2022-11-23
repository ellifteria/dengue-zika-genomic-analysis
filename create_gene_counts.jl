using DataFrames
using CSV

function create_gene_count_file(gene_files, output_file)
    gc_df = DataFrame()
    unique_g_ids = nothing
    for file in gene_files
        file_path, ~ = file
        samp_df = DataFrame(CSV.File(file_path))
        unique_g_ids = union(samp_df[:, "gene_id"], unique_g_ids)
    end
    gc_df.gene_id = unique_g_ids
    unique_g_ids = nothing
    for file in gene_files
        file_path, sample_name = file
        samp_df = DataFrame(CSV.File(file_path))
        unique_g_ids = union(samp_df[:, "gene_id"])
        for gene in unique_g_ids
            g_count = sum(samp_df[samp_df.gene_id .== gene, "expected_count"])
            gc_df[gc_df.gene_id .== gene, sample_name] = g_count
        end
    end
    CSV.write(output_file, gc_df)
end