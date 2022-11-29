using DataFrames
using CSV

function create_gene_counts_file(gene_files, output_file, rounding_digits)
    gc_df = DataFrame()
    unique_g_ids = DataFrame(CSV.File(gene_files[1][2][1])).gene_id
    gc_df.gene_id = unique_g_ids
    unique_g_ids = nothing
    for condition in gene_files
        name, condition_collection = condition
        gene_count_col = nothing
        for file_path in condition_collection
            samp_df = DataFrame(CSV.File(file_path))
            if isnothing(gene_count_col)
                gene_count_col = samp_df.expected_count
            else
                gene_count_col = samp_df.expected_count .+ gene_count_col
            end
        end
        if !(rounding_digits == -1)
            gene_count_col = round.(gene_count_col, digits = rounding_digits)
        end
        gc_df[:, name] = gene_count_col
    end
    CSV.write(output_file, gc_df)
end

function remove_past_dec(input_file, output_file, pos_of_dec)
    df = DataFrame(CSV.File(input_file))
    df.gene_id = SubString.(string.(df.gene_id), 1, pos_of_dec - 1)
    CSV.write(output_file, df)
end

function create_input(ga_file)
    ga_df = DataFrame(CSV.File(ga_file))
    conditions = union(ga_df[:, "donor"] .* "_" .* ga_df[:, "infection_status"])
    num_conds = length(conditions)
    output = []
    for i in 1:num_conds
        cond = conditions[i]
        pivot = findfirst('_', cond)
        donor = cond[1: pivot - 1]
        infection = cond[pivot + 1 : end]
        file_paths = "./expression/" .* ga_df[(ga_df[:, "infection_status"] .== infection) .& (ga_df[:, "donor"] .== donor), :][:, "sample"] .* ".genes.results"
        push!(output, Vector([cond, file_paths]))
    end
    return Vector(output)
end

function drop_duplicates(input_file, output_file, col_name)
    CSV.write(output_file, unique(DataFrame(CSV.File(input_file)), col_name))
end

create_gene_counts_file(create_input("./gene_annotations.csv"), "./pre_processing_gene_counts.csv", 0)

remove_past_dec("./pre_processing_gene_counts.csv", "./in_processing_gene_counts.csv", 16)

drop_duplicates("./in_processing_gene_counts.csv", "./gene_counts.csv", "gene_id")