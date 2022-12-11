using DataFrames
using CSV
using DelimitedFiles
using Statistics

function create_gene_tpm_file(gene_files, output_file, rounding_digits)
    gc_df = DataFrame()
    unique_g_ids = DataFrame(CSV.File(gene_files[1][2][1])).gene_id
    gc_df.gene_id = unique_g_ids
    unique_g_ids = nothing
    for condition in gene_files
        name, condition_collection = condition
        gene_tpm_col = nothing
        for file_path in condition_collection
            samp_df = DataFrame(CSV.File(file_path))
            if isnothing(gene_tpm_col)
                gene_tpm_col = samp_df.TPM
            else
                gene_tpm_col = samp_df.TPM .+ gene_tpm_col
            end
        end
        if !(rounding_digits == -1)
            gene_tpm_col = round.(gene_tpm_col, digits = rounding_digits)
        end
        gc_df[:, name] = gene_tpm_col
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

function norm_log2(input_file, output_file)
    data = DataFrame(CSV.File(input_file))
    numerical_data = data[:, 2:end]
    out_numerical_data = log.(2, numerical_data .+ 1) .- mean.(eachrow(log.(2, numerical_data .+ 1)))
    out_data = copy(data)
    out_data[:, 2:end] = out_numerical_data
    CSV.write(output_file,out_data)
end

function keep_cross_reference(data_input_file, reference_input_file, output_file)
    in_data = DataFrame(CSV.File(data_input_file))
    ref_data = readdlm(reference_input_file)
    out_data = in_data[findall(in(ref_data), in_data.gene_id), :]
    CSV.write(output_file,out_data)
end

create_gene_tpm_file(create_input("./gene_annotations.csv"), "./pre_processing_gene_tpm.csv", 0)

remove_past_dec("./pre_processing_gene_tpm.csv", "./in_processing_gene_tpm_1.csv", 16)

drop_duplicates("./in_processing_gene_tpm_1.csv", "./in_processing_gene_tpm_2.csv", "gene_id")

norm_log2("./in_processing_gene_tpm_2.csv", "./in_processing_gene_tpm_3.csv")

keep_cross_reference("./in_processing_gene_tpm_3.csv", "./gene-ontology/all_sig_genes.txt", "./final_gene_tpm.csv")