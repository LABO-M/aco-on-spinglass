using CSV, DataFrames

# Convert the numeric value to a string format with K, M, G, etc.
function int_to_SI_prefix(value::Int)
    if value % 1_000_000_000 == 0
        return string(value ÷ 1_000_000_000) * "G"
    elseif value % 1_000_000 == 0
        return string(value ÷ 1_000_000) * "M"
    elseif value % 1_000 == 0
        return string(value ÷ 1_000) * "K"
    else
        return string(value)
    end
end

# Function to save Z values to a CSV file with optional downsampling
function save_Z_to_csv(Z_M::Vector{Float64}, filename::String, interval::Int=1000)
    t_values = collect(1:length(Z_M))
    M_downsampled = Z_M[t_values]
    df = DataFrame(t=t_values, M=M_downsampled)
    CSV.write(filename, df)
    println("Saved M values to $filename")
end