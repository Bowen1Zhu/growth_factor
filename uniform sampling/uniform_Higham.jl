using Random, Printf, JLD2

function growth_factor_no_pivot(A)
    A_copy = copy(A)
    N = size(A_copy, 1)
    max_A = maximum(abs.(A_copy))
    max_A_k = max_A
    
    for k in 1:(N-1)
        for i in (k+1):N
            factor = A_copy[i, k] / A_copy[k, k]
            A_copy[i, k:N] .-= factor * A_copy[k, k:N]
            max_val_i = maximum(abs.(A_copy[i, :]))
            if max_val_i > max_A_k
                max_A_k = max_val_i
            end
        end
    end
    ANN = A_copy[N, N]
    return ANN / max_A
end

function uniform_higham(N)
    A = zeros(N, N)
    A[:, N] .= 1
    for k in 1:(N-1)
        local counter = 1
        while true
            for i in 1:k
                A[i, k] = 2*rand() - 1  # uniform in [-1, 1]
            end
            A_ik = -sum(2^(k-i) * A[i, k] for i in 1:k)
            if A_ik >= -1 && A_ik <= 1
                A[(k+1):N, k] .= A_ik
                println(k, "-th column: Success after ", counter, " attempts.") 
                break
            end
            counter += 1
        end
    end
    return A
end

# Parameters
sample_num = 40
N = 25

matrices = Array{Matrix{Float64}, 1}(undef, sample_num)
start_time = time()
for i in 1:sample_num
    println("Generating Matrix ", i, " ...") 
    A = uniform_higham(N)
    @assert growth_factor_no_pivot(A) == 2.0^(N-1) "The growth does not match the expected value."
    matrices[i] = A
    println()
end
time_elapsed = time() - start_time
@printf("Time elapsed: %.2f seconds\n\n", time_elapsed) 

# save data
filename = "higham_matrices_N$(N)_n$(sample_num).jld2"
JLD2.@save filename matrices
println("Samples saved to " * filename) 

# load data by
# data = JLD2.jldopen(filename, "r") do file
#     read(file, "matrices")
# end
