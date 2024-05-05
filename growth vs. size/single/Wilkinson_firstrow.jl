using Random, Printf, LinearAlgebra, Dates

function randn_()
    # return rand(1:10)//rand(11:20)
    denom = BigInt(rand(1:50))
    num = BigInt(rand(1:denom))
    return num // denom
end

function perturbation!(original_matrix, scale, N)
    original_matrix[1, :] .+= scale * randn(N)
end

function wilkinson_rational(N)
    A = zeros(Rational{BigInt}, N, N)
    for i in 1:N, j in 1:N
        if i == j
            A[i, j] = BigInt(1)
        elseif j == N
            A[i, j] = BigInt(1)
        elseif i > j
            A[i, j] = BigInt(-1)
        end
    end
    return A
end

function growth_factor_no_pivot(A)
    A_copy = copy(A)
    N = size(A_copy, 1)
    max_A = maximum(abs.(A_copy))
    max_A_k = max_A

    for k in 1:(N-1)
        # println("k = ", k - 1, ":")
        max_val = abs(A_copy[k, k])

        for i in (k+1):N
            factor = A_copy[i, k] / A_copy[k, k]
            # @printf("factor = %.32f\n", factor)
            A_copy[i, k:N] .-= factor * A_copy[k, k:N]

            max_val_i = maximum(abs.(A_copy[i, :]))
            if max_val_i > max_A_k
                max_A_k = max_val_i
            end
        end
        # @printf("After k = %d: max = %f, Growth = %f\n", k-1, max_A_k, max_A_k/max_A)
    end
    return max_A_k / max_A
end

####################################################

ε = BigInt(1)//BigInt(10000000000) # 1e-10
ε = 1e-8

sample_num = 1

N_values = 10:50

# Store growth factors
growth_factors_perturbed = Float64[]
growth_factors_original = Float64[]

for N in N_values
    A_original = wilkinson_rational(N)
    A_perturbed = copy(A_original)
    perturbation!(A_perturbed, ε, N)

    # Calculate growth factors
    growth_original = growth_factor_no_pivot(A_original)
    growth_perturbed = growth_factor_no_pivot(A_perturbed)

    push!(growth_factors_original, growth_original)
    push!(growth_factors_perturbed, growth_perturbed)
end

# Plotting both sets of growth factors
plot(N_values, growth_factors_original, label="Original Growth Factor", color=:blue, yscale=:log10,
     xlabel="Matrix Size N", ylabel="Growth Factor", title="Growth Factor vs. Matrix Size\nfor Wilkinson Matrix")

plot!(N_values, growth_factors_perturbed, label="Perturbed Growth Factor", color=:red, yscale=:log10)

savefig("wilkinson_" * string(N_values[end]) * "_first_row.png")
