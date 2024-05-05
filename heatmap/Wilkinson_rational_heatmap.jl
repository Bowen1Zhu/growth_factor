using Random, Printf, LinearAlgebra, Dates

function print_matrix(matrix)
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            @printf("%2.12f ", matrix[i, j])
        end
        println()
    end
    println()
end

function randn_()
    # return rand(1:10)//rand(11:20)
    denom = BigInt(rand(1:50))
    num = BigInt(rand(1:denom))
    return num // denom
end

function perturbation!(original_matrix, scale, i, j)
    # for i in eachindex(original_matrix)
    #     original_matrix[i] += scale * randn_()
    # end
    original_matrix[i, j] += scale
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
    ANN = A_copy[N, N]
    return ANN / max_A, ANN
end

####################################################
# Settings
Random.seed!(Int(time_ns()))

ε = BigInt(1)//BigInt(1000)
N = 25
A = Higham(N)
println("A = ")
print_matrix(A)

####################################################
# Now last pivot values instead of ratios
# Check ratios for all pairs of (i, j)
function calculate_observed_and_expected_ratios(i, j, A, ε)
    # Calculate the observed ratio
    A_perturbed = copy(A)
    perturbation!(A_perturbed, ε, i, j)
    # _, growth = growth_factor_no_pivot(A)
    _, growth_perturbed = growth_factor_no_pivot(A_perturbed)
    # observed_ratio = growth_perturbed / growth
    # observed_ratio = growth_perturbed
    # return observed_ratio
    return growth_perturbed
end

observed_ratios = Matrix{Float64}(undef, N, N)

# Iterate over all i, j pairs
for i in 1:N
    for j in 1:N
        observed_ratios[i, j] = calculate_observed_and_expected_ratios(i, j, A, ε)
    end
end

####################################################
# Print the results
function print_ratio_matrix(matrix::Matrix)
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            print(@sprintf("%9.6f ", matrix[i, j]))
        end
        println()
    end
    println()
end

println("Observed Ratios Matrix:")
print_ratio_matrix(observed_ratios)

println("1 - Observed Ratios Matrix:")
print_ratio_matrix(1 .- observed_ratios)

####################################################
# Plot the observed ratios
# import Pkg; Pkg.add("Plots"); Pkg.add("LaTeXStrings")
using Plots, LaTeXStrings
reversed_greys = cgrad(:greys, rev=true)
p = heatmap(log.(10, abs.(observed_ratios)),
        color=reversed_greys,
        # clims=(minimum(log.(10, abs.(observed_ratios))), maximum(log.(10, abs.(observed_ratios)))),
        clims=(0, 8),
        yflip=true,
        xticks=1:N,
        yticks=1:N,
        xlabel="j",
        ylabel="i",
        title="Last pivot perturbing",
        colorbar_title=L"$\log_{10} \left|\!\!p^{(i,j)}_\epsilon\right|$",
        aspect_ratio=:equal,
        grid=false,
        size=(600,623),
        tight_layout=true)

p = annotate!(p, N/2, -2, text("(i, j) element of Wilkinson matrix", 14, :center))

center_x = N * 0.5
center_y = N * 0.5

#
p = annotate!(p, center_x, center_y, "Perturbing the (1, n-1)th entry leads to", font(12, :white))
p = annotate!(p, center_x, center_y+1.5, L"growth $\approx$4000 (independent of n)", font(12, :white))
margin = 1.5
p = plot!(p, [center_x+1, N-margin], [center_y-1.5, margin], legend=false, color="white")
p = plot!(p, [N-margin, N-margin-0.3], [margin, margin], legend=false, color="white")
p = plot!(p, [N-margin, N-margin], [margin, margin+0.3], legend=false, color="white")

#
margin_x = 1.7
margin_y = 1.2
p = annotate!(p, center_x, margin_y, "Smallest last pivot", font(12, :white))
p = plot!(p, [center_x+5, N-margin_x], [margin_y, margin_y], legend=false, color="white")
p = plot!(p, [N-margin_x, N-margin_x-0.25], [margin_y, margin_y+0.25], legend=false, color="white")
p = plot!(p, [N-margin_x, N-margin_x-0.25], [margin_y, margin_y-0.25], legend=false, color="white")

savefig("/mnt/data/code/Plots_final/heatmap_products/heatmap_wilkinson_" * string(N) * ".png")
