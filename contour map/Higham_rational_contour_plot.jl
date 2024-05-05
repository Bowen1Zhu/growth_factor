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

function matrix_multiplication(A, B)
    n = size(A, 1)
    C = zeros(Rational{BigInt}, n, n)
    for i in 1:n
        for j in 1:n
            C[i, j] = BigInt(0)//BigInt(1)
            for k in 1:n
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
    return C
end

function nonsingular_upper_triangular(matrix::Matrix{Rational{BigInt}})
    N = size(matrix, 1)
    lower_bound = -BigInt(10)//BigInt(1)
    upper_bound = BigInt(10)//BigInt(1)

    for i in 1:N-1
        for j in 1:N-1
            if i == j
                matrix[i, j] = randn_() * (upper_bound - lower_bound) + lower_bound
            elseif i < j
                matrix[i, j] = randn_() * (upper_bound - lower_bound) + lower_bound
            else
                matrix[i, j] = BigInt(0)//BigInt(1)  # Set elements below the diagonal to 0
            end
        end
    end

    if det(matrix) == 0//1
        println("The generated matrix is singular. Trying again...")
        nonsingular_upper_triangular(matrix)  # Try again
    end

    return matrix
end

function Higham(N)
    M = [i > j ? BigInt(-1)//BigInt(1) : i == j ? BigInt(1)//BigInt(1) : BigInt(0)//BigInt(1) for i in 1:N, j in 1:N]
    T_theta_d_zero = zeros(Rational{BigInt}, N, N)

    for i in 1:N, j in 1:N
        if i < N && j < N
            T_theta_d_zero[i, j] = (i == j) ? BigInt(1)//BigInt(1) : BigInt(0)//BigInt(1)
        elseif i == N && j != N
            T_theta_d_zero[i, j] = BigInt(0)//BigInt(1)
        elseif j == N
            T_theta_d_zero[i, j] = (BigInt(2)//BigInt(1))^(BigInt((i-1))//BigInt(1))
        end
    end

    U = nonsingular_upper_triangular(T_theta_d_zero)
    A = M * U

    theta = maximum(abs.(A[:, 1:N-1]))
    U[:, N] *= theta

    A = matrix_multiplication(M, U)
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

function compute_ratio_by_case(ε, L_inv, U_inv, l, u, g, i, j)
    N = size(U_inv, 1) + 1
    ratio = 0

    U_inv_u = U_inv * u    # Product U^-1 * u
    lT_L_inv = l' * L_inv  # Product l^T * L^-1

    if i < N && j < N
        # Case 1: i, j < N
        ratio = 1 - (ε / g) * (U_inv_u * lT_L_inv)[j, i] / (1 + ε * (U_inv * L_inv)[j, i])
        println("Case 1")
    elseif i < N && j == N
        # Case 2: i < N, j = N
        ratio = 1 + (ε / g) * lT_L_inv[i]
        println("Case 2")
    elseif i == N && j < N
        # Case 3: i = N, j < N
        ratio = 1 - (ε / g) * U_inv_u[j]
        println("Case 3")
    elseif i == N && j == N
        # Case 4: i = N, j = N
        ratio = 1 + (ε / g)
        println("Case 4")
    end
    ratio
end

####################################################
# Random.seed!(0)

using Plots, LaTeXStrings

# Define ε and N
ε = BigInt(1)//BigInt(1000)
N = 10

# Generate the original Higham matrix
A = Higham(N)

i, j = 1, 1
# Define perturbation matrices P1 and P2
P1 = zeros(Rational{BigInt}, N, N)
P1[i, j] = ε
P2 = zeros(Rational{BigInt}, N, N)
# P2[i, N] = ε

for x in 1:N
    for y in 1:N
        P2[x, y] = randn_() * ε
    end
end
P2[i, j] = 0

# Define the range and step for x and y perturbations
####################################################
## Smallest
x_range = -0.01:0.0005:0.01
y_range = -0.01:0.0005:0.01

# Initialize the matrix to store the ratios
ratios = zeros(length(y_range), length(x_range))

for (x_idx, x) in enumerate(x_range)
    for (y_idx, y) in enumerate(y_range)
        # Perturb the original matrix
        A_perturbed = A + x*P1 + y*P2
        _, ANN = growth_factor_no_pivot(A)
        _, ANN_ = growth_factor_no_pivot(A_perturbed)

        # Calculate the ratio and store it
        ratios[y_idx, x_idx] = ANN_ / ANN
    end
end

# Generate the contour plot
p = contourf(x_range,
        y_range,
        log.(abs.(ratios)),
        levels=collect(-0.125:0.025:0.15),
        xlabel="Perturbation in (1, 1) direction",
        ylabel="Random perturbation without perturbing (1, 1)",
        title="Relative change in last pivot for a Higham matrix\nwhen perturbing (1,1) vs. random perturbation w/o (1,1)",
        colorbar_title=L"$\log_{10} \left|\!\!p_\epsilon / p\ \right|$",
        clabels=true,
        color=:thermal,
        aspect_ratio=:equal,
        size=(700,707))

p = annotate!(p, 0, -0.008, text("Nearly vertical slope implies (1, 1) direction\nis the most significant"))

savefig("/mnt/data/code/Plots/contour_$(N)_$(length(x_range))_points_smallest.png")

####################################################
## Small
x_range = -0.1:0.005:0.1
y_range = -0.1:0.005:0.1

# Initialize the matrix to store the ratios
ratios = zeros(length(y_range), length(x_range))

for (x_idx, x) in enumerate(x_range)
    for (y_idx, y) in enumerate(y_range)
        # Perturb the original matrix
        A_perturbed = A + x*P1 + y*P2
        _, ANN = growth_factor_no_pivot(A)
        _, ANN_ = growth_factor_no_pivot(A_perturbed)

        # Calculate the ratio and store it
        ratios[y_idx, x_idx] = ANN_ / ANN
    end
end

# Generate the contour plot
p = contourf(x_range,
        y_range,
        log.(abs.(ratios)),
        levels=collect(-2:1:8),
        xlabel="Perturbation in (1, 1) direction",
        ylabel="Random perturbation without perturbing (1, 1)",
        title="Relative change in last pivot for a Higham matrix\nwhen perturbing (1,1) vs. random perturbation w/o (1,1)",
        colorbar_title=L"$\log_{10} \left|\!\!p_\epsilon / p\ \right|$",
        clabels=true,
        color=:thermal,
        aspect_ratio=:equal,
        size=(700,707))

savefig("/mnt/data/code/Plots/contour_$(N)_$(length(x_range))_points_small.png")

####################################################
## Larger
x_range = -0.5:0.0125:0.5
y_range = -0.5:0.0125:0.5

# Initialize the matrix to store the ratios
ratios = zeros(length(y_range), length(x_range))

for (x_idx, x) in enumerate(x_range)
    for (y_idx, y) in enumerate(y_range)
        # Perturb the original matrix
        A_perturbed = A + x*P1 + y*P2
        _, ANN = growth_factor_no_pivot(A)
        _, ANN_ = growth_factor_no_pivot(A_perturbed)

        # Calculate the ratio and store it
        ratios[y_idx, x_idx] = ANN_ / ANN
    end
end

p = contourf(x_range,
        y_range,
        log.(abs.(ratios)),
        levels=collect(-2:1:8),
        xlabel="Perturbation in (1, 1) direction",
        ylabel="Random perturbation without perturbing (1, 1)",
        title="Relative change in last pivot for a Higham matrix\nwhen perturbing (1,1) vs. random perturbation w/o (1,1)",
        colorbar_title=L"$\log_{10} \left|\!\!p_\epsilon / p\ \right|$",
        clabels=true,
        color=:thermal,
        aspect_ratio=:equal,
        tight_layout=true,
        size=(700,707))

p = annotate!(p, 0.2, -0.35, text("Belt of the last pivot\nblowing up"))

x_coords0 = [-0.1, 0.1, 0.1, -0.1, -0.1]
y_coords0 = [-0.1, -0.1, 0.1, 0.1, -0.1]
x_coords1 = [-0.01, 0.01, 0.01, -0.01, -0.01]
y_coords1 = [-0.01, -0.01, 0.01, 0.01, -0.01]
x_coords2 = [0.09, 0.11, 0.11, 0.09, 0.09]
y_coords2 = [0, 0, 0.02, 0.02, 0]
plot!(p, x_coords0, y_coords0, line=:solid, color=:blue, linewidth=3, legend=false)
plot!(p, x_coords1, y_coords1, line=:solid, color=:green, linewidth=3, legend=false)
plot!(p, x_coords2, y_coords2, line=:solid, color=:red, linewidth=3, legend=false)

savefig("/mnt/data/code/Plots/contour_$(N)_$(length(x_range))_points_large.png")

####################################################
## Small and focused
x_range = 0.09:0.0005:0.11
y_range = 0:0.0005:0.02

# Initialize the matrix to store the ratios
ratios = zeros(length(y_range), length(x_range))

for (x_idx, x) in enumerate(x_range)
    for (y_idx, y) in enumerate(y_range)
        # Perturb the original matrix
        A_perturbed = A + x*P1 + y*P2
        _, ANN = growth_factor_no_pivot(A)
        _, ANN_ = growth_factor_no_pivot(A_perturbed)

        # Calculate the ratio and store it
        ratios[y_idx, x_idx] = ANN_ / ANN
    end
end

# Generate the contour plot
p = contourf(x_range,
        y_range,
        log.(abs.(ratios)),
        levels=collect(-2:0.5:8.22),
        xlabel="Perturbation in (1, 1) direction",
        ylabel="Random perturbation without perturbing (1, 1)",
        title="Relative change in last pivot for a Higham matrix\nwhen perturbing (1,1) vs. random perturbation w/o (1,1)",
        colorbar_title=L"$\log_{10} \left|\!\!p_\epsilon / p\ \right|$",
        clabels=true,
        color=:thermal,
        aspect_ratio=:equal,
        size=(700,707))

savefig("/mnt/data/code/Plots/contour_$(N)_$(length(x_range))_points_focused.png")

####################################################
## Smallest and focused
x_range = 0.095:0.001:0.105
y_range = 0.003:0.001:0.013

# Initialize the matrix to store the ratios
ratios = zeros(length(y_range), length(x_range))

for (x_idx, x) in enumerate(x_range)
    for (y_idx, y) in enumerate(y_range)
        # Perturb the original matrix
        A_perturbed = A + x*P1 + y*P2
        _, ANN = growth_factor_no_pivot(A)
        _, ANN_ = growth_factor_no_pivot(A_perturbed)

        # Calculate the ratio and store it
        ratios[y_idx, x_idx] = ANN_ / ANN
    end
end

# Generate the contour plot
p = contourf(x_range,
        y_range,
        log.(abs.(ratios)),
        # levels=collect(-2:1:8.22),
        levels=collect(2:0.5:8.22),
        xlabel="Perturbation in (1, 1) direction",
        ylabel="Random perturbation without perturbing (1, 1)",
        title="Relative change in last pivot for a Higham matrix\nwhen perturbing (1,1) vs. random perturbation w/o (1,1)",
        colorbar_title=L"$\log_{10} \left|\!\!p_\epsilon / p\ \right|$",
        clabels=true,
        color=:thermal,
        aspect_ratio=:equal,
        size=(700,707))

savefig("/mnt/data/code/Plots/contour_$(N)_$(length(x_range))_points_focused_smallest.png")
