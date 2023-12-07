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

function print_matrix_g(matrix)
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            @printf("%g ", matrix[i, j])
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

function perturbation!(original_matrix, scale)
    # for i in eachindex(original_matrix)
    #     original_matrix[i] += scale * randn_()
    # end
    original_matrix[1, 1] += scale * randn_()
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

    # Set the random number seed
    Random.seed!(Int(time_ns()))

    for i in 1:N-1
        for j in 1:N-1
            if i <= j
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

# original function to compute the growth factor without pivoting
# function growth_factor_no_pivot(A)
#     N = size(A, 1)
#     max_A = maximum(abs.(A))
#     max_A_k = max_A
    
#     for k in 1:(N-1)
#         println("k = ", k - 1, ":")
#         max_val = abs(A[k, k])
        
#         for i in (k+1):N
#             factor = A[i, k] / A[k, k]
#             @printf("factor = %.16f\n", factor)
#             @printf("value = %.16f\n", abs(factor+1)*(2^52))

#             A[i, k:N] .-= factor * A[k, k:N]
            
#             max_val_i = maximum(abs.(A[i, :]))
#             if max_val_i > max_A_k
#                 max_A_k = max_val_i
#             end
#         end
#         @printf("After k = %d: max = %f, Growth = %f\n", k-1, max_A_k, max_A_k/max_A)
#         # print_matrix(A)
#         println("\n")
#     end
    
#     return max_A_k / max_A
# end

# compute the growth factor and save the matrix at each step
function growth_factor_no_pivot_save(A)
    A_refs = Array{Matrix{Rational{BigInt}},1}()
    push!(A_refs, copy(A))
    N = size(A, 1)
    max_A = maximum(abs.(A))
    max_A_k = max_A
    
    for k in 1:(N-1)
        # println("k = ", k - 1, ":")
        max_val = abs(A[k, k])
        
        for i in (k+1):N
            factor = A[i, k] / A[k, k]
            # @printf("factor = %.16f\n", factor)
            # @printf("value = %.16f\n", abs(factor+1)*(2^52))

            A[i, k:N] .-= factor * A[k, k:N]
            
            max_val_i = maximum(abs.(A[i, :]))
            if max_val_i > max_A_k
                max_A_k = max_val_i
            end
        end
        # @printf("After k = %d: max = %f, Growth = %f\n", k-1, max_A_k, max_A_k/max_A)
        # print_matrix(A)
        # println("\n")
        push!(A_refs, copy(A))
    end
    
    return (max_A_k / max_A), A_refs
end

# compute the growth factor and compare the matrix with A_refs at each step
function growth_factor_no_pivot_compare(A, A_refs)
    @printf("Initial diff:\n")
    print_matrix_g(A - A_refs[1])
    N = size(A, 1)
    max_A = maximum(abs.(A))
    max_A_k = max_A
    
    prev_diff = fill(Rational{BigInt}(0), N)
    for k in 1:(N-1)
        println("k = ", k - 1, ":")
        max_val = abs(A[k, k])
        
        factor = zeros(N,)
        for i in (k+1):N
            factor = A[i, k] / A[k, k]
            @printf("factor = %.16f\n", factor)
            @printf("factor diff = %g\n", factor+1)
            ref_diff = A - A_refs[k]
            if (k == 1)
                # println(factor + 1)
                # println(ref_diff[k, k])
                # println(A[k, k])
                # println(ref_diff[k, k] / A[k, k])
                @assert (A[i, k] + A[k, k]) / A[k, k] == factor - (-1)
                @assert ref_diff[k, k] / A[k, k] == factor - (-1)
            else
                # println(factor + 1)
                # println(ref_diff[k, k])
                # println(A[k, k])
                # println(2 * ref_diff[k, k] / A[k, k])
                @assert (A[i, k] + A[k, k]) / A[k, k] == factor - (-1)
                @assert 2 * ref_diff[k, k] / A[k, k] == factor - (-1)
            end
            # @printf("value = %.16f\n", abs(factor+1)*(2^52))

            A[i, k:N] .-= factor * A[k, k:N]
            
            max_val_i = maximum(abs.(A[i, :]))
            if max_val_i > max_A_k
                max_A_k = max_val_i
            end
        end
        @printf("After k = %d: max = %f, Growth = %f\n", k-1, max_A_k, max_A_k/max_A)
        print_matrix(A)
        @printf("A diff:\n")
        print_matrix_g(A - A_refs[k+1])
        @printf("Predicted A diff:\n")
        
        print_matrix_g((-(factor+1) * A[k, k+1:end] + 2 * prev_diff[2:end])')
        diff = A - A_refs[k+1]
        diff = diff[k+1, k+1:end]
        print(diff == -(factor+1) * A[k, k+1:end] + 2 * prev_diff[2:end])
        @assert diff == -(factor+1) * A[k, k+1:end] + 2 * prev_diff[2:end]
        # diff_part = -(factor+1) * A[k, k+1:end]
        # for i in eachindex(diff_part)
        #     print(sign(diff_part[i]) == sign(prev_diff[i+1]))
        # end
        prev_diff = diff
        
        println("\n")
    end
    
    return max_A_k / max_A
end

function main()
    Random.seed!(Int(round(time())) % 1_000_000_007)

    N = 6
    sample_size = 1
    output_file = "growth_N$(N)_sample$(sample_size).txt"
    
    A = Higham(N)
    println("A =")
    print_matrix(A)

    A_copy = copy(A)
    factor_no_pivot, A_refs = growth_factor_no_pivot_save(A_copy)
    @printf("Original Growth without Pivoting: %.10f\n", factor_no_pivot)

    # Open the file for writing
    open(output_file, "w") do output_file
        # scales = [BigInt(1)//BigInt(10^8), BigInt(1)//BigInt(10^9), BigInt(1)//BigInt(10^10), BigInt(1)//BigInt(10^11), BigInt(1)//BigInt(10^12), BigInt(1)//BigInt(10^13), BigInt(1)//BigInt(10^14), BigInt(1)//BigInt(10^15)]
        scales = [BigInt(1)//BigInt(10^16)]
        
        for scale in scales
            total_growth = BigInt(0) // BigInt(1)
            total_growth_change = BigInt(0) // BigInt(1)
            
            for _ in 1:sample_size
                A_copy = copy(A)
                perturbation!(A_copy, scale)
                factor_ = growth_factor_no_pivot_compare(A_copy, A_refs)
                growth_change = abs(factor_no_pivot - factor_)
                total_growth += factor_
                total_growth_change += growth_change
                
                formatted_string = @sprintf("%e %.16f %.16f\n", float(scale), float(factor_), float(growth_change))
                write(output_file, formatted_string)
            end
            
            average_growth = total_growth / sample_size
            average_growth_change = total_growth_change / sample_size
            @printf("For scale %.10e: Average growth: %.16f, Average growth change: %.16f\n", scale, average_growth, average_growth_change)

            # formatted_string = @sprintf("%e %.16f %.16f\n", scale, float(average_growth), float(average_growth_change))
            # write(output_file, formatted_string)
        end
    end

    time_elapsed = time() - start_time
    @printf("Time elapsed: %f seconds\n", time_elapsed)
end

start_time = time()
main()
