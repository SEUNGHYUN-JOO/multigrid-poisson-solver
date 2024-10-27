using Printf
using LinearAlgebra

# Enumerated types for available smoothers
@enum Smoother JACOBI GAUSS_SEIDEL SOR

# Definition of the finest grid level for code readability
const FINE_MESH = 1

# Set the variable parameters for the simulation here
const NUM_NODES = 1025
const MG_CYCLES = 200
const DISABLE_MG = 0
const NUM_SWEEP = 3
const SMOOTHER = GAUSS_SEIDEL
const VISUALIZE = 0
const FREQUENCY = 1
const TOLERANCE = 20.0
const w_array = [2, 2, 2]
const NoA = length(w_array)

function allocate_arrays(n_nodes)
    # Allocate memory for grid levels
    n_levels = 1
    nodes = n_nodes
    while (nodes - 1) % 2 == 0 && (nodes - 1) / 2 + 1 >= 3
        nodes = (nodes - 1) / 2 + 1
        n_levels += 1
    end
    
    phi = Vector{Matrix{Float64}}(undef, n_levels)
    aux = Vector{Matrix{Float64}}(undef, n_levels)
    f = Vector{Matrix{Float64}}(undef, n_levels)
    x = Vector{Matrix{Float64}}(undef, n_levels)
    y = Vector{Matrix{Float64}}(undef, n_levels)
    phi_exact = Matrix{Float64}(undef, n_nodes, n_nodes)
    
    nodes = n_nodes
    for i_level in 1:n_levels
        phi[i_level] = zeros(nodes, nodes)
        aux[i_level] = zeros(nodes, nodes)
        f[i_level] = zeros(nodes, nodes)
        x[i_level] = zeros(nodes, nodes)
        y[i_level] = zeros(nodes, nodes)
        nodes = (nodes - 1) ÷ 2 + 1
    end
    return phi, phi_exact, f, x, y, aux, n_levels
end

function generate_fine_mesh(x, y, n_nodes, visualize)
    # Initialize the x & y coordinate arrays over [0,1] X [0,1] with uniform spacing
    for i in 1:n_nodes
        for j in 1:n_nodes
            x[i, j] = (i - 1) / (n_nodes - 1)
            y[i, j] = (j - 1) / (n_nodes - 1)
        end
    end
    println("Created fine cartesian grid ($n_nodes x $n_nodes)...")
end

function coarsen_mesh(x, y, n_nodes, n_levels, level)
    # Create coarse MG levels by keeping every other node in the fine mesh
    if level == n_levels
        return
    else
        for i in 1:n_nodes
            for j in 1:n_nodes
                if i % 2 == 1 && j % 2 == 1
                    x[level + 1][(i + 1) ÷ 2, (j + 1) ÷ 2] = x[level][i, j]
                    y[level + 1][(i + 1) ÷ 2, (j + 1) ÷ 2] = y[level][i, j]
                end
            end
        end
        n_coarse = (n_nodes - 1) ÷ 2 + 1
        coarsen_mesh(x, y, n_coarse, n_levels, level + 1)
    end
end

function initialize_solution(phi, phi_exact, f, x, y, n_nodes)
    for i in 1:n_nodes
        for j in 1:n_nodes
            if i == 1 || j == 1 || i == n_nodes || j == n_nodes
                phi[i, j] = exp(x[i, j]) * exp(-2.0 * y[i, j])
            else
                phi[i, j] = 0.0
            end
            phi_exact[i, j] = exp(x[i, j]) * exp(-2.0 * y[i, j])
            f[i, j] = -5.0 * exp(x[i, j]) * exp(-2.0 * y[i, j])
        end
    end
    println("Successfully initialized solution...")
end

function smooth_total(phi, f, aux, n_nodes, n_sweeps, level)
    if SMOOTHER == JACOBI
        smooth_jacobi(phi[level], f[level], aux[level], n_nodes, n_sweeps)
    elseif SMOOTHER == GAUSS_SEIDEL
        smooth_gauss_seidel(phi[level], f[level], aux[level], n_nodes, n_sweeps)
    elseif SMOOTHER == SOR
        smooth_sor(phi[level], f[level], aux[level], n_nodes, n_sweeps)
    else
        println("Unrecognized smoother.")
        exit(1)
    end
end

function multigrid_cycle(phi, f, aux, n_nodes, n_sweeps, n_levels, flevel, w_array)
    # 초기 설정
    base_level = n_levels
    direc = 0
    level = 0
    
    for i in 1:div(NoA, 2)
        if flevel < base_level - 1
            direc = 1
            level = n_levels - (base_level - flevel)
            break
        elseif NoA % 2 == 0 && flevel == base_level - 1 && i == div(NoA, 2)
            level = n_levels - (base_level - flevel)
            break
        elseif flevel < base_level + w_array[i] - 1
            direc = -1
            level = (n_levels - 2) - (flevel - base_level)
            break
        elseif NoA % 2 == 1 && flevel == base_level + w_array[i] - 1 && i + 1 > div(NoA, 2)
            level = (n_levels - 2) - (flevel - base_level)
            break
        end
        base_level += 2 * w_array[i]
    end
    
    # 스무딩 (smoothing)
    smooth_total(phi, f, aux, n_nodes, n_sweeps, level)
    
    # 업/다운 전환 및 재귀 호출
    if direc == 1
        up_down(phi, f, aux, n_nodes, level, n_sweeps, n_levels, flevel, w_array)
    elseif direc == -1
        down_up(phi, f, aux, n_nodes, level, n_sweeps, n_levels, flevel, w_array)
    end
    
    # 스무딩 (마무리)
    smooth_total(phi, f, aux, n_nodes, n_sweeps, level)
end

