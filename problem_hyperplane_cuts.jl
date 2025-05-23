include("constants.jl")
using Random
using StatsBase
using LinearAlgebra

const N = 6  # Dimension of the hypercube
const K = 5  # Number of hyperplanes to find
const max_coeff = 32  # Maximum absolute value for coefficients

# Constants for basis vector representation
const TOTAL_BASIS_DIM = (N+1) * K
const BASIS_TOKENS = ["b$(i)$(sign)" for i in 1:TOTAL_BASIS_DIM for sign in ["+", "-"]]
const RANDOM_TOKENS = ["r$(i)" for i in 1:TOTAL_BASIS_DIM]  # Simplified: reduced number of random tokens
const ALL_TOKENS = vcat(BASIS_TOKENS, RANDOM_TOKENS)

# A hyperplane is represented as a vector of N+1 integers [a_1, a_2, ..., a_N, b]
# It corresponds to the equation: a_1*x_1 + a_2*x_2 + ... + a_N*x_N = b

# Pre-generate all edges of the hypercube
const ALL_HYPERCUBE_EDGES = let
    edges = Vector{Tuple{Vector{Int}, Vector{Int}}}()
    
    # Generate all 2^N vertices of the hypercube
    for i in 0:(2^N-1)
        # Convert i to binary representation as vertex
        vertex = [(i >> j) & 1 for j in 0:(N-1)]
        
        # For each zero coordinate, flip it to get the neighboring vertex
        for j in 1:N
            if vertex[j] == 0
                neighbor = copy(vertex)
                neighbor[j] = 1 
                push!(edges, (vertex, neighbor))
            end
        end
    end
    print("Generated $(length(edges)) edges of the hypercube.\n")
    edges
end

function generate_hypercube_edges()::Vector{Tuple{Vector{Int}, Vector{Int}}}
    """
    Return the precomputed edges of the N-dimensional hypercube.
    """
    return ALL_HYPERCUBE_EDGES
end

function hyperplane_cuts_edge(hyperplane::Vector{Int}, edge::Tuple{Vector{Int}, Vector{Int}})::Bool
    """
    Check if a hyperplane cuts an edge of the hypercube.
    A hyperplane cuts an edge if the endpoints of the edge are on opposite sides.
    """
    v1, v2 = edge
    
    # Extract coefficients and constant term
    coeffs = @view hyperplane[1:N]
    b = hyperplane[N+1]
    
    # Calculate dot products more efficiently
    dot_product1 = dot(coeffs, v1)
    dot_product2 = dot(coeffs, v2)
    
    # Direct comparison - faster than creating intermediate variables
    return (dot_product1 < b && dot_product2 > b) || (dot_product1 > b && dot_product2 < b)
end

function compute_edge_cutting_matrix(hyperplanes::Vector{Vector{Int}}, edges::Vector{Tuple{Vector{Int}, Vector{Int}}})::Matrix{Bool}
    """
    Precompute which hyperplanes cut which edges.
    Returns a matrix where M[h,e] = true if hyperplane h cuts edge e
    """
    num_hyperplanes = length(hyperplanes)
    num_edges = length(edges)
    
    cutting_matrix = falses(num_hyperplanes, num_edges)
    
    for (h_idx, hyperplane) in enumerate(hyperplanes)
        for (e_idx, edge) in enumerate(edges)
            cutting_matrix[h_idx, e_idx] = hyperplane_cuts_edge(hyperplane, edge)
        end
    end
    
    return cutting_matrix
end

function count_cut_edges(hyperplanes::Vector{Vector{Int}}, edges::Vector{Tuple{Vector{Int}, Vector{Int}}})::Int
    """
    Count how many unique edges are cut by at least one of the hyperplanes.
    """
    if isempty(hyperplanes)
        return 0
    end
    
    # Use a precomputed matrix to avoid repeated checks
    cutting_matrix = compute_edge_cutting_matrix(hyperplanes, edges)
    
    # Count edges that are cut by at least one hyperplane
    cut_count = sum(any(cutting_matrix, dims=1))
    
    return cut_count
end

function random_hyperplane()::Vector{Int}
    """
    Generate a random hyperplane using the basis vector approach.
    """
    hyperplane = zeros(Int, N+1)
    
    # Add random basis vectors
    num_vectors = rand(3:8)  # Simplified: reasonable range
    
    for _ in 1:num_vectors
        coef_idx = rand(1:N+1)
        sign = rand([-1, 1])
        hyperplane[coef_idx] += sign
    end
    
    # Ensure at least one non-zero coefficient for the hyperplane equation
    if all(iszero, hyperplane[1:N])
        hyperplane[rand(1:N)] = rand([-1, 1])
    end
    
    return clamp.(hyperplane, -max_coeff, max_coeff)
end

# Fix the token_to_basis_index function to accept AbstractString instead of String
function token_to_basis_index(token::AbstractString)::Tuple{Int, Int, Int}
    """
    Convert a token to basis vector information.
    Returns (hyperplane_index, coefficient_index, value_sign)
    """
    if !startswith(token, "b")
        return (0, 0, 0)
    end
    
    # Find the sign character (last character)
    sign_char = token[end]
    if !(sign_char in ['+', '-'])
        return (0, 0, 0)
    end
    
    # Parse the index
    basis_idx_str = token[2:end-1]
    basis_idx = tryparse(Int, basis_idx_str)
    if basis_idx === nothing || basis_idx < 1 || basis_idx > TOTAL_BASIS_DIM
        return (0, 0, 0)
    end
    
    sign = sign_char == '+' ? 1 : -1
    
    # Map to hyperplane and coefficient indices
    hp_idx = div(basis_idx - 1, (N+1)) + 1
    coef_idx = mod(basis_idx - 1, (N+1)) + 1
    
    return (hp_idx, coef_idx, sign)
end

# Alternative fix would be to convert SubString to String in parse_obj_to_hyperplanes:
function parse_obj_to_hyperplanes(obj::OBJ_TYPE)::Vector{Vector{Int}}
    """
    Parse the token-based representation into hyperplanes.
    """
    hyperplanes = [zeros(Int, N+1) for _ in 1:K]
    
    if isempty(obj)
        return hyperplanes
    end
    
    # Process each token
    for token_substr in split(obj)
        # Convert SubString to String to avoid type issues
        token = String(token_substr)
        hp_idx, coef_idx, sign = token_to_basis_index(token)
        
        if hp_idx >= 1 && hp_idx <= K && coef_idx >= 1 && coef_idx <= (N+1)
            hyperplanes[hp_idx][coef_idx] += sign
        end
    end
    
    # Clamp values to valid range
    for hp in hyperplanes
        for i in 1:length(hp)
            hp[i] = clamp(hp[i], -max_coeff, max_coeff)
        end
    end
    
    return hyperplanes
end

function hyperplanes_to_obj(hyperplanes::Vector{Vector{Int}})::OBJ_TYPE
    """
    Convert hyperplanes to token-based representation.
    """
    tokens = String[]
    
    for (hp_idx, hp) in enumerate(hyperplanes)
        for (coef_idx, value) in enumerate(hp)
            if value == 0
                continue
            end
            
            # Create basis vector tokens for this coefficient
            basis_idx = (hp_idx - 1) * (N+1) + coef_idx
            sign_char = value > 0 ? "+" : "-"
            token = "b$(basis_idx)$(sign_char)"
            
            # Add tokens for the absolute value
            append!(tokens, fill(token, abs(value)))
        end
    end
    
    # Add random tokens if needed for variety
    if isempty(tokens)
        push!(tokens, RANDOM_TOKENS[1])
    elseif length(tokens) < 3
        # Add a few random tokens for short representations
        append!(tokens, sample(RANDOM_TOKENS, 3 - length(tokens), replace=false))
    end
    
    # Shuffle the tokens
    shuffle!(tokens)
    
    return join(tokens, " ")
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Optimized greedy algorithm using token-based representation.
    """
    hyperplanes = parse_obj_to_hyperplanes(obj)
    edges = generate_hypercube_edges()
    total_edges = length(edges)
    
    # Ensure we have K hyperplanes
    while length(hyperplanes) < K
        push!(hyperplanes, random_hyperplane())
    end
    
    # Initial evaluation
    current_cut_count = count_cut_edges(hyperplanes, edges)
    
    # Early exit if we've already cut all edges
    if current_cut_count == total_edges
        return [hyperplanes_to_obj(hyperplanes)]
    end
    
    # Improvement loop
    max_iterations = 10
    for iteration in 1:max_iterations
        improved = false
        
        for i in 1:K
            # Find edges that aren't cut by other hyperplanes
            other_hyperplanes = vcat(hyperplanes[1:i-1], hyperplanes[i+1:end])
            cutting_matrix = compute_edge_cutting_matrix(other_hyperplanes, edges)
            already_cut = [any(cutting_matrix[:, e]) for e in 1:length(edges)]
            uncut_edges = findall(.!already_cut)
            
            if isempty(uncut_edges)
                continue
            end
            
            # Try to find a better hyperplane
            best_cut_count = 0
            best_hyperplane = hyperplanes[i]
            trials = max(100, 200 - iteration * 20)
            
            for _ in 1:trials
                # Generate candidate using basis vector approach
                candidate = random_hyperplane()
                
                # Count how many previously-uncut edges this candidate cuts
                cut_count = sum(edge_idx -> hyperplane_cuts_edge(candidate, edges[edge_idx]), uncut_edges)
                
                if cut_count > best_cut_count
                    best_cut_count = cut_count
                    best_hyperplane = candidate
                    improved = true
                end
            end
            
            hyperplanes[i] = best_hyperplane
        end
        
        # Re-evaluate
        new_cut_count = count_cut_edges(hyperplanes, edges)
        
        # Stop if no improvement or all edges cut
        if new_cut_count == total_edges || new_cut_count <= current_cut_count || !improved
            break
        end
        
        current_cut_count = new_cut_count
    end
    
    return [hyperplanes_to_obj(hyperplanes)]
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Reward is the number of edges cut by the hyperplanes.
    """
    hyperplanes = parse_obj_to_hyperplanes(obj)
    edges = generate_hypercube_edges()
    return count_cut_edges(hyperplanes, edges)
end

function empty_starting_point()::OBJ_TYPE
    """
    Start with a minimal representation.
    """
    num_tokens = min(2, length(RANDOM_TOKENS))
    tokens = sample(RANDOM_TOKENS, num_tokens, replace=false)
    return join(tokens, " ")
end
