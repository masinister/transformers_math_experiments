include("constants.jl")
using Random
using StatsBase
using LinearAlgebra

const N = 6  # Dimension of the hypercube
const K = 5  # Number of hyperplanes to find
const max_coeff = 32  # Maximum absolute value for coefficients

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
    Generate a random hyperplane with integer coefficients.
    Using very small coefficient range to avoid tokenization issues.
    """
    # Coefficients between -max_coeff and max_coeff
    coeffs = rand(-max_coeff:max_coeff, N)
    
    # Avoid generating a zero vector for coefficients
    while all(c -> c == 0, coeffs)
        coeffs = rand(-max_coeff:max_coeff, N)
    end
    
    # Choose the constant term to make the hyperplane intersect the hypercube
    # Even smaller range for b to avoid tokenization issues
    b = rand(0:max_coeff)
    
    return [coeffs..., b]
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Optimized greedy algorithm to find K hyperplanes that cut as many edges as possible.
    """
    # Parse the input object into a list of hyperplanes
    hyperplanes = parse_obj_to_hyperplanes(obj)
    
    # Generate all edges of the hypercube
    edges = generate_hypercube_edges()
    total_edges = length(edges)
    
    # Start with K hyperplanes (add random ones if needed)
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
        
        # Try to improve each hyperplane
        for i in 1:K
            # Remove current hyperplane to see which edges are still cut by others
            other_hyperplanes = vcat(hyperplanes[1:i-1], hyperplanes[i+1:end])
            cutting_matrix = compute_edge_cutting_matrix(other_hyperplanes, edges)
            already_cut = [any(cutting_matrix[:, e]) for e in 1:length(edges)]
            uncut_edges = findall(.!already_cut)
            
            # If all edges are already cut by other hyperplanes, skip this one
            if isempty(uncut_edges)
                continue
            end
            
            # Find a better hyperplane
            best_cut_count = 0
            best_hyperplane = hyperplanes[i]
            trials = max(100, 200 - iteration * 20)
            
            for _ in 1:trials
                candidate = random_hyperplane()
                cut_count = sum(edge_idx -> hyperplane_cuts_edge(candidate, edges[edge_idx]), uncut_edges)
                
                if cut_count > best_cut_count
                    best_cut_count = cut_count
                    best_hyperplane = candidate
                    improved = true
                end
            end
            
            hyperplanes[i] = best_hyperplane
        end
        
        # Re-evaluate overall performance
        new_cut_count = count_cut_edges(hyperplanes, edges)
        
        # Stop if we've cut all edges or if no improvement
        if new_cut_count == total_edges || new_cut_count <= current_cut_count || !improved
            break
        end
        
        current_cut_count = new_cut_count
    end
    
    return [hyperplanes_to_obj(hyperplanes)]
end

function parse_obj_to_hyperplanes(obj::OBJ_TYPE)::Vector{Vector{Int}}
    """
    Parse the string representation into a list of hyperplanes.
    Format: "a1,a2,...,aN,b;a1,a2,...,aN,b;..."
    Handle potential errors gracefully.
    """
    if isempty(obj)
        return []
    end
    
    hyperplanes = []
    for hp_str in split(obj, ";")
        if isempty(hp_str)
            continue
        end
        
        # Parse values with error handling
        try
            parts = split(hp_str, ",")
            if length(parts) == N+1
                coeffs = parse.(Int, parts)
                push!(hyperplanes, coeffs)
            end
        catch e
            # If parsing fails, skip this hyperplane
            continue
        end
    end
    
    # If parsing failed completely, provide a default
    if isempty(hyperplanes) && !isempty(obj)
        for _ in 1:K
            push!(hyperplanes, random_hyperplane())
        end
    end
    
    return hyperplanes
end

function format_two_digits(num::Int)::String
    """
    Format an integer to have exactly 2 digits with appropriate sign.
    Add leading zero for single-digit numbers.
    Add + sign before positive numbers to ensure consistent length.
    """
    if num < 0
        # Negative number
        if num > -10
            return "-0$(abs(num))" # e.g., -5 becomes "-05"
        else
            return "$(num)" # e.g., -15 remains "-15"
        end
    elseif num < 10
        # Single digit positive or zero
        return "+0$(num)" # e.g., 5 becomes "+05", 0 becomes "+00"
    else
        # Double digit positive
        return "+$(num)" # e.g., 15 becomes "+15"
    end
end

function hyperplanes_to_obj(hyperplanes::Vector{Vector{Int}})::OBJ_TYPE
    """
    Convert a list of hyperplanes to string representation.
    Format: "a1,a2,...,aN,b;a1,a2,...,aN,b;..."
    Ensure all values are within a safe range for tokenization.
    Format each number to be exactly 2 digits.
    """
    safe_hyperplanes = []
    for hp in hyperplanes
        safe_hp = [clamp(val, -max_coeff, max_coeff) for val in hp]
        push!(safe_hyperplanes, safe_hp)
    end
    
    # Format each number to be exactly 2 digits
    formatted_hyperplanes = []
    for hp in safe_hyperplanes
        formatted_hp = [format_two_digits(val) for val in hp]
        push!(formatted_hyperplanes, formatted_hp)
    end
    
    # Create string representation
    hp_strs = [join(hp, ",") for hp in formatted_hyperplanes]
    return join(hp_strs, ";")
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Reward is the number of edges cut by the hyperplanes.
    """
    # No need to strip padding as we're using fixed-width number formatting
    hyperplanes = parse_obj_to_hyperplanes(obj)
    edges = generate_hypercube_edges()
    return count_cut_edges(hyperplanes, edges)
end

function empty_starting_point()::OBJ_TYPE
    """
    Start with K hyperplanes with all zero coefficients.
    """
    hyperplanes = [zeros(Int, N+1) for _ in 1:K]
    return hyperplanes_to_obj(hyperplanes)
end
