include("constants.jl")
using StatsBase


const N = 200  # Adjust the size of the set [1,...,N] as needed

function find_all_pythagorean_triples(max_n::Int)::Vector{Tuple{Int, Int, Int}}
    """
    Generate all Pythagorean triples (a, b, c) with 1 ≤ a < b < c ≤ max_n.
    """
    triples = []
    for a in 1:max_n-2
        for b in a+1:max_n-1
            c2 = a^2 + b^2
            c = sqrt(c2)
            if isinteger(c) && c ≤ max_n
                push!(triples, (a, b, c))
            end
        end
    end
    return triples
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Greedy algorithm to construct the largest subset of [1,...,N] avoiding Pythagorean triples.
    """
    # Decode the input object into the subset
    subset = [i for (i, c) in enumerate(obj) if c == '1']
    
    # Find all Pythagorean triples
    triples = find_all_pythagorean_triples(N)

    # Start with the full set and iteratively remove elements to destroy triples
    while true
        # Identify Pythagorean triples in the current subset
        conflicting_triples = filter(t -> all(x -> x in subset, t), triples)
        if isempty(conflicting_triples)
            break  # No more triples to resolve
        end

        # Count frequency of elements in the triples
        frequency = Dict{Int, Int}()
        for (a, b, c) in conflicting_triples
            frequency[a] = get(frequency, a, 0) + 1
            frequency[b] = get(frequency, b, 0) + 1
            frequency[c] = get(frequency, c, 0) + 1
        end

        # Remove a random element with probability weighted by the frequency
        keys_array = collect(keys(frequency))
        weights_array = collect(values(frequency))
        to_remove = sample(keys_array, Weights(weights_array))
        subset = filter(x -> x != to_remove, subset)
    end

    # Convert subset to binary string representation
    result = join([x in subset ? "1" : "0" for x in 1:N])
    return [result]
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Reward is the size of the subset (number of 1s in the binary string).
    """
    return count(isequal('1'), obj)
end

function empty_starting_point()::OBJ_TYPE
    """
    Start with the full set [1,...,N] as the initial binary string.
    """
    return "1" ^ N
end
