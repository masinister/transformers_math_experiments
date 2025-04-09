include("constants.jl")
using StatsBase

const N = 100  # Adjust the size of the set [1,...,N] as needed
const K = 3    # Length of the arithmetic progression to avoid

function find_all_k_aps(subset::Vector{Int}, k::Int)::Vector{Vector{Int}}
    """
    Find all k-term arithmetic progressions in the subset.
    """
    k_aps = []
    for i in subset
        for j in subset
            if i < j
                d = j - i
                k_term_ap = [i + n * d for n in 0:k-1]
                if all(x -> x in subset, k_term_ap)
                    push!(k_aps, k_term_ap)
                end
            end
        end
    end
    return k_aps
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Greedy algorithm to construct the largest subset of [1,...,N] avoiding k-term arithmetic progressions.
    """
    # Decode the input object into the subset
    subset = [i for (i, c) in enumerate(obj) if c == '1']
    
    k_aps = find_all_k_aps(subset, K)

    # Delete the element appearing in the most k-APs until no k-APs are left
    while !isempty(k_aps)
        # Count frequency of each element in k-APs
        element_count = Dict{Int, Int}()
        for ap in k_aps
            for element in ap
                element_count[element] = get(element_count, element, 0) + 1
            end
        end

        # Sample the element to remove, weighted by frequency
        elements = collect(keys(element_count))
        frequencies = collect(values(element_count))
        most_frequent_element = sample(elements, Weights(frequencies))

        # Remove this element from the subset
        subset = setdiff(subset, [most_frequent_element])

        # Update k-APs by removing any that contain the most frequent element
        k_aps = filter(ap -> !(most_frequent_element in ap), k_aps)
    end

    # Now keep adding random integers without creating k-APs, until stuck
    allowed_elements = setdiff(1:N, subset)

    while !isempty(allowed_elements)
        # Randomly select an element to add
        element = allowed_elements[rand(1:length(allowed_elements))]
        subset = push!(subset, element)

        # Check if adding this element creates any k-APs
        k_aps = find_all_k_aps(subset, K)
        if !isempty(k_aps)
            # If it creates k-APs, remove the element and stop
            subset = setdiff(subset, [element])
            break
        end

        # Update allowed elements by removing the added element
        allowed_elements = setdiff(allowed_elements, [element])
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