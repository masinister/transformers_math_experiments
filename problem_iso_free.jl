include("constants.jl")
using StatsBase
using Plots

const N = 30  # Dimension of the grid

function manhattan_distance(x::Vector{Int}, y::Vector{Int})::Int
    return sum(abs.(x .- y))
end

function convert_points_to_string(points::Vector{Int})::String
    return join(points)
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Greedy algorithm to construct the largest subset of the N x N grid avoiding isosceles triangles.
    """
    # Decode the input object into the subset
    subset = [i for (i, c) in enumerate(obj) if c == '1']
    
    # Generate all points in the N x N grid
    points = [[div(i, N), mod(i, N)] for i in 0:(N^2 - 1)]
    
    # Randomly order the subset
    subset = shuffle(subset)
    
    # Start with the full set and iteratively remove points to destroy triangles
    while true
        # Identify the first isosceles triangle in the current subset
        found_triangle = false
        for i in eachindex(subset)
            distance_hash = Dict{Int, Vector{Int}}()
            for j in eachindex(subset)
                if i != j
                    dist = manhattan_distance(points[subset[i]], points[subset[j]])
                    if haskey(distance_hash, dist)
                        for k in distance_hash[dist]
                            subset = filter(x -> x != subset[i], subset)
                            found_triangle = true
                            break
                        end
                        if found_triangle
                            break
                        end
                        push!(distance_hash[dist], subset[j])
                    else
                        distance_hash[dist] = [subset[j]]
                    end
                end
            end
            if found_triangle
                break
            end
        end

        if !found_triangle
            break  # No more triangles to resolve
        end
    end

    # Convert subset to binary string representation
    result = join([x in subset ? "1" : "0" for x in 1:N^2])
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
    Start with the full set of points in the N x N grid as the initial binary string.
    """
    return "1" ^ (N^2)
end

function draw(obj::OBJ_TYPE)
    subset = [i for (i, c) in enumerate(obj) if c == '1']
    points = [[div(i, N), mod(i, N)] for i in 0:(N^2 - 1)]
    x_coords = [points[i][1] for i in subset]
    y_coords = [points[i][2] for i in subset]
    
    # Draw the N x N grid with gridlines
    p = plot(seriestype=:scatter, markershape=:circle, xlims=(-1, N), ylims=(-1, N), aspect_ratio=:equal, xaxis=false, yaxis=false)
    for i in 0:N-1
        plot!([i, i], [0, N-1], color=:gray, linestyle=:dash, label=false)
        plot!([0, N-1], [i, i], color=:gray, linestyle=:dash, label=false)
    end
    
    # Highlight the subset points
    scatter!(p, x_coords, y_coords, markershape=:circle, color=:red, label=false)
end