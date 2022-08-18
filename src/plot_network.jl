
PT = GeometryBasics.Point{2,Float64}
const DEFAULT_LON = "37.8270"
const DEFAULT_LAT = "-122.8230"

function color_nodes!(
    g,
    sys,
    node_colors::Vector{RGB{Colors.N0f8}} = map(
        x -> RGB{Colors.N0f8}(1.0, 0.0, 0.0),
        1:nv(g),
    ),
)
    alpha = Vector()
    for (ix, b) in enumerate(get_components(Bus, sys))
        a = (haskey(get_ext(b), "lat") & haskey(get_ext(b), "lon")) ? 1.0 : 0.1
        push!(alpha, a)
    end
    set_prop!(g, :nodecolor, node_colors)
    set_prop!(g, :alpha, alpha)
end

function color_nodes!(g, sys, color_by::Type{T}) where {T<:AggregationTopology}
    # Generate n maximally distinguishable colors in LCHab space.
    areas = get_components(color_by, sys)
    buses = get_components(Bus, sys)
    area_colors =
        Dict(zip(get_name.(areas), distinguishable_colors(length(areas), colorant"blue")))
    node_colors = getindex.(Ref(area_colors), get_name.(get_area.(buses)))
    color_nodes!(g, sys, node_colors)
    set_prop!(g, :group, get_name.(get_area.(buses)))
end

# construct a graph from a PowerSystems.System
function make_graph(sys::PowerSystems.System; kwargs...)
    @info "creating graph from System"
    buses = get_components(Bus, sys)
    g = MetaGraph(length(buses))
    for (i, b) in enumerate(buses)
        set_props!(
            g,
            i,
            Dict(
                :name => get_name(b),
                :number => get_number(b),
                :area => get_name(get_area(b)),
            ),
        )# :nodecolor = node_colors[i], :x = network[i][1], :y = network[i][2]))
    end
    set_indexing_prop!(g, :name)

    for a in get_components(Arc, sys)
        from = get_name(get_from(a))
        to = get_name(get_to(a))
        add_edge!(g, g[from, :name], g[to, :name])
    end

    # color nodes
    #fixed_colors = Dict(zip(unique(fixed), distinguishable_colors(2, colorant"blue")))
    #node_colors = getindex.(Ref(fixed_colors), fixed)
    color_nodes!(g, sys, Area)

    # layout nodes
    a = adjacency_matrix(g) # generates a sparse adjacency matrix
    #a = Adjacency(sys, check_connectivity = false).data
    fixed = .!isempty.(get_ext.(buses))
    initial_positions = map(x -> PT(0.0, 0.0), 1:size(a, 1))
    for (ix, b) in enumerate(buses)
        if !isempty(get_ext(b))
            ext = get_ext(b)
            lat = tryparse(Float64, "$(get(ext, "latitude", get(ext, "lat", nothing)))")
            lon = tryparse(Float64, "$(get(ext, "longitude", get(ext, "lon", nothing)))")
            if !isnothing(lat) && !isnothing(lon)
                initial_positions[ix] = PT(lon, lat)
            else
                fixed[ix] = false
            end
        end
    end
    nz_pos = initial_positions[initial_positions.!=PT(0.0, 0.0)]
    max_pos = PT(maximum(first.(nz_pos)), maximum(last.(nz_pos)))
    min_pos = PT(minimum(first.(nz_pos)), minimum(last.(nz_pos)))
    lat_range = first(max_pos) - first(min_pos)
    lon_range = last(max_pos) - last(min_pos)

    initial_positions[initial_positions.==PT(0.0, 0.0)] = map(
        x -> (lat_range, lon_range) .* rand(PT) .+ (first(min_pos), last(min_pos)),
        1:(length(initial_positions)-length(nz_pos)),
    )

    @info "calculating node locations with SFDP"
    K = get(kwargs, :K, 0.1)
    network = NetworkLayout.SFDP.layout(
        a,#(a .- 1) .* -1,
        NetworkLayout.SFDP.Point2f0,
        tol = 1.0,
        C = 0.0002,
        K = K,
        iterations = 100,
        fixed = fixed,
        startpositions = initial_positions,
    ) # generate 2D layout

    set_prop!(g, :x, first.(network))
    set_prop!(g, :y, last.(network))
    name_accessor = get(kwargs, :name_accessor, get_name)
    set_prop!(g, :name, name_accessor.(buses))

    return g
end

# plot
#draw(PDF("net.pdf"), gplot(g, first.(network), last.(network), nodefillc=node_colors))

function plot_net(g; kwargs...)
    p = plot()
    plot_net!(p, g; kwargs)
end

function plot_net!(p::Plots.Plot, g; kwargs...)
    x = get_prop(g, :x)
    y = get_prop(g, :y)
    m = lonlat_to_webmercator.(x, y)
    lines = get(kwargs, :lines, false)
    linealpha = get(kwargs, :linealpha, 0.2)
    nodesize = get(kwargs, :nodesize, 2.0)
    linewidth = get(kwargs, :linewidth, 1.0)
    linecolor = get(kwargs, :linecolor, "black")
    buffer = get(kwargs, :buffer, 0.75e5)
    size = get(kwargs, :size, (600,800))
    xlim = get(kwargs, :xlim, (minimum(first.(m)) - buffer, maximum(first.(m)) + buffer))
    ylim = get(kwargs, :ylim, (minimum(last.(m)) - buffer, maximum(last.(m)) + buffer))

    if lines
        @info "plotting lines"
        for e in edges(g)
            p = plot!(
                p,
                [first(m[e.src]), first(m[e.dst])],
                [last(m[e.src]), last(m[e.dst])],
                color = linecolor,
                label = "",
                alpha = linealpha,
                linewidth = linewidth,
            )
        end
    end

    @info "plotting nodes"
    shownodelegend = get(kwargs, :shownodelegend, true)
    group = shownodelegend ? get_prop(g, :group) : ["node" for g in get_prop(g, :group)]
    p = scatter!(
        p,
        first.(m),
        last.(m),
        markercolor = get(kwargs, :nodecolor, get_prop(g, :nodecolor)),
        markeralpha = get(kwargs, :nodealpha, get_prop(g, :alpha)),
        markerstrokecolor = get(kwargs, :nodecolor, get_prop(g, :nodecolor)),
        markersize = nodesize,
        group = group,
        hover = get_prop(g, :name),
        legend = shownodelegend,
        legend_ontcolor = :red,
        xlim = xlim,
        ylim = ylim,
        size = size
    )
    return p
end

# borrowed from
function lonlat_to_webmercator(xLon, yLat)

    # Check coordinates are in range
    abs(xLon) <= 180 || throw("Maximum longitude is 180.")
    abs(yLat) < 85.051129 || throw(
        "Web Mercator maximum lattitude is 85.051129. This is the lattitude at which the full map becomes a square.",
    )

    # Ellipsoid semi-major axis for WGS84 (metres)
    # This is the equatorial radius - the Polar radius is 6356752.0
    a = 6378137.0

    # Convert to radians
    λ = xLon * 0.017453292519943295    # λ = xLon * π / 180
    ϕ = yLat * 0.017453292519943295    # ϕ = yLat * π / 180

    # Convert to Web Mercator
    # Note that:
    # atanh(sin(ϕ)) = log(tan(π/4 + ϕ/2)) = 1/2 * log((1 + sin(ϕ)) / (1 - sin(ϕ)))
    x = a * λ
    y = a * atanh(sin(ϕ))

    return x, y
end

function lonlat_to_webmercator(xy::Tuple)
    return lonlat_to_webmercator(first(xy), last(xy))
end

function lonlat_to_webmercator(xy::Shapefile.Point)
    return lonlat_to_webmercator(xy.x, xy.y)
end

function lonlat_to_webmercator(shp::Vector{Union{Missing, Shapefile.Polygon}})
    return [lonlat_to_webmercator(polygon) for polygon in shp]
end

function lonlat_to_webmercator(poly::Shapefile.Polygon)
    (left, bottom) = lonlat_to_webmercator(poly.MBR.left, poly.MBR.bottom)
    (right, top) = lonlat_to_webmercator(poly.MBR.right, poly.MBR.top)
    return Shapefile.Polygon(
        Shapefile.Rect(left, bottom, right, top),
        poly.parts,
        [Shapefile.Point(x, y) for (x, y) in lonlat_to_webmercator.(poly.points)]
    )
end

function plot_components!(p, components, color, markersize, label)
    labels = get_name.(components)
    lat_lons = get_ext.(get_bus.(components))
    xy =
        x -> lonlat_to_webmercator(
            parse.(Float64, (get(x, "longitude", get(x, "lon", DEFAULT_LON)), get(x, "latitde", get(x, "lat", DEFAULT_LAT)))),
        )
    scatter!(
        p,
        first.(xy.(lat_lons)),
        last.(xy.(lat_lons)),
        color = color,
        markersize = markersize,
        hover = labels,
        label = label,
    )
end

component_locs(components) = collect(
    zip(
        [parse(Float64, get(x, "latitude", get(x, "lat", DEFAULT_LAT))) for x in get_ext.(get_bus.(components))],
        [parse(Float64, get(x, "longitude", get(x, "lon", DEFAULT_LON))) for x in get_ext.(get_bus.(components))]
    ),
)
component_locs(components::PowerSystems.IS.FlattenIteratorWrapper{Bus}) = collect(
    zip(
        [parse(Float64, get(x, "latitude", get(x, "lat", DEFAULT_LAT))) for x in get_ext.(components)],
        [parse(Float64, get(x, "longitude", get(x, "lon", DEFAULT_LON))) for x in get_ext.(components)]
    ),
)
