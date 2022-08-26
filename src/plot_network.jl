
PT = GeometryBasics.Point{2, Float64}
const DEFAULT_LON = "-122.8230"
const DEFAULT_LAT = "37.8270"

function color_nodes!(
    g,
    sys,
    node_colors::Vector{Colors.RGB{Colors.N0f8}} = map(
        x -> Colors.RGB{Colors.N0f8}(1.0, 0.0, 0.0),
        1:nv(g),
    ),
)
    alpha = Vector()
    for (ix, b) in enumerate(get_components(Bus, sys))
        ext = get_ext(b)
        a =
            (haskey(ext, "lat") & haskey(ext, "lon")) ||
            (haskey(ext, "latitude") & haskey(ext, "longitude")) ? 1.0 : 0.1
        push!(alpha, a)
    end
    set_prop!(g, :nodecolor, node_colors)
    set_prop!(g, :alpha, alpha)
end

function color_nodes!(g, sys, color_by::Type{T}) where {T <: AggregationTopology}
    # Generate n maximally distinguishable colors in LCHab space.
    areas = get_components(color_by, sys)
    buses = get_components(Bus, sys)
    area_colors = Dict(
        zip(
            get_name.(areas),
            Colors.distinguishable_colors(length(areas), Colors.colorant"blue"),
        ),
    )
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
    initial_positions = Vector{PT}()
    for (ix, b) in enumerate(buses)
        if !isempty(get_ext(b))
            ext = get_ext(b)
            lat = tryparse(Float64, "$(get(ext, "latitude", get(ext, "lat", nothing)))")
            lon = tryparse(Float64, "$(get(ext, "longitude", get(ext, "lon", nothing)))")
            if !isnothing(lat) && !isnothing(lon)
                push!(initial_positions, PT(lon, lat))
            else
                fixed[ix] = false
            end
        end
    end
    sort_ids = vcat(findall(fixed), findall(.!fixed))
    orig_ids = sortperm(sort_ids)
    sorted_a = a[sort_ids, sort_ids]

    @info "calculating node locations with SFDP_Fixed"
    K = get(kwargs, :K, 0.1)
    network = sfdp_fixed(
        sorted_a,
        tol = 1.0,
        C = 0.0002,
        K = K,
        iterations = 100,
        fixed = true,
        initialpos = initial_positions,
    )[orig_ids] # generate 2D layout and sort back to order of a

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
    size = get(kwargs, :size, (600, 800))
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
        size = size,
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
        [Shapefile.Point(x, y) for (x, y) in lonlat_to_webmercator.(poly.points)],
    )
end

function plot_components!(p, components, color, markersize, label)
    labels = get_name.(components)
    lat_lons = get_ext.(get_bus.(components))
    plot_components!(p, labels, lat_lons, color, markersize, label)
end

function plot_components!(
    p,
    components::PowerSystems.FlattenIteratorWrapper{Bus},
    color,
    markersize,
    label,
)
    labels = get_name.(components)
    lat_lons = get_ext.(components)
    plot_components!(p, labels, lat_lons, color, markersize, label)
end

function get_longitude(ext)
    lon = get(ext, "longitude", get(ext, "lon", DEFAULT_LON))
    return parse(Float64, "$lon")
end
function get_latitude(ext)
    lat = get(ext, "latitude", get(ext, "lat", DEFAULT_LAT))
    return parse(Float64, "$lat")
end

function plot_components!(p, labels, lat_lons, color, markersize, label)
    xy = [lonlat_to_webmercator((get_longitude(l), get_latitude(l))) for l in lat_lons]
    x = first.(xy)
    y = last.(xy)
    scatter!(p, x, y, color = color, markersize = markersize, hover = labels, label = label)
end

component_locs(components) = collect(
    zip(
        [get_latitude(x) for x in get_ext.(get_bus.(components))],
        [get_longitude(x) for x in get_ext.(get_bus.(components))],
    ),
)
component_locs(components::PowerSystems.IS.FlattenIteratorWrapper{Bus}) = collect(
    zip(
        [get_latitude(x) for x in get_ext.(components)],
        [get_longitude(x) for x in get_ext.(components)],
    ),
)
