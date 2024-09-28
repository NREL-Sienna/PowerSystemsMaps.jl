
sys = make_test_sys()
@testset "test maps" begin
    g = make_graph(sys; K = 0.01)
    @test typeof(g) <: PSM.MetaGraphsNext.MetaGraph
    @test length(PSM.get_prop(g, :x)) == 200
    shapefile = joinpath(TEST_DIR, "test_data", "IL_BNDY_County", "IL_BNDY_County_Py.shp")
    shp = PSM.Shapefile.shapes(PSM.Shapefile.Table(shapefile))
    shp = PSM.lonlat_to_webmercator(shp) #adjust coordinates
    @test length(shp) == 102

    # plot a map from shapefile
    p = plot(
        shp;
        fillcolor = "grey",
        background_color = "white",
        linecolor = "darkgrey",
        axis = nothing,
        border = :none,
        label = "",
        legend_font_color = :red,
    )

    # plot the network on the map
    p = plot_net!(
        p,
        g;
        nodesize = 3.0,
        linecolor = "blue",
        linewidth = 0.6,
        lines = true,
        #nodecolor = "red",
        nodealpha = 1.0,
        shownodelegend = true,
        size = (1500, 800),
        buffer = 0.4e4,
    )
    @test typeof(p) == PSM.Plots.Plot{PSM.Plots.GRBackend}

    p = plot_map(
        sys,
        joinpath(TEST_DIR, "test_data", "IL_BNDY_County", "IL_BNDY_County_Py.shp");
        map_fillcolor = "grey",
        map_background_color = "white",
        map_linecolor = "darkgrey",
        map_axis = nothing,
        map_border = :none,
        map_label = "",
        map_legend_font_color = :red,
        nodesize = 3.0,
        linecolor = "blue",
        linewidth = 0.6,
        lines = true,
        nodealpha = 1.0,
        shownodelegend = true,
        size = (1500, 800),
        buffer = 0.4e4,
    )
    @test typeof(p) == PSM.Plots.Plot{PSM.Plots.GRBackend}
end
