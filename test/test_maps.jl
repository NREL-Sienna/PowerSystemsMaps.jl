

sys = make_test_sys()

@testset "test maps" begin
    g = make_graph(sys, K = 0.01)
    @test typeof(g) == PSM.MetaGraphs.MetaGraph{Int64, Float64}
    @test length(PSM.get_prop(g, :x)) == 200
    shp = PSM.Shapefile.shapes(PSM.Shapefile.Table(joinpath(TEST_DIR, "test_data", "IL_BNDY_County", "IL_BNDY_County_Py.shp")))
    shp = PSM.lonlat_to_webmercator(shp) #adjust coordinates
    @test length(shp) == 102

    # plot a map from shapefile
    p = plot(
        shp,
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
        g,
        nodesize = 3.0,
        linecolor = "blue",
        linewidth = 0.6,
        lines = true,
        #nodecolor = "red",
        nodealpha = 1.0,
        shownodelegend = true,
        size = (1500,800),
        buffer = 0.4e4
    )
    @test typeof(p) == PSM.Plots.Plot{PSM.Plots.GRBackend}
end