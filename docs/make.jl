using Documenter, PowerSystemsMaps

makedocs(
    modules = [PowerSystemsMaps],
    format = Documenter.HTML(),
    sitename = "PowerSystemsMaps.jl",
    authors = "Clayton Barrows",
    pages = ["Home" => "index.md", "Function Index" => "api.md"],
)

Documenter.deploydocs(
    repo = "github.com/NREL-SIIP/PowerSystemsMaps.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"],
)
