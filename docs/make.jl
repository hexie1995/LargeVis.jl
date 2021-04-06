using Documenter, LargeVis
push!(LOAD_PATH,"../src/")

makedocs(sitename="LargeVis.jl Documentation")
deploydocs(
    repo = "https://github.com/hexie1995/LargeVis.jl.git",
    devbranch = "main"
)
