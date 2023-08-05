push!(LOAD_PATH,"../src/")

using Documenter
using ForceAnalysis

makedocs(
    sitename = "ForceAnalysis.jl",
    doctest = false,
    authors = "Oliver Lindemann",
    pages = [
        "index.md",
    ],
)

#deploydocs(; repo = "github.com/lindemann09/ForceAnalysis.jl", push_preview = true)