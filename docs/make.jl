using NonSmoothDynamics
using Documenter

DocMeta.setdocmeta!(
  NonSmoothDynamics,
  :DocTestSetup,
  :(using NonSmoothDynamics);
  recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
  file for file in readdir(joinpath(@__DIR__, "src")) if
  file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
  modules = [NonSmoothDynamics],
  authors = "Tangi Migot <tangi.migot@gmail.com>",
  repo = "https://github.com/tmigot/NonSmoothDynamics.jl/blob/{commit}{path}#{line}",
  sitename = "NonSmoothDynamics.jl",
  format = Documenter.HTML(; canonical = "https://tmigot.github.io/NonSmoothDynamics.jl"),
  pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/tmigot/NonSmoothDynamics.jl")
