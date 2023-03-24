using LineBroadenings
using Documenter

DocMeta.setdocmeta!(LineBroadenings, :DocTestSetup, :(using LineBroadenings);
                    recursive = true)

makedocs(;
         modules = [LineBroadenings],
         authors = "Northwemko <yinyangxiangsheng@gmail.com> and contributors",
         repo = "https://github.com/lihua-cat/LineBroadenings.jl/blob/{commit}{path}#{line}",
         sitename = "LineBroadenings.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://lihua-cat.github.io/LineBroadenings.jl",
                                  edit_link = "main",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "line profile functions" => "profiles.md",
             "spectral line broadenings" => "spectra.md",
         ])

deploydocs(;
           repo = "github.com/lihua-cat/LineBroadenings.jl",
           devbranch = "main")
