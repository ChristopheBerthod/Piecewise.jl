using Documenter
using Piecewise, PiecewiseHilbert, PiecewiseLorentz

makedocs(
    repo = Documenter.Remotes.GitHub("ChristopheBerthod", "Piecewise"),
    sitename = "Piecewise",
    format = Documenter.HTML(prettyurls = false,
        edit_link = "main",
        inventory_version = "v0.1.0"
        ),
    modules = [Piecewise, PiecewiseHilbert, PiecewiseLorentz],
    pages = [
        "index.md",
        "hilbert.md",
        "lorentz.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/ChristopheBerthod/Piecewise.git",
    devbranch = "main"
)
