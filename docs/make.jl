using Documenter

push!(LOAD_PATH, "../src")
using FiniteStateAutomata

makedocs(
    sitename = "FAST",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Usage" => [
           "Creating WFSTs" => "create.md",
           "Operations" => "operations.md"
        ],
        "Library" => "library.md",
        "Contributing" => "contrib.md"
    ]
)
