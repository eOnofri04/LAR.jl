push!(LOAD_PATH,"../src/")

using Documenter, LAR


makedocs(
	format = :html,
	sitename = "LAR.jl",
	assets = ["assets/lar.css", "assets/logo.png"],
	pages = [
		"Home" => "index.md",
	]
)
