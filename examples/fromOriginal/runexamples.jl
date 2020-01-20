include("./../../src/Lar.jl")

include("./randomarrangement2d.jl")
include("./randomshapes.jl")
include("./stresstest2d.jl")
include("./svg2lar.jl")

randomarrangement2d(100)
randomshapes(100)
stresstest2d()

#svgarrangement("./test/svg/curved.svg")
svgarrangement("./test/svg/holes.svg")
#svgarrangement("./test/svg/howto.svg")
#svgarrangement("./test/svg/howto2.svg")
#svgarrangement("./test/svg/Lar.svg")
#svgarrangement("./test/svg/new.svg")
#svgarrangement("./test/svg/paths.svg")
#svgarrangement("./test/svg/twopaths.svg")
