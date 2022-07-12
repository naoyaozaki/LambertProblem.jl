using LambertProblem
using Test
using Downloads

import SPICE

# Download generic kernels
genker_naif0012 = Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls")
genker_gm_de431 = Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc")
genker_de440 = Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets//de440.bsp")

# Load generic kernels
SPICE.furnsh(genker_naif0012) # Leap seconds kernel
SPICE.furnsh(genker_gm_de431) # Gravity Constant
SPICE.furnsh(genker_de440) # Planetary ephemeris kernel

@testset "LambertProblem.jl" begin
    # Gravity Constant
    μ = SPICE.bodvrd("SUN", "GM", 1)[1];

    # Flight Time
    et1 = SPICE.str2et("2031/03/01 00:00:00 UTC") # 地球出発日時, UTC
    et2 = SPICE.str2et("2032/01/01 00:00:00 UTC") # 火星到着日時, UTC
    tof = et2 - et1 # 飛行時間

    #　Earth and Mars Orbital Element
    state_earth, _ = SPICE.spkez(399, et1, "ECLIPJ2000", "NONE", 10)
    state_mars, _ = SPICE.spkez(4, et2, "ECLIPJ2000", "NONE", 10)

    r1 = state_earth[1:3]
    r2 = state_mars[1:3]

    # Solve Lambert's Problem
    v1, v2, num_sol = lambert_problem(r1,r2,tof,μ,2)

    println(v1,v2,num_sol)

end
