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
    @testset "Case 1" begin
        # Gravity Constant
        μ = SPICE.bodvrd("SUN", "GM", 1)[1]

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

        @test num_sol == 1
        @test v1[:,1] ≈ [-12.146901903784691, -30.249196666164899, 1.732823081867305]
        @test v2[:,1] ≈ [-5.062557143526494, 22.773747568054521, -1.043526507750940]

    end

    @testset "Case 2" begin
        # Example used in pykep Tutorial
        μ = 1.32712440018e+11
        r1 = [-25216645.728283768, 144924279.08132498, -38.276915766745136]
        r2 = [177909722.13822687, -105168473.55967535, -6575244.5882079229]
        tof = 55296000

        # Solve Lambert's Problem
        v1, v2, num_sol = lambert_problem(r1,r2,tof,μ,2)

        @test num_sol == 3
        @test v1[:,1] ≈ [-35.499818966294406, 0.40742940300237478, 1.4595203390221081]
        @test v2[:,1] ≈ [3.0376984971968541, 27.064521805431243, -0.31914600319656404]
        @test v1[:,2] ≈ [-27.81451513546858, -15.91479886802767, 1.2599207890195353]
        @test v2[:,2] ≈ [16.828356734556462, 14.965476740707376, -0.80053299320499104]
        @test v1[:,3] ≈ [-26.997843155205635, -17.749606456280773, 1.2394297881268258]
        @test v2[:,3] ≈ [18.383111848123743, 13.641214500991762, -0.85508959337278111]

    end

end
