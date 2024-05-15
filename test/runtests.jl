using ApproxMasterEqs
using Test

N = 10
alpha = 2
K = 3
seed = 1

eta = 1.0
rargs = [eta]

answ_CME = CME_KSAT(rate_FMS_KSAT, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, seed_g=seed)
answ_CDA = CDA_KSAT(rate_FMS_KSAT, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, seed_g=seed)

@testset "ApproxMasterEqs.jl" begin
    
    @test answ_CME.u[end][end] > 0 && answ_CME.u[end][end] < 1
    @test answ_CDA.u[end][end] > 0 && answ_CDA.u[end][end] < 1
end
