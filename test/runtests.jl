using ApproxMasEq
using Test

N = 10
alpha = 2
K = 3
c = K * alpha
p0 = 0.5
seed = 1

eta = 1.0
alg_str = "FMS"
rf = rate_FMS_KSAT
rargs = [eta]

answ_CME = CME_KSAT(rf, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha)
answ_CDA = CDA_KSAT(rf, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha)

@testset "ApproxMasEq.jl" begin
    
    @test answ_CME.u[end][end] > 0 && answ_CME.u[end][end] < 1
    @test answ_CDA.u[end][end] > 0 && answ_CDA.u[end][end] < 1
end
