
module hamiltonian

export bilinear_biquadratic_hamiltonian, optimize_finite_groundstate, optimize_infinite_groundstate

using MPSKitModels, MPSKit, TensorKit, TensorOperations

function bilinear_biquadratic_hamiltonian(lattice=InfiniteChain(1); spin=1, J=1.0, θ=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    return @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * SS{i,j} * SS{i,j}) for (i, j) in nearest_neighbours(lattice))
end

function optimize_finite_groundstate(H; N=2, spin=1, max_bond::Int64=10)
    @assert N > 1

    Ψ = FiniteMPS(randn, ComplexF64, N, ℂ^Int(2 * spin + 1), ℂ^max_bond)
    Ψ, envs, δ = find_groundstate(Ψ, H, DMRG())
    return Ψ
end

function optimize_infinite_groundstate(H; spin=1, max_bond::Int64=10)
    @assert spin >= 1

    Ψ = InfiniteMPS(2 * spin + 1, max_bond)
    Ψ, envs, δ = find_groundstate(Ψ, H, VUMPS())
    return Ψ
end

end
