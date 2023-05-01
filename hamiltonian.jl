
function bilinear_biquadratic_hamiltonian(lattice=InfiniteChain(1); spin=1, J=1.0, θ=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    return @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * SS{i,j} * SS{i,j}) for (i, j) in nearest_neighbours(lattice))
end

function optimize_groundstate(; N=2, spin=1, J=0.5, θ=0.0, bond::Int64=10, maxiter::Int64=100)
    @assert N > 1

    H = bilinear_biquadratic_hamiltonian(spin=spin, J=J, θ=θ)

    if (N == Inf)
        Ψ = InfiniteMPS(ℂ^Int(2 * spin + 1), ℂ^(bond))

        algorithm = VUMPS(maxiter=maxiter)

        Ψ, envs, δ = find_groundstate(Ψ, H, algorithm)
    else
        Ψ = FiniteMPS(randn, ComplexF64, N, ℂ^Int(2 * spin + 1), ℂ^bond)
        
        algorithm = DMRG(maxiter=maxiter)

        Ψ, envs, δ = find_groundstate(Ψ, H, algorithm)
    end

    return Ψ
end
