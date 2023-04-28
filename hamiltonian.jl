
module hamiltonian

export bilinear_biquadratic_hamiltonian, optimize_groundstate

using MPSKitModels, MPSKit, TensorKit, TensorOperations

function bilinear_biquadratic_hamiltonian(lattice=InfiniteChain(1); spin=1, J=1.0, θ=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    return @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * SS{i,j} * SS{i,j}) for (i, j) in nearest_neighbours(lattice))
end

function optimize_groundstate(; N=2, spin=1, J=0.5, θ=0.0, max_bond::Int64=10, maxiter::Int64=nothing)
    @assert N > 1

    H = bilinear_biquadratic_hamiltonian(spin=spin, J=J, θ=θ)

    if (N == Inf) #could do this with a different function too
        Ψ = InfiniteMPS(ℂ^Int(2 * spin + 1), ℂ^(max_bond))

        if (isnothing(maxiter))
            algorithm = VUMPS()
        else
            algorithm = VUMPS(maxiter=maxiter)
        end

        Ψ, envs, δ = find_groundstate(Ψ, H, algorithm)
    else
        Ψ = FiniteMPS(randn, ComplexF64, N, ℂ^Int(2 * spin + 1), ℂ^max_bond)
        
        if (isnothing(maxiter))
            algorithm = DMRG()
        else
            algorithm = DMRG(maxiter=maxiter)
        end

        Ψ, envs, δ = find_groundstate(Ψ, H, algorithm)
    end

    return Ψ
end

end
