
struct Hamiltonian
    H
    spin
end

function bilinear_biquadratic_hamiltonian(lattice=InfiniteChain(1); spin=1, J=1.0, θ=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    quad = SS * SS
    H = @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * quad{i,j}) for (i, j) in nearest_neighbours(lattice))

    return Hamiltonian(H, spin)
end

function optimize_groundstate(ham; bond::Int64=10, maxiter::Int64=100)
    Ψ = InfiniteMPS(ℂ^Int(2 * ham.spin + 1), ℂ^(bond))
    Ψ, envs, δ = find_groundstate(Ψ, ham.H, VUMPS(maxiter=maxiter))
    return Ψ
end
