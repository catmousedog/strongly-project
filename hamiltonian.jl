using MPSKit, TensorKit, TensorOperations, MPSKitModels


function hamiltonian(eltype=ComplexF64, symmetry=ℤ{1}, lattice=InfiniteChain(1); J=1.0, θ=0.0, spin=1 // 2)
    XX = sigma_xx(eltype, symmetry; spin=spin)
    YY = sigma_yy(eltype, symmetry; spin=spin)
    ZZ = sigma_zz(eltype, symmetry; spin=spin)
    return @mpoham sum(J * (cos(θ) * (XX{i,j} + YY{i,j} + ZZ{i,j}) + sin(θ) * (XX{i,j} * XX{i,j} + YY{i,j} * YY{i,j} + ZZ{i,j} * ZZ{i,j}))
                       for (i, j) in nearest_neighbours(lattice))
end

N=4

H = hamiltonian(J=1.0, θ=0.0)
Ψ = FiniteMPS(randn, ComplexF64, N, ℂ^2, ℂ^10)
Ψ, envs, δ = find_groundstate(Ψ, H, DMRG())
E = sum(expectation_value(Ψ, H))
println(E)