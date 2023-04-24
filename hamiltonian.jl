
module hamiltonian

export AKLT_hamiltonian

using MPSKitModels, MPSKit, TensorKit, TensorOperations

function AKLT_hamiltonian(eltype=ComplexF64, symmetry=ℤ{1}, lattice=InfiniteChain(1); J=1.0, θ=0.0, spin=1)
    SS = sigma_exchange(eltype, symmetry; spin=spin)
    return @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * SS{i,j} * SS{i,j}) for (i, j) in nearest_neighbours(lattice))
end

end