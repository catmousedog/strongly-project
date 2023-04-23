
module hamiltonian

export AKLT_hamiltonian

using MPSKitModels, MPSKit, TensorKit, TensorOperations

function AKLT_hamiltonian(eltype=ComplexF64, symmetry=ℤ{1}, lattice=InfiniteChain(1); J=1.0, θ=0.0, spin=1)
    XX = sigma_xx(eltype, symmetry; spin=spin)
    YY = sigma_yy(eltype, symmetry; spin=spin)
    ZZ = sigma_zz(eltype, symmetry; spin=spin)
    return @mpoham sum(J * (cos(θ) * (XX{i,j} + YY{i,j} + ZZ{i,j}) + sin(θ) * (XX{i,j} * XX{i,j} + YY{i,j} * YY{i,j} + ZZ{i,j} * ZZ{i,j}))
                       for (i, j) in nearest_neighbours(lattice))
end

end