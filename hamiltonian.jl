
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

"
Breaks all symmetries
g → +∞: no phase transition (stays gapped)
"
function HAFM_staggered(lattice=InfiniteChain(1); spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    G = sigma_z(spin=spin) # still need to add staggered factor!

    H = @mpoham sum(J * SS{i,j} + g * G{i} for (i, j) in nearest_neighbours(lattice))

    return Hamiltonian(H, spin)
end

"
Doesn't break ℤ2 × ℤ2 and ℤ2_Time.
g → +∞: phase transition
"
function HAFM_zz(lattice=InfiniteChain(1); spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    G = sigma_z(spin=spin) * sigma_z(spin=spin)

    H = @mpoham sum(J * SS{i,j} + g * G{i} for (i, j) in nearest_neighbours(lattice))

    return Hamiltonian(H, spin)
end
"
Doesn't break ℤ2 × ℤ2
g → +∞: phase transition
"
function HAFM_xy(lattice=InfiniteChain(1); spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    G = sigma_x(spin=spin) * sigma_y(spin=spin)

    H = @mpoham sum(J * SS{i,j} + g * G{i} for (i, j) in nearest_neighbours(lattice))

    return Hamiltonian(H, spin)
end
