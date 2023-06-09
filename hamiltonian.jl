
"
Objects of this struct are denoted by 'ham' whilst mpoham objects are denoted by 'H'
"
struct Hamiltonian
    H
    spin
    unit_cells
end

function bilinear_biquadratic_hamiltonian(; spin=1, J=1.0, θ=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    quad = SS * SS
    H = @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * quad{i,j}) for (i, j) in nearest_neighbours(InfiniteChain(1)))

    return Hamiltonian(H, spin, 1)
end

"
Breaks all symmetries
g → +∞: no phase transition (stays gapped)
"
function HAFM_staggered(; spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    quad = SS * SS # weird bug
    Sz = sigma_z(spin=spin)

    H = @mpoham sum(J * SS{i,j} + 0.0 * quad{i,j} for (i, j) in nearest_neighbours(InfiniteChain(2)))
    Sz_staggered = @mpoham sum((-1)^(linearize_index(i)) * Sz{i} for i in vertices(InfiniteChain(2)))
    H = H + g*Sz_staggered

    return Hamiltonian(H, spin, 2)
end

"
Doesn't break ℤ2 × ℤ2, ℤ2_T and ℤ2_P.
g → +∞: product state |000...>
g → +∞: phase transition
"
function HAFM_zz(; spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    G = sigma_z(spin=spin) * sigma_z(spin=spin)

    H = @mpoham sum(J * SS{i,j} + g * G{i} for (i, j) in nearest_neighbours(InfiniteChain(1)))

    return Hamiltonian(H, spin, 1)
end

"
Doesn't break ℤ2 × ℤ2, ℤ2_T and ℤ2_P.
g → +∞: phase transition
"
function HAFM_xy(; spin=1, J=1.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    G = sigma_x(spin=spin) * sigma_y(spin=spin)

    H = @mpoham sum(J * SS{i,j} + g * G{i} for (i, j) in nearest_neighbours(InfiniteChain(1)))

    return Hamiltonian(H, spin, 1)
end
