
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
Attempt to break all symmetries and return to a trivial state (non-SPT)
"
function bilinear_biquadratic_hamiltonian_perturbed(lattice=InfiniteChain(1); spin=1, J=1.0, θ=0.0, g=0.0)
    SS = sigma_exchange(ComplexF64, ℤ{1}; spin=spin)
    quad = SS * SS
    G = sigma_x(spin=spin) * sigma_y(spin=spin)

    H = @mpoham sum(J * (cos(θ) * SS{i,j} + sin(θ) * quad{i,j}) + g * G{i} for (i, j) in nearest_neighbours(lattice))

    return Hamiltonian(H, spin)
end


