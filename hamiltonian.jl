using MPSKit, TensorKit, TensorOperations, LinearAlgebra

J = 1.0
θ = 0.0
c = cos(θ)
s = sin(θ)
N = 3
@assert N > 1

σx = TensorMap(zeros, ComplexF64, ℂ^2 ← ℂ^2)
σx[2, 1] = 1
σx[1, 2] = 1
σy = TensorMap(zeros, ComplexF64, ℂ^2 ← ℂ^2)
σy[1, 2] = -im
σy[2, 1] = im
σz = TensorMap(zeros, ComplexF64, ℂ^2 ← ℂ^2)
σz[1, 1] = 1
σz[2, 2] = -1

# ℂ^2 ← ℂ^2 operator that acts on site i
function on_site(operator, i::Int64)
    return mapfoldl(⊗, 1:N) do j
        j == i ? operator : id(ℂ^2)
    end
end

function local_hamiltonian(i::Int64)
    # periodic boundary conditions
    i_next = (i % N) + 1
    # interaction term
    h_int = on_site(σx, i) * on_site(σx, i_next) + on_site(σy, i) * on_site(σy, i_next) + on_site(σz, i) * on_site(σz, i_next)
    return J * h_int
end

function hamiltonian()
    H = mapfoldl(+, 1:N) do i
        local_hamiltonian(i)
    end
    return H
end

H = hamiltonian()
Ψ = FiniteMPS(randn, ComplexF64, N, (ℂ^2)^N, ℂ^10)
Ψ, envs, δ = find_groundstate(Ψ, H, Dmrg())
E = sum(expectation_value(Ψ, H))
println(E)