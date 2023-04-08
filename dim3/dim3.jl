# using OhMyREPL

d = 3

include("../lib.jl")


u(t) = [1; λ(t); λ(t)]
v(t) = [α(t); β(t); β(t)]

α(t) = sqrt( (5+4t+4*√(1+7t-8t^2)) / 81 )

λ(t) = (-1+√(1+4*ρ(t))) / 2

ρ(t) = (1+2t-√(1+7t-8t^2)) / (1-4t)

β(t) = -α(t)*(λ(t)+1)

x(t) = √3 * norm(v(t)) * u(t)
y(t) = v(t) / norm(v(t))

function Φ(M1,M2)
    ret = zeros(ComplexF64, 9, 9)
    for i=1:d, j=1:d
        M1ij = T(i,j) * M1 * T(i,j)'
        M2ij = T(i,j) * M2 * T(i,j)'
        ret += M1ij ⊗ M2ij
    end
    return ret/d^2
end 

function check(t)
    Pt = Φ(x(t)*x(t)', y(t)*y(t)')
    return Pt - 2t/3*Psym()-(1-4t)/9*I(9) |> norm
end
