using LinearAlgebra
using SparseArrays
using Memoize

import Base.|>

const ⊗ = kron

# bra ket 0-based 
@memoize Dict{Tuple{Int64,Int64}, Any} (
ket(dim, i) = sparsevec([1+mod(i,dim)],[1.0+0*im],dim)
)
ket(i) = ket(d, i)
bra(dim, i) = ket(dim, i)'
bra(i) = ket(i)'
    
# Weyl-Heisenberg group
@memoize Dict{Tuple{Int64,Int64,Int64}, Any} (
function T(dim,i,j::Int64)::SparseMatrixCSC{ComplexF64,Int64}
    ω = exp(2pi*im/dim)
    τ = -exp(pi*im/dim) # w = tau^2 
    ⊕ = (x,y) -> mod(x+y-1,dim)+1

    ret = spzeros(ComplexF64, dim, dim)

    i = mod(i,dim)
    j = mod(j,dim)

    if i==0 && j==0 
        return ret + I
    end

    if j==0 # shift
        for k=1:dim
            ret[k ⊕ i, k] = 1
        end

        return ret      
    end

    if i==0  # clock
        for k=1:dim
            ret[k, k] = ω^((k-1)*j)
        end

        return ret
    end

    return τ^(i*j) * T(dim,i,0) * T(dim,0,j)
end
)  
T(i,j) = T(d,i,j)
  
@memoize Dict{Tuple{Int64}, Any} (
function Psym(d)
    ret = spzeros(ComplexF64,d^2,d^2)

    for i=1:d, j=1:i
        if i==j
        b = sparsevec( [1 + (d+1)*(i-1)] , [1], d^2 )
        ret += b*b'
        continue
        end
        b = sparsevec( [d*(i-1)+j, d*(j-1)+i] , [1,1], d^2 )          
        ret += b*b'/2
    end
    return ret
end
)
Psym() = Psym(d)

# out = index of subsystem to trace out (1 or 2)
function partial_trace(M, out::Int64)
    # todo rectangular case
    # todo block tracing?
    d = size(M)[1] |> sqrt |> round |> Int 
    res = zeros(ComplexF64,d,d)
    if out == 1 
        for i=1:d
            res += (bra(d,i)⊗I(d)) * M * (ket(d,i)⊗I(d))
        end
    else # out == 2
        for i=1:d
            res += (I(d)⊗bra(d,i)) * M * (I(d)⊗ket(d,i))
        end
    end
    return res    
end