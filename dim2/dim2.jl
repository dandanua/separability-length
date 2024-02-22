# using OhMyREPL

d = 2

include("../lib.jl")

x = Dict()
y = Dict()

# x[0] = [complex(-0.010141, 0)
#     complex(0.202484, -0.979233)]
# y[0] = [complex(-0.00906845, 0)
#         complex(0.414726, -0.909901)]

# x[1] = [complex(-0.820867, 0)
#         complex(0.481841, -0.306607)]
# y[1] = [complex(0.820585, 0)
#         complex(-0.480642, 0.30923)]

# x[2] = [complex(0.816714, 0)
#         complex(0.515464, 0.259376)]
# y[2] = [complex(0.818038, 0)
#         complex(0.513417, 0.259261)]

# x[3] = [complex(0.811821, 0)
#         complex(-0.0288314, -0.583194)]
# y[3] = [complex(-0.810784, 0)
#         complex(0.0269192, 0.584726)]

# x[0] =
# [1 + 0 * im;
# 0 + 0 * im]
# y[0] =
# [0.866025 + 0 * im;
# 0.5 + 0 * im]

# x[1] =
# [0.774208 + 0 * im;
# -0.198412 + -0.601028 * im]
# y[1] =
# [0.645497 + 0 * im;
# -0.736991 + -0.200443 * im]

# x[2] =
# [0.607436 + 0 * im;
# -0.0338794 + 0.793646 * im]
# y[2] =
# [0.645497 + 0 * im;
# -0.685221 + 0.33735 * im]

# x[3] =
# [0.177832 + 0 * im;
# 0.979533 + -0.094293 * im]
# y[3] =
# [0.645497 + 0 * im;
# 0.751392 + -0.136907 * im]

x[0] =
    [1
        0]
y[0] =
    [0.866025
        0.5]

x[1] =
    [0.472821
        0.134042 + 0.870904 * im]
y[1] =
    [0.645497
        -0.597623 + 0.475584 * im]

x[2] =
    [0.862153
        -0.274544 + -0.425815 * im]
y[2] =
    [0.645497
        -0.753041 + -0.127524 * im]

x[3] =
    [0.182024
        0.952186 + -0.245374 * im]
y[3] =
    [0.645497
        0.679844 + -0.34806 * im]

function check(t)
    Pt = zeros(ComplexF64, d^2, d^2)
    for i = 0:3
        Pt += (x[i] * x[i]') ⊗ (y[i] * y[i]')
    end

    return Pt / d^2 - 2t / d * Psym() - (1 - t * (d + 1)) / d^2 * I(d^2) |> norm
end

function G(x)
    G = zeros(ComplexF64, 4, 4)
    for i = 1:4, j = 1:4
        G[i, j] = x[i-1]' * x[j-1]
    end
    return G
end

function solG()
    x = Dict()
    # x[0] = 1/sqrt(6)*[sqrt(3+sqrt(3)); exp(im*pi/4)*sqrt(3-sqrt(3))]
    x[0] = 1 / sqrt(6) * [-sqrt(3 - sqrt(3)); exp(im * pi / 4) * sqrt(3 + sqrt(3))]
    for i = 0:1, j = 0:1
        x[i+2j] = T(i, j) * x[0]
    end
    G = zeros(ComplexF64, 4, 4)
    for i = 1:4, j = 1:4
        G[i, j] = x[i-1]' * x[j-1]
    end
    @show G
end


