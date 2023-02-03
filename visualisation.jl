import automaton.jl
import parametrization.jl


p = 3
q = 3
r = 4

Triangle = buildFundTriangle(p,q,r)


θ = π
α = π

My_Rep = parametrization(θ,α,p,q,r)






function stereographic!(point3R,point3C)

    c = (-real(point3C[1]) + 1)

    point3R[1] = imag(point3C[1]) / c
    point3R[2] = real(point3C[2]) / c
    point3R[3] = imag(point3C[2]) / c
end

function siegel!(point2R,point3C)
    x =  point3C[2] / (-point3C[1] + 1)

    point2R[1] = real(x)
    point2R[2] = imag(x)
end