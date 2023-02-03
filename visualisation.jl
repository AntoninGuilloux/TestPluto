function stereographic!(point3R,point3C)
    #quaternion_div!(x_point)

    #a = point3C[1] + point3C[2]
    #b = point3C[1] - point3C[2]

    #x =  a + b*1im
    #x = x_point[1] * conj(b[1]) + x_point[2] * conj(b[2])
    #y = -b - a*1im
    #y = - x_point[1] * b[2] + x_point[2] * b[1]

    #c = -real(x) + 1
    #point3C /= point3C[3]
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