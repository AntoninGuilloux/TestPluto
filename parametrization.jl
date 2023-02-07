using LinearAlgebra

function fixedPoints(chosenCycles,Rep)
    # On cherche les points fixes des cycles choisis. A optimiser:
    # gérer les points fixes paraboliques
    # ne pas tout diagonaliser (utiliser la dynamique contractante)
    fPoints = Dict()
    for s in chosenCycles
        for i in 1:length(s)#On doit considérer toutes les permutations circulaires
            ss = join(circshift(collect(s),i))
            s2 = ss*ss #On double pour avoir des mots pairs
            Rep_s = wordToMatrix(s2,Rep)
            eigen_s = eigen(Rep_s)
            for (i, eig) in enumerate(eigen_s.values)
                if norm(eig)>1+10^-5
                    point_s = eigen_s.vectors[:,i]
                    fPoints[ss] = point_s/point_s[3]
                    break
                end
            end
        end
    end
    return fPoints
end



function wordToPoint(w,Rep,fixedPoints_computed)
    (p,s) = w
    Rep_p = wordToMatrix(p,Rep)
    point_s = fixedPoints_computed[s]
    point_w = Rep_p * point_s
    return point_w/point_w[3]
end

function wordxToPoint(w,Rep,fixedPoints_computed)
    (p,pp,s) = w
    p *= pp
    Rep_p = wordToMatrix(p,Rep)
    point_s = fixedPoints_computed[s]
    point_w = Rep_p * point_s
    return point_w/point_w[3]
end

function wordToMatrix(w,Rep)

    if (length(w)%2)==1
        
        throw(DomainError("The word w = $w must have even length"))
    else
        l = Int(length(w)/2)
        if (l%2)==1 #The length is not even anymore
            l=l-1
        end
        if w in keys(Rep)
            return Rep[w]
        elseif l == 0
            throw(DomainError("The word w = $w has not been found"))
        else
            w1 = w[1:l]
            w2 = w[l+1:end]
            M1 = wordToMatrix(w1,Rep)
            M2 = wordToMatrix(w2,Rep)
            M = M1*M2
            Rep[w]=M
            return M
        end
    end
end

function representation(θ,α,p,q,r)
    R1,R2,R3 = Symmetries(θ,α,p,q,r)
    Rep = Dict("" => Diagonal([1,1,1]),
                "11" => Diagonal([1,1,1]),
                "22" => Diagonal([1,1,1]),
                "33" => Diagonal([1,1,1]),
                "12" => R1*inv(R2),
                "23" => R2*inv(R3),
                "31" => R3*inv(R1),
                "21" => R2*inv(R1),
                "32" => R3*inv(R2),
                "13" => R1*inv(R3),
                )
    return Rep
end

function Symmetries(θ,α,p,q,r)
    sp, sp2, sq, sq2, sr, sr2 = precompute(p,q,r)
    θ_min, θ_max = θ_range(p)
    if θ <θ_min || θ >θ_max
        throw(DomainError("θ is not acceptable for this p"))
    else
        α_min, α_max = α_range(θ,sp2,sq2,sr2)
        if α <α_min || α >α_max
            throw(DomainError("α is not acceptable for this p, q, r, θ"))
        else
            sθ = sin(θ/2)
            ϕ_a = (asin( sp/sθ))
            ϕ_b = (asin( sq/sθ))
            ϕ_c = (asin( sr/sθ))
            
            sphib = sin(ϕ_b)
            sphic = sin(ϕ_c)
            
            tmp = cos(ϕ_a)*exp(-im*α)-cos(ϕ_b)*cos(ϕ_c)
            ρ_a = (abs(tmp)/(sphib*sphic))
            ν = angle(tmp)
            
            R3 = τ(θ)
            R2 = M(ϕ_b)τ(θ)*inv(M(ϕ_b))
            R1 = σ(-1*ν)T(ρ_a)M(ϕ_c)τ(θ)*inv(M(ϕ_c))*inv(T(ρ_a))σ(ν)
            return (R1,R2,R3)
        end
    end
end


function θ_range(p)
    return (2*π/p, 2*π*(1-1/p))
end

function α_range(θ,sp2,sq2,sr2)
	cos_max = cos_α_max(θ,sp2,sq2,sr2)
	if cos_max>1
		return 0
	elseif cos_max<-1
		return π
	else
		return (acos(cos_max),2*π-acos(cos_max))
	end
end

function cos_α_max(θ,sp2,sq2,sr2)
	sθ = sin(θ/2)
	sθ2 = sθ^2
	cos_α_m = sθ * (2 * sθ2 - (sp2+sq2+sr2))/(2*sqrt(π*(sθ2-sp2))) 
	return cos_α_m
end

function τ(x)
    return Diagonal([1,exp(im*x),1])
end
function σ(x)
    return Diagonal([exp(im*x),1,1])
end
function M(x)
    return [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1]
end
function T(x)
    return [x 0 sqrt(x^2-1); 0 1 0; sqrt(x^2-1) 0 x]
end

function precompute(p,q,r)
    sp = (sin(π/p))
    sp2 = sp*sp
    sq = (sin(π/q))
    sq2 = sq*sq
    sr = (sin(π/r))
    sr2 = sr*sr;
    return sp, sp2, sq, sq2, sr, sr2
end