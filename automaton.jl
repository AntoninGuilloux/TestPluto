using IntervalArithmetic, LinearAlgebra, Symbolics

function buildFundTriangle(p,q,r)
    if p == Inf
        return [[im,-1],[-1,1],[1,im]]#outside of the triangle
    else
        P = [-1,1]
        Q = [exp(im*(2*pi-pi/p)),exp(im*(pi-pi/p))]
        #building of the 3rd geodesic
        cp = cos(π/p)
        cq = cos(π/q)
        cr = cos(π/r)
        
        r_a0 = (cp*tan(pi/p))/(cq+cp*cr)
        d = r_a0*cr
        
        
        α = sqrt(d^2+1)/r_a0
        
        r_a = 1/sqrt(α^2-1)

        dist = α*r_a
        dist2 = dist^2
        
        
        x = -1*dist*cos(atan(d))
        y = dist*sin(atan(d))
        
        tmp0 = 1-r_a^2
        
        tmp1 = (1/2+tmp0/(2*dist2))*(x+im*y)
        tmp2 = 1/2*sqrt(1-(tmp0^2/dist2^2))*(y-im*x)
        
        
        R = [tmp1+tmp2,tmp1-tmp2]
        return [P,Q,R]
    end
end
	
function checkCircle(x,y,r,α)
    print("\n")
    print("\n")
    #orthogonal to unit
    print(x^2+y^2-(1+r^2)) 
    print("\n")
    # angle with (-1,1)
    print(pi/acos(y/r))
    print("\n")
    # angle with the other: y=-x/cos(α)
    print(pi/acos(cos(α)*(-1*x*tan(α)-y)/r))
    print("\n")
    print(abs((x+y*cos(α))/sqrt(1+(cos(α))^2))-r)
end


#given 2 points and a third, compute the image of the third by the homography fixing the others
#homographies preserves the cross ratios
function getImage(a,b,c)
    return (c*(a+b)-2*a*b)/(2*c-(a+b))
end

dummyAuto = [[2,3,4],
                [0,3,4],
                [5,0,4],
                [6,7,0],
                [0,0,4],
                [0,3,8],
                [9,0,0],
                [0,7,0],
                [0,0,10],
                [11,7,0],
                [0,12,0],
                [5,0,13],
                [14,0,0],
                [0,5,8]
        ];

function findQuarter(min,max)
    if 0<=min && min<pi && 0<=max && max<pi 
        return 1
    elseif pi/2<=min && min<pi*3/2 && pi/2<=max && max<pi*3/2 
        return 2
    elseif pi<=min && min<2*pi && pi<=max && max<2*pi 
        return 3
    elseif pi*3/2<=min || min<pi/2 && pi*3/2<=max || max<pi/2 
        return 4
    end
end

function between(a,b,c)
    aa = mod(angle(a),2*pi)
    ab = mod(angle(b),2*pi)
    if ab<aa
        ab+=2*pi
    end
    ac = mod(angle(c),2*pi)
    if ac<aa
        ac+=2*pi
    end
    return aa<ac && ac<ab   
end

function betweenI(a,b,c,d)
    aa = mod(angle(a),2*pi)
    ab = mod(angle(b),2*pi)
    if ab<aa
        ab+=2*pi
    end
    ac = mod(angle(c),2*pi)
    if ac<aa
        ac+=2*pi
    end
    ad = mod(angle(d),2*pi)
    while ad<ac
        ad+=2*pi
    end
    return (aa<=ac && ac<ab)# && ab-ac>(ad-ac)/10) || (aa=<ad && ad=<ab && ad-aa>(ad-ac)/10)
end

function intervalSize(a,b)
    aa = mod(angle(a),2*pi)
    ab = mod(angle(b),2*pi)
    if ab<aa
        ab+=2*pi
    end
    return (ab-aa)
end

function toAngles(a,b)
    aa = mod(angle(a),2*pi)
    ab = mod(angle(b),2*pi)
    if ab<aa
        ab+=2*pi
    end
    return [aa/(2*pi),ab/(2*pi)]
end

function allAngles(intervals)
    l = [[0.0],[0.0],[0.0]]
    for i in [1,2,3]
        l[i] = toAngles(intervals[i][1],intervals[i][2])
    end
    return l
end

function finished(min,max,interval)        
    between(interval[1],interval[2],min) && between(interval[1],interval[2],max)
end

function updateIntevals(symmetryAxisIndex,intervals)
    fix1 = intervals[symmetryAxisIndex][1]
    fix2 = intervals[symmetryAxisIndex][2]
    for i in [1,2,3] 
        if i != symmetryAxisIndex
            tmp = getImage(fix1,fix2,intervals[i][1])
            intervals[i][1] = getImage(fix1,fix2,intervals[i][2])
            intervals[i][2] = tmp
        end
    end
    intervals[symmetryAxisIndex][1] = fix2
    intervals[symmetryAxisIndex][2] = fix1
end
            
function allSize(intervals)
    l = [0.0,0.0,0.0]
    for i in [1,2,3]
        l[i] = intervalSize(intervals[i][1],intervals[i][2])
    end
    return l
end

function findWord(min,max,auto,baseTri)
    currentState = 1
    currentIntervals = copy.(baseTri)
    emin = exp(im*(min))
    emax = exp(im*(max))
    size = intervalSize(emin,emax)
    word = ""
    while true
        b = false
        for i in [1,2,3]
            if auto[currentState][i]!=0
                if betweenI(currentIntervals[i][1],currentIntervals[i][2],emin,emax)
                    b = true
                    word*=string(i)
                    if intervalSize(currentIntervals[i][1],currentIntervals[i][2])< size
                        return word
                    else
                        updateIntevals(i,currentIntervals)
                        currentState = auto[currentState][i]
                    end
                    break
                end
            end
        end
        if !b
            @assert(false)
        end
    end
end