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

#(a*a - 2*a*b+b*b)/(4(c-(a+b)/2)) + (a+b)/2



function betterGetImage(tri)
    function partialFractionDecomp(l)
        a = l[1]
        b = l[2]
        return [-(a+b)/2,(a*a-2*a*b+b*b)/4, (a+b)/2]
    end
    l1 = partialFractionDecomp(tri[1])
    l2 = partialFractionDecomp(tri[2])
    l3 = partialFractionDecomp(tri[3])
    return function f(i,c)
        if i == 1
            return l1[2]/(c+l1[1])+l1[3]
        elseif i == 2
            return l2[2]/(c+l2[1])+l2[3]
        elseif i==3
            return l3[2]/(c+l3[1])+l3[3]
        else
            @assert(false)
        end
    end
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

function updateInterval(symmetryAxis,interval)
    min = interval[1]
    max = interval[2]
    interval[1] = getImage(symmetryAxis[1],symmetryAxis[2],max)
    interval[2] = getImage(symmetryAxis[1],symmetryAxis[2],min)
end

function updateInterval2(symmetryAxisIndex,interval)
    min = interval[1]
    max = interval[2]
    interval[1] = getImage2(symmetryAxisIndex,max)
    interval[2] = getImage2(symmetryAxisIndex,min)
end

function updateIntervals(symmetryAxisIndex,intervals)
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


function intersection(a,b,c,d)
    m=c
    M=d
    if between(c,d,b)
        m=b
    end
    if between(c,d,a)
        M=d
    end
    return [m,M]
end


function xfindWord(min,max,auto,baseTri)
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
                        updateIntervals(i,currentIntervals)
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


function locate(baseTri,target,localisation)
    localisation[1] = locateElem(baseTri,target[1])
    localisation[2] = locateElem(baseTri,target[2])
    
    if localisation[1]>3
        if localisation[1]==4
            if real(target[1])<0
                localisation[1] = 1
            else
                localisation[1] = 2
            end
        else
            localisation[1] = 3
        end
    end
    
    if localisation[2]>3
        if localisation[2]==4
            if real(target[2])<0
                localisation[2] = 3
            else
                localisation[2] = 1
            end
        else
            localisation[2] = 2
        end
    end
    
    #check the orientation
    if localisation[1] == localisation[2]
        if (localisation[1] == 1 && real(target[1])<real(target[2])) || 
            (localisation[1] !=1 && real(target[2])<real(target[1]))
            localisation[3] = 1
        else
            localisation[3] = 0
        end
    end
end



    

function locateElem(baseTri,target)
    if imag(target)<-1e-16
        return 1
    elseif imag(target)>1e-16 && abs(abs(real(target))-abs(real(baseTri[2][2])))<1e-16
        return 5
    elseif abs(imag(target))<1e-16
        return 4
    elseif real(target)>real(baseTri[2][2])
        return 2
    else
        return 3
    end
end

function completeWord(baseTri,auto,target,currentState,parity)
    localisation = [0,0,0]
    word = ""
    
    while true
#         if length(word)>200
#             print(word)
#             @assert(false)
#         end
        locate(baseTri,target,localisation)
        if localisation[1] == localisation[2] && localisation[3] == 1
            updateInterval(baseTri[localisation[1]],target)
            parity = mod(parity+1,2)
            word*=string(localisation[1])
#             if (auto[currentState][localisation[1]]==0)
#                 @assert(1==4)
#             end
            currentState = auto[currentState][localisation[1]]
        else
#             if mod(length(word),2) != w2
#                 loc = localisation[parity+1]
#                 word*=string(loc)
#                 if (auto[currentState][loc]==0)
#                     @assert(1==5)
#                 end
#                 currentState = auto[currentState][loc]
#             end
            return (word,currentState)
        end
    end
end




function findFirstPrefix(baseTri,auto,min,max)
    currentState = 1
    login = [min,max]
    emin = exp(im*(min))
    emax = exp(im*(max))
    target = [emin,emax]
    localisation = [0,0,0]
    word = ""
    
    while true
#         if length(word)>200
#             print(word)
#             @assert(false)
#         end
        locate(baseTri,target,localisation)
        if localisation[1] == localisation[2] && localisation[3] == 1
            updateInterval(baseTri[localisation[1]],target)
#             updateInterval2(localisation[1],target)
            word*=string(localisation[1])
            
#             if (auto[currentState][localisation[1]]==0)
#                 @assert(1==2)
#             end
            currentState = auto[currentState][localisation[1]]
        elseif mod(localisation[1]+1,3) == mod(localisation[2],3) #split
#             w2 = mod(length(word),2)
            
            split = 0
#             if localisation[1] == 1
#                 split = 0
#             elseif localisation[1] == 2
#                 split = imag(log(baseTri[2][2]))
#             else
#                 split = π
#             end
            
            if localisation[1] == 1
                split = 1
            elseif localisation[1] == 2
                split = baseTri[2][2]
            else
                split = -1
            end
            
#             min0 = imag(log(target[1]))
#             max0 = imag(log(target[2]))
            
            wordc=""
            
            if abs(split-min)>abs(max-split)
                target[2] = split
                (wordc,currentState) = completeWord(baseTri,auto,target,currentState,0)
            else
                target[1] = split
                (wordc,currentState) = completeWord(baseTri,auto,target,currentState,1)
            end
            return (word*wordc,currentState)
        else #finished (make the work pair)
#             if mod(length(word),2) == 1
#                 loc = mod(localisation[1],3)+1
#                 word*=string(loc)
#                 currentState = loc
#             end
            return (word,currentState)
        end
    end
end
    

function cutBaseInterval(min,max,triangle)
    
end

function explore(graph,distances,paths,distance,path,current,base)
    for i in [1,2,3]
        tmpPath =path*string(i)
        if graph[current][i] != 0
            if distances[base,graph[current][i]] == -1 || distances[base,graph[current][i]]>distance
                distances[base,graph[current][i]]=distance
                paths[base,graph[current][i]]=tmpPath
                explore(graph,distances,paths,distance+1,tmpPath,graph[current][i],base)
            end
        end
    end
end




function dijkstra(auto)
    chosenPrefixes = Array{String}(undef,length(auto))
    chosenCycles = Array{String}(undef,length(auto))
    distances = Array{Int64}(undef,length(auto),length(auto))
    paths = Array{String}(undef,length(auto),length(auto))
    
    for vertex in 1:length(auto)
        for i in 1:length(auto)
            distances[vertex,i]=-1
            paths[vertex,i]=""
        end
        explore(auto,distances,paths,1,"",vertex,vertex)
    end
    
    for vertex in 1:length(auto)
        index = -1
        mn = 0
        if distances[vertex,vertex]!=-1
            index = vertex
            mn = distances[vertex,vertex]
        end
        for i in 1:length(auto)
            if distances[vertex,i] != -1
                if distances[vertex,i]!=-1 && distances[i,i]!=-1 && (index==-1 || distances[i,i]<mn || 
                        distances[i,i]<=mn && distances[vertex,i]<=distances[vertex,index])
                    index = i
                    mn = distances[i,i]
                end
            end
        end
        if index != vertex
            chosenPrefixes[vertex] = paths[vertex,index]
            chosenCycles[vertex] = paths[index,index]
        else
            chosenPrefixes[vertex] = ""
            chosenCycles[vertex] = paths[index,index]
        end
    end
    
    
    return (chosenPrefixes,chosenCycles)
end

function permute(str)
    return str[2:length(str)]*str[1]
end



function findWord(baseTri,auto,chosenPrefixes,chosenCycles,min,max)
    (prefix,currentState) = findFirstPrefix(baseTri,auto,min,max)
    prefix*=chosenPrefixes[currentState]
    if mod(length(prefix),2) == 1
        prefix*=chosenCycles[currentState][1]
        return (prefix,permute(chosenCycles[currentState]))
    else
        return (prefix,chosenCycles[currentState])
    end
end