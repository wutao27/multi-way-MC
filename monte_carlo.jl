# one random walk simulation for the MC method
# x = Hx + b
using Distributions

function read_matrix(matrixPath, skipLen, flag)
    f = readdlm(matrixPath,skipstart=skipLen)
    H = sparse(round(Int64,f[:,1]),round(Int64,f[:,2]),f[:,3])
    d = diag(H)
    if flag
        H = H + H' - spdiagm(d)
    end

    n = size(H,1)
    H = speye(n)-spdiagm(1./d)*H
    validInd = sum(abs(H),2).!=0
    while sum(validInd)!=n
        H = H[validInd,validInd]
        n = sum(validInd)
        validInd = sum(abs(H),2).!=0
    end
    return H
end



function tran_tensor(HT, m)
    # HT is the transpose of the matrix H, we use HT instead of H is beacuse
    # Julia's sparse mtrix is CSC
    # each of the return trnsition matrix is coliumn stochastic
    HTplus = abs(HT)
    n = size(H,1)
    w = ones(1,n)
    P = []; p = []; indMap = []
    for k = m:-1:1
        yeta = w*HTplus
        tempP = (spdiagm(vec(w))*HTplus)*spdiagm(1./vec(yeta))
        push!(P, tempP)
        tempp = []
        for i = 1:n
            push!(tempp, Categorical(tempP[:,i].nzval))
        end
        push!(p, tempp)
        w = yeta
    end
    for i = 1:n
        push!(indMap, HT[:,i].rowval)
    end
    return P, p, indMap
end

function compute_var(H, h, b, P, p, indMap, chainLen, numTrial)
    n = size(H,1); m = length(P)
    ht = h'
    x = \(speye(n) - H, b)
    ex = (ht*x)[1]
    h = h/ex; ht = h'
    answerVec = zeros(chainLen); res = zeros(chainLen)
    Hb = H*b; total = b
    for i = 1:chainLen
        total += Hb
        answerVec[i] = (ht*total)[1]
        Hb = H*Hb
    end
    # println(answerVec)
    p0 = abs(h)/sum(abs(h))
    P0 = Categorical(p0)
    for i = 1:numTrial
        ind0 = rand(P0)
        W = h[ind0]/p0[ind0]
        X = W*b[ind0]
        for j = 1:chainLen
            k = m-(j-1)%m
            ind1 = indMap[ind0][rand(p[k][ind0])]
            W = W*H[ind0,ind1]/P[k][ind1,ind0]
            X += W*b[ind1]
            res[j] += (X - answerVec[j])^2
            ind0 = ind1
        end
    end
    return res/numTrial
end

function multi_mc(H, h, b, P, p, indMap, chainLen, numTrial, stepLen, fileName)
    n = size(H,1); m = length(P)
    ht = h'
    x = \(speye(n) - H, b)
    ex = (ht*x)[1]
    h = h/ex; ht = h'
    res = []
    p0 = abs(h)/sum(abs(h))
    P0 = Categorical(p0)
    total = 0; count = 1
    for i = 1:numTrial
        ind0 = rand(P0)
        W = h[ind0]/p0[ind0]
        X = W*b[ind0]
        for j = 1:chainLen
            k = m-(j-1)%m
            ind1 = indMap[ind0][rand(p[k][ind0])]
            W = W*H[ind0,ind1]/P[k][ind1,ind0]
            X += W*b[ind1]
            ind0 = ind1
        end
        total += X
        if i == stepLen[count]
            count += 1
            # push!(res, abs(total/i-1))
            f = open(fileName, "a")
            write(f, string(abs(total/i-1)));write(f, '\n')
            close(f)
            println("i = $(i) ------ current error is $(abs(total/i-1))")
        end
    end
    return res
end