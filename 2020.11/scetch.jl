




Base.copy(c::Chain) = Chain(c.layers)


# (c::Chain)(x) =
#     ifelse(
#         ((c.bottleneck == nothing) |
#         (c.bottleneck = c.body(x)))
#         ,
#         [h(c.body(x)) for h in c.heads],
#         [h(c.bottleneck)) for h in c.heads]
#     )

"""
NEED TO SOLVE BOTTLENECK HEAD CONNECTION
"""
#
# (c::Chain)(x) =
#     (
#         for h in c.heads
#             if (c.bottleneck == nothing)
#                 c.bottleneck = c.body(x)
#             end
#             h(c.bottleneck)
#         end
#     )




# (h::Head)(x, c::Chain) = h(c.body(x))

# (c::Chain)(x, y, lossF, h) = lossF(c.heads[h](x, c), y)

# (c::Chain)(x, y) = [
#     (multinomial_nll(c.heads[h1](x, c), y),
#     mse(c.heads[h2](x, c), y))
#     for (h1, h2) in zip(1:2:c.n, 2:2:c.n)
#
# ]



############################

############################################
#
# x = rand(Float32, (1000, 4,1,1))
# outputs = bpnet(x)
#
# ytrue = rand(Float32, (1000, 1, 2, 1));
# ypred = outputs[1]

# ypred = log.(ypred)

# probabilites = exp.(scores) ./ sum(exp.(scores),dims=1)


function mytrain!(
    c::Chain,
    trainData::Data, testData::Data,
    period::Int=10, iters::Int=500
    )
    # Your code here
    trainLoss, testLoss = [], []

    for i=1:period:iters
        push!(trainLoss, @diff c(trainData))
        push!(testLoss, @diff c(testData))
        for w in params(c)
            g = grad(loss,w)



        for (h1, h2) in zip(1:2:c.n, 2:2:c.n)
            pLoss = adam(c(..., multinomial_nll, h1), take(cycle(trainData),period))
            cLoss = adam(c(..., mse, h2), take(cycle(trainData),period))
        end
    end
    push!(trainLoss, model(trainData))
    push!(testLoss, model(testData))
    return 0:period:iters, trainLoss, testLoss
end

########################################

using Knet: minibatch

batchsize = 20
shuffle = false


for y in [y1, y2, y1, y2]
    push!(data, minibatch(x, y, batchsize; shuffle))
end

for d in data[1]
    println(typeof(d))
    bpnet(d)
end

bpnet(x, [y1,y2,y1,y2])


#############################

function eltype(d::Data)
    Data
end

function HasEltype(d::Data)
    Data
end



(c::Chain)(d::Data) = (args)
    body
end


]




for data in DataBatch
    adam(multinomial_nll, data)



(factorial.(sum(ytrue, dims=1)) ./ prod(factorial.(ytrue), dims=1)) * prod(pV .^ ytrue)

total = sum(ypred, dims=1)


yp = reshape(ypred, (699,2))[:,1]

probabilites = exp.(ypred) ./ sum(exp.(ypred),dims=1)

-(log(probabilites[labels[1],1]) + log(probabilites[labels[2],2]))/2



Multinomial(699, yp)
Distributions._logpdf()

h(bpnet.body(x))







[
    h.lossF(x,y1);cH.lossF(x,y2)
    for (pH, cH) in zip(c.profileHeads, y)
]


[
    i;j
    for (i,j) in zip(1:3,3:5)
]














f = ceil.(abs.(rand(Float16, 699, 1) * 10))
r = ceil.(abs.(rand(Float16, 699, 1) * 10))

nF = sum(f)
nR = sum(r)


bpnet
update!()


for batch in minibatches
    currentLoss = bpnet(batch)



















bottleneck = copy(bpnet)



for task in tasks




###########################

x = cat(x,x, dims=3)

o = bpnet(x)

d = rand(25,1,2,64)
dC = deconv4(d, o, padding=(12,0))


pool

pool(o; window=(699,1))


size(x)

w's shape is (W1,W2,...,I,O)
    - Where Ws are the kernel dimentions. (Spatial)
    - I is the input channels.
    - O is the output channels.
w = rand(25,4,1,64)

x's shape is (X1,X2,...,I,N)
    - Where Xs are input of one example's  shape. (Spatial)
    - I is the input channels.
    - N is the number of instances
x = rand(1000,4,1,100)

y's shape is (Y1,Y2,...,O,N)
    - Where Ys are output dimentions. (Spatial)
    - O is the output channels.
    - N is the number of instances
conv4(w, x, padding=(12,0), stride=(1,1))


size(w)

size(o)











a = Conv(25,25,1,64, padding=12)

b = conv4(25,1,1,64, padding=(12,0),  stride=(1,0))

w = rand(25,4, 1,64)
conv4(w,x,padding=(12,0),stride=(1,0))


o = a(x)

f = rand(Float32,25,1,1,64)

deconv4(f,o)

b = Conv(3,3,64,64, dilation=2^i)



r = reshape(o, 4,699,1,64,1)

l = rand(Float32, 1,25,1,2)
deconv4(l,o)



(Int(max((21 - 1) * 1 + 5 - 21, 0)/2),
 Int(max((51 - 1) * 1 + 5 - 51, 0)/2))







x2 = rand(4,101,1,1)


w = rand(4,25,1,64)

conv4(w, dropout(x,0), padding=(0,12), stride=1)





P = ((0-1)*W-S+F)/2











(c::Conv)(x) = c.f.(conv4(c.w, dropout(x,c.pDrop), padding=) .+ c.b))

Conv(w1::Int,w2::Int,cx::Int,cy::Int,f=relu;pdrop=0) = Conv(param(w1,w2,cx,cy), param0(1,1,cy,1), f, pdrop)









Yi=1+floor((Xi+2*padding[i]-Wi)/stride[i])

paddinf[i] = (Yi-1) * stride[i]

cy-1

1 + floor((25+2*12-25)/1)

F = 25
S = 0


((W - F * 2P)/S) + 1

((101 - 25 *2))



cy=699; s1=1; w1=25; cx=699
pad1 = Int(max((cy - 1) * s1 + w1 - cx, 0)/2)





######################################


function Base.iterate(d::iterN, state=(0)) # here the start point of state is 0
    if state == 0 && d.shuffle
        d.idx = randperm(d.n)
    end

    if state >= d.n
        return nothing
    end

    i = state + 1 # here the beginning of the current slice
    j = i + min(d.batchsize, length(d.idx)) - 1 # here the end of the current slice

    xbatch = d.x[:,i:j]
    ybatch = d.y[i:j]

    return ((xbatch, ybatch), j) # here it returns the batch and the next state of the iteration
end
