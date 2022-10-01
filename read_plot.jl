using PyPlot
using LinearAlgebra
include("quat.jl")


function readData(dataFilePath::String, n::Int, m::Int)

    data = zeros(n, m)

    file = open(dataFilePath, "r")

    for i=1:n
        line = readline(file)
        numberString = split(line, ',');
        for j = 1:m
            data[i,j] = parse(Float64, numberString[j])
        end
    end

    return data

end

function readData(dataFilePath::String, m::Int)

    data = Array{Float64}(undef, 0, m)
    rowInsert = zeros(1,m)

    file = open(dataFilePath, "r")
    
    while !eof(file)
        line = readline(file)
        numberString = split(line, ',')
        for j = 1:m
            rowInsert[1,j] = parse(Float64, numberString[j])
        end
        data = vcat(data,rowInsert)
    end

    return data

end

function plotx3t(flightHistory::Matrix{Float64})

    pygui(true)
    return plot(flightHistory[:,1], flightHistory[:,4])

end

function get3DQuiverPlot_py(flightHistory::Matrix{Float64})
    #z: input matrix of state vectors
    #gives quiverplot with scatter overlay
    
    z = flightHistory[:,2:end]
    #max array of rocket pointing vectors
    b3_I = zeros(size(z)[1], 3)
    for i = 1:size(z)[1]
        b3_I[i,:] = rotateFrame([0,0,1.0], quatInv(z[i, 7:10]))
    end

    #shorten array of b3_I vectors (not as many)
    numToPlot = 30
    array = floor.(Int, collect(LinRange(1, size(z)[1], numToPlot)))
    x = zeros(size(array))
    y = zeros(size(array))
    k = zeros(size(array))
    shortb3_I = zeros(size(array)[1], 3)

    shortIndex = 1
    for index in array
        x[shortIndex] = z[index, 1]
        y[shortIndex] = z[index, 2]
        k[shortIndex] = z[index, 3]
        shortb3_I[shortIndex, :] = b3_I[index, :]
        shortIndex = shortIndex + 1
    end


    #plot  all positions
    plot3D(z[:,1], z[:,2], z[:,3])
    #overlay dedensified pointing arrows
    quiver(x, y, k, shortb3_I[:,1], shortb3_I[:,2], shortb3_I[:,3], length = 10, colors = "red")
    autoscale()
    xlabel("e1 (m)")
    ylabel("e2 (m)")
    zlabel("e3 (m)")

end


let 

    n = 1000
    m = 14
    data = readData("export_data.txt", m)

    plotx3t(data)
    
end