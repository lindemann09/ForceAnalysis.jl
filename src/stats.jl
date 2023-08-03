
function mean(fe::ForceEpochs;
        rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = mean(fe.dat, dims=1)
        bsl = mean(fe.baseline)
    else
        dat = mean(fe.dat[rows, :], dims=1)
        bsl = mean(fe.baseline[rows, :])
    end
    return ForceEpochs(dat, fe.sr, DataFrame(), [bsl], fe.zero_sample)
end

function median(fe::ForceEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = median(fe.dat, dims=1)
        bsl = median(fe.baseline)
    else
        dat = median(fe.dat[rows, :], dims=1)
        bsl = median(fe.baseline[rows, :])
    end
    return ForceEpochs(dat, fe.sr, DataFrame(), [bsl], fe.zero_sample)
end

function var(fe::ForceEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = var(fe.dat, dims=1)
        bsl = var(fe.baseline)
    else
        dat = var(fe.dat[rows, :], dims=1)
        bsl = var(fe.baseline[rows, :])
    end
    return ForceEpochs(dat, fe.sr, DataFrame(), [bsl], fe.zero_sample)
end

function std(fe::ForceEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = std(fe.dat, dims=1)
        bsl = std(fe.baseline)
    else
        dat = std(fe.dat[rows, :], dims=1)
        bsl = std(fe.baseline[rows, :])
    end
    return ForceEpochs(dat, fe.sr, DataFrame(), [bsl], fe.zero_sample)
end


function diff(fe::ForceEpochs{T}; dims::Integer) where {T}
    mtx = fe.dat
    if dims == 1
        z = zeros(T, 1, size(mtx, 2))
        dat = vcat(z, diff(mtx; dims=1))
    elseif dims == 2
        z = zeros(T, size(mtx, 1), 1)
        dat = hcat(z, diff(mtx; dims=2))
    else
        throw(ArgumentError("dims has to be 1 or 2 and not $dims"))
    end
    return ForceEpochs(dat, fe.sr, fe.design, fe.baseline, fe.zero_sample)
end;


"""
    minimum(fe:ForceEpochs)

Minimum of each epoch.
"""
minimum(fe::ForceEpochs) = return vec(minimum(fe.dat, dims=2))

"""
    maximum(fe:ForceEpochs)

Minimum of each epoch.
"""
maximum(fe::ForceEpochs) = return vec(maximum(fe.dat, dims=2))
