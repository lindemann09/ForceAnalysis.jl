
"""
    mean(fe::BeForEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)

Computes the mean of the epochs, i.e., column-wise means for each sample.
If `rows` is defined, only the selected rows are considered.
"""
function Statistics.mean(fe::BeForEpochs;
        rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = mean(fe.dat, dims=1)
        bsl = mean(fe.baseline)
    else
        dat = mean(fe.dat[rows, :], dims=1)
        bsl = mean(fe.baseline[rows, :])
    end

    return BeForEpochs(dat, fe.sampling_rate;
                        baseline = [bsl], zero_sample = fe.zero_sample)
end

"""
    median(fe::BeForEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)

Computes the median of the epochs, i.e., column-wise medians for each sample.
If `rows` is defined, only the selected rows are considered.
"""
function Statistics.median(fe::BeForEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = median(fe.dat, dims=1)
        bsl = median(fe.baseline)
    else
        if length(rows) != 0
            dat = median(fe.dat[rows, :], dims=1)
            bsl = median(fe.baseline[rows, :])
        else
            dat = Matrix{eltype(fe.dat)}(undef, 1, size(fe.dat,2))
            fill!(dat, NaN)
            bsl = NaN
        end
    end
    return BeForEpochs(dat, fe.sampling_rate;
                        baseline = [bsl], zero_sample = fe.zero_sample)

end

"""
    var(fe::BeForEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)

Computes the variance of the epochs, i.e., column-wise variances for each sample.
If `rows` is defined, only the selected rows are considered.
"""
function Statistics.var(fe::BeForEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = var(fe.dat, dims=1)
        bsl = var(fe.baseline)
    else
        dat = var(fe.dat[rows, :], dims=1)
        bsl = var(fe.baseline[rows, :])
    end
    return BeForEpochs(dat, fe.sampling_rate;
                        baseline = [bsl], zero_sample = fe.zero_sample)

end

"""
    std(fe::BeForEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)

Computes the standard deviation of the epochs, i.e., column-wise standard
deviations for each sample. If `rows` is defined, only the selected rows are considered.
"""
function Statistics.std(fe::BeForEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = std(fe.dat, dims=1)
        bsl = std(fe.baseline)
    else
        dat = std(fe.dat[rows, :], dims=1)
        bsl = std(fe.baseline[rows, :])
    end
    return BeForEpochs(dat, fe.sampling_rate;
                        baseline = [bsl], zero_sample = fe.zero_sample)
end


"""
    sderr(fe::BeForEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)

Computes the standard error of the epochs, i.e., column-wise standard
deviations for each sample. If `rows` is defined, only the selected rows are considered.
"""
function sderr(fe::BeForEpochs;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        n = size(fe.dat, 1)
        dat = std(fe.dat, dims=1)
        bsl = std(fe.baseline)
    else
        n = size(fe.dat[rows, :], 1)
        dat = std(fe.dat[rows, :], dims=1)
        bsl = std(fe.baseline[rows, :])
    end
    return BeForEpochs(dat/sqrt(n), fe.sampling_rate;
                        baseline = [bsl], zero_sample = fe.zero_sample)
end


function Base.diff(fe::BeForEpochs; dims::Integer)
    mtx = fe.dat
    if dims == 1
        z = zeros(Float64, 1, size(mtx, 2))
        dat = vcat(z, diff(mtx; dims=1))
    elseif dims == 2
        z = zeros(Float64, size(mtx, 1), 1)
        dat = hcat(z, diff(mtx; dims=2))
    else
        throw(ArgumentError("dims has to be 1 or 2 and not $dims"))
    end
    return BeForEpochs(dat, fe.sampling_rate, fe.design, fe.baseline, fe.zero_sample)
end;


"""
    minimum(fe:BeForEpochs)

Minimum of each epoch.
"""
Base.minimum(fe::BeForEpochs) = return vec(minimum(fe.dat, dims=2))

"""
    maximum(fe:BeForEpochs)

Minimum of each epoch.
"""
Base.maximum(fe::BeForEpochs) = return vec(maximum(fe.dat, dims=2))
