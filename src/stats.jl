
function mean(fp::ForceProfiles;
        rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = mean(fp.dat, dims=1)
        bsl = mean(fp.baseline)
    else
        dat = mean(fp.dat[rows, :], dims=1)
        bsl = mean(fp.baseline[rows, :])
    end
    @info typeof(bsl), typeof(dat)
    return ForceProfiles(dat, fp.sr, DataFrame(), [bsl], fp.zero_sample)
end

function median(fp::ForceProfiles;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = median(fp.dat, dims=1)
        bsl = median(fp.baseline)
    else
        dat = median(fp.dat[rows, :], dims=1)
        bsl = median(fp.baseline[rows, :])
    end
    return ForceProfiles(dat, fp.sr, DataFrame(), [bsl], fp.zero_sample)
end

function var(fp::ForceProfiles;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = var(fp.dat, dims=1)
        bsl = var(fp.baseline)
    else
        dat = var(fp.dat[rows, :], dims=1)
        bsl = var(fp.baseline[rows, :])
    end
    return ForceProfiles(dat, fp.sr, DataFrame(), [bsl], fp.zero_sample)
end

function std(fp::ForceProfiles;
    rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = std(fp.dat, dims=1)
        bsl = std(fp.baseline)
    else
        dat = std(fp.dat[rows, :], dims=1)
        bsl = std(fp.baseline[rows, :])
    end
    return ForceProfiles(dat, fp.sr, DataFrame(), [bsl], fp.zero_sample)
end

function z_transform(fp::ForceProfiles; corrected::Bool=true)
    ## FIXME not tested
    m = mean(fp.dat, dims=1)
    sd = std(fp.dat, dims=1, corrected=corrected)
    dat = (fp.dat .- m)./sd
    bsl = (fp.baseline .- m)./sd
    return ForceProfiles(dat, fp.sr, DataFrame(), bsl, fp.zero_sample)
end

"""
    minimum(fp:ForceProfiles)

Minimum of each profile.
"""
minimum(fp::ForceProfiles) = return vec(minimum(fp.dat, dims=2))

"""
    maximum(fp:ForceProfiles)

Minimum of each profile.
"""
maximum(fp::ForceProfiles) = return vec(maximum(fp.dat, dims=2))
