
function mean(fp::ForceProfiles;
        rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)
    if isnothing(rows)
        dat = mean(fp.dat, dims=1)
        bsl = mean(fp.baseline)
    else
        dat = mean(fp.dat[rows, :], dims=1)
        bsl = mean(fp.baseline[rows, :])
    end
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


function diff(fp::ForceProfiles{T}; dims::Integer) where {T}
    mtx = fp.dat
    if dims == 1
        z = zeros(T, 1, size(mtx, 2))
        dat = vcat(z, diff(mtx; dims=1))
    elseif dims == 2
        z = zeros(T, size(mtx, 1), 1)
        dat = hcat(z, diff(mtx; dims=2))
    else
        throw(ArgumentError("dims has to be 1 or 2 and not $dims"))
    end
    return ForceProfiles(dat, fp.sr, fp.design, fp.baseline, fp.zero_sample)
end;


# function z_transform(fp::ForceProfiles; corrected::Bool=true)
#     m = mean(fp.dat, dims=2)
#     sd = std(fp.dat, dims=2, corrected=corrected)
#     dat = (fp.dat .- m)./sd
#     bsl = (fp.baseline .- vec(m))./vec(sd)
#     return ForceProfiles(dat, fp.sr, fp.design, bsl, fp.zero_sample)
# end

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
