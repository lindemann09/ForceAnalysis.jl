
FloatOrMissing = Union{Missing, Float64}

## helper functions
function randFloatOrMissing(size; percent_missings::AbstractFloat=0.1)
    rtn = convert(Matrix{FloatOrMissing},
                rand(Float64, size))
    n_missings = ceil(Int64, length(rtn)*percent_missings)
    ids = randperm(length(rtn))[1:n_missings]
    rtn[ids] .= missing
    return rtn
end;

function z_transform(vector::AbstractVector{T};
            corrected::Bool=true) where T<:FloatOrMissing
    # handles missings
    v = skipmissing(vector)
    return((vector .- mean(v))./std(v; corrected))
end

function eachcol_select_rows(x::Matrix{T}; rows=Int64[])  where T<:FloatOrMissing
    ## eachrow
    if length(rows)>0
        return eachcol(x[rows, :])
    else
        return eachcol(x)
    end;
end;

function eachrow_select_cols(x::Matrix{T}; cols=Int64[])  where T<:FloatOrMissing
    ## eachrow
    if length(cols)>0
        return eachrow(x[:, cols])
    else
        return eachrow(x)
    end;
end;


## column and row wise statistics of matrices (Float64) with missings


function column_mean(x::Matrix{T}; rows=Int64[])  where T<:FloatOrMissing
    return map(x-> mean(skipmissing(x)), eachcol_select_rows(x; rows))
end

function column_median(x::Matrix{T}; rows=Int64[])  where T<:FloatOrMissing
    return map(x-> median(skipmissing(x)), eachcol_select_rows(x; rows))
end


function column_var(x::Matrix{T};
                    corrected::Bool=true, rows=Int64[]) where T<:FloatOrMissing
    return map(x-> var(skipmissing(x);corrected), eachcol_select_rows(x; rows))
end


function column_std(x::Matrix{T};
                    corrected::Bool=true, rows=Int64[]) where T<:FloatOrMissing
    return map(x-> std(skipmissing(x);corrected), eachcol_select_rows(x; rows))
end


function column_sum(x::Matrix{T}; rows=Int64[])  where T<:FloatOrMissing
    return map(x-> sum(skipmissing(x)), eachcol_select_rows(x; rows))
end


function column_minmax(x::Matrix{T}; rows=Int64[])  where T<:FloatOrMissing
    minx = map(x-> minimum(skipmissing(x)), eachcol_select_rows(x; rows))
    maxx = map(x-> maximum(skipmissing(x)), eachcol_select_rows(x; rows))
    return (minx, maxx)
end

# rows wise

function row_mean(x::Matrix{T}; cols=Int64[])  where T<:FloatOrMissing
    return map(x-> mean(skipmissing(x)), eachrow_select_cols(x; cols))
end

function row_median(x::Matrix{T}; cols=Int64[])  where T<:FloatOrMissing
    return map(x-> median(skipmissing(x)), eachrow_select_cols(x; cols))
end

function row_var(x::Matrix{T};
                corrected::Bool=true, cols=Int64[])  where T<:FloatOrMissing
    return map(x-> var(skipmissing(x);corrected), eachrow_select_cols(x; cols))
end


function row_std(x::Matrix{T};
                corrected::Bool=true, cols=Int64[])  where T<:FloatOrMissing
    return map(x-> std(skipmissing(x);corrected), eachrow_select_cols(x; cols))
end


function row_sum(x::Matrix{T}; cols=Int64[])  where T<:FloatOrMissing
    return map(x-> sum(skipmissing(x)), eachrow_select_cols(x; cols))
end

function row_minmax(x::Matrix{T}; cols=Int64[])  where T<:FloatOrMissing
    minx = map(x-> minimum(skipmissing(x)), eachrow_select_cols(x; cols))
    maxx = map(x-> maximum(skipmissing(x)), eachrow_select_cols(x; cols))
    return (minx, maxx)
end