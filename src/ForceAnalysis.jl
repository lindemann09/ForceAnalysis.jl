module ForceAnalysis

using Random
using Statistics
using DataFrames
using DSP # signal processing, filtering


import Base.copy
import DataAPI.ncol
import DataFrames: subset, aggregate

export ForceExpData,
        ForceProfiles,
        copy,
        ncol,
        # preprocessing
        peak_difference,
        baselines,
        lowpass_filter,
        force_profile_matrix,
        force_data_preprocess,
        profile_parameter,
        aggregate,
        subset,
        #stats missings
        FloatOrMissing,
        randFloatOrMissing,
        z_transform,
        eachcol_select_rows,
        eachrow_select_cols,
        column_mean,
        column_minmax,
        column_sum,
        column_var,
        column_std,
        row_mean,
        row_std,
        row_minmax,
        row_sum,
        row_var


include("stats_missings.jl")
include("data_struct.jl")
include("preprocessing.jl")


end
