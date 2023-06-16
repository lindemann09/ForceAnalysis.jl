module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using Statistics
using DataFrames

import Base.copy
import DataFrames: subset, aggregate

export  FloatOrMissing,
    ForceData, # force data
    MultiForceData,
    ForceProfiles,
    force,
    timestamps,
    sampling_rate,
    n_samples,
    n_profiles,
    copy,
    # preprocessing
    peak_difference,
    lowpass_filter,
    force_profile_matrix,
    force_data_preprocess,
    profile_parameter,
    aggregate,
    subset,
    detect_force_onset,
    #stats missings,
    randFloatOrMissing,
    z_transform,
    eachcol_select_rows,
    eachrow_select_cols,
    column_mean,
    column_minmax,
    column_sum,
    column_var,
    column_std,
    column_diff,
    row_mean,
    row_std,
    row_minmax,
    row_sum,
    row_var


FloatOrMissing = Union{Missing, Float64}

include("stats_missings.jl")
include("data_struct.jl")
include("preprocessing.jl")

end
