module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using Statistics
using DataFrames

import Base.copy
import DataFrames: subset, aggregate

export FloatOrMissing,
    ForceData, # force data
    MultiForceData,
    ForceProfiles,
    force,
    baseline,
    timestamps,
    sampling_rate,
    n_samples,
    n_profiles,
    copy,
    duration,
    # preprocessing
    lowpass_filter,
    force_profile_matrix,
    force_data_preprocess,
    # processing
    peak_difference,
    profile_parameter,
    aggregate,
    subset,
    # response detection
    OnsetCriterion,
    ResponseBounds,
    detect_responses,
    response_onset,
    response_offset,
    peak_force,
    impulse_size,
    latency,
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

FloatOrMissing = Union{Missing,Float64}

include("stats_missings.jl")
include("data_struct.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

end
