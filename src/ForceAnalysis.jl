module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using DataFrames
using JLD2, CodecZlib

import Base: copy, minimum, maximum
import DataFrames: subset, aggregate
import CategoricalArrays: unique
import Statistics: mean, median, var, std

export ForceData, # force data
    MultiForceData,
    ForceProfiles,
    load_force_profiles,
    load_force_data,
    save_force,
    force,
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
    ForceResponse,
    response_detection,
    response_onset,
    response_offset,
    peak_force,
    impulse_size,
    latency,
    extract_response,
    minimum,
    maximum,
    mean,
    std,
    var,
    median,
    z_transform

include("data_struct.jl")
include("stats.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

end
