module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using DataFrames
using JLD2, CodecZlib

import Base: minimum, maximum, diff
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
    duration,
    scale_force!,
    # preprocessing
    lowpass_filter!,
    extract_force_profiles,
    adjust_baseline!,
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
    diff

include("data_struct.jl")
include("stats.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

end
