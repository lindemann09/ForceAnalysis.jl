module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using DataFrames
using JLD2, CodecZlib
using CSV, JSON
using ZipArchives

import Base: minimum, maximum, diff
import DataFrames: subset, aggregate
import CategoricalArrays: unique
import Statistics: mean, median, var, std
import FileIO: save, load

export ForceData, # force data
    MultiForceData,
    ForceEpochs,
    force,
    duration,
    scale_force!,
    save,
    load,
    # preprocessing
    lowpass_filter!,
    extract_force_epochs,
    adjust_baseline!,
    # processing
    peak_difference,
    epoch_parameter,
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
include("io.jl")
include("stats.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

end
