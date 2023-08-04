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
    # io
    save,
    load,
    # preprocessing
    scale_force!,
    lowpass_filter!,
    epochs,
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
    peak_force,
    impulse_size,
    duration,
    latency,
    # force statistics
    minimum,
    maximum,
    mean,
    std,
    var,
    median,
    diff

include("data_structs.jl")
include("io.jl")
include("stats.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

end
