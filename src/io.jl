fl_formats = Dict(:jld2 => ".jld2", :csv => ".csv.zip")

function save(filename::String, forces::ForceData)
	return jldsave(_check_suffix(filename, :jdl2), true; forces)
end

function load(::Type{ForceData}, filename::String)
	return convert(ForceData, load_object(filename))
end;

## Force profiles

function save(filename::String, format::Symbol, fp::ForceProfiles; compress = true)
	filename = _check_suffix(filename, format)
	if format == :jld2
		jldsave(filename, compress; fp)

	elseif format == :csv
		dataname = _dataname(filename)
		ZipWriter(filename) do zipfl
			zip_newfile(zipfl, dataname * ".forces.csv"; compress)
			CSV.write(zipfl, CSV.Tables.table(fp.dat); writeheader = false)

			zip_newfile(zipfl, dataname * ".design.csv"; compress)
			CSV.write(zipfl, fp.design)

			zip_newfile(zipfl, dataname * ".baseline.csv"; compress)
			CSV.write(zipfl, CSV.Tables.table(fp.baseline); writeheader = false)

			zip_newfile(zipfl, dataname * ".meta.json"; compress)
			JSON.print(zipfl, Dict(:sampling_rate => fp.sampling_rate,
					:zero_sample => fp.zero_sample), 1)
		end
	end
	return nothing
end

# jdl2 as default save method
function save(filename::String, fp::ForceProfiles; compress = true)
	save(filename, :jdl2, fp; compress)
end

function load(::Type{ForceProfiles}, format::Symbol, filename::String)
	if format == :jld2
		return convert(ForceProfiles, load_object(filename))
	elseif format == :csv
		dataname = _dataname(filename)

		zipfl = zip_open_filereader(filename)
        txt = zip_readentry(zipfl, dataname * ".meta.json", String)
        meta = JSON.parse(txt)
        fl = zip_openentry(zipfl, dataname * ".forces.csv")
        dat = CSV.read(fl, CSV.Tables.matrix; header = false)
        fl = zip_openentry(zipfl, dataname * ".baseline.csv")
        baseline = CSV.read(fl, CSV.Tables.matrix; header = false)
        fl = zip_openentry(zipfl, dataname * ".design.csv")
        design = CSV.read(fl, DataFrame)
        close(zipfl)

        return ForceProfiles(dat, meta["sampling_rate"], design,
			vec(baseline), meta["zero_sample"])
	else
		_format_error(format)
	end
end;

# guess format, if not defined
function load(::Type{ForceProfiles}, filename::String)
	format = :unknown
	for (k, v) in fl_formats
		if endswith(filename, v)
			format = k
		end
	end
	if format == :unknown
		throw(DomainError(filename, "Please specify file format. Can't guess it from filename."))
	end
	return load(ForceProfiles, format, filename)
end


## helper
_format_error(format) = throw(DomainError(format,
	"Unkown fileformat. Possible formats are $(keys(fl_formats))"))

_dataname(flname) = first(split(last(splitdir(flname)), "."))

function _check_suffix(flname::AbstractString, format::Symbol)
	# ensure correct suffix based on format
	if haskey(fl_formats, format)
		suffix = fl_formats[format]
		return endswith(flname, suffix) ? flname : flname * suffix
	else
		_format_error(format)
	end
end
