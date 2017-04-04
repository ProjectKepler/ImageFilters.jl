include("boxfiltermono.jl")
include("boxfiltertriple.jl")
include("datatypes.jl")

function boxfilter(imagedata::SinglePrecisionImages, filtersize::Tuple{Int64,
    Int64})
  if typeof(imagedata) == Array{Float32, 2}
    filteredimage = boxfiltermono(imagedata, filtersize)
  elseif typeof(imagedata) == Array{Float32, 3}
    filteredimage = boxfiltertriple(imagedata, filtersize)
  end
  return filteredimage
end

function boxfilter!(filteredimage::SinglePrecisionImages,
    buffer::SinglePrecisionImages, imagedata::SinglePrecisionImages,
    filtersize::Tuple{Int64, Int64})
  if typeof(imagedata) == Array{Float32, 2}
    boxfiltermono!(filteredimage, buffer, imagedata, filtersize)
  elseif typeof(imagedata) == Array{Float32, 3}
    boxfiltertriple!(filteredimage, buffer, imagedata, filtersize)
  end
end
