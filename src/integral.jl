function getgrayintegraldata!(imagedata::Array{Float64,2},
    integraldata::Array{Float64,2})
  imagewidth, imageheight = size(imagedata)
  result = 0.0
  @simd for pixelx = 1:imagewidth
    @inbounds result += imagedata[pixelx, 1]
    @inbounds integraldata[pixelx, 1] = result
  end

  @simd for pixely = 2:imageheight
    result = 0.0
    @simd for pixelx = 1:imagewidth
      @inbounds result += imagedata[pixelx, pixely]
      @inbounds integraldata[pixelx, pixely] = result + integraldata[pixelx,
        pixely-1]
    end
  end
  return integraldata
end
