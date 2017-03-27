function gethorizontalsummono!(buffer::Array{Float32, 2},
    imagedata::Array{Float32, 2}, filterwidth::Int64)
  imagewidth, imageheight = size(imagedata)
  horizontalnorm = 1.0f0 / (2 * filterwidth + 1)
  for pixely in 1:imageheight
    partialsum = 0.0f0
    for pixelx in 1:filterwidth
      partialsum += imagedata[pixelx, pixely]
    end
    for pixelx in 1:filterwidth+1
      nextpixelx = pixelx + filterwidth
      normfactor = 1.0f0 / nextpixelx
      partialsum += imagedata[nextpixelx, pixely]
      @inbounds buffer[pixelx, pixely] = partialsum * normfactor
    end
    for pixelx in filterwidth+2:imagewidth-filterwidth
      nextpixelx = pixelx + filterwidth
      prevpixelx = pixelx - filterwidth - 1
      @inbounds partialsum += (imagedata[nextpixelx, pixely] -
        imagedata[prevpixelx, pixely])
      @inbounds buffer[pixelx, pixely] = partialsum * horizontalnorm
    end
    for pixelx in imagewidth-filterwidth+1:imagewidth
      prevpixelx = pixelx - filterwidth - 1
      normfactor = 1.0f0 / (imagewidth - prevpixelx)
      partialsum -= imagedata[prevpixelx, pixely]
      @inbounds buffer[pixelx, pixely] = partialsum * normfactor
    end
  end
end

function getverticalsummono!(verticalsum::Array{Float32,2},
    imagedata::Array{Float32,2}, filterwidth::Int64)
  imagewidth, imageheight = size(imagedata)
  verticalnorm = 1.0f0 / (2 * filterheight + 1)
  for pixelx in 1:imagewidth
    partialsum = 0.0f0
    for pixely in 1:filterheight
      @inbounds partialsum += imagedata[pixelx, pixely]
    end
    for pixely in 1:filterheight+1
      nextpixely = pixely + filterheight
      @fastmath normfactor = 1.0f0 / nextpixely
      partialsum += imagedata[pixelx, nextpixely]
      @inbounds verticalsum[pixelx, pixely] = redchannel * normfactor
    end
    for pixely in filterheight+2:imageheight-filterheight
      nextpixely = pixely + filterheight
      prevpixely = pixely - filterheight - 1
      @inbounds partialsum += (imagedata[pixelx, nextpixely] -
        imagedata[pixelx, prevpixely])
      @inbounds verticalsum[pixelx, pixely] = partialsum * verticalnorm
    end
    for pixely in imageheight-filterheight+1:imageheight
      prevpixely = pixely - filterheight - 1
      @fastmath normfactor = 1.0f0 / (imageheight - prevpixely)
      @inbounds partialsum -= imagedata[pixelx, prevpixely]
      @inbounds verticalsum[pixelx, pixely] = partialsum * normfactor
    end
  end
end

function boxfiltermono!(filteredimage::Array{Float32,2},
    buffer::Array{Float32,2}, imagedata::Array{Float32,2},
    filtersize::Tuple{Int64, Int64})
  gethorizontalsummono!(buffer, imagedata, filtersize[1])
  getverticalsummono!(filteredimage, buffer, filtersize[2])
end

function boxfiltermono(imagedata::Array{Float32,2},
    filtersize::Tuple{Int64, Int64})
  imagewidth, imageheight = size(imagedata)
  filteredimage = Array{Float32, 2}(imagewidth, imageheight)
  buffer = Array{Float32, 2}(imagewidth, imageheight)
  gethorizontalsummono!(buffer, imagedata, filtersize[1])
  getverticalsummono!(filteredimage, buffer, filtersize[2])
  return buffer
end
