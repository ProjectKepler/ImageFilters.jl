function gethorizontalsumtriple!(horizontalsum::Array{Float32,3},
    imagedata::Array{Float32,3}, filterwidth::Int64)
  channels, imagewidth, imageheight = size(imagedata)
  horizontalnorm = 1.0f0 / (2 * filterwidth + 1)
  for pixely in 1:imageheight
    redchannel = 0.0f0
    greenchannel = 0.0f0
    bluechannel = 0.0f0
    for pixelx in 1:filterwidth
      @inbounds redchannel += imagedata[1, pixelx, pixely]
      @inbounds greenchannel += imagedata[2, pixelx, pixely]
      @inbounds bluechannel += imagedata[3, pixelx, pixely]
    end
    for pixelx in 1:filterwidth+1
      nextpixelx = pixelx + filterwidth
      normfactor = 1.0f0 / nextpixelx
      @inbounds redchannel += imagedata[1, nextpixelx, pixely]
      @inbounds greenchannel += imagedata[2, nextpixelx, pixely]
      @inbounds bluechannel += imagedata[3, nextpixelx, pixely]
      @inbounds horizontalsum[1, pixelx, pixely] = redchannel * normfactor
      @inbounds horizontalsum[2, pixelx, pixely] = greenchannel * normfactor
      @inbounds horizontalsum[3, pixelx, pixely] = bluechannel * normfactor
    end
    for pixelx in filterwidth+2:imagewidth-filterwidth
      nextpixelx = pixelx + filterwidth
      prevpixelx = pixelx - filterwidth - 1
      @inbounds redchannel += (imagedata[1, nextpixelx, pixely] -
        imagedata[1, prevpixelx, pixely])
      @inbounds greenchannel += (imagedata[2, nextpixelx, pixely] -
        imagedata[2, prevpixelx, pixely])
      @inbounds bluechannel += (imagedata[3, nextpixelx, pixely] -
        imagedata[3, prevpixelx, pixely])
      @inbounds horizontalsum[1, pixelx, pixely] = redchannel * horizontalnorm
      @inbounds horizontalsum[2, pixelx, pixely] = greenchannel * horizontalnorm
      @inbounds horizontalsum[3, pixelx, pixely] = bluechannel * horizontalnorm
    end
    for pixelx in imagewidth-filterwidth+1:imagewidth
      prevpixelx = pixelx - filterwidth - 1
      normfactor = 1.0f0 / (imagewidth - prevpixelx)
      @inbounds redchannel -= imagedata[1, prevpixelx, pixely]
      @inbounds greenchannel -= imagedata[2, prevpixelx, pixely]
      @inbounds bluechannel -= imagedata[3, prevpixelx, pixely]
      @inbounds [1, pixelx, pixely] = redchannel * normfactor
      @inbounds horizontalsum[2, pixelx, pixely] = greenchannel * normfactor
      @inbounds horizontalsum[3, pixelx, pixely] = bluechannel * normfactor
    end
  end
end

function getverticalsumtriple!(verticalsum::Array{Float32,3},
    imagedata::Array{Float32,3}, filterheight::Int64)
  channels, imagewidth, imageheight = size(imagedata)
  verticalnorm = 1.0f0 / (2 * filterheight + 1)
  for pixelx in 1:imagewidth
    redchannel = 0.0f0
    greenchannel = 0.0f0
    bluechannel = 0.0f0
    for pixely in 1:filterheight
      @inbounds redchannel += imagedata[1, pixelx, pixely]
      @inbounds greenchannel += imagedata[2, pixelx, pixely]
      @inbounds bluechannel += imagedata[3, pixelx, pixely]
    end
    for pixely in 1:filterheight+1
      nextpixely = pixely + filterheight
      normfactor = 1.0f0 / nextpixely
      @inbounds redchannel += imagedata[1, pixelx, nextpixely]
      @inbounds greenchannel += imagedata[2, pixelx, nextpixely]
      @inbounds bluechannel += imagedata[3,pixelx, nextpixely]
      @inbounds verticalsum[1, pixelx, pixely] = redchannel * normfactor
      @inbounds verticalsum[2, pixelx, pixely] = greenchannel * normfactor
      @inbounds verticalsum[3, pixelx, pixely] = bluechannel * normfactor
    end
    for pixely in filterheight+2:imageheight-filterheight
      nextpixely = pixely + filterheight
      prevpixely = pixely - filterheight - 1
      @inbounds redchannel += (imagedata[1, pixelx, nextpixely] -
        imagedata[1, pixelx, prevpixely])
      @inbounds greenchannel += (imagedata[2, pixelx, nextpixely] -
        imagedata[2, pixelx, prevpixely])
      @inbounds bluechannel += (imagedata[3, pixelx, nextpixely] -
        imagedata[3, pixelx, prevpixely])
      @inbounds verticalsum[1, pixelx, pixely] = redchannel * verticalnorm
      @inbounds verticalsum[2, pixelx, pixely] = greenchannel * verticalnorm
      @inbounds verticalsum[3, pixelx, pixely] = bluechannel * verticalnorm
    end
    for pixely in imageheight-filterheight+1:imageheight
      prevpixely = pixely - filterheight - 1
      normfactor = 1.0f0 / (imageheight - prevpixely)
      @inbounds redchannel -= imagedata[1, pixelx, prevpixely]
      @inbounds greenchannel -= imagedata[2, pixelx, prevpixely]
      @inbounds bluechannel -= imagedata[3, pixelx, prevpixely]
      @inbounds verticalsum[1, pixelx, pixely] = redchannel * normfactor
      @inbounds verticalsum[2, pixelx, pixely] = greenchannel * normfactor
      @inbounds verticalsum[3, pixelx, pixely] = bluechannel * normfactor
    end
  end
end

function boxfiltertriple!(filteredimage::Array{Float32,3},
    buffer::Array{Float32,3}, imagedata::Array{Float32,3},
    filtersize::Tuple{Int64, Int64})
  gethorizontalsumtriple!(buffer, imagedata, filtersize[1])
  getverticalsumtriple!(filteredimage, buffer, filtersize[2])
end

function boxfiltertriple(imagedata::Array{Float32,3},
    filtersize::Tuple{Int64, Int64})
  channels, imagewidth, imageheight = size(imagedata)
  filteredimage = Array{Float32, 3}(channels, imagewidth, imageheight)
  buffer = Array{Float32, 3}(channels, imagewidth, imageheight)
  gethorizontalsumtriple!(buffer, imagedata, filtersize[1])
  getverticalsumtriple!(filteredimage, buffer, filtersize[2])
  return filteredimage
end
