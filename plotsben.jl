"""
    carcaj(x0, y0, z0, Fx, Fy, Fz; samppx=2^3, samppy=2^3, samppz=2^3)

Computes the final coordinates (xf,yf,zf) for a given three-dimensional vector (Fx,Fy,Fz) and a starting point (x0,y0,z0) """
function carcaj(Xs, Ys, Zs, Fx, Fy, Fz; samppx=2^3,samppy=2^3,samppz=2^3)
    # CAREFUL!!!  Xs, Ys, Zs, Fx, Fy, Fz are three dimensional arrays!

    # Samples from the matrices that define the ellipses
    pointsx = size(Xs,1)
    pointsy = size(Ys,2)
    pointsz = size(Zs,3)
    spacex = div(pointsx,samppx)
    spacey = div(pointsy,samppy)
    spacez = div(pointsz,samppz)

    Xsamp = Xs[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]
    Ysamp = Ys[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]
    Zsamp = Zs[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]

    Fxsamp = Fx[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]
    Fysamp = Fy[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]
    Fzsamp = Fz[spacex+1:spacex:end, spacey+1:spacey:end, spacez+1:spacez:end]
    Fsamp = sqrt.(Fxsamp.^2 .+ Fysamp.^2 .+ Fzsamp.^2)

    minl = minimum(Fsamp)
    maxl = maximum(Fsamp)

    # Gets the size of the sampled arrays
    sampointsx = size(Fxsamp,1)
    sampointsy = size(Fxsamp,2)
    sampointsz = size(Fxsamp,3)

    # Scaling factor
    t = 0.3

    # Final coordinates
    Xf = Xsamp .+ log10.(Fsamp .+ 1.0).*Fxsamp*t./Fsamp
    Yf = Ysamp .+ log10.(Fsamp .+ 1.0).*Fysamp*t./Fsamp
    Zf = Zsamp .+ log10.(Fsamp .+ 1.0).*Fzsamp*t./Fsamp

    # Arrow head
    pointxx = 10
    xxmax = 0.25
    tp = collect(range(0.0,length=pointxx,stop=xxmax))
    ss = cat(50tp,4,dims=1)

    # Plots the lines given [x0, xf],[y0, yf],[z0, zf]
    for zind = 1:sampointsz
        for xind = 1:sampointsx
            for yind = 1:sampointsy
                xi = Xsamp[:,:,zind][xind,yind]
                xif = Xf[:,:,zind][xind,yind] .+ log10.(Fsamp[:,:,zind][xind,yind] .+ 1.0).*Fxsamp[:,:,zind][xind,yind].*tp./Fsamp[:,:,zind][xind,yind]
                yi = Ysamp[:,:,zind][xind,yind]
                yif = Yf[:,:,zind][xind,yind] .+ log10.(Fsamp[:,:,zind][xind,yind] .+ 1.0).*Fysamp[:,:,zind][xind,yind].*tp./Fsamp[:,:,zind][xind,yind]
                zi = Zsamp[:,:,zind][xind,yind]
                zif = Zf[:,:,zind][xind,yind] .+ log10.(Fsamp[:,:,zind][xind,yind] .+ 1.0).*Fzsamp[:,:,zind][xind,yind].*tp./Fsamp[:,:,zind][xind,yind]

                xp = cat(xi, xif, dims=1)
                yp = cat(yi, yif, dims=1)
                zp = cat(zi, zif, dims=1)

                # Condition to avoid
                if abs(xi+xif[1]+yi+yif[1]+zi+zif[1])>0
                    plot3d!(xp, yp, zp, linewidth=ss[end:-1:1],
                                line_z = ones(pointxx).*Fsamp[:,:,zind][xind,yind],
                                clims = (minl,maxl),
                                color = :cividis)
                end

            end
        end
    end
    plot3d!(leg=false, xlabel = "x", ylabel = "y", zlabel = "z", colorbar = false)
end

"""
    flecha(x0, y0, z0, Fx, Fy, Fz; tamano=0.3, ffirst=1)

Computes the final coordinates (xf,yf,zf) for a given three-dimensional vector (Fx,Fy,Fz) and a starting point (x0,y0,z0) """
function flecha(x, y, z, fx, fy, fz; tamano=0.3, cl = :white)

    # Location of the field
    xi=x
    yi=y
    zi=z

    # Scaling factor
    t = tamano

    # Magnitude of the field
    F = sqrt(fx^2 + fy^2 + fz^2)

    # Final coordinates
    xf = x + log10(F + 1.0)*fx*t/F
    yf = y + log10(F + 1.0)*fy*t/F
    zf = z + log10(F + 1.0)*fz*t/F

    # Head arrows
    pointxx = 10
    xxmax = 0.25
    tp = collect(range(0.0,length=pointxx,stop=xxmax))
    ss = cat(50tp,4,dims=1)

    xif = xf .+ log10(F + 1.0)*fx.*tp./F
    yif = yf .+ log10(F + 1.0)*fy.*tp./F
    zif = zf .+ log10(F + 1.0)*fz.*tp./F

    xp = cat(xi, xif, dims=1)
    yp = cat(yi, yif, dims=1)
    zp = cat(zi, zif, dims=1)

    # Condition to avoid
    if abs(xi+xif[1]+yi+yif[1]+zi+zif[1])>0
        plot3d!(xp, yp, zp, linewidth=ss[end:-1:1],
                color = cl)
    end
end

"""
    ejes(winsize; cl=:white)

Generates axis of the give size and color """
function ejes(winsize; cl = :white)
    # Draw axis to overcome weird scaling behaviour
    pointxx = 10
    xxmax = 0.1
    tp = collect(range(0.0,length=pointxx,stop=xxmax))
    ss = cat(winsize*35*tp,1,1,dims=1)
    a1 = cat(-winsize, winsize, winsize.+tp, dims=1)
    a2 = zeros(pointxx+2)
    a3 = zeros(pointxx+2)
    plot3d!(a1, a2, a3, color=cl, linewidth=ss[end:-1:1])  # X-axis
    plot3d!(a2, a1, a3, color=cl, linewidth=ss[end:-1:1])  # Y-axis
    plot3d!(a2, a3, a1, color=cl, linewidth=ss[end:-1:1])  # Z-axis
end


"""
    volumen(Xs, Ys, Zs; nslice=16)

Plots a three dimensional scalar field """
function volumen(Xs, Ys, Zs, F; nslice=16)

    # Calculates spacing according to the desired number of slices
    pointsx, pointsy, pointsz = size(Xs)
    space = div(pointsz,nslice)

    xs = Xs[:,1,1]
    ys = Ys[1,:,1]
    zs = Zs[1,1,:]
    fillz = F
    minl = minimum(fillz)
    maxl = maximum(fillz)
    ml = maximum([abs(minl) abs(maxl)])
    fillz[1,1,:] .= ml
    fillz[pointsx-1,1,:] .= -ml

    for ii=1:space:pointsz-1
    surface!(xs, ys, zs[ii] .+ 0.0.*Xs[:,:,1], color=:lisbon, alpha=0.3,
                  xlabel = "x", ylabel = "y", zlabel = "z",
                  fill_z=fillz[:,:,ii])  # Flipped Xs and Ys
    end
end
