using Pkg; Pkg.activate(@__DIR__)
using StructArrays
import JuKeBOX as Jk
using GLMakie
using Makie

using JuKeBOX
using ForwardDiff
using ComradeBase

x = range(-15.0, 15.0, length=128)
y = range(-15.0, 15.0, length=128)

ix,iy = Tuple(ind)

g, acc = bam(2, spin, 1.0, 6.0, 5.0, 0.9, π/2, π/2)
o = Observer(1.0, inc)

s = Jk.SimpleModel(acc, g, o)


function f(spin, inc, x, y)
    g, acc = bam(1, spin, 1.0, 6.0, 5.0, 0.9, π/2, π/2)
    o = Jk.Observer(1.0, inc)
    first(raytrace(x, y, g, o, acc))
end

function img(spin, inc)
    p = [spin, inc]
    x = range(-15.0, 15.0, length=128)
    y = range(-15.0, 15.0, length=128)
    println(x)
    broadcast(x,y') do xx,yy
        g(x) = f(x[1], x[2], xx, yy)
        ForwardDiff.gradient(g, p)
    end
end

pimg = f.(spin, inc, x, y')
heatmap(pimg)

spin = rand()
inc = rand()*π/2

dimg = img(spin, inc)
incderiv = last.(dimg)
ind = sum(isnan.(incderiv))
heatmap(incderiv)
aderiv = first.(dimg)
heatmap(aderiv)

fov = 30.0
alpha = range(-15.0, 15.0, length=124)
beta  = range(-15.0, 15.0, length=124)
polim = StructArray((I=zeros(length(alpha), length(beta)), Q = zeros(length(alpha), length(beta)), U = zeros(length(alpha), length(beta))))

inc = Observable(0.01)
spin = Observable(0.5)
g, acc = Jk.bam(2, 0.25, 1.0, 6.0, 5.0, 0.9, π/2, π/2)
l1 = on(inc) do val
    obs = Jk.Observer(1.0, val)
    g = Jk.Kerr(spin.val)
    t = @elapsed Jk.traceimg!(polim, alpha, beta, g, obs, acc)
    println("I took $(t*1000.0) ms to run")
    io[] = polim.I'
end

l2 = on(spin) do val
    obs = Jk.Observer(1.0, inc.val)
    g = Jk.Kerr(val)
    t = @elapsed Jk.traceimg!(polim, alpha, beta, g, obs, acc)
    println("I took $(t*1000.0) ms to run")
    println("Using spin $val and inc $(inc.val)")
    io[] = polim.I'
end

io = Observable(polim.U')

fig, ax = image(alpha, beta, io, colormap=:afmhot,
                axis=(aspect=1, xlabel="RA (μas)", ylabel="DEC (μas)"),
                #title="Inclination = $inc"
                )
Makie.deactivate_interaction!(ax, :rectanglezoom)
spoint = select_point(ax.scene)

on(spoint) do z
    x,y = z
    inc[] =  (fov/2+y)/fov*π/2
    spin[] = x/(fov/2)
end


function test(a, θ)
    g, acc = bam(1, a, 1.5, 5.0, 1.0, 0.9, π/2, π/3)
    o = Observer(1.0, θ)
    #s = Jk.SimpleModel(acc, g, o)
    #sI = ComradeBase.SingleStokes(s, :I)
    v = raytrace(5.0, 2.0, g, o, acc)
    return norm(StokesVector(v[1],v[2],v[3],v[4]))
end


function foo(u)
    x = u[1]
    y = u[2]
    v = ComradeBase.StokesVector((x^2, y^2, (x-y)^2, x^2+y^2))
    return norm(v)
end
