abstract type TraitSymmetry end

struct Axisymmetric <: TraitSymmetry end
struct NoSymmetry <: TraitSymmetry end


abstract type RayTraceDirection end

struct Forward <: RayTraceDirection end
struct Reverse <: RayTraceDirection end


abstract type AccretionExtent end

struct Planar <: AccretionExtent end
struct Full <: AccretionExtent end
