
struct Composite{A,B,C}
    a::A
    b::B
    c::B
end

struct Ax
    x::Float64
end
struct Bx
    x::Float64
end

struct Cx
    x::Float64
end



abstract type Symmetry end

struct HasSym <: Symmetry end
struct NoSym  <: Symmetry end
struct AllSym <: Symmetry end

symmetry(::Type{Ax}) = HasSym()
symmetry(::Type{Bx}) = NoSym()
symmetry(::Type{Cx}) = AllSym()

getsym(a::T) where {T}  = getsym(symmetry(T), a)

getsym(::HasSym, x) = x+1
getsym(::NoSym, x) = x-1


getsym(A, 0.0)
