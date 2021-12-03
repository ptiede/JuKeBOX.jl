using LinearAlgebra
using PyCall
@pyimport ehtim as eh



function amp_debias(amp, σ, force_nonzero=False)
    amp = abs(amp)
    σ = abs(σ)
    deb² = amp^2 - σ^2
    deb²[amp.<σ].=0 
    if force_nonzero
        deb²[amp.<σ] = σ^2
    end
    return √deb²
end


"""
    Stealing this from eht-imaging, as long as
    we use the quadrangles from ehtim then
    the definitions of numerators and denominators
    will work out.
"""
function make_log_closure_amplitude(n1amp, n2amp, d1amp, d2amp, n1err, n2err, d1err, d2err, debias=true)
    p1 = abs(n1amp)
    p2 = abs(n2amp)
    p3 = abs(d1amp)
    p4 = abs(d2amp)
    
    if debias
        p1 = amp_debias(p1, n1err)
        p2 = amp_debias(p2, n2err)
        p3 = amp_debias(p3, d1err)
        p4 = amp_debias(p4, d2err)
    end
    snr1 = p1/n1err
    snr2 = p2/n2err
    snr3 = p3/d1err
    snr4 = p4/d2err
    
    logcamp = log(p1)+log(p2)-log(p3)-log(p4)
    logcamp_err = √(1/snr1^2 + 1/snr2^2 + 1/snr3^2 + 1/snr4^2)
    if debias
        #WIP
    end
end

function amp_add_syserr(amp, amp_error, fractional=0, additive=0)
    σ = √(amp_error^2+(fractional*amp)^2+additive^2)
    return σ
end

