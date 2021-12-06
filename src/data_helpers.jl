using LinearAlgebra
using PyCall
@pyimport ehtim as eh



function amp_debias(amp, σ; force_nonzero=false)
    amp = abs(amp)
    σ = abs(σ)
    if amp < σ
        if force_nonzero
            deb² = σ^2
        else
            deb² = 0
        end
    else
        deb² = amp^2 - σ^2
    end
    return √deb²
end

function logcamp_debias(log_camp, snr1, snr2, snr3, snr4)
    log_camp_debias = log_camp + 0.5*(1/snr1^2 + 1/snr2^2 - 1/snr3^2 -1/snr4^2)
    return log_camp_debias
end


"""
    Stealing this from eht-imaging, as long as
    we use the quadrangles from ehtim then
    the definitions of numerators and denominators
    will work out.
"""
function make_log_closure_amplitude(n1amp, n2amp, d1amp, d2amp, n1err, n2err, d1err, d2err; debias=true)
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
        logcamp = logcamp_debias(logcamp, snr1, snr2, snr3, snr4)
    end
    return logcamp, logcamp_err
end

function amp_add_syserr(amp, amp_error; fractional=0, additive=0)
    σ = √(amp_error^2+(fractional*amp)^2+additive^2)
    return amp, σ
end

function logcamp_add_syserr(n1amp, n2amp, d1amp, d2amp, n1err, n2err, d1err, d2err; fractional=0, additive = 0, debias=true)
    n1amp, n1err = amp_add_syserr(n1amp, n1err, fractional=fractional, additive=additive)
    n2amp, n2err = amp_add_syserr(n2amp, n2err, fractional=fractional, additive=additive)
    d1amp, d1err = amp_add_syserr(d1amp, d1err, fractional=fractional, additive=additive)
    d2amp, d2err = amp_add_syserr(d2amp, d2err, fractional=fractional, additive=additive)
    return make_log_closure_amplitude(n1amp, n2amp, d1amp, d2amp, n1err, n2err, d1err, d2err, debias=debias)
end

