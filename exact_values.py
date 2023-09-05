#!/usr/bin/env python3
from exact_fp import exact_parafermions_matrix, make_pf_spectrum2, rhoX, einf, gap_inf, einf2
from exact_fp import spectrum_matrix, pfs_matrix, fix_pfs, pfs_matrix_limited, make_m1
from fp_polynomial import find_eps_poly
from polytest import sympy_eps2
import numpy
from scipy.special import hyp2f1
from scipy.special import ellipeinc
import math

def gamma_new(g):
    # return (1.+g) / g
    return g / (1.+g)

L_limit = 500
def fix_parameters(s):
    if s['model_type'] == "fp":
        if not s['lambda_real'] == [] and not s['lambda_imag'] == []:
        # if not s['lambda_real'] is None and not s['lambda_imag'] is None:
            a = s['lambda_real']
            b = s['lambda_imag']
            c = a + b * 1.j
            s.config['lambda'] = abs(c)
            u = numpy.angle(c) * s['N'] / 2. / numpy.pi
            s.config['phi'] = u
            s.config['lambda_full'] = s['lambda'] * numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            s.config['lambda_angle'] = numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            l = s['lambda']
            s.config['gamma'] = l / (1.+l)
            # s['lambda'] = l
        # elif s['gamma'] is None and not s['lambda'] is None:
        elif s['gamma'] == [] and not s['lambda'] == []:
            l = s['lambda']
            s.config['gamma'] = l / (1.+l)
            s.config['lambda_full'] = s['lambda'] * numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            s.config['lambda_angle'] = numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            s.config['lambda_real'] = s.config['lambda_full'].real
            s.config['lambda_imag'] = s.config['lambda_full'].imag
        else:
            gamma = s['gamma']
            if abs(gamma) < 1E-60:
                l = 1E60
            # l = (1.0-gamma)/gamma
            if abs(1.-gamma) < 1E-16:
                l = 1E16
            else:
                l = gamma / (1.-gamma)
            if s['flipZ']:
                l = -l
            s.config['lambda'] = l
            s.config['lambda_full'] = s['lambda'] * numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            s.config['lambda_angle'] = numpy.exp(2.j * numpy.pi / s['N'] * s['phi'])
            s.config['lambda_real'] = s.config['lambda_full'].real
            s.config['lambda_imag'] = s.config['lambda_full'].imag

        s['omega'] = numpy.exp(2.j * numpy.pi / s['N'])
        gamma = s['gamma']
        if abs(gamma - 1) < 1E-60:
            s['lambda_scale'] = 1E60
        else:
            s['lambda_scale'] = 1./(1.-gamma)
    elif s['model_type'] == 'xyh':
        s['omega'] = numpy.exp(2.j * numpy.pi / s['N'])
        if not s['lambda_real'] == [] and not s['lambda_imag'] == []:
            a = s['lambda_real']
            b = s['lambda_imag']
            c = a + b * 1.j
            s.config['lambda'] = abs(c)
            g = (1-c) / (1+c)
            s['gamma'] = g
            # print(a,b)
            # exit()
        if not s['gamma_real'] == [] and not s['gamma_imag'] == []:
            a = s['gamma_real']
            b = s['gamma_imag']
            c = a + b * 1.j
            s.config['gamma'] = abs(c)
            u = numpy.angle(c) / 2. / numpy.pi * s["N"]
            s.config['phi'] = u
            s.config['gamma_full'] = c
        elif not s['lambda'] == [] and not s['phi'] == []:
            lf = s['lambda'] * numpy.exp(2.j*numpy.pi * s['phi']/s['N'])
            s.config['lambda_full'] = lf
            s.config['lambda_real'] = lf.real
            s.config['lambda_imag'] = lf.imag
            print(lf)
            c = lf
            g = (1-c) / (1+c)
            s.config['gamma'] = abs(g)
            s.config['gamma_full'] = g
        else:
            gf = s['gamma'] * numpy.exp(2.j*numpy.pi * s['phi']/s['N'])
            s.config['gamma_full'] = gf
            s.config['gamma_real'] = gf.real
            s.config['gamma_imag'] = gf.imag
        # s.config['lambda_real'] = s.config['gamma_real']
        # s.config['lambda_imag'] = s.config['gamma_imag']

    elif s['model_type'] == 'xyh2':
        s['omega'] = numpy.exp(2.j * numpy.pi / s['N'])
        if not s['eta_real'] == [] and not s['eta_imag'] == []:
            a = s['eta_real']
            b = s['eta_imag']
            c = a + b * 1.j
            s.config['eta_abs'] = abs(c)
            u = numpy.angle(c) / 2. / numpy.pi * s["N"]
            s.config['phi'] = u
            s.config['eta_full'] = c
            s.config['gamma_full'] = (1. - c) / (1. + c)
        elif not s['eta'] == [] and not s['phi'] == []:
            lf = s['eta'] * numpy.exp(2.j*numpy.pi * s['phi']/s['N'])
            s.config['eta_full'] = lf
            s.config['eta_real'] = lf.real
            s.config['eta_imag'] = lf.imag
            c = lf
            g = (1-c) / (1+c)
            s.config['gamma'] = abs(g)
            s.config['gamma_full'] = g

def make_exact_values_post_model(s):
    if s['skip_exact']:
        if not s['gamma_real'] == [] and not s['gamma_imag'] == []:
            1

    L = s['L']
    N = s['N']
    if L > L_limit:
        L = 500
    model_type = s.get('model_type', 'fp')
    gamma = s['gamma']

    sc = 1.
    mu = s['mu']
    phi = s['phi']
    if s['scale'] == 'gamma':
        sc = 1./(1.-gamma)
    sc *= numpy.exp(mu * phi * 2.j * numpy.pi / N)
    if model_type == 'fpXYW':
        L = s['L']
        if L > 15:
            return
        # pfs = find_eps_poly(s)
        # pfs = [p / L for p in pfs]
        pfs = sympy_eps2(s)
        s.pfs = pfs
        s['pfs'] = pfs
        e = 0
        for p in pfs:
            e -= p
        # e /= s['L']
        s['e_exact'] = e*L
        # s.spectrum = make_pf_spectrum2(e, pfs, s['N'])
        s['energies_pf'] = [x * L for x in s.spectrum]
        s['es_pf'] = [x for x in s.spectrum]
        # 1
def make_exact_values(s):
    if s['skip_exact']:
        return
    L = s['L']
    N = s['N']
    if L > L_limit:
        L = 500
    model_type = s.get('model_type', 'fp')
    gamma = s['gamma']

    sc = 1.
    mu = s['mu']
    phi = s['phi']
    if s['scale'] == 'gamma':
        sc = 1./(1.-gamma)
    sc *= numpy.exp(mu * phi * 2.j * numpy.pi / N)

    if model_type == 'fp':
        s['e_inf'] = einf(s) * sc
        s['e_inf2'] = einf2(s) * sc
        s['gap_inf'] = gap_inf(s) * sc

    if model_type == 'fp' and s.get('bc', 0) == 1:
        s['rhoX_exact'] = rhoX(s['N'], s['lambda'])
        s['rhoZ_exact'] = rhoX(s['N'], 1./s['lambda'])
    if model_type == 'fp' and s.get('bc', 0) == 0:
        e = None
        m = s['method']
        if m == 'poly' or 'matrix' in m: return
        if s['skip_pfs']:
            return

        # if abs(s['phi']) > 1E-16:
        #     pfs, e0 = pfs_matrix(s)
        #     s['e_exact'] = e0
        #     s['e_exact_lambda'] = e0 * s['lambda_scale']
        #     s.pfs = pfs
        # else:
        #     e, x, pfs = exact_parafermions_matrix(L, s['N'], s['lambda'], s['flipZ'])
        #     s['e_exact_lambda'] = e
        #     s['e_exact'] = e / s['lambda_scale']
        #     s.pfs = [p / s['lambda_scale'] for p in pfs]
        pfs, e0 = pfs_matrix(s)
        s['e_exact'] = e0 / s['L']
        # s['e_exact_lambda'] = e0 * s['lambda_scale']
        s.pfs = pfs
        s["pfs"] = pfs
        s['rhoX_exact'] = rhoX(s['N'], s['lambda'])
        s['rhoZ_exact'] = rhoX(s['N'], 1./s['lambda'])

        smallest_pf = None
        for p in pfs:
            if abs(p) > 1E-6:
                smallest_pf = p
                break
        spf = pfs[-1] * L
        lpf = smallest_pf * L
        if s['method'] == 'exact':
            s.spectrum = make_pf_spectrum2(s['e_exact'], s.pfs, s['N'])
        else:
            pass

    elif model_type == 'fpXYW':
        pass

    elif model_type == 'sxx':

        res = s

        h = s.get('h', 0)
        k = s['lambda']
        kk = k
        if (abs(kk-1) < 1e-10):
            kk = 1+1E-6

        gs1 = -math.sqrt((kk-1.)**2) * hyp2f1(-0.5, 0.5, 1, -4.*kk/(kk-1.)/(kk-1.))
        gs1_new = -math.sqrt((kk+1.)**2) * hyp2f1(-0.5, 0.5, 1, 4.*kk/(kk+1.)/(kk+1.))

        gs3 = -abs(h)

        q = (1.+k)**2 - h**2
        q /= 4.*k
        if q < 0:
            q = 0
        q = numpy.sqrt(q)
        if q > 1:
            gs2 = 0
        else:
            q = numpy.arcsin(q)
            gs2 = -2./math.pi*abs(1.+k) * ellipeinc(q, 4.*k/(1.+k)/(1.+k))
            gs2 += -2.*abs(h)/math.pi*(math.pi/2.-q)
            if numpy.isnan(gs2):
                gs2 = 0

        q1 = 1. - h**2 / 4.
        np = numpy
        q1 = np.arcsin(np.sqrt(q1))
        gs_k1 = -4./math.pi * ellipeinc(q1, 1.)
        gs_k1 += -2. * abs(h) / math.pi * (math.pi/2.-q1)
        res['gs_k1'] = gs_k1
        res['gs2'] = gs2
        res['gs1'] = gs1
        res['gs1_new'] = gs1_new
        res['gs3'] = gs3

        k = s['lambda']
        if (h <= abs(1-k)):
            res['gap_exact'] = abs(2 * abs(1-k) - 2*h)
            res['e_exact'] = gs1
        elif (h > abs(1+k)):
            res['gap_exact'] = abs(2 * abs(1+k) - 2*h)
            res['e_exact'] = gs3
        else:
            res['gap_exact'] = 0
            res['e_exact'] = gs2
        delta = abs(h-abs(1+k))-abs(h+abs(1+k))
        delta = abs(delta)
        res['ge2'] = delta
