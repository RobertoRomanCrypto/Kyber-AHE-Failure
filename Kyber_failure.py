import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from proba_util import *

def p2_cyclotomic_final_error_distribution_homomorphic_2(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    # e,s
    chis = build_centered_binomial_law(ps.ks)
    # e1,e2
    chie = build_centered_binomial_law(ps.ke_ct)
    # r
    chie_pk = build_centered_binomial_law(ps.ke)
    # cu
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)
    # cv
    R2 = build_mod_switching_error_law(ps.q, ps.rq2)

    # r1 + r2
    chir1r2 = law_convolution(chie_pk, chie_pk)
    # e * (r1 + r2)
    B1 = law_product(chis, chir1r2)
    C1 = iter_law_convolution(B1, ps.m * ps.n)

    # e1,1 + e1,2
    chie11e12 = law_convolution(chie, chie)
    # cu,1 + cu,2
    chicu1cu2 = law_convolution(Rc, Rc)
    # cu,1 + cu,2 + cu
    chicu = law_convolution(chicu1cu2, Rc)
    # (e1,1 + e1,2) + (cu,1 + cu,2 + cu)
    chie1cu = law_convolution(chie11e12, chicu)
    # s * ((e1,1 + e1,2) + (cu,1 + cu,2 + cu))
    B2 = law_product(chis, chie1cu)
    C2 = iter_law_convolution(B2, ps.m * ps.n)

    # e2,1 + e2,2
    C3 = law_convolution(chie, chie)

    # cv,1 + cv,2
    chicv1cv2 = law_convolution(R2, R2)
    # cv,1 + cv,2 + cv
    C4 = law_convolution(chicv1cv2, R2)

    # C1 + C2
    C5 = law_convolution(C1, C2)
    # C3 + C4
    C6 = law_convolution(C3, C4)

    # C1 + C2 + C3 + C4
    C7 = law_convolution(C5, C6)

    return C7

def p2_cyclotomic_final_error_distribution_homomorphic_3(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    # -- Cómputo distribuciones --
    # e,s
    chis = build_centered_binomial_law(ps.ks)
    # e1,e2
    chie = build_centered_binomial_law(ps.ke_ct)
    # r
    chie_pk = build_centered_binomial_law(ps.ke)
    # cu
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)
    # cv
    R2 = build_mod_switching_error_law(ps.q, ps.rq2)

    # r1 + r2
    chir1r2 = law_convolution(chie_pk, chie_pk)
    # r1 + r2 + r3
    chir1r2r3 = law_convolution(chir1r2, chie_pk)

    # e * (r1 + r2 + r3)
    B1 = law_product(chis, chir1r2r3)
    C1 = iter_law_convolution(B1, ps.m * ps.n)

    # e1,1 + e1,2
    chie11e12 = law_convolution(chie, chie)
    # e1,1 + e1,2 + e1,3
    chie11e12e13 = law_convolution(chie11e12, chie)

    # cu,1 + cu,2
    chicu1cu2 = law_convolution(Rc, Rc)
    # cu,1 + cu,2 + cu,3
    chicu1cu2cu3 = law_convolution(chicu1cu2, Rc)

    # cu,1 + cu,2 + cu,3 + cu
    chicu = law_convolution(chicu1cu2cu3, Rc)

    # (e1,1 + e1,2 + e1,3) + (cu,1 + cu,2 + cu,3 + cu)
    chie1cu = law_convolution(chie11e12e13, chicu)

    # s * ((e1,1 + e1,2) + (cu,1 + cu,2 + cu))
    B2 = law_product(chis, chie1cu)
    C2 = iter_law_convolution(B2, ps.m * ps.n)

    # e2,1 + e2,2
    C3temp = law_convolution(chie, chie)
    # e2,1 + e2,2 + e2,3
    C3 = law_convolution(C3temp, chie)

    # cv,1 + cv,2
    chicv1cv2 = law_convolution(R2, R2)
    # cv,1 + cv,2 + cv,3
    chicv1cv2cv3 = law_convolution(chicv1cv2, R2)

    # cv,1 + cv,2 + cv,3 + cv
    C4 = law_convolution(chicv1cv2cv3, R2)

    # C1 + C2
    C5 = law_convolution(C1, C2)
    # C3 + C4
    C6 = law_convolution(C3, C4)

    # C1 + C2 + C3 + C4
    C7 = law_convolution(C5, C6)

    return C7

def p2_cyclotomic_error_probability_homomorphic_2(ps):
    F = p2_cyclotomic_final_error_distribution_homomorphic_2(ps)
    proba = tail_probability(F, ps.q/4)
    return F, ps.n*proba

def p2_cyclotomic_error_probability_homomorphic_3(ps):
    F = p2_cyclotomic_final_error_distribution_homomorphic_3(ps)
    proba = tail_probability(F, ps.q/4)
    return F, ps.n*proba