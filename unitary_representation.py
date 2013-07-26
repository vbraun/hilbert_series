"""
Unitary Representations from Characters

Representations of the unitary group can be encoded in (Laurent) polynomials:

  * U(1) characters are Laurent monomials 

  * U(n) characters are symmetric polynomials

  * SU(n) characters are symmetric polynomials with the relation that
    the product of all (the determinant) equals one.


EXAMPLES:

Single SU(3) instanton::

    sage: load('monomial_ideal.pyx')
    sage: load('unitary_representation.py')
    sage: R.<Q0, Q1, Q2, Qt0, Qt1, Qt2, psiA, psiB1> = QQ[]
    sage: H = R.ideal(Q0*Qt0 + Q1*Qt1 + Q2*Qt2)
    sage: Z.<z, u0, u1, u2> = LaurentPolynomialRing(QQ)
    sage: action = [u0*z^-1, u1*z^-1, u2*z^-1, u1*u2*z, u0*u2*z, u0*u1*z, 1, 1]
    sage: ideal = leading_terms_ideal(H.groebner_basis())
    sage: hs = ideal.Hilbert_series(refine=action, deg_bound=10)
    sage: t = hs.parent().gen(0)
    sage: hs = hs * (1-t)^2
    sage: for i, sf in hs.dict().iteritems():
    ....:     print i, UnitaryRepresentation(sf).powerof(z,0).yng(u0, u1, u2)
    0 s[]
    1 0
    2 s[2, 1]
    3 0
    4 s[4, 2]
    5 0
    6 s[6, 3]
    7 0
    8 s[8, 4]
    9 0
    10 s[10, 5]

`N=2` gauge theory with `N_c=3` and `N_f=6` ::

    sage: R.<Q00, Q01, Q02, Q03, Q04, Q05, Q10, Q11, Q12, Q13, Q14, Q15, Q20, Q21, Q22, Q23, Q24, Q25, Qt00, Qt01, Qt02, Qt10, Qt11, Qt12, Qt20, Qt21, Qt22, Qt30, Qt31, Qt32, Qt40, Qt41, Qt42, Qt50, Qt51, Qt52> = QQ[]
    sage: H = R.ideal(Q00*Qt00 - Q20*Qt02 + Q01*Qt10 - Q21*Qt12 + Q02*Qt20 - Q22*Qt22 
    ....:     + Q03*Qt30 - Q23*Qt32 + Q04*Qt40 - Q24*Qt42 + Q05*Qt50 - Q25*Qt52, 
    ....:     Q10*Qt00 + Q11*Qt10 + Q12*Qt20 + Q13*Qt30 + Q14*Qt40 + Q15*Qt50, 
    ....:     Q20*Qt00 + Q21*Qt10 + Q22*Qt20 + Q23*Qt30 + Q24*Qt40 + Q25*Qt50, 
    ....:     Q00*Qt01 + Q01*Qt11 + Q02*Qt21 + Q03*Qt31 + Q04*Qt41 + Q05*Qt51,
    ....:     Q10*Qt01 - Q20*Qt02 + Q11*Qt11 - Q21*Qt12 + Q12*Qt21 - Q22*Qt22 
    ....:     + Q13*Qt31 - Q23*Qt32 + Q14*Qt41 - Q24*Qt42 + Q15*Qt51 - Q25*Qt52,
    ....:     Q20*Qt01 + Q21*Qt11 + Q22*Qt21 + Q23*Qt31 + Q24*Qt41 + Q25*Qt51, 
    ....:     Q00*Qt02 + Q01*Qt12 + Q02*Qt22 + Q03*Qt32 + Q04*Qt42 + Q05*Qt52, 
    ....:     Q10*Qt02 + Q11*Qt12 + Q12*Qt22 + Q13*Qt32 + Q14*Qt42 + Q15*Qt52)
    sage: Z.<z0,z1,z2, x0,x1,x2,x3,x4,x5> = QQ[]
    sage: action_full = [x0*z1*z2,
    ....:     x1*z1*z2,
    ....:     x2*z1*z2,
    ....:     x3*z1*z2,
    ....:     x4*z1*z2,
    ....:     x5*z1*z2,
    ....:     x0*z0*z2,
    ....:     x1*z0*z2,
    ....:     x2*z0*z2,
    ....:     x3*z0*z2,
    ....:     x4*z0*z2,
    ....:     x5*z0*z2,
    ....:     x0*z0*z1,
    ....:     x1*z0*z1,
    ....:     x2*z0*z1,
    ....:     x3*z0*z1,
    ....:     x4*z0*z1,
    ....:     x5*z0*z1,
    ....:     x1*x2*x3*x4*x5*z0,
    ....:     x1*x2*x3*x4*x5*z1,
    ....:     x1*x2*x3*x4*x5*z2,
    ....:     x0*x2*x3*x4*x5*z0,
    ....:     x0*x2*x3*x4*x5*z1,
    ....:     x0*x2*x3*x4*x5*z2,
    ....:     x0*x1*x3*x4*x5*z0,
    ....:     x0*x1*x3*x4*x5*z1,
    ....:     x0*x1*x3*x4*x5*z2,
    ....:     x0*x1*x2*x4*x5*z0,
    ....:     x0*x1*x2*x4*x5*z1,
    ....:     x0*x1*x2*x4*x5*z2,
    ....:     x0*x1*x2*x3*x5*z0,
    ....:     x0*x1*x2*x3*x5*z1,
    ....:     x0*x1*x2*x3*x5*z2,
    ....:     x0*x1*x2*x3*x4*z0,
    ....:     x0*x1*x2*x3*x4*z1,
    ....:     x0*x1*x2*x3*x4*z2]
    sage: action = [z1*z2, z1*z2, z1*z2, z1*z2, z1*z2, z1*z2, z0*z2,
    ....:     z0*z2, z0*z2, z0*z2, z0*z2, z0*z2, z0*z1, z0*z1, z0*z1, z0*z1, z0*z1,
    ....:     z0*z1, z0, z1, z2, z0, z1, z2, z0, z1, z2, z0, z1, z2, z0, z1, z2,
    ....:     z0, z1, z2]
    sage: ideal = leading_terms_ideal(H.groebner_basis())
    sage: hs = ideal.Hilbert_series(refine=action, deg_bound=5)
    sage: for i, sf in hs.dict().iteritems():
    ....:     print i, UnitaryRepresentation(sf).special_invariants(z0, z1, z2)
    0 1
    1 0
    2 36
    3 40
    4 630
    5 1120
""" 
from sage.rings.all import QQ
from sage.structure.sage_object import SageObject



class UnitaryRepresentation(SageObject):

    def __init__(self, symmetric_function):
        self._sf = symmetric_function
        self._ring = symmetric_function.parent()
        self._sym = SymmetricFunctions(QQ).schur()

    def _repr_(self):
        return repr(self._sf)

    def collected_iter(self, *variables):
        """
        Collect terms by given ``variables``

        EXAMPLES::

            sage: R.<u,v,x,y> = QQ[]
            sage: r = UnitaryRepresentation((1+u+v)*(1+x+y) + u*v)
            sage: i = r.collected_iter(u,v)
            sage: i.next()
            (x + y, u + v + 1)
            sage: i.next()
            (1, u*v + u + v + 1)
            sage: i.next()
            Traceback (most recent call last):
            ...
            StopIteration
        """
        R = self._ring
        position = [R.gens().index(x) for x in variables]
        # first, collect {remaining_monomial:variables_polynomial}
        coll = dict()
        for c, m in self._sf:
            e = list(m.exponents()[0])
            e_pos = [0] * R.ngens()
            for i in position:
                e_pos[i] = e[i]
                e[i] = 0
            e = tuple(e)
            e_pos = tuple(e_pos)
            d = coll.get(e, dict())
            d[e_pos] = c
            coll[e] = d
        # second: combine monomials with same variables_polynomial
        # TODO: normalize first
        polynomials = dict()
        for e, var_poly in coll.iteritems():
            var_poly = R(var_poly)
            d = polynomials.get(var_poly, dict())
            d[e] = R.base_ring().one()
            polynomials[var_poly] = d
        # return result
        for var_poly, remaining in polynomials.iteritems():
            yield R(remaining), UnitaryRepresentation(var_poly)
        # check
        total = R.zero()
        for var_poly, remaining in polynomials.iteritems():
            total += R(remaining) * var_poly
        assert total == self._sf

    def yng(self, *variables):
        """
        Convert symmetric function into Young tableaux

        EXAMPLES::

            sage: R.<u,v,w> = QQ[]
            sage: r = UnitaryRepresentation(u^3 + v^3 + w^3)
            sage: r.yng()
            s[1, 1, 1] - s[2, 1] + s[3]
            sage: r.yng(u, v, w)
            s[1, 1, 1] - s[2, 1] + s[3]
            sage: t = UnitaryRepresentation(u^3*w + v^3*w)
            sage: t.yng(u, v)
            -s[2, 1] + s[3]
        """
        if len(variables) == 0:
            variables = self._ring.gens()
        n = len(variables)
        FF = self._ring.base_ring()
        R = PolynomialRing(FF, 'z', n)
        subs = [R.one()] * self._ring.ngens()
        for i, x in enumerate(variables):
            subs[self._ring.gens().index(x)] = R.gen(i)
        sf = R(self._sf(subs))
        return self._sym.from_polynomial(sf).restrict_partition_lengths(n, exact=False)

    def invariants(self, *variables):
        """
        Return the `U(n)` invariants

        EXAMPLES::

            sage: R.<u,v,x,y> = QQ[]
            sage: r = UnitaryRepresentation((1+u+v)*(1+x) + (1+u^2*v^2)*y)
            sage: r.invariants(u,v)
            x + y + 1
            sage: r.special_invariants(u,v)
            x + 2*y + 1
        """
        if len(variables) == 0:
            variables = self._ring.gens()
        inv = self._ring.zero()
        for remaining, rep in self.collected_iter(*variables):
            sf = rep.yng(*variables)
            for partition, coeff in sf:
                if partition == []:
                    inv += remaining*coeff
        return UnitaryRepresentation(inv)

    def special_invariants(self, *variables):
        """
        Return the `SU(n)` invariants

        See also :meth:`invariants`.

        EXAMPLES::

            sage: R.<u,v,x,y> = QQ[]
            sage: r = UnitaryRepresentation((1+u+v)*(1+x) + (1+u^2*v^2)*y)
            sage: r.invariants(u,v)
            x + y + 1
            sage: r.special_invariants(u,v)
            x + 2*y + 1
        """
        if len(variables) == 0:
            variables = self._ring.gens()
        inv = self._ring.zero()
        for remaining, rep in self.collected_iter(*variables):
            sf = rep.yng(*variables)
            for partition, coeff in sf:
                if partition == [] or partition == [partition[0]]*len(variables):
                    inv += remaining*coeff
        return UnitaryRepresentation(inv)

    def powerof(self, variable, exponent):
        """
        Return the coefficient of ``pow(variable, eponent)``

        EXAMPLES::

            sage: Z.<z, u0, u1, u2> = LaurentPolynomialRing(QQ)
            sage: r = UnitaryRepresentation(2*z^2+z^-1+u0+u1*u2)
            sage: r.powerof(z,-1)
            1
            sage: r.powerof(z,0)
            u1*u2 + u0
            sage: r.powerof(z,1)
            0
            sage: r.powerof(z,2)
            2
        """
        R = self._ring
        pos = R.gens().index(variable)
        coeff = dict()
        for c, m in self._sf:
            e = list(m.exponents()[0])
            if e[pos] != exponent:
                continue
            e[pos] = 0
            e = tuple(e)
            coeff[e] = c
        return UnitaryRepresentation(R(coeff))
    
    def dimension(self):
        """
        Return the dimension
        """
        return sum(self._sf.coefficients())
