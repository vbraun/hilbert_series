"""
Gauge Theories

    sage: X = gauge_theories.threeSU(3)
    sage: X.betti(10)
    [1, 36, 658, 8422, 82232, 657996, 4496356, 26978272, 145008501, 708930508]
    sage: X.Higgs_branch().hilbert_series().factor()
    (t - 1)^-28 * (t + 1)^8


    sage: deg_bound = 5
    sage: X = gauge_theories.E6theory().Higgs_branch()
    sage: gb = X.groebner_basis(deg_bound=deg_bound)
    sage: i = leading_terms_ideal(gb)
    sage: i.Hilbert_series(deg_bound=deg_bound)
"""
from sage.all_cmdline import *

from sage.misc.cachefunc import cached_method, cached_function
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.combinat.cartesian_product import CartesianProduct
from sage.matrix.constructor import matrix





#####################################################

@cached_function
def _Hilbert_series_recursion(R, ideal_gens, T, t):
    r"""
    Recursively compute the Hilbert series

    OUPUT:

    The numerator $P(t)$ of $HS(t) = P(t) / (1-t)^n$, where $n$ is the
    number of variables in the polynomial ring.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: T.<t> = QQ[]
        sage: _Hilbert_series_recursion(R, (x**3*y, x**2*y**2, x*y**7), T, t)
        t^9 - t^8 + t^5 - 2*t^4 + 1
    """
    if len(ideal_gens) == 0:
        return T.one()
    f = ideal_gens[0]
    f_deg = f.total_degree()
    ideal_gens = tuple(ideal_gens[1:])
    if len(ideal_gens) == 0:
        return T.one() - t ** f_deg
    HS_remaining = _Hilbert_series_recursion(R, ideal_gens, T, t)
    J = [i // f.gcd(i) for i in ideal_gens]
    # J = R.ideal(J).interreduced_basis()
    J = R.ideal(J, coerce=False)
    from sage.rings.polynomial.multi_polynomial_ideal_libsingular import interred_libsingular
    J = interred_libsingular(J)
    
    ### We could wrap things up now, but we will implement the step A below
    # HS_J = _Hilbert_series_recursion(R, tuple(J), T, t)
    # return HS_remaining - t ** f_deg * HS_J

    J = sorted(J)
    J_linear = []
    J_nonlinear = []
    for j in J:
        if j.total_degree() == 1:
            J_linear.append(j)
        else:
            J_nonlinear.append(j)
    HS_J = _Hilbert_series_recursion(R, tuple(J_nonlinear), T, t)
    return HS_remaining - t ** f_deg * HS_J * (1 - t) ** len(J_linear)


    
def Hilbert_series(self, deg_bound=5, name='t', output='powerseries'):
    """
    Compute the Hilbert series of $R/I$

    ALGORITHM:
    
    Uses the exact sequence

    .. math::

        0 \to R/J \to R/I' \to R/I \to 0

    where $I'$ is $I$ with one generator, say, $f$, removed. The
    kernel is
    
    .. math::

        J = \{i/gcd(i,f) | i \in I\}

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: I = R.ideal(x**3*y, x**2*y**2, x*y**7)
        sage: Hilbert_series(I, 10)
        1 + 2*t + 3*t^2 + 4*t^3 + 3*t^4 + 3*t^5 + 3*t^6 + 3*t^7 + 2*t^8 + 2*t^9 + 2*t^10 + O(t^11)
    """
    n = self.ring().ngens()
    gb = self.groebner_basis(deg_bound=deg_bound)
    # Theorem: HS(LM(I)) = HS(I)
    # In the monomial ideal all coefficients are 1, so use GF(2) for speed
    from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
    R = MPolynomialRing_libsingular(GF(2), n, reversed(self.ring().gens()))
    ideal_gens = tuple(sorted([R(gen.lm()) for gen in gb]))
    T = PowerSeriesRing(QQ, name, default_prec=deg_bound+1, order='lex')
    t = T.gen(0)
    HS = _Hilbert_series_recursion(R, ideal_gens, T, t)
    #return HS
    return HS / (1-t)**n







#####################################################
def cartesian_product_01(n):
    return CartesianProduct(*([range(0,2)]*n))


#####################################################
epsilon = matrix(ZZ,[[0,1],[-1,0]])

#####################################################
class NodeUnitary(SageObject):
    def __init__(self, rank, index=0, external=False, special=False):
        self.rank = rank
        self.external = external
        self.index = index
        self.special = special

    def _repr_(self):
        s = ''
        if self.external:
            s += 'global '
        s += 'U_'+str(self.index)+'('+str(self.rank)+')'
        return s
    
    def character_str(self):
        return ['z{0}_{1}'.format(str(self.index), str(i)) 
                for i in range(self.rank)]

    def characters(self):
        R = self.theory.character_ring()
        return map(R, self.character_str())

    

#####################################################
class FieldBifundamental(SageObject):
    def __init__(self, src_node, dst_node, name):
        self.src = src_node
        self.dst = dst_node
        assert self.src.index != self.dst.index
        self.name = name
    
    def _repr_(self):
        s = self.name + \
            '({0},{1})'.format(str(self.src.index), str(self.dst.index))
        return s
    
    def components_str(self):
        result = []
        base = self.name + '{0}{1}_'.format(str(self.src.index), str(self.dst.index))
        for i in range(self.src.rank):
            for j in range(self.dst.rank):
                result.append(base + '{0}{1}'.format(i,j))
        return result

    def components(self):
        R = self.theory.ring()
        return map(R, self.components_str())

    def matrix(self):
        R = self.theory.ring()
        return matrix(R, self.src.rank, self.dst.rank, self.components())
    
    def group_action(self):
        """
        The group action on the component fields

        EXAMPLES::

            sage: X = gauge_theories.Bifundametal(3,2)
            sage: Q = X.fields()[0]
            sage: Q.group_action()
            {Q01_21: z0_0*z0_1*z1_1,
             Q01_20: z0_0*z0_1*z1_0,
             Q01_11: z0_0*z0_2*z1_1,
             Q01_10: z0_0*z0_2*z1_0,
             Q01_01: z0_1*z0_2*z1_1,
             Q01_00: z0_1*z0_2*z1_0}
        """
        src_chars = self.src.characters()
        dst_chars = self.dst.characters()
        dual = prod(src_chars)
        src_rep = [dual // x for x in src_chars]
        dst_rep = dst_chars
        action = dict()
        m = self.matrix()
        for i in range(self.src.rank):
            for j in range(self.dst.rank):
                action[m[i,j]] = src_rep[i] * dst_rep[j]
        return action


class FieldAdjoint(SageObject):
    def __init__(self, node, name):
        self.node = node
        self.name = name
        self.exclude_trace = node.special
    
    def _repr_(self):
        s = self.name + '({0})'.format(str(self.node.index))
        return s
    
    def components_str(self):
        result = []
        base = self.name + '{0}_'.format(str(self.node.index))
        n = self.node.rank
        for i in range(n):
            for j in range(n):
                if not (self.exclude_trace and (i,j) == (n-1,n-1 )):
                    result.append(base + '{0}{1}'.format(i,j))
        return result

    def components(self):
        R = self.theory.ring()
        return map(R, self.components_str())

    def matrix(self):
        R = self.theory.ring()
        n = self.node.rank
        if self.exclude_trace:
            components = self.components() + [R.zero()]
        else:
            components = self.components()
        m = matrix(R, n, n, components)
        if self.exclude_trace:
            tr = m.trace()
            m[n-1, n-1] = -m.trace()
        return m
        
    def group_action(self):
        """
        The group action on the component fields

        EXAMPLES::

            sage: X = gauge_theories.Adjoint(3)
            sage: phi = X.fields()[0]
            sage: phi.group_action()
            {phi0_21: z0_0*z0_1^2,
             phi0_20: z0_0^2*z0_1,
             phi0_12: z0_0*z0_2^2,
             phi0_11: z0_0*z0_1*z0_2,
             phi0_10: z0_0^2*z0_2,
             phi0_02: z0_1*z0_2^2,
             phi0_01: z0_1^2*z0_2,
             phi0_00: z0_0*z0_1*z0_2}
        """
        src_chars = self.node.characters()
        dst_chars = self.node.characters()
        dual = prod(src_chars)
        src_rep = [dual // x for x in src_chars]
        dst_rep = dst_chars
        action = dict()
        n = self.node.rank
        m = self.matrix()
        for i in range(n):
            for j in range(n):
                if not (self.exclude_trace and (i,j) == (n-1,n-1)):
                    action[m[i,j]] = src_rep[i] * dst_rep[j]
        return action

        


#####################################################
class GaugeTheory(SageObject):
    """
    N=1 Gauge Theory
    """

    def __init__(self):
        super(GaugeTheory, self).__init__()
        self._gauge_groups = []
        self._fields = []
        self._from_vector = []
        # self._ring_extra_args = {'order':'lex'}
        self._ring_extra_args = {}

    def ring(self):
        """
        The polynomial ring of field components
        """
        components = []
        for field in self.fields():
            components += field.components_str()
        R = PolynomialRing(QQ, components, **self._ring_extra_args)
        return R

    def _Higgs_branch_ring(self):
        """
        The polynomial ring of hypermultiplet field components.
        """
        components = []
        for field in self.fields():
            if field not in self._from_vector:
                components += field.components_str()
        R = PolynomialRing(QQ, components, **self._ring_extra_args)
        return R

    def character_ring(self):
        characters = []
        for G in self.gauge_groups():
            characters += G.character_str()
        R = PolynomialRing(QQ, characters, **self._ring_extra_args)
        return R

    def gauge_groups(self):
        return tuple(self._gauge_groups)

    def get_gauge_group(self, G):
        if G in self.gauge_groups():
            return G
        return self.gauge_groups()[G]

    def fields(self):
        return tuple(self._fields)

    def add_gauge_group(self, group, rank, external):
        n = len(self._gauge_groups)
        if group == 'SU':
            node = NodeUnitary(rank, index=n, external=external, special=True)
        elif group == 'U':
            node = NodeUnitary(rank, index=n, external=external, special=False)
        else:
            raise NotImplementedError('unknown gauge group')
        node.theory = self
        self._gauge_groups.append(node)
        return node

    def add_SU(self, rank, external=False):
        return self.add_gauge_group('SU', rank, external=external)

    def add_U(self, rank, external=False):
        return self.add_gauge_group('U', rank, external=external)
        
    def add_bifundamental(self, src, dst, name):
        src = self.get_gauge_group(src)
        dst = self.get_gauge_group(dst)
        field = FieldBifundamental(src, dst, name)
        field.theory = self
        self._fields.append(field)
        return field

    def add_adjoint_hyper(self, node, name):
        """
        Add an N=1 adjoint field from an N=2 hypermultiplet
        """
        node = self.get_gauge_group(node)
        field = FieldAdjoint(node, name)
        field.theory = self
        self._fields.append(field)
        return field

    def add_adjoint_vector(self, node, name):
        """
        Add an N=1 adjoint field coming from an N=2 hypermultiplet

        The distinction to :meth:`add_adjoint_hyper` is only to
        remember what is set to zero in the :meth:`Higgs_branch`
        """
        field = self.add_adjoint_hyper(node, name)
        self._from_vector.append(field)
        return field

    def set_superpotential(self, W):
        self.W = W

    def Fterms(self, field):
        return [self.W.derivative(x) for x in field.components()]
        
    def Fideal(self, *fields):
        relations = []
        for field in fields:
            relations.extend(self.Fterms(field))
        return self.ring().ideal(relations)

    def Higgs_branch(self):
        ideal = self.Fideal(*self._from_vector)
        R = self._Higgs_branch_ring()
        return R.ideal(ideal)

    def plot(self, **kwds):
        """
        EXAMPLES::
    
            sage: X = gauge_theories.threeSU(3)
            sage: print X.plot().description()
        """
        g = DiGraph(loops=True, sparse=True, multiedges=True)
        for G in self._gauge_groups:
            g.add_vertex(G)
        for field in self._fields:
            if isinstance(field, FieldBifundamental):
                g.add_edge(field.src, field.dst, field)
            if isinstance(field, FieldAdjoint):
                g.add_edge(field.node, field.node, field)
        return g.plot(vertex_labels=True, edge_labels=True, graph_border=True)

    def Hilbert_series(self, limit=4):
        """
        EXAMPLES::

            sage: gauge_theories.threeSU(3).Hilbert_series()
        """
        return Hilbert_series(self.Higgs_branch(), deg_bound=limit)

    def Fflat_betti(self, limit=4):
        """
        Hilbert function of the Higgs branch F-flat space

        EXAMPLES::

            gauge_theories.threeSU(3).n_Fflat(10)
            [1, 36, 658, 8422, 82232, 657996, 4496356, 26978272, 145008501, 708930508]
        """
        F = self.Higgs_branch()
        ngens = F.ring().ngens()
        F_gb = F.groebner_basis(deg_bound=limit)
        F_lm = set([f.lm() for f in F_gb])
        F_lm_count = [
            sum(1 for f in F_lm if f.total_degree()==i)
            for i in range(limit) ]
        betti = [binomial(ngens+i-1, i) - F_lm_count[i] for i in range(limit)]
        return betti

    def ngens_groebner_basis(self, limit=3):
        """
        Number of the generators in the groebner basis of the Higgs branch F-flat space

        EXAMPLES::

            sage: gauge_theories.threeSU(3).ngens_groebner_basis(10)
            [0, 0, 8, 14, 19, 12, 32, 56, 12, 0]
        """
        F = self.Higgs_branch()
        F_gb = F.groebner_basis(deg_bound=limit)
        F_lm = set([f.lm() for f in F_gb])
        F_lm_count = [
            sum(1 for f in F_lm if f.total_degree()==i)
            for i in range(limit) ]
        return F_lm_count
    
    def _group_action(self):
        """
        Return the group action (characters) on the field components

        EXAMPLES::

            sage: X = gauge_theories.threeSU(3)
            sage: phi = X.fields()[-1];  phi
            phi(1)
            sage: F = X.Fterms(phi)
            sage: act = X._group_action()
            sage: for f in F:
            ....:     print set(m.subs(act).exponents()[0] for m in f.monomials())
            set([(1, 1, 1, 1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1, 1, 1, 1)])
            set([(1, 1, 1, 2, 0, 1, 0, 0, 0), (0, 0, 0, 2, 0, 1, 1, 1, 1)])
            set([(1, 1, 1, 2, 1, 0, 0, 0, 0), (0, 0, 0, 2, 1, 0, 1, 1, 1)])
            set([(1, 1, 1, 0, 2, 1, 0, 0, 0), (0, 0, 0, 0, 2, 1, 1, 1, 1)])
            set([(1, 1, 1, 1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1, 1, 1, 1)])
            set([(0, 0, 0, 1, 2, 0, 1, 1, 1), (1, 1, 1, 1, 2, 0, 0, 0, 0)])
            set([(0, 0, 0, 0, 1, 2, 1, 1, 1), (1, 1, 1, 0, 1, 2, 0, 0, 0)])
            set([(0, 0, 0, 1, 0, 2, 1, 1, 1), (1, 1, 1, 1, 0, 2, 0, 0, 0)])
        """
        action = dict()
        for field in self.fields():
            action.update(field.group_action())
        return action

    def _group_action_Higgs_branch(self):
        R = self._Higgs_branch_ring()
        action_all = self._group_action()
        action = dict()
        for monomial, character in action_all.iteritems():
            try:
                action[R(monomial)] = character
            except TypeError:
                pass
        return action

    def monomials_characters(self, degree):
        """
        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()

            sage: X = GaugeTheory()
            sage: G1 = X.add_SU(2, external=True)
            sage: G3 = X.add_SU(2, external=True)
            sage: G2 = X.add_SU(3)
            sage: Q1tilde = X.add_bifundamental(G1, G2, 'Q1t')
            sage: Q1 = X.add_bifundamental(G2, G1, 'Q1')
            sage: Q2tilde = X.add_bifundamental(G2, G3, 'Q2t')
            sage: Q2 = X.add_bifundamental(G3, G2, 'Q2')
            sage: phi = X.add_adjoint_vector(G2, 'phi')
            sage: X.set_superpotential((Q1tilde.matrix() * phi.matrix() * Q1.matrix()).trace() + (Q2.matrix() * phi.matrix() * Q2tilde.matrix()).trace() )
            
            sage: #chi = X.monomials_characters(4)
            sage: chi = X.relation_characters(3)
            sage: R.<u,v,w> = QQ[]
            sage: chi = chi.subs(z0_0=1, z0_1=1, z1_0=1, z1_1=1, z2_0=u, z2_1=v, z2_2=w)
            sage: s.from_polynomial(chi)



            sage: X = GaugeTheory()
            sage: G1 = X.add_SU(1, external=True)
            sage: G2 = X.add_SU(3)
            sage: Qtilde = X.add_bifundamental(G1, G2, 'Qt')
            sage: Q = X.add_bifundamental(G2, G1, 'Q')
            sage: phi = X.add_adjoint_vector(G2, 'phi')
            sage: X.set_superpotential((Qtilde.matrix() * phi.matrix() * Q.matrix()).trace())
            sage: #chi = X.monomials_characters(4)
            sage: chi = X.relation_characters(3)
            sage: R.<u,v,w> = QQ[]
            sage: chi = chi.subs(z0_0=1, z1_0=u, z1_1=v, z1_2=w)
            sage: s.from_polynomial(chi)


        
            sage: X = gauge_theories.threeU(3)

            sage: chi = X.monomials_characters(1)
            sage: s = SymmetricFunctions(QQ).schur()
            sage: R.<u,v,w> = QQ[]
            sage: chi = chi.subs(z0_0=1, z0_1=1, z0_2=1, z1_0=1, z1_1=1, z1_2=1, z2_0=u, z2_1=v, z2_2=w)
            sage: s.from_polynomial(chi)

            sage: chi = X.relation_characters(1)
            sage: chi = chi.subs(z0_0=1, z0_1=1, z0_2=1, z1_0=1, z1_1=1, z1_2=1, z2_0=u, z2_1=v, z2_2=w)
            sage: s.from_polynomial(chi)
        """
        F = self.Higgs_branch()
        R = F.ring()
        M = R.ideal(R.gens())
        M = M**degree
        M = set(M.gens())
        action = self._group_action_Higgs_branch()
        return sum([m.subs(action) for m in M])
        

    def relation_characters(self, degree):
        """
        Return the relations in the given degree as group characters

        EXAMPLES::
        
            sage: X = gauge_theories.threeSU(3)
            sage: X = gauge_theories.threeU(3)
            sage: chi = X.relation_characters(3)
            sage: zip(chi.coefficients(), chi.exponents())

        """
        F = self.Higgs_branch()
        R = F.ring()
        F_gb = F.groebner_basis(deg_bound=degree)
        F_lm = set([f.lm() for f in F_gb if f.total_degree() == degree])
        action_all = self._group_action()
        action = dict()
        for monomial, character in action_all.iteritems():
            try:
                action[R(monomial)] = character
            except TypeError:
                pass
        return sum([f.subs(action) for f in F_lm])
    


#####################################################
class GaugeTheoryLibrary:

    def Bifundametal(self, N1, N2):
        """
        A single bifundamental field from G1 to G2

        EXAMPLES::

            sage: X = gauge_theories.Bifundametal(2,3)
            sage: X.gauge_groups()
            (SU_0(2), SU_1(3))
            sage: X.fields()
            (Q(0,1),)
            sage: X.fields()[0].matrix()
            [Q01_00 Q01_01 Q01_02]
            [Q01_10 Q01_11 Q01_12]
        """
        X = QuiverGaugeTheory()
        G1 = X.add_SU(N1)
        G2 = X.add_SU(N2)
        Q = X.add_bifundamental(G1, G2, 'Q')
        return X

    def Adjoint(self, N):
        """
        A single $SU(N)$ adjoint field

        EXAMPLES::

            sage: X = gauge_theories.Adjoint(3)
            sage: X.gauge_groups()
            (SU_0(3),)
            sage: X.fields()
            (Q(0,1),)
            sage: X.fields()[0].matrix()
            [           phi0_00            phi0_01            phi0_02]
            [           phi0_10            phi0_11            phi0_12]
            [           phi0_20            phi0_21 -phi0_00 - phi0_11]
        """
        X = QuiverGaugeTheory()
        G = X.add_SU(N)
        phi = X.add_adjoint_hyper(G, 'phi')
        return X
            
    def InstantonNequals2(self, N_global, N_local):
        X = GaugeTheory()
        G_global = X.add_SU(N_global, external=True)
        G_local = X.add_SU(N_local)
        Q = X.add_bifundamental(G_local, G_global, 'Q')
        Qtilde = X.add_bifundamental(G_global, G_local, 'Qt')
        phi = X.add_adjoint_vector(G_local, 'phi')
        psiA = X.add_adjoint_hyper(G_local, 'psiA')
        psiB = X.add_adjoint_hyper(G_local, 'psiB')
        X.set_superpotential(
            (Qtilde.matrix() * phi.matrix() * Q.matrix()).trace() + 
            (psiA.matrix() * phi.matrix() * psiB.matrix() -
             psiB.matrix() * phi.matrix() * psiA.matrix()   ).trace())
        return X

    def SingleInstanton(self, N):
        """
        EXAMPLES::

            sage: X = gauge_theories.SingleInstanton(2)
        """
        X = GaugeTheory()
        G_global = X.add_SU(N, external=True)
        G_local = X.add_U(1)
        Q = X.add_bifundamental(G_local, G_global, 'Q')
        Qtilde = X.add_bifundamental(G_global, G_local, 'Qt')
        phi = X.add_adjoint_vector(G_local, 'phi')
        psiA = X.add_adjoint_hyper(G_local, 'psiA')
        psiB = X.add_adjoint_hyper(G_local, 'psiB')
        X.set_superpotential(
            (Qtilde.matrix() * phi.matrix() * Q.matrix()).trace() + 
            (psiA.matrix() * phi.matrix() * psiB.matrix() -
             psiB.matrix() * phi.matrix() * psiA.matrix()   ).trace())
        return X

    def FlavorNequals2(self, N_c, N_f):
        """
        EXAMPLES::

            sage: X = gauge_theories.FlavorNequals2(3, 6)
        """
        X = GaugeTheory()
        G_global = X.add_U(N_f, external=True)
        G_local = X.add_SU(N_c)
        Q = X.add_bifundamental(G_local, G_global, 'Q')
        Qtilde = X.add_bifundamental(G_global, G_local, 'Qt')
        phi = X.add_adjoint_vector(G_local, 'phi')
        X.set_superpotential(
            (Qtilde.matrix() * phi.matrix() * Q.matrix()).trace())
        return X
        

    def E6theory(self):
        """
        N=2 theory defined by the N=2 quiver

        SU(3)_global --q1-- SU(3) --q2-- SU(3) --q3-- SU(3)_global
        """
        X = GaugeTheory()
        G1 = X.add_SU(3, external=True)
        G2 = X.add_SU(3)
        G3 = X.add_SU(3)
        G4 = X.add_SU(3, external=True)
        phi2 = X.add_adjoint_vector(G2, 'phi2')
        phi3 = X.add_adjoint_vector(G3, 'phi3')
        Q1tilde = X.add_bifundamental(G1, G2, 'Q1t')
        Q1      = X.add_bifundamental(G2, G1, 'Q1')
        Q2tilde = X.add_bifundamental(G2, G3, 'Q2t')
        Q2      = X.add_bifundamental(G3, G2, 'Q2')
        Q3tilde = X.add_bifundamental(G3, G4, 'Q3t')
        Q3      = X.add_bifundamental(G4, G3, 'Q3')
        X.set_superpotential(
            (Q1tilde.matrix() * phi2.matrix() * Q1.matrix()).trace() + 
            (Q2.matrix() * phi2.matrix() * Q2tilde.matrix()).trace() + 
            (Q2tilde.matrix() * phi3.matrix() * Q2.matrix()).trace() + 
            (Q3.matrix() * phi3.matrix() * Q3tilde.matrix()).trace() )
        return X

    def threeSU(self, n):
        """
        N=2 theory defined by the N=2 quiver

        SU(n)_global --q1-- SU(n) --q2-- SU(n)_global
        """
        X = GaugeTheory()
        G1 = X.add_SU(n, external=True)
        G3 = X.add_SU(n, external=True)
        G2 = X.add_SU(n)
        Q1tilde = X.add_bifundamental(G1, G2, 'Q1t')
        Q1 = X.add_bifundamental(G2, G1, 'Q1')
        Q2tilde = X.add_bifundamental(G2, G3, 'Q2t')
        Q2 = X.add_bifundamental(G3, G2, 'Q2')
        phi = X.add_adjoint_vector(G2, 'phi')
        X.set_superpotential(
            (Q1tilde.matrix() * phi.matrix() * Q1.matrix()).trace() + 
            (Q2.matrix() * phi.matrix() * Q2tilde.matrix()).trace() )
        return X
        

    def threeU(self, n):
        """
        N=2 theory defined by the N=2 quiver

        U(n)_global --q1-- U(n) --q2-- U(n)_global
        """
        X = GaugeTheory()
        G1 = X.add_U(n, external=True)
        G3 = X.add_U(n, external=True)
        G2 = X.add_U(n)
        Q1tilde = X.add_bifundamental(G1, G2, 'Q1t')
        Q1 = X.add_bifundamental(G2, G1, 'Q1')
        Q2tilde = X.add_bifundamental(G2, G3, 'Q2t')
        Q2 = X.add_bifundamental(G3, G2, 'Q2')
        phi = X.add_adjoint_vector(G2, 'phi')
        X.set_superpotential(
            (Q1tilde.matrix() * phi.matrix() * Q1.matrix()).trace() + 
            (Q2.matrix() * phi.matrix() * Q2tilde.matrix()).trace() )
        return X
        


gauge_theories = GaugeTheoryLibrary()



X = gauge_theories.threeSU(3)
#E6 = WeylCharacterRing('E6', style='coroots')




#X = gauge_theories.InstantonNequals2(3,2)
#print X.Higgs_branch().hilbert_series().factor()
#print ascii_art([Q.matrix(), Qtilde.matrix(), phi.matrix(), psiA.matrix(), psiB.matrix()])
#print X.Fideal(phi).hilbert_series()

#for i in []:# range(5):
#    X = gauge_theories.InstantonNequals2(i+1,2)
#    print i+1, X.Higgs_branch().hilbert_series().factor()




