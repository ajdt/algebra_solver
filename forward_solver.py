import pdb  # used for debugging
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import random # for random eqn generation


# symbol for sympy polynomials 
x_symb = sp.symbols('x')

# Utility Functions
def simplifyPolyTerms(term_list, default, constructor):
	if len(term_list) == 0 :
		return default
	elif len(term_list) == 1 :
		return term_list[0]
	else:
		return constructor(term_list)
def mergeLists(list_of_lists):
	return reduce(lambda x, y: x + y, list_of_lists, [])

class CorePoly(object):
	"""basic polynomial, mostly contains code to be overridden"""
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################
	def __init__(self):
		self.is_zero = False
		self.is_linear = self.degree() < 2
		self.is_one = False
    # these operations don't do any simplification
	def addition(self, poly) : 
		if isinstance(poly, SumPoly):
			return poly.addition(self)
		else:
			return SumPoly([self, poly])
	def subtract(self, poly) : return self.addition(poly.neg())
	def mult(self, poly) : # TODO: apply same fix to additon/subtraction etc.
		if isinstance(poly, ProdPoly):
			return poly.mult(self)
		else:
			return ProdPoly([self, poly])
	def divide(self, poly) : return RatPoly(num=self, denom=poly)
	def neg(self): raise NotImplementedError

	# UTILITY FUNCTIONS
	def degree(self) : raise NotImplementedError
	def copy(self): raise NotImplementedError
	def __str__(self): raise NotImplementedError
	def __eq__(self, other): raise NotImplementedError
	def __ne__(self, other): raise NotImplementedError

	# BOOLEANS
	def hasFractions(self): raise NotImplementedError
	def isFactored(self): return False
	def isConstTerm(self): return self.degree() == 0 # NOTE: this works for all  poly types

	############################## TERMS AND FACTORS ##############################
	def isConstantTerm(self): return False

class SumPoly(CorePoly):
	def __init__(self, poly_list):
		self.subpoly = poly_list
		CorePoly.__init__(self)

	# OVERRIDE OPERATION POLYNOMIALS
	def addition(self, poly) :
		if isinstance(poly, SumPoly):
			return SumPoly(self.subpoly + poly.subpoly)
		else:
			return SumPoly(self.subpoly + [poly])
	def subtract(self, poly) : return self.addition(poly.neg())
	def neg(self): return SumPoly([p.neg() for p in self.subpoly])

	# HELPERS
	def degree(self) : return max([p.degree() for p in self.subpoly])
	def copy(self): return SumPoly( [p.copy() for p in self.subpoly] )
	def __str__(self): return  ' + '.join([str(p) for p in self.subpoly])

	def __eq__(self, other):
		# same class types and same number of terms
		if not (self.__class__ == other.__class__ and len(self.subpoly) == len(other.subpoly)):
			return False
		# every subpoly occurs in 'other' as well
		return all([ p in other.subpoly for p in self.subpoly ])

	def __ne__(self, other): return not (self == other)

	# BOOLEANS
	def hasFractions(self): return any( [isinstance(p, RatPoly) for p in self.subpoly] )
	def isFactored(self): return all([p.isFactored() for p in self.subpoly])
	def isSameTerm(self, other): return self == other

	# MISC HELPERS
	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
	def coeffOf(self, base): # TODO: assumes only one monomial of each base type exists
		# only returns the first subpoly of that base
		for p in self.subpoly:
			if isinstance(p, StdPoly) and p.degree() == base:
				return p.coef()
		return None

	############################## terms and factors ##############################
	# TODO: implement all of these
	def __add__(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		return SumPoly(self.subpoly + other.subpoly).sumCommonTerms()
	def isConstantTerm(self): return self.degree() == 0
	def getNonConstTerms(self): return mergeLists([ p.getNonConstTerms() for p in self.subpoly if isinstance(p, StdPoly)])
	def getConstTerms(self): return mergeLists([ p.getConstTerms() for p in self.subpoly if isinstance(p, StdPoly)])
	def distribute(multiplier): return SumPoly( [p.mult(multiplier) for p in self.subpoly] )
	def hasCommonTerms(self):
		for idx, poly in enumerate(self.subpoly):
			for jidx, other in enumerate(self.subpoly):
				if idx == jidx:
					continue
				if poly.isSameTerm(other):
					return True
		return False

	def sumCommonTerms(self):
		ls = self.subpoly
		condensed = [] # common terms after being added together
		while len(ls) > 0:
			# in each pass, grab all subpoly that have sameTerms and add them
			result = reduce(lambda x,y: x + y, filter(lambda x: ls[0].isSameTerm(x), ls) )
			condensed.append(result)
			ls = filter(lambda x: not ls[0].isSameTerm(x), ls)

		# return new poly
		return simplifyPolyTerms(condensed, StdPoly.zero(), SumPoly)

class ProdPoly(CorePoly):
	def __init__(self, poly_list):
		self.subpoly = poly_list
		CorePoly.__init__(self)

	# OVERRIDE MATH OPERATIONS
	def mult(self,poly):
		if isinstance(poly, ProdPoly):
			return ProdPoly(self.subpoly + poly.subpoly)
		elif isinstance(poly, RatPoly):
			return RatPoly(self.mult(poly.num), poly.denom)
		else:
			return ProdPoly(self.subpoly + [poly])

	# TODO: should we make copies here?
	def neg(self): return ProdPoly( [self.subpoly[0].neg()] + self.subpoly[1:] )

	# UTILITY FUNCTIONS
	def degree(self) : return sum([p.degree() for p in self.subpoly])
	def copy(self): return ProdPoly([p.copy() for p in self.subpoly])
	def __str__(self): return ' * '.join(['('+str(p) +')' for p in self.subpoly])
	def __eq__(self, other):
		# same class types and same number of terms
		if not (self.__class__ == other.__class__ and len(self.subpoly) == len(other.subpoly)):
			return False
		# every subpoly occurs in 'other' as well
		return all([ p in other.subpoly for p in self.subpoly ])

	def __ne__(self, other): return not (self == other)

	# BOOLEANS
	def hasFractions(self): any([ isinstance(p, RatPoly) for p in self.subpoly])
	def isFactored(self): return all([p.isFactored() for p in self.subpoly])
	def isSameTerm(self, other): return self == other
	def haveCommonFactors(self, other): # used by RatPoly.numDenomHaveCommonFactors
		if not isinstance(other, ProdPoly):
			return False
		for p in self.subpoly:
			if p in other.subpoly:
				return True
		return False


	# MISC HELPERS
	@staticmethod
	def foil(poly1, poly2):
		""" multiply every term between both polynomials together """
		# TODO: input validation, what if given product polys or rational polys?
		terms1 =  poly1.subpoly if isinstance(poly1, SumPoly)  else [poly1]
		terms2 =  poly2.subpoly if isinstance(poly2, SumPoly)  else [poly2]

		# multiply each pair of terms
		new_terms = []
		for p in terms1:
			for q in terms2:
				if isinstance(p, StdPoly) and isinstance(p, StdPoly):
					new_terms.append(p * q)
				else:
					new_terms.append(p.mult(q))

		# return sumpoly as a result
		return simplifyPolyTerms(new_terms, StdPoly.zero(), SumPoly)

	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
	def getNonConstTerms(self): return [p for p in self.subpoly if not p.isConstTerm()]
	def getConstTerms(self): return [p for p in self.subpoly if p.isConstTerm()]
	############################## TERMS AND FACTORS ##############################
	def __add__(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		if not self.isSameTerm(other):
			return self.addition(other)
		else:
			return ProdPoly(self.subpoly + [StdPoly(2, x_symb)])
	def commonFactors(self): raise NotImplementedError
	def commonFactors(self, other): raise NotImplementedError
	def isConstantTerm(self): return False
	def cancelFactors(self, factors): # used by RatPoly
		new_subpoly = [ p for p in self.subpoly if p not in factors]
		return simplifyPolyTerms(new_subpoly, StdPoly.one(), ProdPoly)

class RatPoly(CorePoly):
	def __init__(self, num, denom):
		self.num, self.denom = num, denom
		CorePoly.__init__(self)
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################
	def mult(self, other) :
		if isinstance(other, RatPoly):
			return RatPoly(self.num.mult(other.num), self.denom.mult(other.denom))
		elif isinstance(other, ProdPoly):
			return RatPoly(self.num.mult(other), self.denom)
		else:
			return ProdPoly([self, other])
	def divide(self, other): return RatPoly(self.num, self.denom.mult(other))
	def neg(self): return RatPoly(self.num.neg(), self.denom)

	# UTILITY FUNCTIONS
	def degree(self) : return max([self.num.degree(), self.denom.degree()])
	def copy(self) : return RatPoly(self.num.copy(), self.denom.copy())
	def __str__(self): return '(' + str(self.num) + ')/(' + str(self.denom) + ')'
	def __eq__(self, other):
		if not (self.__class__ == other.__class__) :
			return False
		return self.num == other.num and self.denom == other.denom
	def __ne__(self, other): return not (self == other)

	# BOOLEANS
	def hasFractions(self): return True # true if you have or are a fraction
	def isFactored(self): return self.num.isFactored() and self.denom.isFactored()
	def isSameTerm(self, other): return self.__class__ == other.__class__ and self.denom == other.denom


	# MISC HELPERS
	def getFractions(self): return [self]
	def getNonConstTerms(self): return [self]
	def getConstTerms(self): return [] # TODO: change this, in case of const fraction?
	############################## MISC OPERATIONS ##############################
	def __add__(self, other):
		if not self.isSameTerm(other):
			return self.addition(other)
		else:
			return RatPoly(self.num + other.num, self.denom.copy())
	def cancelNumDenom(self) : raise NotImplementedError# if num and denom have common terms cancel them, #TODO: maybe keep as a function?
	def multReduce(self, other): raise NotImplementedError# if other occurs in denom, then cancel otherwise just multiply
	def reciprocal(self): return RatPoly(self.denom.copy(), self.num.copy())
	def numDenomShareFactors(self):
		if isinstance(self.num, ProdPoly) and isinstance(self.denom, ProdPoly):
			return self.num.haveCommonFactors(self.denom)
		elif isinstance(self.num, ProdPoly) and self.denom in self.num.subpoly:
			return True
		elif isinstance(self.denom, ProdPoly) and self.num in self.denom.subpoly:
			return True
		elif self.num == self.denom:
			return True
		return False

	def cancelCommonFactors(self): # used by simp5
		if isinstance(self.num, ProdPoly) and isinstance(self.denom, ProdPoly):
			common = [p for p in self.num.subpoly if p in self.denom.subpoly]
			new_num = self.num.cancelFactors(common)
			new_denom = self.denom.cancelFactors(common)
			if new_denom.is_one:
				return new_num
			else:
				return RatPoly(new_num, new_denom) # handle more cases here
		elif isinstance(self.num, ProdPoly) and self.denom in self.num.subpoly:
			new_num = self.num.cancelFactors([self.denom])
			return new_num # no longer a fraction
		elif isinstance(self.denom, ProdPoly) and self.num in self.denom.subpoly:
			new_denom = self.denom.cancelFactors([self.num])
			return RatPoly(StdPoly.one(), new_denom)
		elif self.num == self.denom:
			return StdPoly.one()
		return self

class Bases: # used as an enum
	CONST, X, X2, X3 = range(4)

class StdPoly(sp.Poly, CorePoly):
	def __init__(self, *args):
		sp.Poly.__init__(self, *args)

	# UTILITY FUNCTIONS
	def coef(self): return self.coeffs()[0] # TODO: remove this, temporary while refactoring takes place
	def __str__(self): return str(self.as_expr()).replace('x_symb', 'x') # don't like how sp.Poly prints polynomials
	def copy(self): return StdPoly(self.as_expr(), x_symb) # TODO: only way to have constants make copies of themselves since I can't pass generator through Poly.copy() method

	# BOOLEANS
	def hasFractions(self): return False
	def isFactored(self): return self.degree() < 2
	def isConstTerm(self): return self.degree() == 0
	def isSameTerm(self, other): 
		if not self.__class__ == other.__class__:
			return False
		my_coeffs, other_coeffs = list(reversed(self.all_coeffs())), list(reversed(other.all_coeffs()))
		same_terms = lambda x, y: x is not None and y is not None and x != 0 and y != 0
		return  any(map(same_terms, my_coeffs, other_coeffs))

	def coeffOf(self, base): return self.coeffs()[0] if base == self.degree() else None # TODO: revise this to return self.coeff_monomial(x_symb**base), causes bugs to crop up

	# MISC HELPERS
	@staticmethod
	def zero(): return StdPoly(0, x_symb)
	@staticmethod
	def one(): return StdPoly(1, x_symb)
	def getFractions(self): return [] # TODO: try to remove this
	def getNonConstTerms(self): return [ StdPoly((self - self.all_coeffs()[-1]).as_expr() ) ] if self.degree() > 0 else []
	def getConstTerms(self): return [ StdPoly(self.all_coeffs()[-1], x_symb) ]
	############################## TERMS AND FACTORS ##############################
class Eqn:

	# random generation globals
	rand = random.Random()
	ADD, MUL, FRAK, STD = range(4)
	CONST, LIN, QUAD, CUBIC = range(4)
	node_types = [(ADD, 4), (MUL, 2), (FRAK,3), (STD,16)]
	std_types = [(CONST,4), (LIN,8), (QUAD,10), (CUBIC,1)]

	# lists to select from
	node_list, std_list = [], []
	for node, weight in node_types:
		node_list += [node]*weight
	for std, weight in std_types:
		std_list += [std]*weight

	def __init__(self, eqn_string):
		sides = eqn_string.split('=') # use x_symb variable, and split into two sides
		l, r  = sp.sympify(sides[0], evaluate=False), sp.sympify(sides[1], evaluate=False)
		self.left, self.right = Eqn.convertToPolyTree(l), Eqn.convertToPolyTree(r)

	def __str__(self):
		return str(self.left) + '=' + str(self.right)
	def copy(self): 
		new_eqn = Eqn('3*x = 2')
		new_eqn.left  = self.left.copy()
		new_eqn.right = self.right.copy()
		return new_eqn
	def degree(self): return max([self.left.degree(), self.right.degree()])

	@staticmethod
	def strToPolyTree(string):
		return  Eqn.convertToPolyTree( sp.sympify(string, evaluate=False) ) 
	@staticmethod
	def convertToPolyTree(sympoly):
		""" convert a sympy.add.Add or sympy.mul.Mul to the tree structure used in rest of program
		NOTE: when we refactor to use sympy Add and Mul classes, this will no longer be necessary
		"""
		# TODO: revisit to avoid automatic simplification
		if sympoly.is_polynomial() and all([sp.Poly(p, x_symb).is_monomial for p in sympoly.args]):
			return StdPoly(sympoly, x_symb) 
			#TODO: this is inefficient, figure out way to avoid polynomial evaluation!
			#return simplifyPolyTerms([StdPoly(p, x_symb) for p in sympoly.args], StdPoly.zero(), SumPoly)
		elif sympoly.is_Mul : 
			if any([p.is_Pow and (-1 in p.args) for p in sympoly.args ]): # rational poly
				denom_factors = [Eqn.convertToPolyTree(p.args[0]) for p in sympoly.args if p.is_Pow and (-1 in p.args) ]
				num_factors = [Eqn.convertToPolyTree(p) for p in sympoly.args if not (p.is_Pow and (-1 in p.args)) ]

				denom_poly = simplifyPolyTerms(denom_factors, StdPoly.one(), ProdPoly)
				num_poly = simplifyPolyTerms(num_factors, StdPoly.one(), ProdPoly)
				return RatPoly(num_poly, denom_poly)
			else :
				return ProdPoly( [ Eqn.convertToPolyTree(p) for p in sympoly.args ] )
			return RatPoly(Eqn.convertToPolyTree(sympoly.args[1]), Eqn.convertToPolyTree(sympoly.args[0].args[0]) )
		elif sympoly.is_Add :
			return SumPoly( [ Eqn.convertToPolyTree(p) for p in sympoly.args ] )
		elif sympoly.is_Pow :
			return RatPoly(StdPoly.one(), Eqn.convertToPolyTree(sympoly.args[0]))
		else :
			raise TypeError

	@staticmethod
	def genRandTree():
		node_type = random.sample(Eqn.node_list,1)[0]
		if node_type == Eqn.ADD:
			num_terms = Eqn.rand.randint(2,3)
			terms = [Eqn.genRandTree() for i in range(num_terms)]
			return SumPoly(terms)
		elif node_type == Eqn.MUL:
			num_terms = Eqn.rand.randint(2,2)
			terms = [Eqn.genRandTree() for i in range(num_terms)]
			return ProdPoly(terms)
		elif node_type == Eqn.FRAK:
			num = Eqn.genRandTree()
			denom = Eqn.genRandTree()
			return RatPoly(num,denom)
		elif node_type == Eqn.STD:
			degree = random.sample(Eqn.std_list, 1)[0]
			coeff = [random.randint(0,20) for i in range(4)]
			coeff[degree] = random.randint(1,20)
			# ensure the coeff of higher degree monomials is zero
			for i in range(degree+1, len(coeff)):
				coeff[i] = 0
			d,c,b,a = coeff # coeffs are listed in ascending order
			return StdPoly(a*x_symb**3 + b*x_symb**2 + c*x_symb + d, x_symb)
		return StdPoly.one()

	@staticmethod
	def genRandEqn():
		eqn = Eqn('x = 3')
		eqn.right = Eqn.genRandTree()
		eqn.left = Eqn.genRandTree()
		return eqn

class WorkingMem:
	SET_RHS_ZERO = 'set rhs zero'
	def __init__(self, steps_in_latex=False):
		self.steps_in_latex = steps_in_latex
		self.backtrack_stack = []
		self.steps = []
		self.goals = []
	def copy(self):
		new_wm = WorkingMem(self.steps_in_latex)
		new_wm.backtrack_stack = list(self.backtrack_stack)
		new_wm.steps = list(self.steps)
		new_wm.goals = list(self.goals)
		return new_wm

	def hasGoal(self, goal): return goal in self.goals
	def addGoal(self, goal): self.goals.append(goal)
	def removeGoal(self, goal): self.goals.remove(goal)
	def addStep(self, eqn, rule) : 
		if not self.steps_in_latex:
			self.steps.append(str(eqn) + ': ' + rule.__name__)
		else :
			# TODO: this will be simplified when we switch completely to Mult and Add classes in sympy
			sympy_eqn = sp.Eq( sp.sympify(str(eqn.left), evaluate=False), sp.sympify(str(eqn.right), evaluate=False))
			self.steps.append('$$' + sp.latex(sympy_eqn) + '\t\t: \\text{ ' + rule.__name__ + '}$$')
	def btpeek(self): return self.backtrack_stack[-1]
	def btpop(self) : return self.backtrack_stack.pop()
	def btpush(self, eqn) : return self.backtrack_stack.append(eqn)

class EqnRule(object):
	"""A rule that applies at the level of a single equation, and uses working mem"""
	def __init__(self, cond, action, desc="", name=""):
		self._cond, self._action, self._desc, self.__name__ = cond, action, desc, name
	def checkCondition(self, eqn, working_mem):
		""" @return: true if condition applies"""
		return self._cond(eqn, working_mem)
	def applyAction(self, eqn, working_mem):
		""" Apply given action to eqn. NOTE: Assumed eqn is mutable; nothing is returned."""
		self._action(eqn, working_mem)

class PolyRule(EqnRule):
	"""A rule that applies to polynomials. Accepts a single equation and traverses the eqn tree."""
	def __init__(self, *args): EqnRule.__init__(self, *args)
	def checkCondition(self, eqn, working_mem):
		return self._checkPoly(eqn.left) or self._checkPoly(eqn.right)
	def applyAction(self, eqn, working_mem):
		""" apply the rules action to given equation, and return the equation"""
		eqn.left, changed = self._applyActionRecursive(eqn.left)
		if not changed:
			eqn.right, changed = self._applyActionRecursive(eqn.right)
		return eqn

	def _checkPoly(self, poly):
		""" recursively check if rule applies to polynomial"""
		if self._cond(poly):
			return True
		else:
			if isinstance(poly, StdPoly):
				return False
			elif isinstance(poly, RatPoly):
				return self._checkPoly(poly.num) or self._checkPoly(poly.denom)
			else : # SumPoly or ProdPoly
				return any( [self._checkPoly(p) for p in poly.subpoly] )

	def _applyActionRecursive(self, poly):
		""" Find a polynomial for which rule applies, apply rule. Boolean indicates rule was applied
		@return: (new_polynomial, True/False) 
		"""
		if self._cond(poly):
			return (self._action(poly), True)
		else:
			if isinstance(poly, StdPoly):
				return (poly, False)
			elif isinstance(poly, RatPoly):
				new_num, changed = self._applyActionRecursive(poly.num)
				if changed:
					return (RatPoly(new_num, poly.denom), changed)
				else:
					new_denom, changed = self._applyActionRecursive(poly.denom)
					return (RatPoly(poly.num, new_denom), changed)
			else : # SumPoly or ProdPoly
				ls = [self._applyActionRecursive(p) for p in poly.subpoly]
				terms, bools = map(lambda x: x[0], ls), map(lambda x: x[1], ls)
				result = SumPoly(terms) if isinstance(poly, SumPoly) else ProdPoly(terms)
				return (result, any(bools))
		
class RuleHelper:
	""" contains mostly static methods to help rules operate on polys"""
	@staticmethod
	def computeLCM( poly_list):
		"""
		compute the lcm of a list of fractions
		@return: list of polynomials in the lcm
		"""
		# TODO: overly complex, simplify
		terms = []
		for p in poly_list:
			if isinstance(p, ProdPoly):
				for subpoly in p.subpoly:
					num_p = len([x for x in p.subpoly if x == subpoly])
					num_terms = len([x for x in terms if x == subpoly])
					# if subpoly occurs in p more times than already
					# accounted for, then include it until we have
					# enough multiples
					if num_p > num_terms:
						for i in range(num_p - num_terms):
							terms.append(subpoly)
			elif p not in terms:
				terms.append(p)
		return terms


	@staticmethod
	def completeSquare(std_poly):
		if not isinstance(std_poly, StdPoly):
			raise TypeError
		# TODO: for now assumes a = 1, to avoid fractions
		d = std_poly.coeff_monomial(x_symb)**2/4 # take coeff attached to x term
		c = std_poly.coeff_monomial(x_symb**0) # coef of constant term
		poly = (std_poly - c + d).factor()
		if isinstance(poly, sp.Pow): # factoring was successful
			factor = StdPoly(poly.args[0])
			return (SumPoly([ProdPoly([factor, factor.copy()]), StdPoly(c - d, x_symb)]), True)
		else:
			return (std_poly, False)

	@staticmethod
	def factorCubic(std_poly):
		"""
		factors a sumpoly that is in standard form
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		if not isinstance(std_poly, StdPoly):
			raise TypeError

		poly = std_poly.factor()
		if isinstance(poly, sp.Mul): # factoring was successful
			return (ProdPoly([StdPoly(p) for p in poly.args]), True)
		else:
			return (std_poly, False)

	@staticmethod
	def factor(std_poly):
		"""
		factors a standard poly
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		# TODO: remove boolean return value
		if not isinstance(std_poly, StdPoly):
			raise TypeError

		poly = std_poly.factor()
		if isinstance(poly, sp.Mul): # factoring was successful
			return (ProdPoly([ StdPoly(p, x_symb) for p in poly.args ]), True)
		else:
			return (std_poly, False)

	@staticmethod
	def _removeZeroes(sum_poly):
		no_zeroes = [p for p in sum_poly.subpoly if not p.is_zero ]
		return simplifyPolyTerms(no_zeroes, StdPoly.zero(), SumPoly)
	@staticmethod
	def removeFactorsFrom(factor, remove_from):
		"""
		return a new list with the given polynomial removed
		"""
		new_list = list(remove_from)
		if isinstance(factor, ProdPoly):
			for p in factor.subpoly:
				new_list.remove(p)
		else :
			new_list.remove(factor)
		return new_list

	@staticmethod
	def mult4Helper(sum_poly):
		lcm = Solver.computeLCM( [ i.denom for i in sum_poly.getFractions() ])
		ls = []
		# multiply all fractions to get a common denominator
		for poly in sum_poly.subpoly:
			if isinstance(poly, RatPoly): # compute the correct multiplier to get common denominator
				multiplier = simplifyPolyTerms( Solver.removeFactorsFrom(poly.denom, lcm), StdPoly.one(), ProdPoly )
				mult_frac = RatPoly(multiplier, multiplier.copy())
				ls.append(poly.mult(mult_frac))
			else:
				ls.append(poly)
		return SumPoly(ls)
	@staticmethod
	def mult5Helper(prod_poly):
		new_terms = [ProdPoly.foil(prod_poly.subpoly[0], prod_poly.subpoly[1]) ] + prod_poly.subpoly[2:]
		return simplifyPolyTerms(new_terms, StdPoly.zero(), ProdPoly)
	@staticmethod
	def getRHSNonConstLHSConst(eqn):
		right, left = eqn.right, eqn.left
		to_remove = right.getNonConstTerms() + left.getConstTerms()
		to_remove = [p for p in to_remove if not p.is_zero] # remove all zero terms
		return to_remove
	@staticmethod
	def moveConstRHSNonConstLHS(eqn):
		to_remove = RuleHelper.getRHSNonConstLHSConst(eqn)
		remove_poly =  simplifyPolyTerms(to_remove, StdPoly.zero(), SumPoly)
		# subtract the terms from both sides
		eqn.left, eqn.right  = eqn.left.subtract(remove_poly), eqn.right.subtract(remove_poly)
	@staticmethod
	def subtractFromEqn(eqn, poly):
		eqn.right = eqn.right.subtract(poly)
		eqn.left = eqn.left.subtract(poly)
	@staticmethod
	def setLHSToNum(eqn):
		eqn.left = eqn.left.num
	@staticmethod
	def simp8Helper(eqn):
		divisor = eqn.left.coeffOf(Bases.X)
		eqn.left = StdPoly(x_symb)
		eqn.right = eqn.right.divide(StdPoly(divisor,x_symb))
	@staticmethod
	def mult2Helper(eqn):
		""" if both sides of eqn have fractions, then multiply each side by the lcm over all fractions.  """
		# TODO: remove right.is_zero condition, after fixing the condition that requires moving everything to lhs
		# TODO: remove the code for checking if rhs is zero, and adjusting likewise
		# TODO: SIMPLIFY THIS CODE!!
		# TODO: make atomic, right now does distributive and multiplicative step.
		# get list of denom from both sides
		left_denom = [i.denom for i in eqn.left.getFractions()]
		right_denom = [] if eqn.right.is_zero else [i.denom for i in eqn.right.getFractions()]
		# compute lcm and multiply
		lcm = simplifyPolyTerms(RuleHelper.computeLCM(left_denom + right_denom), StdPoly.one(), ProdPoly)
		left = SumPoly([p.mult(lcm) for p in eqn.left.subpoly]) if isinstance(eqn.left, SumPoly) else eqn.left.mult(lcm)
		eqn.left = left
		eqn.right = eqn.right if eqn.right.is_zero else eqn.right.mult(lcm)

############################## RULES ##############################	
# condition, action, description, name
SIMP0 =	PolyRule(	lambda x : isinstance(x, SumPoly) and any([p.is_zero for p in x.subpoly]), 
										RuleHelper._removeZeroes,
										""" simp0: if zeroes exist as additive terms, then remove them """,
										'simp0'
					)
SIMP1 =	EqnRule(	lambda eq, wm : eq.degree() >= 2 and not wm.hasGoal(WorkingMem.SET_RHS_ZERO),
					lambda eq, wm : wm.addGoal(WorkingMem.SET_RHS_ZERO),
					""" simp1: if degree is >= 2, then set working mem goal to make rhs zero """,
					'simp1'
					)
SIMP2 =	PolyRule(	lambda x : isinstance(x, SumPoly) and SumPoly.hasCommonTerms(x), 
										SumPoly.sumCommonTerms,
										""" simp2: if sumpoly has common terms, then add them together """,
										'simp2'
					)
SIMP3 =	EqnRule(	lambda eq, wm : not wm.hasGoal(WorkingMem.SET_RHS_ZERO) and not eq.right.isConstTerm() and len(RuleHelper.getRHSNonConstLHSConst(eq)) > 0,
					lambda eq, wm : RuleHelper.moveConstRHSNonConstLHS(eq),
					""" if solving a linear eqn, cancel all constant terms on the lhs and all non-constant terms on the rhs """,
					'simp3'
					)
SIMP4 =	EqnRule(	lambda eq, wm : wm.hasGoal(WorkingMem.SET_RHS_ZERO) and not eq.right.is_zero,
					lambda eq, wm : RuleHelper.subtractFromEqn(eq, eq.right),
					""" if our goal is to set rhs to zero, then subtract all rhs terms from lhs""",
					'simp4'
					)
SIMP5 =	PolyRule(	lambda p: isinstance(p, RatPoly) and RatPoly.numDenomShareFactors(p), 
										RatPoly.cancelCommonFactors,
										""" simp5: if num and denom of a rational polynomial have common factors, then cancel these factors """,
										'simp5'
					)
SIMP6 =	PolyRule(	lambda x : isinstance(x, RatPoly) and x.num.is_zero , 
										lambda x : StdPoly.zero(),
										""" simp6: if zero exists in a numerator, remove the fraction involved """,
										'simp6'
					)
SIMP7 =	EqnRule(	lambda eq, wm : isinstance(eq.left, RatPoly) and eq.right.is_zero, 
										lambda eq,wm : RuleHelper.setLHSToNum(eq),
										""" if lhs is a rational polynomial, and rhs is zero, solve for numerator """,
										'simp7'
					)
SIMP8 =	EqnRule(	lambda eq, wm : eq.right.isConstTerm() and eq.left.is_linear and eq.left.coeffOf(Bases.X) != 1 and eq.left.coeffOf(Bases.CONST) is None, 
										lambda eq,wm : RuleHelper.simp8Helper(eq),
										""" if equation has form ax = b, divide by a """,
										'simp8'
					)
SIMP9 =	EqnRule(	lambda eq, wm : eq.degree() < 2 and wm.hasGoal(WorkingMem.SET_RHS_ZERO), 
										lambda eq,wm : wm.removeGoal(WorkingMem.SET_RHS_ZERO),
										""" if SET_RHS_ZERO is a goal and we've reduced problem to linear eqn, then remove this goal""",
										'simp9'
					)
MULT1 =	PolyRule(	lambda p: isinstance(p, RatPoly) and isinstance(p.denom, RatPoly), 
										lambda p: ProdPoly([p.num, p.denom.reciprocal()]),
										""" if denom of rational poly is a fraction, the multiply by its reciprocal """,
										'mult1'
					)
MULT2 =	EqnRule(	lambda eq, wm : eq.left.hasFractions() and  (eq.right.hasFractions() or eq.right.is_zero), 
										lambda eq,wm : RuleHelper.mult2Helper(eq),
										""" if both sides of eqn have fractions, then multiply each side by the lcm over all fractions.  """,
										'mult2'
					)
MULT4 =	PolyRule( lambda p: isinstance(p, SumPoly) and p.hasFractions() and len(p.getFractions()) > 1,
										RuleHelper.mult4Helper,
										""" if a polynomial is a sum over rational polynomials, then multiply every polynomial by lcm/lcm""",
										'mult4'
									)
MULT5 =	PolyRule(lambda p: isinstance(p, ProdPoly),
									RuleHelper.mult5Helper,
									""" if a there is a product polynomial, then foil the first two factors""",
									'mult5'
									)

HEUR1 =	PolyRule(lambda p:  isinstance(p, StdPoly) and p.degree() == 2 and isinstance(p.factor(x_symb), sp.Mul) # TODO: look for is_factorable() method
									,lambda p : RuleHelper.factor(p)[0]
									,""" if a 2nd degree polynomial occurs anywhere, then attempt to factor it """,
									'heur1'
									)

HEUR2 =	PolyRule(lambda p: isinstance(p, StdPoly) and p.degree() == 2 and RuleHelper.completeSquare(p)[1]
									,lambda p : RuleHelper.completeSquare(p)[0]
									,""" if a 2nd degree polynomial occurs anywhere, then factor it by completing the square """,
									'heur2'
									)

HEUR3 =	PolyRule(lambda p: p.degree() == 3 and isinstance(p, StdPoly) and RuleHelper.factorCubic(p)[1]
									,lambda p : RuleHelper.factorCubic(p)[0]
									,""" if a 3rd degree polynomial occurs anywhere, then attempt to factor it """,
									'heur3'
				)
class Solver:
	def __init__(self, eqn, rule_ord=lambda rule: Solver.RULE_ORDERING.index(rule)): 
		self.eqn, self.working_mem	= eqn, WorkingMem()
		self.rule_ordering			= rule_ord

	def __str__(self): return '\n'.join(self.working_mem.steps)
	def copy(self):
		new_solver = Solver(self.eqn.copy())
		new_solver.working_mem = self.working_mem.copy()
		return new_solver
	############################## win conditions  ##############################
	def win1(self): # """ case a = b"""
		""" contradiction:  reduced to a = b, but a =/= b """
		right, left = self.eqn.right, self.eqn.left
		return (left.isConstTerm() and right.isConstTerm() and left != right)
	def win2(self):
		""" win condition: ax = b problem is solved """
		right, left = self.eqn.right, self.eqn.left
		return left.is_linear and right.isConstTerm() and left.coeffOf(Bases.X) == 1 and left.coeffOf(Bases.CONST) is None
	def win3(self):
		""" win condition: lhs is completely factored, rhs is zero """
		right, left = self.eqn.right, self.eqn.left
		# TODO: revisit this rule, it's gotten complex
		return ( self.eqn.degree() >= 2 and left.isFactored() and right.is_zero and not isinstance(left, SumPoly) and not isinstance(left, RatPoly))

	############################## checking and solving ##############################
	## list of rules and precedences they take
	SIMP_RULES		= [SIMP6,SIMP7, SIMP8, SIMP0, SIMP1, SIMP2, SIMP3, SIMP4, SIMP5, SIMP9 ]
	WIN_RULES		= [win1, win2, win3]
	MULT_RULES		= [MULT1, MULT2, MULT4, MULT5]
	MISC_RULES		= []
	HEURISTICS		= [HEUR1, HEUR2, HEUR3]
	ALL_RULES 		= [SIMP_RULES, HEURISTICS, MULT_RULES, MISC_RULES ]
	RULE_ORDERING 	= SIMP_RULES + HEURISTICS + MULT_RULES + MISC_RULES 

	## solve the problem
	def solve(self):
		"""solve the equation given, return steps to the solution"""
		self.working_mem.steps.append(str(self.eqn) + ': ' + self.solve.__name__)
		while not self.checkWinCond():
			triggered_rules = self.getTriggeredRules()
			if len(triggered_rules) == 0 : # stuck, there's no more todo
				break

			# select a rule and apply it
			applied_rule = self.selectRule(triggered_rules)
			applied_rule.applyAction(self.eqn, self.working_mem)
			self.working_mem.addStep(self.eqn, applied_rule) 
		return self.working_mem.steps

	def checkWinCond(self):
		""" check win conditions """
		# if any win condition is active, then we've solved the equation
		for rule in self.WIN_RULES:
			if rule(self) :
				#print "win condition " + rule.__name__ + " applies"
				self.working_mem.addStep( self.eqn, rule)
				return True
		return False

	def getTriggeredRules(self):
		"""@return: list of all rules triggered by current game state"""
		return [ rule for rule in Solver.RULE_ORDERING if rule.checkCondition(self.eqn, self.working_mem) ]

	def selectRule(self, triggered_rules): return min(triggered_rules, key=self.rule_ordering)

class SuperSolver():
	def __init__(self, eqn_str):
		self.eqn = Eqn(eqn_str)

	def allSolns(self):
		solutions = []
		solvers = [Solver(self.eqn)]
		pdb.set_trace()
		while len(solvers) > 0:
			# take the top solver
			soln = solvers.pop()
			# if finished solving, add it's solution
			if soln.checkWinCond():
				solutions.append(soln.working_mem.steps)
				continue
			# ... otherwise generate a solver for each triggered rule and apply the rule
			triggered_rules = soln.getTriggeredRules()
			if len(triggered_rules) == 0: # deadend, no solution down this path
				continue
			for rule in triggered_rules:
				new_solver = soln.copy()
				rule.applyAction(new_solver.eqn, new_solver.working_mem)
				new_solver.working_mem.addStep(new_solver.eqn, rule)
				solvers.append(new_solver)
		return solutions

# Notes:
#	all/any for reducing a list of booleans

# certifier for trace on a solution
# enumerate all valid traces
# test traces
# borrow ideas from drools  jess
# quick check
# use model checker?
# cmbc
# strict hierarchy, or limit to tree depth
#
# python to run specific tests: python -m unittest test_fs.SampleCasesTest.test_solve2
