import pdb # used for debugging

# sympy functionality
import sympy as sp
x_symb = sp.symbols('x')
poly = sp.Poly(3*x_symb)

# Utility Functions
def simplifyPolyTerms(term_list, default, constructor):
	if len(term_list) == 0 :
		return default
	elif len(term_list) == 1 :
		return term_list[0]
	else:
		return constructor(term_list)

class CorePoly(object):
	"""basic polynomial, mostly contains code to be overridden"""
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################
	def add(self, poly) : return SumPoly([self, poly])
	def sub(self, poly) : return SumPoly([self, poly.negate() ])
	def mult(self, poly) : return ProdPoly([self, poly])
	def divide(self, poly) : return RatPoly(num=self, denom=poly)
	def negate(self): raise NotImplementedError

	# UTILITY FUNCTIONS
	def order(self) : raise NotImplementedError
	def copy(self): raise NotImplementedError
	def __str__(self): raise NotImplementedError
	def __eq__(self, other): raise NotImplementedError
	def __ne__(self, other): raise NotImplementedError

	# BOOLEANS
	def hasFractions(self): raise NotImplementedError
	def isFactored(self): return False
	def isZero(self): return False
	def isOne(self): return False
	def isLinearStdPoly(self): return False  # override for sum and monomial classes
	def isConstTerm(self): return self.order() == 0 # NOTE: this works for all  poly types

	############################## TERMS AND FACTORS ##############################
	def isConstantTerm(self): return False

class SumPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE OPERATION POLYNOMIALS
	def add(self, poly) :
		if isinstance(poly, SumPoly):
			return SumPoly(self.subpoly + poly.subpoly)
		else:
			return SumPoly(self.subpoly + [poly])
	def sub(self, poly) : return self.add(poly.negate())
	def negate(self): return SumPoly([p.negate() for p in self.subpoly])

	# HELPERS
	def order(self) : return max([p.order() for p in self.subpoly])
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
	def isLinearStdPoly(self): return all([isinstance(p, StdPoly) for p in self.subpoly]) and self.order() < 2
	def isFactored(self): self.order() < 2
	def isSameTerm(self, other): return self == other

	# MISC HELPERS
	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
	def coeffOf(self, base): # TODO: assumes only one monomial of each base type exists
		# only returns the first subpoly of that base
		for p in self.subpoly:
			if isinstance(p, StdPoly) and p.order() == base:
				return p.coef()
		return None

	############################## terms and factors ##############################
	# TODO: implement all of these
	def sumSameTerm(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		return sumCommonTerms(SumPoly(self.subpoly + other.subpoly))
	def isConstantTerm(self): return self.order() == 0
	def getNonConstTerms(self): return [ p for p in self.subpoly if not p.isConstTerm() ]
	def getConstTerms(self): return [ p for p in self.subpoly if p.isConstTerm() ]
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
			result = reduce(lambda x,y: x.sumSameTerm(y), filter(lambda x: ls[0].isSameTerm(x), ls) )
			condensed.append(result)
			ls = filter(lambda x: not ls[0].isSameTerm(x), ls)

		# return new poly
		return simplifyPolyTerms(condensed, StdPoly.zero(), SumPoly)

class ProdPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE MATH OPERATIONS
	def mult(self,poly):
		if isinstance(poly, ProdPoly):
			return ProdPoly(self.subpoly + poly.subpoly)
		else:
			return ProdPoly(self.subpoly + [poly])

	# TODO: should we make copies here?
	def negate(self): return ProdPoly( [self.subpoly[0].negate()] + self.subpoly[1:] )

	# UTILITY FUNCTIONS
	def order(self) : return sum([p.order() for p in self.subpoly])
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
	def isFactored(self): return max([p.order() for p in self.subpoly]) <= 1
	def isLinearStdPoly(self): return self.order < 2
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
					new_terms.append(StdPoly( (p.poly * q.poly).as_expr() ))
				else:
					new_terms.append(p.mult(q))

		# return sumpoly as a result
		return simplifyPolyTerms(new_terms, StdPoly.zero(), SumPoly)

	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
	def getNonConstTerms(self): return [p for p in self.subpoly if not p.isConstTerm()]
	def getConstTerms(self): return [p for p in self.subpoly if p.isConstTerm()]
	############################## TERMS AND FACTORS ##############################
	def sumSameTerm(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return ProdPoly(self.subpoly + [StdPoly(2)])
	def commonFactors(self): raise NotImplementedError
	def commonFactors(self, other): raise NotImplementedError
	def isConstantTerm(self): return False
	def cancelFactors(self, factors): # used by RatPoly
		new_subpoly = [ p for p in self.subpoly if p not in factors]
		return simplifyPolyTerms(new_subpoly, StdPoly.one(), ProdPoly)

class RatPoly(CorePoly):
	def __init__(self, num, denom): self.num, self.denom = num, denom
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################
	def mult(self, other) :
		if isinstance(other, RatPoly):
			return RatPoly(self.num.mult(other.num), self.denom.mult(other.denom))
		else:
			return ProdPoly([self, other])
	def divide(self, other): return RatPoly(self.num, self.denom.mult(other))
	def negate(self): return RatPoly(self.num.negate(), self.denom)

	# UTILITY FUNCTIONS
	def order(self) : return max([self.num.order(), self.denom.order()])
	def copy(self) : return RatPoly(self.num.copy(), self.denom.copy())
	def __str__(self): return '(' + str(self.num) + ')/(' + str(self.denom) + ')'
	def __eq__(self, other):
		if not (self.__class__ == other.__class__) :
			return False
		return self.num == other.num and self.denom == other.denom
	def __ne__(self, other): return not (self == other)

	# BOOLEANS
	def hasFractions(self): return True # true if you have or are a fraction
	def isFactored(self): return max([self.num.order(), self.denom.order()]) < 2
	def isZero(self): return self.num.isZero()
	def isLinearStdPoly(self): return False  # override for sum and monomial classes
	def isSameTerm(self, other): return self.__class__ == other.__class__ and self.denom == other.denom


	# MISC HELPERS
	def getFractions(self): return [self]
	def getNonConstTerms(self): return [self]
	def getConstTerms(self): return [] # TODO: change this, in case of const fraction?
	############################## MISC OPERATIONS ##############################
	def sumSameTerm(self, other):
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return RatPoly(self.num.add(other.num), self.denom.copy())
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
			if new_denom.isOne():
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


class StdPoly(CorePoly):
	def __init__(self, expr):
		"""
		@param expr: a sympy expression
		"""
		if isinstance(expr, int) or isinstance(expr, sp.Integer):
			self.poly = sp.Poly(expr + x_symb) - sp.Poly(x_symb) # TODO: slight work around since I can't initialize a constant monomial easily
		else:
			self.poly = sp.Poly(expr)

	# Overriden operations
	def negate(self):
		return StdPoly( (self.poly*-1).as_expr() ) # TODO: revise when changed to standard poly

	# UTILITY FUNCTIONS
	def coef(self): return self.poly.coeffs()[0] # TODO: remove this, temporary while refactoring takes place
	def order(self) : return self.poly.degree()
	def copy(self):
		return StdPoly(self.poly.as_expr()) # TODO: revise when changed to standard poly
	def __str__(self): return str(self.poly.as_expr())
	def __eq__(self, other): return (self.__class__ == other.__class__) and self.poly == other.poly
	def __ne__(self, other): return not (self == other)

	# BOOLEANS
	def hasFractions(self): return False
	def isFactored(self): return self.poly.degree() < 2
	def isZero(self): return self.poly.is_zero
	def isOne(self): return self.poly.coeffs()[0] == 1 and self.poly.degree() == 0
	def isLinearStdPoly(self): return True
	def isConstTerm(self): return self.order() == 0
	def isSameTerm(self, other): return self.__class__ == other.__class__ and self.poly.degree() == other.poly.degree()
	def coeffOf(self, base): return self.poly.coeffs()[0] if base == self.poly.degree() else None

	# MISC HELPERS
	@staticmethod
	def zero(): return StdPoly(0)
	@staticmethod
	def one(): return StdPoly(1)
	def getFractions(self): return []
	def getNonConstTerms(self): return [self] if self.order() > 0  else []
	def getConstTerms(self): return [self] if self.order() == 0  else []
	############################## TERMS AND FACTORS ##############################
	def sumSameTerm(self, other):
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return StdPoly( (self.poly + other.poly).as_expr() )

class Eqn:
	def __init__(self, left, right):
		self.left, self.right = left, right
	def __str__(self):
		return str(self.left) + '=' + str(self.right)
	def copy(self): return Eqn(self.left.copy(), self.right.copy())
	def order(self): return max([self.left.order(), self.right.order()])

class WorkingMem:
	SET_RHS_ZERO = 'set rhs zero'
	def __init__(self):
		self.backtrack_stack = []
		self.steps = []
		self.goals = []
	def hasGoal(self, goal): return goal in self.goals
	def addGoal(self, goal): self.goals.append(goal)
	def removeGoal(self, goal): self.goals.remove(goal)
	def addStep(self, step) : self.steps.append(step)
	def btpeek(self): return self.backtrack_stack[-1]
	def btpop(self) : return self.backtrack_stack.pop()
	def btpush(self, eqn) : return self.backtrack_stack.append(eqn)

class Solver:
	def __init__(self, eqn): self.eqn, self.working_mem	 = eqn, WorkingMem()

	############################## win conditions  ##############################
	def win1(self): # """ case a = b"""
		right, left = self.eqn.right, self.eqn.left
		return (left.isConstTerm() and right.isConstTerm() and left != right)
	def win2(self):
		""" case ax = b"""
		# TODO: break this up into another rule
		right, left = self.eqn.right, self.eqn.left
		return left.isLinearStdPoly() and right.isConstTerm() and left.coeffOf(Bases.X) != 1 and left.coeffOf(Bases.CONST) is None
	def win3(self):
		right, left = self.eqn.right, self.eqn.left
		return ( self.eqn.order() >= 2 and left.isFactored() and right.isZero() )

	############################## rules ##############################
	# structure of recursive rules:
	def polyRule(self, poly, condition, action):
		# will return new polynomial and T/F depending on whether action was performed
		if condition(poly):
			return (action(poly), True)
		else:
			if isinstance(poly, StdPoly):
				return (poly, False)
			elif isinstance(poly, RatPoly):
				new_num, changed = self.polyRule(poly.num, condition, action)
				if changed:
					return (RatPoly(new_num, poly.denom), changed)
				else:
					new_denom, changed = self.polyRule(poly.denom, condition, action)
					return (RatPoly(poly.num, new_denom), changed)
			else : # SumPoly or ProdPoly
				ls = [self.polyRule(p, condition, action) for p in poly.subpoly]
				terms, bools = map(lambda x: x[0], ls), map(lambda x: x[1], ls)
				result = SumPoly(terms) if isinstance(poly, SumPoly) else ProdPoly(terms)
				return (result, any(bools))

	def checkEqnForRule(self, cond, action):
		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed

	@staticmethod
	def _removeZeroes(sum_poly):
		no_zeroes = [p for p in sum_poly.subpoly if not p.isZero() ]
		return simplifyPolyTerms(no_zeroes, StdPoly.zero(), SumPoly)

	def simp0(self):
		""" if sumpoly has zeroes remove them """
		cond	= lambda x : isinstance(x, SumPoly) and any([p.isZero() for p in x.subpoly])
		action	= Solver._removeZeroes
		return self.checkEqnForRule(cond, action)

	def simp1(self):
		""" set working mem goal, if order is >= 2 """
		if self.eqn.order() >= 2 and not self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO):
			self.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)
			return True
		return False

	def simp2(self):
		""" if sum poly has common terms, then add them together """
		cond	= lambda x : isinstance(x, SumPoly) and SumPoly.hasCommonTerms(x)
		action	= SumPoly.sumCommonTerms
		return self.checkEqnForRule(cond, action)

	def simp3(self):
		""" if solving linear eqn, move everything except constant terms to lhs """
		left, right = self.eqn.left, self.eqn.right
		# TODO: may have to fix when this rule fires, what if we have x+3 = 2, can't proceed
		if not self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO) and not right.isConstTerm():
			to_remove = right.getNonConstTerms() + left.getConstTerms()
			if len(to_remove) == 0:
				return False
			else:
				remove_poly =  simplifyPolyTerms(to_remove, StdPoly.zero(), SumPoly)
			# subtract the terms from both sides
			self.eqn.left, self.eqn.right  = self.eqn.left.sub(remove_poly), self.eqn.right.sub(remove_poly)
			return True
		return False

	def simp4(self):
		""" if equation is higher than 1st degree set rhs to zero"""
		if self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO) and not self.eqn.right.isZero():
			self.eqn.left = self.eqn.left.sub(self.eqn.right)
			self.eqn.right = self.eqn.right.sub(self.eqn.right)
			return True
		return False

	def simp5(self):
		""" if num and denom of a rational polynomial have common factors, remove them"""
		cond = lambda p: isinstance(p, RatPoly) and RatPoly.numDenomShareFactors(p)
		action = RatPoly.cancelCommonFactors
		return self.checkEqnForRule(cond, action)

	def mult1(self):
		""" if denom of rational poly is a fraction, multiply by its reciprocal """
		cond	= lambda p: isinstance(p, RatPoly) and isinstance(p.denom, RatPoly)
		action	= lambda p: ProdPoly([p.num, p.denom.reciprocal()])
		return self.checkEqnForRule(cond, action)

	def mult2(self):
		""" if both sides of eqn have fractions, then multiply each side by lcm
			over all fractions.
		"""
		if self.eqn.left.hasFractions() and  self.eqn.right.hasFractions():
			# get list of denom from both sides
			left_denom = [i.denom for i in self.eqn.left.getFractions()]
			right_denom = [i.denom for i in self.eqn.right.getFractions()]
			# compute lcm and multiply
			lcm = Solver.computeLCM(left_denom + right_denom)
			self.eqn.right = self.eqn.right.mult(lcm)
			self.eqn.left = self.eqn.left.mult(lcm)
			return True
		return False

	@staticmethod
	def mult4Helper(sum_poly):
		lcm = Solver.computeLCM( [ i.denom for i in sum_poly.getFractions() ])
		rp = RatPoly(lcm, lcm.copy())
		ls = []
		for i in range(len(p.subpoly)):
			if isinstance(p.subpoly[i], RatPoly):
				ls.append(p.subpoly[i].mult(rp))
			else:
				ls.append(p.subpoly[i])
		return SumPoly(ls)

	def mult4(self):
		cond = lambda p: isinstance(p, SumPoly) and p.hasFractions() and len(p.getFractions()) > 1
		action = Solver.mult4Helper
		return self.checkEqnForRule(cond, action)

	@staticmethod
	def mult5Helper(prod_poly):
		new_terms = [ProdPoly.foil(prod_poly.subpoly[0], prod_poly.subpoly[1]) ] + prod_poly.subpoly[2:]
		return simplifyPolyTerms(new_terms, StdPoly.zero(), ProdPoly)

	def mult5(self):
		cond = lambda p: isinstance(p, ProdPoly)
		action = Solver.mult5Helper
		return self.checkEqnForRule(cond, action)

	def heur1(self):
		""" attempt to factor a polynomial of degree  == 2 """
		cond = lambda p: p.order() == 2 and isinstance(p, SumPoly) and Solver.factor(p)[1]
		action = lambda p : Solver.factor(p)[0]
		return self.checkEqnForRule(cond, action)

	def heur2(self):
		""" factor by completing the square """
		cond = lambda p: p.order() == 2 and isinstance(p, SumPoly) and Solver.completeSquare(p)[1]
		action = lambda p : Solver.completeSquare(p)[0]
		return self.checkEqnForRule(cond, action)

	def heur3(self):
		""" attempt to factor a polynomial of degree  == 3 """
		cond = lambda p: p.order() == 3 and isinstance(p, SumPoly) and Solver.factorCubic(p)[1]
		action = lambda p : Solver.factorCubic(p)[0]
		return self.checkEqnForRule(cond, action)

	@staticmethod
	def computeLCM( poly_list):
		""" compute the lcm of a list of fractions"""
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
		return simplifyPolyTerms(terms, StdPoly.one(), ProdPoly)


	@staticmethod
	def completeSquare(sum_poly):
		if not isinstance(sum_poly, SumPoly):
			raise TypeError
		a, b, c, = sum_poly.coeffOf(Bases.X2), sum_poly.coeffOf(Bases.X), sum_poly.coeffOf(Bases.CONST)
		# TODO: for now assumes a = 1, to avoid fractions
		d = b**2/4

		poly = sp.Poly(a*x_symb**2 + b*x_symb + d).factor()
		if isinstance(poly, sp.Pow): # factoring was successful
			left = SumPoly([StdPoly(poly.args[0].coeff(x_symb)**x_symb), StdPoly(sp.Poly(poly.args[0]).coeffs()[-1])])
			right = left.copy()
			return (SumPoly([ProdPoly([left, right]), StdPoly(c - d)]), True)
		else:
			return (sum_poly, False)



	@staticmethod
	def factorCubic(sum_poly):
		"""
		factors a sumpoly that is in standard form
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		if not isinstance(sum_poly, SumPoly):
			raise TypeError
		# TODO: switch to using sympy in the future
		a, b, c, d = sum_poly.coeffOf(Bases.X3), sum_poly.coeffOf(Bases.X2), sum_poly.coeffOf(Bases.X), sum_poly.coeffOf(Bases.CONST)

		poly = sp.Poly(a*x_symb**3 + b*x_symb**2 + c*x_symb + d).factor()
		if isinstance(poly, sp.Mul): # factoring was successful
			if len(poly.args) == 3: # factored to linear terms
				left = SumPoly([StdPoly(poly.args[0].coeff(x_symb)*x_symb), StdPoly(sp.Poly(poly.args[0]).coeffs()[-1])])
				middle = SumPoly([StdPoly(poly.args[1].coeff(x_symb)*x_symb), StdPoly(sp.Poly(poly.args[1]).coeffs()[-1])])
				right = SumPoly([StdPoly(poly.args[2].coeff(x_symb)*x_symb), StdPoly(sp.Poly(poly.args[2]).coeffs()[-1])])
				return (ProdPoly([left, middle, right]), True)

			else: # factored to one linear and one square term
				linear, quad = poly.args if poly.args[0].coeff(x_symb**2) == 0 else tuple(reversed(poly.args))
				left = SumPoly([StdPoly(linear.coeff(x_symb)*x_symb), StdPoly(sp.Poly(linear).coeffs()[-1])])
				# TODO: Uuuuugly, rework the line below. It looks awful
				right = Solver._removeZeroes(SumPoly([StdPoly(quad.coeff(x_symb**2)*x_symb**2), StdPoly(quad.coeff(x_symb)*x_symb), StdPoly(sp.Poly(quad).coeffs()[-1])]) )
			return (ProdPoly([left, right]), True)
		else:
			return (sum_poly, False)

	@staticmethod
	def factor(sum_poly):
		"""
		factors a sumpoly that is in standard form
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		if not isinstance(sum_poly, SumPoly):
			raise TypeError
		# TODO: switch to using sympy in the future
		a, b, c, = sum_poly.coeffOf(Bases.X2), sum_poly.coeffOf(Bases.X), sum_poly.coeffOf(Bases.CONST)

		poly = sp.Poly(a*x_symb**2 + b*x_symb + c).factor()
		if isinstance(poly, sp.Mul): # factoring was successful
			left = SumPoly([StdPoly(poly.args[0].coeff(x_symb)*x_symb), StdPoly(sp.Poly(poly.args[0]).coeffs()[-1])])
			right = SumPoly([StdPoly(poly.args[1].coeff(x_symb)*x_symb), StdPoly(sp.Poly(poly.args[1]).coeffs()[-1])])
			return (ProdPoly([left, right]), True)
		else:
			return (sum_poly, False)

	############################## checking and solving ##############################
	## list of rules and precedences they take
	SIMP_RULES		= [simp0, simp1, simp2, simp3, simp4, simp5]
	WIN_RULES		= [win1, win2, win3]
	MULT_RULES		= [mult1, mult2]
	MISC_RULES		= []
	HEURISTICS		= []

	## solve the problem
	def solve(self):
		"""solve the equation given"""
		while not self.checkWinCond():
			#print str(self.eqn)
			if self.checkRuleSet(self.SIMP_RULES):
				continue
			elif self.checkRuleSet(self.MULT_RULES):
				continue
			elif self.checkRuleSet(self.MISC_RULES):
				continue
			elif self.checkRuleSet(self.HEURISTICS):
				continue
			else:
				# TODO: change this later
				#print "no rules apply"
				#print str(self.eqn)
				return str(self.eqn)
		#print str(self.eqn)
		return str(self.eqn)

	def checkWinCond(self):
		""" check win conditions """
		# if any win condition is active, then we've solved the equation
		for rule in self.WIN_RULES:
			if rule(self) :
				#print "win condition " + rule.__name__ + " applies"
				return True
		return False

	def checkRuleSet(self, ruleset):
		""" check an argument ruleset"""
		for rule in ruleset:
			if rule(self) :
				#print 'applied: ' + rule.__name__
				return True
		return False

# Notes:
#	all/any for reducing a list of booleans
