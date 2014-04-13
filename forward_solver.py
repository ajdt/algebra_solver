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
	def isLinearStdPoly(self): return False  # override for sum and monomial classes
	def isConstTerm(self): return self.order() == 0 # NOTE: this works for all  poly types
	def shareFactors(self, other): return False

	############################## TERMS AND FACTORS ##############################	
	def isConstantTerm(self): return False

class SumPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE OPERATION POLYNOMIALS
	def add(self, poly) : 
		self.subpoly.append(poly) # TODO: these don't return a new polynomial!
		return self
	def sub(self, poly) : 
		self.subpoly.append(poly.negate()) 
		return self
	def negate(self, poly): return SumPoly([p.negate() for p in self.subpoly])

	# IMPLEMENT HELPERS
	def order(self) : return max([p.order() for p in self.subpoly])
	def copy(self): return SumPoly( [p.copy() for p in self.subpoly] )
	def __str__(self): ' + '.join([str(p) for p in self.subpoly])

	def __eq__(self, other): 
		# same class types and same number of terms
		if not (self.__class__ == other.__class__ and len(self.subpoly) == len(other.subpoly)):
			return False
		# every subpoly occurs in 'other' as well
		return all([ p in other.subpoly for p in self.subpoly ])

	def __ne__(self, other): return not (self == other)

	# BOOLEANS 
	def hasFractions(self): return any( [isinstance(p, RatPoly) for p in self.subpoly] )
	def isLinearStdPoly(self): return all([isinstance(p, Monomial) for p in self.subpoly]) and self.order() < 2
	def isFactored(self): self.order() < 2
	def isSameTerm(self, other): return self == other
	def shareFactors(self, other): return self == other

	# MISC HELPERS
	def coeffOf(self, base): 
		# only returns the first subpoly of that base
		for p in self.subpoly:
			if p.base == base:
				return p.coef
		return None

	############################## terms and factors ##############################	
	# TODO: implement all of these
	def sumSameTerm(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		return sumCommonTerms(SumPoly(self.subpoly + other.subpoly))
	def isConstantTerm(self): return self.order() == 0
	def getNonConstTerms(self): return [ p for p in self.subpoly if not p.isConstTerm() ] 
	def getConstantTerms(self): return [ p for p in self.subpoly if p.isConstTerm() ] 
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
		if len(condensed) > 1:
			return SumPoly(condensed)
		else:
			return condensed[0]

class ProdPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE MATH OPERATIONS
	def mult(self,poly): 
		self.subpoly.append(poly) 
		return self
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


	############################## TERMS AND FACTORS ##############################	
	def sumSameTerm(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return ProdPoly(self.subpoly + [Monomial(coef=2, base=Bases.CONST)])
	def foil(self): raise NotImplementedError #return result of multiplying terms together # foil specific terms?
	def commonFactors(self): raise NotImplementedError
	def commonFactors(self, other): raise NotImplementedError
	def isConstantTerm(self): return False
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors

class RatPoly(CorePoly):
	def __init__(self, num, denom): self.num, self.denom = num, denom
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################	
	def divide(self, other): 
		self.num = self.num.mult(other) 
		return self
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
	def hasFractions(self): return False 
	def isFactored(self): return max([self.num.order(), self.denom.order()]) < 2
	def isZero(self): return self.num.isZero()
	def isLinearStdPoly(self): return False  # override for sum and monomial classes
	def isSameTerm(self, other): return self.denom() == other.denom()


	############################## MISC OPERATIONS ##############################	
	def sumSameTerm(self, other):
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return RatPoly(self.num.add(other.num), self.denom.copy())
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors
	def cancelNumDenom(self) : raise NotImplementedError# if num and denom have common terms cancel them, #TODO: maybe keep as a function?
	def multReduce(self, other): raise NotImplementedError# if other occurs in denom, then cancel otherwise just multiply
	def cancelCommonFactors(self): raise NotImplementedError# look @ num and denom and cancel any common factors they have
	def reciprocal(self): return RatPoly(self.denom.deepCopy(), self.num.deepCopy())

class Bases: # used as an enum
	CONST, X, X2, X3 = range(4)

class Monomial(CorePoly):
	def __init__(self, coef, base): self.coef, self.base = coef, base
	
	# Overriden operations
	def negate(self): return Monomial(self.coef*-1, self.base)

	# UTILITY FUNCTIONS
	def order(self) : return self.base
	def copy(self): return Monomial(self.coef, self.base)
	def __str__(self): 
		if self.coef == 0:
			return 0
		elif self.base == Bases.CONST:
			return str(self.coef)
		elif self.base == Bases.X:
			return str(self.coef) +'x'
		else:
			return str(self.coef) + 'x' + str(self.base)
	def __eq__(self, other): 
		if not (self.__class__ == other.__class__) :
			return False
		return self.coef == other.coef and self.base == other.base
	def __ne__(self, other): return not (self == other)

	# BOOLEANS 
	def hasFractions(self): return False
	def isFactored(self): return self.order() < 2 
	def isZero(self): return self.coef == 0
	def isLinearStdPoly(self): return True  
	def isConstTerm(self): return self.base == Bases.CONST
	def isSameTerm(self, other): return self.base == other.base
	def coeffOf(self, base): return self.coef if base == self.base else None

	############################## TERMS AND FACTORS ##############################	
	def sumSameTerm(self, other):
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return Monomial(self.coef + other.coef, self.base)
	def shareFactors(self, other): raise NotImplementedError
class Eqn:
	def __init__(self, left, right):
		self.left, self.right = left, right
	def __str__(self):
		return str(self.left) + '=' + str(self.right)
	def copy(self): return Eqn(self.left.copy(), self.right.copy())
	def order(self): return max([self.left.order(), self.right.order()])

class WorkingMem:
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
		return all([left.isLinearStdPoly(),right.isConstTerm(),left.coeffOf(Bases.X) != 1])
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
			if isinstance(poly, BasicPoly):
				return (poly, False)
			elif isinstance(poly, RatPoly):
				new_num, changed = func(poly.num)
				if changed:
					return (RatPoly(new_num, poly.denom), changed)
				else:
					new_denom, changed = fun(poly.denom)
					return (RatPoly(poly.num, new_denom), changed)
			else : # SumPoly or ProdPoly
				ls = [self.polyRule(p, condition, action) for p in poly.subpoly]
				terms = map(lambda x: x[0], ls)
				bools = map(lambda x: x[1], ls)
				result = SumPoly(terms) if isinstance(poly, SumPoly) else ProdPoly(terms)
				return (result, any(bools))

	def simp1(self): return False
	def simp2(self): return False
	def simp3(self): return False
	def simp4(self): return False
	def simp5(self): return False
	def simp6(self): return False

	############################## checking and solving ##############################	
	## list of rules and precedences they take
	SIMP_RULES		= [simp1, simp2, simp3, simp4, simp5, simp6]
	WIN_RULES		= [win1, win2, win3]
	MULT_RULES		= []
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
				print "no rules apply"
				return
		#print str(self.eqn)

	def checkWinCond(self):
		""" check win conditions """
		# if any win condition is active, then we've solved the equation
		for rule in self.WIN_RULES:
			if rule(self) :
				return True
		return False

	def checkRuleSet(self, ruleset):
		""" check an argument ruleset"""
		for rule in ruleset:
			if rule(self) :
				return True
		return False

# testing code		
left = SumPoly([Monomial(10, Bases.X), Monomial(3, Bases.CONST), Monomial(-3, Bases.CONST)])
right = Monomial(10, Bases.CONST)

solver = Solver(Eqn(left, right))
solver.solve()
# Notes:
#	all/any for reducing a list of booleans
