import pdb # used for debugging
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
	def shareFactors(self, other): return False

	############################## TERMS AND FACTORS ##############################	
	def isConstantTerm(self): return False

class SumPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE OPERATION POLYNOMIALS
	def add(self, poly) : 
		if isinstance(poly, SumPoly):
			return SumPoly(self.subpoly + poly.subpoly)
		else:
			self.subpoly.append(poly) # TODO: these don't return a new polynomial!
			return self
	def sub(self, poly) : return self.add(poly.negate())
	def negate(self): return SumPoly([p.negate() for p in self.subpoly])

	# IMPLEMENT HELPERS
	def order(self) : return max([p.order() for p in self.subpoly])
	def copy(self): return SumPoly( [p.copy() for p in self.subpoly] )
	def __str__(self): 
		return  ' + '.join([str(p) for p in self.subpoly])

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
	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
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
	def haveCommonFactors(self, other): # used by RatPoly.numDenomHaveCommonFactors
		if not isinstance(other, ProdPoly): 
			return False
		for p in self.subpoly:
			if p in other.subpoly:
				return True
		return False


	# MISC HELPERS
	def getFractions(self): return [ p for p in self.subpoly if isinstance(p, RatPoly) ]
	def getNonConstTerms(self): return [p for p in self.subpoly if not p.isConstTerm()]
	def getConstTerms(self): return [p for p in self.subpoly if p.isConstTerm()]
	############################## TERMS AND FACTORS ##############################	
	def sumSameTerm(self, other): # TODO: may need to fix this, for when we want to add 3(x+2)(x-3) + (x +2)(x-3)
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return ProdPoly(self.subpoly + [Monomial(coef=2, base=Bases.CONST)])
	def foil(self): raise NotImplementedError
		# does two actions, foil and simplify foiled terms
	def commonFactors(self): raise NotImplementedError
	def commonFactors(self, other): raise NotImplementedError
	def isConstantTerm(self): return False
	def shareFactors(self, other): return any([p in other.subpoly for p in self.subpoly])
	def cancelFactors(self, factors): # used by RatPoly
		new_subpoly = [ p for p in self.subpoly if p not in factors]
		if len(new_subpoly) == 0:
			return Monomial(1, Bases.CONST)
		elif len(new_subpoly) == 1:
			return new_subpoly[0]
		else :
			return ProdPoly(new_subpoly)

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
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors
	def cancelNumDenom(self) : raise NotImplementedError# if num and denom have common terms cancel them, #TODO: maybe keep as a function?
	def multReduce(self, other): raise NotImplementedError# if other occurs in denom, then cancel otherwise just multiply
	def reciprocal(self): return RatPoly(self.denom.deepCopy(), self.num.deepCopy())
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
			return RatPoly(Monomial(1, Bases.CONST), new_denom)
		elif self.num == self.denom:
			return Monomial(1, Bases.CONST)
		return self

		

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
		if self.base == Bases.CONST:
			return str(self.coef)
		elif self.base == Bases.X:
			if self.coef == 1:
				return 'x'
			else:
				return str(self.coef) +'x'
		else:
			if self.coef == 1:
				return 'x^' + str(self.base)
			else :
				return str(self.coef) + 'x^' + str(self.base)
	def __eq__(self, other): 
		if not (self.__class__ == other.__class__) :
			return False
		return self.coef == other.coef and self.base == other.base
	def __ne__(self, other): return not (self == other)

	# BOOLEANS 
	def hasFractions(self): return False
	def isFactored(self): return self.order() < 2 
	def isZero(self): return self.coef == 0
	def isOne(self): return self.coef == 1 and self.base == Bases.CONST
	def isLinearStdPoly(self): return True  
	def isConstTerm(self): return self.base == Bases.CONST
	def isSameTerm(self, other): return self.__class__ == other.__class__ and self.base == other.base
	def coeffOf(self, base): return self.coef if base == self.base else None

	# MISC HELPERS
	def getFractions(self): return []
	def getNonConstTerms(self): return [self] if self.order() > 0  else []
	def getConstTerms(self): return [self] if self.order() == 0  else []
	############################## TERMS AND FACTORS ##############################	
	def sumSameTerm(self, other):
		if not self.isSameTerm(other):
			return self.add(other)
		else:
			return Monomial(self.coef + other.coef, self.base)
	def shareFactors(self, other): return self == other
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
			if isinstance(poly, Monomial):
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
				terms = map(lambda x: x[0], ls)
				bools = map(lambda x: x[1], ls)
				result = SumPoly(terms) if isinstance(poly, SumPoly) else ProdPoly(terms)
				return (result, any(bools))

	# TODO: move these helper methods elsewhere
	@staticmethod
	def _removeZeroes(sum_poly):
		no_zeroes = [p for p in sum_poly.subpoly if not p.isZero() ]
		if len(no_zeroes) == 0:
			return Monomial(0, Bases.CONST)
		elif len(no_zeroes) == 1:
			return no_zeroes[0]
		else:
			return SumPoly(no_zeroes)

	def simp0(self): 
		""" if sumpoly has zeroes remove them """
		cond	= lambda x : isinstance(x, SumPoly) and any([p.isZero() for p in x.subpoly]) 
		action	= Solver._removeZeroes
		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed
	def simp1(self): 
		if self.eqn.order() > 2 and not self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO):
			self.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)
			return True
		return False
	def simp2(self): 
		cond	= lambda x : isinstance(x, SumPoly) and SumPoly.hasCommonTerms(x)
		action	= SumPoly.sumCommonTerms
		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed
	def simp3(self): 
		left, right = self.eqn.left, self.eqn.right
		# TODO: may have to fix when this rule fires, what if we have x+3 = 2, can't proceed
		if not self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO) and not right.isConstTerm():
			to_remove = right.getNonConstTerms() + left.getConstTerms()
			if len(to_remove) == 0:
				return False
			if len(to_remove) == 1:
				remove_poly = to_remove[0]
			else :
				remove_poly = SumPoly(to_remove)
			# subtract the terms from both sides
			self.eqn.left = self.eqn.left.sub(remove_poly)
			self.eqn.right = self.eqn.right.sub(remove_poly)
			return True
		return False

	def simp4(self): 
		if self.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO) and not self.eqn.right.isZero():
			self.eqn.left = self.eqn.left.sub(self.eqn.right)
			self.eqn.right = self.eqn.right.sub(self.eqn.right)
			return True
		return False

	def simp5(self): 
		cond = lambda p: isinstance(p, RatPoly) and RatPoly.numDenomShareFactors(p)
		action = RatPoly.cancelCommonFactors

		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed

	def mult1(self):
		cond	= lambda p: isinstance(p, RatPoly) and isinstance(p.denom, RatPoly)
		action	= lambda p: ProdPoly([p.num, p.denom.reciprocal()])
		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed

	def mult2(self):
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
		self.eqn.left, changed = self.polyRule(self.eqn.left, cond, action)
		if not changed:
			self.eqn.right, changed = self.polyRule(self.eqn.right, cond, action)
		return changed
	#def mult5(self):
	#	cond = lambda p: isinstance(p, ProdPoly)
	#	action = ProdPoly.foil()

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
		if len(terms) == 1:
			return terms[0]
		else :
			return ProdPoly(terms)


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
			print str(self.eqn)
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
				print str(self.eqn)
				return str(self.eqn)
		print str(self.eqn)
		return str(self.eqn)

	def checkWinCond(self):
		""" check win conditions """
		# if any win condition is active, then we've solved the equation
		for rule in self.WIN_RULES:
			if rule(self) :
				print "win condition " + rule.__name__ + " applies"
				return True
		return False

	def checkRuleSet(self, ruleset):
		""" check an argument ruleset"""
		for rule in ruleset:
			if rule(self) :
				print 'applied: ' + rule.__name__
				return True
		return False

# TODO: remove later
# testing code		
def testSolve1():
	left = SumPoly([Monomial(10, Bases.X), Monomial(3, Bases.CONST), Monomial(-3, Bases.CONST)])
	right = Monomial(10, Bases.CONST)

	solver = Solver(Eqn(left, right))
	solver.solve()
def testSolve2():
	left = SumPoly([Monomial(10, Bases.X), Monomial(3, Bases.CONST)])
	right = SumPoly([Monomial(10, Bases.CONST), Monomial(3, Bases.X)])

	solver = Solver(Eqn(left, right))
	solver.solve()
def testSolve3():
	left = SumPoly([Monomial(3, Bases.X), Monomial(1, Bases.X2)])
	right = SumPoly([Monomial(3, Bases.CONST), Monomial(1, Bases.X2)])

	solver = Solver(Eqn(left, right))
	solver.solve()
def testSolve4():
	sp1 = SumPoly([Monomial(1, Bases.X), Monomial(1, Bases.CONST)])
	sp2 = SumPoly([Monomial(1, Bases.X), Monomial(3, Bases.CONST)])

	num = ProdPoly([Monomial(1, Bases.X2), sp1, sp2])
	denom  = ProdPoly([sp1, sp2])

	left = SumPoly([RatPoly(num,denom), Monomial(3, Bases.X)])
	right = SumPoly([Monomial(3, Bases.CONST), Monomial(1, Bases.X2)])

	solver = Solver(Eqn(left, right))
	#pdb.set_trace()
	solver.solve()
	print "\n"

def testSolve5():
	sp1 = SumPoly([Monomial(1, Bases.X), Monomial(1, Bases.CONST)])
	sp2 = SumPoly([Monomial(1, Bases.X), Monomial(3, Bases.CONST)])

	left = RatPoly(Monomial(1, Bases.CONST), sp1)
	right = RatPoly(Monomial(1, Bases.CONST), sp2)
	solver = Solver(Eqn(left, right))
	pdb.set_trace()
	solver.solve()
def main():
	testSolve1()
	testSolve2()
	testSolve3()
	testSolve4()
	testSolve5()
#main()
# Notes:
#	all/any for reducing a list of booleans
