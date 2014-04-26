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

#class PolyUtils:
	#""" contains static methods that help add polys and do other such things """

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
		self.left, self.right  = sp.sympify(sides[0], evaluate=False), sp.sympify(sides[1], evaluate=False)

	def __str__(self):
		return str(self.left) + '=' + str(self.right)
	def copy(self): 
		new_eqn = Eqn('3*x = 2')
		new_eqn.left  = self.left.copy()
		new_eqn.right = self.right.copy()
		return new_eqn
	def degree(self): return max([sp.degree(self.left, gens=x_symb), sp.degree(self.right, gens=x_symb)])

	@staticmethod # TODO: unnecessary!
	def strToPolyTree(string):
		return  sp.sympify(string, evaluate=False)

	## TODO: change to generate rand string!
	#@staticmethod
	#def genRandTree():
		## TODO: rewrite to generate rand string, that is then parsed!
		#node_type = random.sample(Eqn.node_list,1)[0]
		#if node_type == Eqn.ADD:
			#num_terms = Eqn.rand.randint(2,3)
			#terms = [Eqn.genRandTree() for i in range(num_terms)]
			#return SumPoly(terms)
		#elif node_type == Eqn.MUL:
			#num_terms = Eqn.rand.randint(2,2)
			#terms = [Eqn.genRandTree() for i in range(num_terms)]
			#return ProdPoly(terms)
		#elif node_type == Eqn.FRAK:
			#num = Eqn.genRandTree()
			#denom = Eqn.genRandTree()
			#return RatPoly(num,denom)
		#elif node_type == Eqn.STD:
			#degree = random.sample(Eqn.std_list, 1)[0]
			#coeff = [random.randint(0,20) for i in range(4)]
			#coeff[degree] = random.randint(1,20)
			## ensure the coeff of higher degree monomials is zero
			#for i in range(degree+1, len(coeff)):
				#coeff[i] = 0
			#d,c,b,a = coeff # coeffs are listed in ascending order
			#return StdPoly(a*x_symb**3 + b*x_symb**2 + c*x_symb + d, x_symb)
		#return StdPoly.one()

	#@staticmethod
	#def genRandEqn():
		#eqn = Eqn('x = 3')
		#eqn.right = Eqn.genRandTree()
		#eqn.left = Eqn.genRandTree()
		#return eqn

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
			if not poly.is_polynomial(): # rational polynomial
				numer, denom = poly.as_numer_denom()
				return self._checkPoly(numer) or self._checkPoly(denom)
			elif poly.is_Add or poly.is_Mul:
				return any( [self._checkPoly(p) for p in poly.args] )
			else : 
				return False

	def _applyActionRecursive(self, poly):
		""" Find a polynomial for which rule applies, apply rule. Boolean indicates rule was applied
		@return: (new_polynomial, True/False) 
		"""
		# TODO: how to deal with sp.power.Pow polynomials??
		if self._cond(poly):
			return (self._action(poly), True)
		elif not poly.is_polynomial(): # rational poly
			numer,denom = poly.as_numer_denom()
			new_num, changed = self._applyActionRecursive(numer)
			if changed:
				return (new_num/denom, True)
			else :
				new_denom, changed = self._applyActionRecursive(denom)
				return (new_num/new_denom, changed)
		elif poly.is_Add :
			# recursively try every subpoly, generate a new Add poly from their result, and return true if any of the subpolys had the rule applied to it
			recursive_apply = [self._applyActionRecursive(p) for p in poly.args]
			new_polys, bools = [res[0] for res in recursive_apply], [res[1] for res in recursive_apply]
			return (sp.add.Add.fromiter(new_polys), any(bools))
		elif poly.is_Mul :
			# recursively try every subpoly, generate a new Add poly from their result, and return true if any of the subpolys had the rule applied to it
			new_polys, bools = map(list, zip(*sorted( self._applyActionRecursive(p) for p in poly.args )))
			return (sp.mul.Mul.fromiter(new_polys), any(bools))
		else:
			return (poly, False)
		
class RuleHelper:
	""" contains mostly static methods to help rules operate on polys"""
	@staticmethod
	def computeLCM( poly_list):
		"""
		compute the lcm of a list of fractions
		@return: list of polynomials in the lcm
		"""
		# TODO: dealing with finding the lcm of a list of 
		return reduce(sp.lcm, poly_list, 0*x_symb)

	@staticmethod
	def completeSquare(poly):
		if not poly.is_Add:
			raise TypeError
		# TODO: for now assumes a = 1, to avoid fractions
		x_coef, const_coef = poly.coeff(x_symb)**2/4, poly.coeff(x_symb**0)

		factored_poly = (std_poly - const_coef + x_coef).factor()
		if factored_poly.is_Pow: # factoring was successful
			return (factored_poly + (const_coef - x_coef), True)
		else:
			return (poly, False)

	@staticmethod
	def factor(std_poly):
		"""
		factors a standard poly
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		# TODO: remove boolean return value
		if not std_poly.is_Add:
			raise TypeError

		factored_poly = std_poly.factor()
		if isinstance(poly, sp.Mul): # factoring was successful
			return (factored_poly, True)
		else:
			return (std_poly, False)

	@staticmethod
	def _removeZeroes(sum_poly):
		return sp.add.Add.fromiter( filter( lambda p : not p.is_zero, sum_poly.args) )

	@staticmethod
	def mult4Helper(sum_poly):
		# find the lcm over all denominators
		denom_list = map(lambda (num,deno): deno, [p.as_numer_denom() for p in sum_poly.args] )
		lcm = RuleHelper.computeLCM(denom_list)
		ls = []
		# multiply every fraction, so that it's denom is the lcm
		for poly in sum_poly.args:
			if not poly.is_polynomial: # a fraction
				numer, denom = poly.as_numer_denom()
				multiplier = (lcm/numer).cancel() # compute the correct multiplier to get common denominator
				ls.append( sp.mul.Mul.fromiter([poly, multiplier/multiplier]) )
			else:
				ls.append(poly)
		return sp.add.Add.fromiter(ls)

	@staticmethod
	def mult5Helper(prod_poly):
		""" foil all terms in the product poly """
		if not prod_poly.is_Mul or not prod_poly.is_polynomial: #  ensure this is not a fraction also
			raise TypeError
		# multiply all polynoms together, then convert to an expression
		return reduce(__mul__, [sp.Poly(p) for p in prod_poly.args]).as_expr() 

	@staticmethod
	def getRHSNonConstLHSConst(eqn):
		""" @return: a polynomial of const terms on lhs and non_const terms on rhs or zero if there are no terms"""
		const_terms	= filter(lambda x: x.is_Number, eqn.left.as_ordered_terms()) 
		var_terms	= filter(lambda x: not x.is_Number, eqn.right.as_ordered_terms()) 
		return sp.add.Add.fromiter(const_terms + var_terms)

	@staticmethod
	def subtractPoly(poly, to_subtract):
		""" subtract one polynomial from the other without simplification"""
		return sp.sympify(str(poly) + str(-1*to_subtract), evaluate=False)
	@staticmethod
	def moveConstRHSNonConstLHS(eqn):
		""" subtract all lhs constant terms and rhs nonconstant terms from both sides """
		remove_poly = RuleHelper.getRHSNonConstLHSConst(eqn)
		RuleHelper.subtractFromEqn(eqn, remove_poly) 

	@staticmethod
	def subtractFromEqn(eqn, poly):
		eqn.left, eqn.right  = RuleHelper.subtractPoly(eqn.left, poly), RuleHelper.subtractPoly(eqn.right, poly)

	@staticmethod
	def setLHSToNum(eqn):
		eqn.left, _ = eqn.left.as_numer_denom() # returns (num,denom) tuple

	@staticmethod
	def simp8Helper(eqn):
		divisor = eqn.left.coeff(x_symb)
		eqn.left, eqn.right = eqn.left/divisor, eqn.right/divisor 

	@staticmethod
	def mult2Helper(eqn):
		""" if both sides of eqn have fractions, then multiply each side by the lcm over all fractions.  """
		# only do this if both sides have addition as top level operation
		if not ( eqn.left.is_Add and eqn.right.is_Add):
			return
		left_denom = map(lambda x: x.as_numer_denom()[1], left.args)
		left_denom = map(lambda x: x.as_numer_denom()[1], right.args)
		# compute lcm and multiply
		lcm = RuleHelper.computeLCM( left_denom + right_denom )
		eqn.left	= sp.mul.Mul.fromiter([eqn.left, lcm])
		eqn.right 	= sp.mul.Mul.fromiter([eqn.right, lcm])

	@staticmethod
	def heur4Helper(std_poly):
		"""
		factors a poly of the form a*x**2 - b
		@return: (factored_poly, True) otherwise (original_poly, False)
		XXX: assumes poly is in standard form
		"""
		# TODO: refactor conversion to ProdPoly (repeated in factor() function
		if not isinstance(std_poly, StdPoly):
			raise TypeError

		quad_coeff, const_coef = std_poly.coeff(x_symb**2), std_poly.coeff(x_symb**0)
		factored_poly = std_poly.factor(extension=sp.sqrt(abs(const_coeff))/sp.sqrt(quad_coeff))
		if isinstance(factored_poly, sp.Mul): # factoring was successful
			return (factored_poly, True)
		else:
			return (std_poly, False)

	@staticmethod
	def isFactored(poly):
		all([sp.degree(p, gens=x_symb) < 2 for p in poly.args ])

############################## RULES ##############################	
# condition, action, description, name
SIMP0 =	PolyRule(	lambda x : x.is_Add and any([p.is_zero for p in x.args]), 
					RuleHelper._removeZeroes,
					""" simp0: if zeroes exist as additive terms, then remove them """,
					'simp0'
					)
SIMP1 =	EqnRule(	lambda eq, wm : eq.degree() >= 2 and not wm.hasGoal(WorkingMem.SET_RHS_ZERO), # TODO: remove this rule and add check to other simp rules instead!
					lambda eq, wm : wm.addGoal(WorkingMem.SET_RHS_ZERO),
					""" simp1: if degree is >= 2, then set working mem goal to make rhs zero """,
					'simp1'
					)
SIMP2 =	PolyRule(	lambda x : x.is_Add and not (x.collect(x_symb) == x),  # collect adds like terms
					lambda x: x.collect(x_symb),
					""" simp2: if sumpoly has common terms, then add them together """,
					'simp2'
					)
SIMP3 =	EqnRule(	lambda eq, wm : eq.degree() == 1 and not RuleHelper.getRHSNonConstLHSConst(eq).is_zero,
					lambda eq, wm : RuleHelper.moveConstRHSNonConstLHS(eq),
					""" if solving a linear eqn, cancel all constant terms on the lhs and all non-constant terms on the rhs """,
					'simp3'
					)
SIMP4 =	EqnRule(	lambda eq, wm : wm.hasGoal(WorkingMem.SET_RHS_ZERO) and not eq.right.is_zero,
					lambda eq, wm : RuleHelper.subtractFromEqn(eq, eq.right),
					""" if our goal is to set rhs to zero, then subtract all rhs terms from lhs""",
					'simp4'
					)
SIMP5 =	PolyRule(	lambda p: not p.is_polynomial()  and p.cancel() != p,  # TODO: better way to check if have common factors
					sp.polys.cancel,
					""" simp5: if num and denom of a rational polynomial have common factors, then cancel these factors """,
					'simp5'
					)
#SIMP6 =	PolyRule(	lambda x : isinstance(x, RatPoly) and x.num.is_zero , 
										#lambda x : StdPoly.zero(),
										#""" simp6: if zero exists in a numerator, remove the fraction involved """,
										#'simp6'
					#)
#SIMP7 =	EqnRule(	lambda eq, wm : isinstance(eq.left, RatPoly) and eq.right.is_zero, 
										#lambda eq,wm : RuleHelper.setLHSToNum(eq),
										#""" if lhs is a rational polynomial, and rhs is zero, solve for numerator """,
										#'simp7'
					#)
#SIMP8 =	EqnRule(	lambda eq, wm : eq.right.isConstTerm() and eq.left.is_linear and eq.left.coeffOf(Bases.X) != 1 and eq.left.coeffOf(Bases.CONST) is None, 
										#lambda eq,wm : RuleHelper.simp8Helper(eq),
										#""" if equation has form ax = b, divide by a """,
										#'simp8'
					#)
## TODO: simp9 will become obsolete if we remove simp1
#SIMP9 =	EqnRule(	lambda eq, wm : eq.degree() < 2 and wm.hasGoal(WorkingMem.SET_RHS_ZERO), 
										#lambda eq,wm : wm.removeGoal(WorkingMem.SET_RHS_ZERO),
										#""" if SET_RHS_ZERO is a goal and we've reduced problem to linear eqn, then remove this goal""",
										#'simp9'
					#)
#MULT1 =	PolyRule(	lambda p: isinstance(p, RatPoly) and isinstance(p.denom, RatPoly), 
										#lambda p: ProdPoly([p.num, p.denom.reciprocal()]),
										#""" if denom of rational poly is a fraction, the multiply by its reciprocal """,
										#'mult1'
					#)
#MULT2 =	EqnRule(	lambda eq, wm : eq.left.hasFractions() and  (eq.right.hasFractions() or eq.right.is_zero), 
										#lambda eq,wm : RuleHelper.mult2Helper(eq),
										#""" if both sides of eqn have fractions, then multiply each side by the lcm over all fractions.  """,
										#'mult2'
					#)
#MULT4 =	PolyRule( lambda p: isinstance(p, SumPoly) and p.hasFractions() and len(p.getFractions()) > 1,
										#RuleHelper.mult4Helper,
										#""" if a polynomial is a sum over rational polynomials, then multiply every polynomial by lcm/lcm""",
										#'mult4'
									#)
#MULT5 =	PolyRule(lambda p: isinstance(p, ProdPoly) and len([poly for poly in p.subpoly if not poly.isConstTerm()]) > 1, # there must be at least two non constant terms for this to make sense
									#RuleHelper.mult5Helper,
									#""" if a there is a product polynomial, then foil the first two factors""",
									#'mult5'
									#)

#HEUR1 =	PolyRule(lambda p:  isinstance(p, StdPoly) and p.degree() == 2 and isinstance(p.factor(x_symb), sp.Mul) # TODO: look for is_factorable() method
									#,lambda p : RuleHelper.factor(p)[0]
									#,""" if a 2nd degree polynomial occurs anywhere, then attempt to factor it """,
									#'heur1'
									#)

#HEUR2 =	PolyRule(lambda p: isinstance(p, StdPoly) and p.degree() == 2 and RuleHelper.completeSquare(p)[1]
									#,lambda p : RuleHelper.completeSquare(p)[0]
									#,""" if a 2nd degree polynomial occurs anywhere, then factor it by completing the square """,
									#'heur2'
									#)

#HEUR3 =	PolyRule(lambda p: p.degree() == 3 and isinstance(p, StdPoly) and RuleHelper.factorCubic(p)[1]
									#,lambda p : RuleHelper.factorCubic(p)[0]
									#,""" if a 3rd degree polynomial occurs anywhere, then attempt to factor it """,
									#'heur3'
				#)

## TODO: see if this rule can be written with other fraction rules
#HEUR4 =	PolyRule(lambda p:  isinstance(p, StdPoly) and p.degree() == 2 and p.coeff_monomial(x_symb)==0 and p.coeff_monomial(x_symb**0) < 0	 # TODO: look for is_factorable() method
									#,lambda p : RuleHelper.heur4Helper(p)[0]
									#,""" if a polynomial of the form ax**2 -b occurs anywhere, then factor it as (x + sqrt(a/b)) (x - sqrt(a/b)) """,
									#'heur4'
									#)
class Solver:
	def __init__(self, eqn, rule_ord=lambda rule: Solver.RULE_ORDERING.index(rule)): 
		self.eqn, self.working_mem	= eqn, WorkingMem()
		self.rule_ordering			= rule_ord

	def __str__(self): return '\n'.join(self.working_mem.steps)
	def copy(self):
		new_solver = Solver(self.eqn.copy(), self.rule_ordering)
		new_solver.working_mem = self.working_mem.copy()
		return new_solver
	############################### win conditions  ##############################
	def win1(self): # """ case a = b"""
		""" contradiction:  reduced to a = b, but a =/= b """
		right, left = self.eqn.right, self.eqn.left
		return (left.is_Number and right.is_Number and left != right)
	def win2(self):
		""" win condition: ax = b problem is solved """
		right, left = self.eqn.right, self.eqn.left
		return sp.degree(left,gens=x_symb) == 1 and right.is_Number and left.coeff(x_symb) == 1 and left.coeff(x_symb**0) == 0
	def win3(self):
		""" win condition: lhs is completely factored, rhs is zero """
		right, left = self.eqn.right, self.eqn.left
		# TODO: revisit this rule, it's gotten complex
		return ( self.eqn.degree() >= 2 and RuleHelper.isFactored(left) and right.is_zero and not left.is_Add and not left.is_polynomial) # not a rational polynomial

	############################### checking and solving ##############################
	### list of rules and precedences they take
	#SIMP_RULES		= [SIMP6,SIMP7, SIMP8, SIMP0, SIMP1, SIMP2, SIMP3, SIMP4, SIMP5, SIMP9 ]
	#WIN_RULES		= [win1, win2, win3]
	#MULT_RULES		= [MULT1, MULT2, MULT4, MULT5]
	#MISC_RULES		= []
	#HEURISTICS		= [HEUR1, HEUR2, HEUR3, HEUR4]
	#ALL_RULES 		= [SIMP_RULES, HEURISTICS, MULT_RULES, MISC_RULES ]
	#RULE_ORDERING 	= SIMP_RULES + HEURISTICS + MULT_RULES + MISC_RULES 

	### solve the problem
	#def solve(self):
		#"""solve the equation given, return steps to the solution"""
		#self.working_mem.steps.append(str(self.eqn) + ': ' + self.solve.__name__)
		#while not self.checkWinCond():
			#triggered_rules = self.getTriggeredRules()
			#if len(triggered_rules) == 0 : # stuck, there's no more todo
				#break

			## select a rule and apply it
			#applied_rule = self.selectRule(triggered_rules)
			#applied_rule.applyAction(self.eqn, self.working_mem)
			#self.working_mem.addStep(self.eqn, applied_rule) 
		#return self.working_mem.steps

	#def checkWinCond(self):
		#""" check win conditions """
		## if any win condition is active, then we've solved the equation
		#for rule in self.WIN_RULES:
			#if rule(self) :
				##print "win condition " + rule.__name__ + " applies"
				#self.working_mem.addStep( self.eqn, rule)
				#return True
		#return False

	#def getTriggeredRules(self):
		#"""@return: list of all rules triggered by current game state"""
		#return [ rule for rule in Solver.RULE_ORDERING if rule.checkCondition(self.eqn, self.working_mem) ]

	#def selectRule(self, triggered_rules): return min(triggered_rules, key=self.rule_ordering)

#class SuperSolver():
	#def __init__(self, eqn_str):
		#self.eqn = Eqn(eqn_str)
		#self.steps_tried = set()

	#def allSolns(self):
		#solutions = []
		#solvers = [Solver(self.eqn)]
		##pdb.set_trace()
		#while len(solvers) > 0:
			## take the top solver
			#soln = solvers.pop()
			#for s in soln.working_mem.steps:
				#print s
			#print "##############################"
			## if finished solving, add it's solution
			#if soln.checkWinCond():
				#solutions.append(soln.working_mem.steps)
				#continue
			## ... otherwise generate a solver for each triggered rule and apply the rule
			#triggered_rules = soln.getTriggeredRules()
			#if len(triggered_rules) == 0: # deadend, no solution down this path
				#continue
			#for rule in triggered_rules:
				#new_solver = soln.copy()
				#if (str(new_solver.eqn), rule.__name__) not in self.steps_tried:
					#self.steps_tried.add((str(new_solver.eqn), rule.__name__)) # remember what you've tried to do
					#rule.applyAction(new_solver.eqn, new_solver.working_mem)
					#new_solver.working_mem.addStep(new_solver.eqn, rule)
					#solvers.append(new_solver)
		#return solutions

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

# causes infinite loop: 3*x**2/(x+1) = 0

# issues: zero not returning as a constant term
#solver = Solver(Eqn('3*x**2/(x+1) = 0'))
#pdb.set_trace()
#applied = MULT5.checkCondition(solver.eqn, solver.working_mem)
#MULT5.applyAction(solver.eqn, solver.working_mem)

#print 'hot dog'
#soln_gen = SuperSolver('3*x**2/(x+1) = 0')
#for soln in soln_gen.allSolns():
#	print soln

print 'puppies'
