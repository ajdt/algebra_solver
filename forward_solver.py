class CorePoly(object):
	"""basic polynomial, mostly contains code to be overridden"""
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################	
	def add(self, poly) : return SumPoly([self, poly])
	def sub(self, poly) : return SumPoly([self, poly.negate() ])
	def mult(self, poly) : return ProdPoly([self, poly])
	def divide(self, poly) : return RatPoly(num=self, denom=poly)
	def negate(self): raise NotImplementedError

	# UTILITY FUNCTIONS
	def order() : raise NotImplementedError
	def copy(self): raise NotImplementedError
	def __str__(self): raise NotImplementedError
	def __eq__(self, other): raise NotImplementedError
	def __ne__(self, other): raise NotImplementedError

	# BOOLEANS 
	def hasFractions(self): raise NotImplementedError 
	def isFactored(self): return False 
	def isZero(): return False
	def isLinearStdPoly(self): return False  # override for sum and monomial classes

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
	def isLinearStdPoly(self): return self.order() == 1
	def isFactored(self): self.order() < 2

	############################## terms and factors ##############################	
	# TODO: implement all of these
	def distribute(multiplier): return SumPoly( [p.mult(multiplier) for p in self.subpoly] )
	def isConstantTerm(self): return self.order() == 0
	def hasCommonTerms(self): return False # override for SumPoly!!
	def getNonConstTerms(self): return [] # override for SumPoly and BasicPoly
	def getConstantTerms(self): return [] # override for SumPoly and BasicPoly, don't recurse on this
	def sumCommonTerms(self): raise NotImplementedError # add together common terms

class ProdPoly(CorePoly):
	def __init__(self, poly_list): self.subpoly = poly_list

	# OVERRIDE MATH OPERATIONS
	def mult(self,poly): 
		self.subpoly.append(poly) 
		return self
	def negate(self): return ProdPoly( [self.subpoly[0].negate()] + self.subpoly[1:] )

	# UTILITY FUNCTIONS
	def order() : return sum([p.order() for p in self.subpoly])
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
	def isZero(): return False
	def isLinearStdPoly(self): return self.order < 2  


	############################## TERMS AND FACTORS ##############################	
	def foil(self): raise NotImplementedError #return result of multiplying terms together # foil specific terms?
	def commonFactors(self): raise NotImplementedError
	def commonFactors(self, other): raise NotImplementedError
	def isConstantTerm(self): return False
	def hasCommonTerms(self): return False # override for SumPoly!!
	def nonConstantTerms(self): return [] # override for SumPoly and BasicPoly
	def constantTerms(self): return [] # override for SumPoly and BasicPoly, don't recurse on this
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors

class RatPoly(CorePoly):
	def __init__(self, num, denom): self.num, self.denom = num, denom
	############################## OPERATIONS ON POLYNOMIAL (NON-REDUCED VERSIONS)  ##############################	
	def divide(self, other): 
		self.num = self.num.mult(other) 
		return self
	def negate(self): return RatPoly(self.num.negate(), self.denom)

	# UTILITY FUNCTIONS
	def order() : return max([self.num.order(), self.denom.order()])
	def copy(self) : return RatPoly(self.num.copy(), self.denom.copy())
	def __str__(self): return '(' + str(self.num) + ')/(' + str(self.denom) + ')'
	def __eq__(self, other): 
		if not (self.__class__ == other.__class__) :
			return False
		return self.num == other.num and self.denom == other.denom
	def __ne__(self, other): return not (self == other)

	# BOOLEANS 
	def hasFractions(self): return False 
	def isFactored(self): return False 
	def isZero(): return False
	def isLinearStdPoly(self): return False  # override for sum and monomial classes


	############################## MISC OPERATIONS ##############################	
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
	def order() : return self.base
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
	def isFactored(self): return True 
	def isZero(): return self.coef == 0
	def isLinearStdPoly(self): return True  
	def isConstantTerm(self): return self.base == Bases.CONST


	############################## TERMS AND FACTORS ##############################	
	def nonConstantTerms(self): return ([self] if self.base != Bases.CONST  else [] )
	def constantTerms(self): return ([self] if self.base == Bases.CONST  else  [] )
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors
# testing code		
p = CorePoly()
print "ok"

# Notes:
#	all/any for reducing a list of booleans
