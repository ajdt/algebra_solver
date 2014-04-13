class CorePoly(object):
	"""basic polynomial, mostly contains code to be overridden"""
	############################## operations on polynomial (non-reduced versions)  ##############################	
	def add(self, poly) : return SumPoly([self, poly])
	def sub(self, poly) : return SumPoly([self, poly.negate() ])
	def mult(self, poly) : return ProdPoly([self, poly])
	def divide(self, poly) : return RatPoly(num=self, denom=poly)

	# utility functions
	def order() : raise NotImplementedError
	def deepCopy(self): raise NotImplementedError

	# booleans 
	def hasFractions(self): raise NotImplementedError 
	def isFactored(self): return False 
	def isZero(): return False
	def isLinearStdPoly(self): return False  # override for sum and monomial classes


	############################## terms and factors ##############################	
	def isConstantTerm(self): return False
	def hasCommonTerms(self): return False # override for SumPoly!!
	def nonConstantTerms(self): return [] # override for SumPoly and BasicPoly
	def constantTerms(self): return [] # override for SumPoly and BasicPoly, don't recurse on this
	def shareFactors(self, other): raise NotImplementedError # override, return true if these polys share factors

# testing code		
p = CorePoly()
print "ok"
