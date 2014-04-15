from forward_solver import *
import unittest 
import pdb

class PolyTest(unittest.TestCase):
	def setUp(self):
		# monomials
		self.x_term = Monomial(coef=1,base=Bases.X)
		self.const = Monomial(coef=3,base=Bases.CONST)
		self.x2_term = Monomial(coef=5, base=Bases.X2)
		self.x_term2 = Monomial(coef=2, base=Bases.X)

		self.prod_poly = ProdPoly([self.x_term, self.const, self.x2_term])
		self.new_term = self.prod_poly.divide(ProdPoly([self.x_term2, self.const]) )

	def test_eq(self):
		self.assertTrue( self.x_term == self.x_term ) 				# equal to self
		self.assertTrue( self.x_term  == Monomial(coef=1,base=Bases.X) )	# equal to unique copy of self
		self.assertFalse(self.x_term == self.const)				

		self.assertFalse(self.x_term == self.prod_poly)
		self.assertFalse(self.new_term == self.prod_poly)
	def test_ne(self):
		self.assertTrue( self.x_term != self.const)
		self.assertTrue( self.x_term != self.prod_poly)
		self.assertTrue( self.new_term != self.prod_poly)
	def test_commonTerms(self):
		sum_poly = SumPoly([self.x_term, self.const, self.x_term2])
		self.assertTrue(sum_poly.hasCommonTerms())

		# now add the common terms, and test for equality
		result = SumPoly([self.const, Monomial(self.x_term.coef + self.x_term2.coef, self.x_term.base)])
		self.assertTrue(sum_poly.sumCommonTerms() == result)
class SolverTest(unittest.TestCase):
	def setUp(self):
		# monomials
		self.x_term = Monomial(coef=1,base=Bases.X)
		self.const = Monomial(coef=3,base=Bases.CONST)
		self.const2 = Monomial(coef=1, base=Bases.CONST)
		self.x2_term = Monomial(coef=5, base=Bases.X2)
		self.x_term2 = Monomial(coef=2, base=Bases.X)

		# prod poly
		self.prod_poly = ProdPoly([self.x_term, self.const, self.x2_term])
		self.new_term = self.prod_poly.divide(ProdPoly([self.x_term2, self.const]) )

	def test_win1(self):
		solver = Solver(Eqn(self.const, self.const2))
		self.assertTrue(solver.win1())

		# set the two sides equal
		solver.eqn = Eqn(self.const, self.const)
		self.assertFalse(solver.win1())

		# set non-constant terms
		solver.eqn = Eqn(self.x2_term, self.prod_poly)
		self.assertFalse(solver.win1())
	def test_win2(self):
		solver = Solver(Eqn(self.x2_term, self.const2))
		self.assertTrue(solver.win2())
	def test_win3(self):
		solver = Solver(Eqn(ProdPoly([self.x_term2, self.x_term2]), Monomial(0, Bases.CONST)))
		self.assertTrue(solver.win3())

		# set to a higher degree poly, should be false
		solver.eqn = Eqn(self.x2_term, Monomial(0, Bases.CONST))
		print 'is factored ' + str(solver.eqn.left.isFactored()) + str(solver.eqn.left.order())
		self.assertFalse(solver.win3())
	#def
class SampleCasesTest(unittest.TestCase):
	def test_solve1(self):
		left = SumPoly([Monomial(10, Bases.X), Monomial(3, Bases.CONST), Monomial(-3, Bases.CONST)])
		right = Monomial(10, Bases.CONST)

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '10x=10')
	def test_solve2(self):
		left = SumPoly([Monomial(10, Bases.X), Monomial(3, Bases.CONST)])
		right = SumPoly([Monomial(10, Bases.CONST), Monomial(3, Bases.X)])

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '7x=7')
	def test_solve3(self):
		left = SumPoly([Monomial(3, Bases.X), Monomial(1, Bases.X2)])
		right = SumPoly([Monomial(3, Bases.CONST), Monomial(1, Bases.X2)])

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '3x=3')
	def test_solve4(self):
		sp1 = SumPoly([Monomial(1, Bases.X), Monomial(1, Bases.CONST)])
		sp2 = SumPoly([Monomial(1, Bases.X), Monomial(3, Bases.CONST)])

		num = ProdPoly([Monomial(1, Bases.X2), sp1, sp2])
		denom  = ProdPoly([sp1, sp2])

		left = SumPoly([RatPoly(num,denom), Monomial(3, Bases.X)])
		right = SumPoly([Monomial(3, Bases.CONST), Monomial(1, Bases.X2)])

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '3x + -3=0')

class SolverRuleTest(unittest.TestCase):
	def setUp(self):
		self.sp1 = SumPoly([Monomial(1, Bases.X), Monomial(1, Bases.CONST)]) # (x + 1)
		self.sp2 = SumPoly([Monomial(1, Bases.X), Monomial(3, Bases.CONST)]) # (x + 3)
		self.sp3 = SumPoly([Monomial(3, Bases.X), Monomial(3, Bases.CONST)]) # (3x + 3)
		self.solver = Solver(Eqn(self.sp1, self.sp2))

	def test_simp0(self):
		rp = RatPoly(self.sp1.add(Monomial(0, Bases.X)), self.sp2) # (x+1 + 0)
		self.solver.eqn = Eqn(rp, self.sp2)
		self.assertTrue(self.solver.simp0())
		#print "simp0 " + str(self.solver.eqn.left)
		#print "simp0 " + str(RatPoly(self.sp1, self.sp2))
		self.assertEqual(self.solver.eqn.left, RatPoly(self.sp1, self.sp2))
	
	def test_simp1(self):
		self.solver.eqn = Eqn(ProdPoly([self.sp1, self.sp2]), self.sp2.copy())
		self.assertTrue(self.solver.simp1())
		self.assertTrue(self.solver.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO))
	def test_simp2(self):
		self.solver.eqn = Eqn(self.sp1.add(self.sp2), self.sp1)
		self.assertTrue(self.solver.simp2())
		self.assertEqual(self.solver.eqn.left, SumPoly([Monomial(2, Bases.X), Monomial(4, Bases.CONST)]) )
	def test_simp3(self):
		self.solver.eqn = Eqn(self.sp3, self.sp1) # (3x + 3 = x +1 )
		self.assertTrue(self.solver.simp3())
		left = self.sp3.add( SumPoly([Monomial(-1, Bases.X), Monomial(-3, Bases.CONST)]) )
		right = self.sp1.add(SumPoly([Monomial(-3, Bases.CONST), Monomial(-1, Bases.X)]))
		self.assertEqual(self.solver.eqn.left, left) # (3x + 2 - x -2 )
		self.assertEqual(self.solver.eqn.right, right) # (x + 1 -2 -x)
	def test_simp4(self):
		self.solver.working_mem = WorkingMem() # ensure goal is in wm
		self.solver.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)

		# set eqn to (x+1)(x+3) = (3x + 3)
		prod_poly = ProdPoly([self.sp1, self.sp2])
		self.solver.eqn = Eqn(prod_poly, self.sp3)
		self.assertTrue( self.solver.simp4())
		self.assertEqual(self.solver.eqn.left, prod_poly.sub(self.sp3))
		self.assertTrue(self.sp3.sub(self.sp3.copy()))
	def test_simp5(self):
		rp = RatPoly(self.sp3, ProdPoly([self.sp1, self.sp3]))
		self.solver.eqn = Eqn(rp, self.sp2)
		self.assertTrue(self.solver.simp5())
		self.assertEqual(self.solver.eqn.left, RatPoly(Monomial(1, Bases.CONST), self.sp1) )
	
	def test_mult1(self):
		rp = RatPoly(self.sp3, ProdPoly([self.sp1, self.sp3]))
		rp2 = RatPoly(self.sp2, rp) # (x+3) / (3x+3)/((x+1)(3x+3))
		self.solver.eqn = Eqn(rp2, self.sp2)
		self.assertTrue(self.solver.mult1())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp2, rp.reciprocal()]))

	def test_mult2(self):
		sp1, sp2 = self.sp1, self.sp2
		left = RatPoly(Monomial(1, Bases.CONST), sp1)
		right = RatPoly(Monomial(1, Bases.CONST), sp2)
		solver = Solver(Eqn(left, right))
		solver.mult2()
		# result should be...
		left = ProdPoly([left, ProdPoly([sp2, sp1]) ])
		right = ProdPoly([right, ProdPoly([sp2, sp1]) ])
		self.assertEqual(solver.eqn.left, left)
		self.assertEqual(solver.eqn.right, right)

	# TODO: test for mult4
	def test_mult5(self):
		self.solver.eqn = Eqn(ProdPoly([self.sp1, self.sp2]), self.sp3)
		self.assertTrue(self.solver.mult5())
		self.assertEqual(self.solver.eqn.left, SumPoly([ Monomial(1, Bases.X2), Monomial(3, Bases.X), Monomial(1, Bases.X), Monomial(3, Bases.CONST)]) ) # x^2 + 3x + x + 3

class SolverUtilTest(unittest.TestCase):
	def test_computeLCM(self):
		m1, m2, m3 = Monomial(1, Bases.X), Monomial(3, Bases.CONST), Monomial(2, Bases.CONST)
		sp1, sp2 = SumPoly([m1, m2]), SumPoly([m1, m3])
		solver = Solver(Eqn(sp1, sp2))
		self.assertEqual(solver.computeLCM([sp1, sp2]), ProdPoly([sp1, sp2]))
		print str(solver.computeLCM([sp1, sp1.copy()]))
		self.assertEqual(solver.computeLCM([sp1, sp1.copy()]), sp1)

#pdb.set_trace()

# TODO: untested code
# mult1
