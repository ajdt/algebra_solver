from forward_solver import *
import unittest
import pdb

class PolyTest(unittest.TestCase):
    # TODO: revise these tests to test other polynomial types
	def setUp(self):
		# monomials
		self.x_term = StdPoly(x_symb)
		self.const = StdPoly(3, x_symb)
		self.x2_term = StdPoly(5*x_symb**2)
		self.x_term2 = StdPoly(2*x_symb)

		self.prod_poly = ProdPoly([self.x_term, self.const, self.x2_term])
		self.new_term = self.prod_poly.divide(ProdPoly([self.x_term2, self.const]) )

	def test_eq(self):
		self.assertTrue( self.x_term == self.x_term ) 				# equal to self
		self.assertTrue( self.x_term  == StdPoly(x_symb) )	# equal to unique copy of self
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
		result = SumPoly([self.const, StdPoly(3*x_symb)])
		self.assertTrue(sum_poly.sumCommonTerms() == result)
class SolverTest(unittest.TestCase):
    # TODO: rewrite tests (more cases)
	def setUp(self):
		# monomials
		self.x_term = StdPoly(x_symb)
		self.const = StdPoly(3, x_symb)
		self.const2 = StdPoly(1, x_symb)
		self.x2_term = StdPoly(5*x_symb**2)
		self.x_term2 = StdPoly(2*x_symb)

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
		solver = Solver(Eqn(self.x_term, self.const2))
		self.assertTrue(solver.win2())
	def test_win3(self):
		solver = Solver(Eqn(ProdPoly([self.x_term2, self.x_term2]), StdPoly(0, x_symb)))
		self.assertTrue(solver.win3())

		# set to a higher degree poly, should be false
		solver.eqn = Eqn(self.x2_term, StdPoly(0, x_symb))
		self.assertFalse(solver.win3())

class SampleCasesTest(unittest.TestCase):
	"""
	def test_solve1(self):
		left = SumPoly([StdPoly(10*x_symb + 3), StdPoly(-3)])
		right = StdPoly(10)

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '10x=10')
	def test_solve2(self):
		left = StdPoly( 10*x_symb + 3 )
		right = StdPoly(10 + 3*x_symb )

		solver = Solver(Eqn(left, right))
		self.assertEqual(solver.solve(), '7x=7')
		"""
	def test_solve3(self):
		left = StdPoly(3*x_symb + x_symb**2)
		right = StdPoly(3 + x_symb**2)

		solver = Solver(Eqn(left, right))
		expected_steps = 	[	'x**2 + 3*x=x**2 + 3:initial state',
								'x**2 + 3*x=x**2 + 3: set working mem goal, if degree is >= 2 '
							] 
		self.assertEqual(solver.solve(), expected_steps)

	def test_solve4(self):
		sp1 = StdPoly(x_symb + 1)
		sp2 = StdPoly(x_symb + 3)

		num = ProdPoly([StdPoly(x_symb**2), sp1, sp2])
		denom  = ProdPoly([sp1, sp2])

		left = SumPoly([RatPoly(num,denom), StdPoly(3*x_symb)])
		right = StdPoly(3 + x_symb**2)

		solver = Solver(Eqn(left, right))
		expected_steps = 	[	'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3:initial state', 
								'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3: set working mem goal, if degree is >= 2 '
							]
		self.assertEqual(solver.solve(), expected_steps)
		self.assertEqual(str(solver.eqn), '3x + -3=0')


class SolverRuleTest(unittest.TestCase):
	def setUp(self):
		self.sp1 = StdPoly(x_symb + 1) # (x + 1)
		self.sp2 = StdPoly(x_symb + 3) # (x + 3)
		self.sp3 = StdPoly(3*x_symb + 3) # (3x + 3)
		self.sp4 = StdPoly(x_symb**2 + 4*x_symb + 3) # (x^2 + 4x + 3)
		self.sp5 = StdPoly(x_symb**2 + 2*x_symb + 3) # (x^2 + 2x + 3) , for completing the square test
		self.sp6 = StdPoly(x_symb**3 + x_symb**2 + -9*x_symb + -9) # (x^3 + x^2 -9x -9 )
		self.sp7 = StdPoly(1*x_symb + -3) # (x - 3)
		self.sp8 = StdPoly(1*x_symb**3 + 5*x_symb**2 + x_symb + 5) # (x^3 + 5x^2 +x +5 )
		self.sp9 = StdPoly(x_symb + 5) # (x + 3)
		self.sp10 = StdPoly(x_symb**2 + 1) # (x^2 + 1)
		self.solver = Solver(Eqn(self.sp1, self.sp2))

	def test_simp0(self):
		rp = RatPoly(SumPoly([self.sp1, StdPoly(0, x_symb)]), self.sp2) # (x+1 + 0)
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
		self.solver.eqn = Eqn(SumPoly([self.sp1,self.sp2]), self.sp1)
		self.assertTrue(self.solver.simp2())
		self.assertEqual(self.solver.eqn.left, StdPoly(2*x_symb + 4) )
	def test_simp3(self):
		self.solver.eqn = Eqn(self.sp3, self.sp1) # (3x + 3 = x +1 )
		self.assertTrue(self.solver.simp3())
		left = SumPoly([self.sp3, StdPoly(-1*x_symb-3)])
		right = SumPoly([self.sp1, StdPoly(-3-1*x_symb)])
		self.assertEqual(str(self.solver.eqn.left), str(left)) # (3x + 2 - x -2 )
		self.assertEqual(self.solver.eqn.right, right) # (x + 1 -2 -x)
	def test_simp4(self):
		self.solver.working_mem = WorkingMem() # ensure goal is in wm
		self.solver.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)

		# set eqn to (x+1)(x+3) = (3x + 3)
		prod_poly = ProdPoly([self.sp1, self.sp2])
		self.solver.eqn = Eqn(prod_poly, self.sp3)
		self.assertTrue( self.solver.simp4())
		self.assertEqual(self.solver.eqn.left, prod_poly.subtract(self.sp3))
		self.assertTrue(self.sp3.subtract(self.sp3.copy()))
	def test_simp5(self):
		rp = RatPoly(self.sp3, ProdPoly([self.sp1, self.sp3]))
		self.solver.eqn = Eqn(rp, self.sp2)
		self.assertTrue(self.solver.simp5())
		self.assertEqual(self.solver.eqn.left, RatPoly(StdPoly(1, x_symb), self.sp1) )

	def test_mult1(self):
		rp = RatPoly(self.sp3, ProdPoly([self.sp1, self.sp3]))
		rp2 = RatPoly(self.sp2, rp) # (x+3) / (3x+3)/((x+1)(3x+3))
		self.solver.eqn = Eqn(rp2, self.sp2)
		self.assertTrue(self.solver.mult1())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp2, rp.reciprocal()]))

	def test_mult2(self):
		sp1, sp2 = self.sp1, self.sp2
		left = RatPoly(StdPoly(1, x_symb), sp1)
		right = RatPoly(StdPoly(1, x_symb), sp2)
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
		self.assertEqual(self.solver.eqn.left, StdPoly(1*x_symb**2 + 3*x_symb + 1*x_symb + 3) ) # x^2 + 3x + x + 3

	def test_heur1(self):
		self.solver.eqn = Eqn(self.sp4, self.sp3)
		self.assertTrue(self.solver.heur1())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp1, self.sp2]))

	def test_heur2(self):
		self.solver.eqn = Eqn(self.sp5, self.sp3)
		self.assertTrue(self.solver.heur2())
		self.assertEqual(str(self.solver.eqn.left), str(SumPoly([ProdPoly([self.sp1, self.sp1]), StdPoly(2, x_symb)])))

	def test_heur3(self):
		# factor down to linear terms
		self.solver.eqn = Eqn(self.sp6, self.sp3)
		self.assertTrue(self.solver.heur3())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp1, self.sp2, self.sp7]))

		# factor when can't reduce a quadratic term
		self.solver.eqn = Eqn(self.sp8, self.sp3)
		self.assertTrue(self.solver.heur3())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp9, self.sp10]))

class SolverUtilTest(unittest.TestCase):
	def test_computeLCM(self):
		sp1, sp2 = StdPoly(x_symb + 3), StdPoly(x_symb + 2)
		solver = Solver(Eqn(sp1, sp2))
		self.assertEqual(solver.computeLCM([sp1, sp2]), ProdPoly([sp1, sp2]))
		print str(solver.computeLCM([sp1, sp1.copy()]))
		self.assertEqual(solver.computeLCM([sp1, sp1.copy()]), sp1)

#pdb.set_trace()

# TODO: untested code
# mult1
# sample test cases!
