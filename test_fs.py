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
	def test_win1(self):
		solver = Solver(Eqn('3=1'))
		self.assertTrue(solver.win1())

		# set the two sides equal
		solver.eqn = Eqn('3=3')
		self.assertFalse(solver.win1())

		# set non-constant terms
		solver.eqn = Eqn('5*x**2=3*(x)*(5*x**2)')
		self.assertFalse(solver.win1())
	def test_win2(self):
		solver = Solver(Eqn('x = 1'))
		self.assertTrue(solver.win2())
	def test_win3(self):
		solver = Solver(Eqn('(2+x)*(2+x)=0'))
		self.assertTrue(solver.win3())

		# set to a higher degree poly, should be false
		solver.eqn = Eqn('5*x**2=0')
		self.assertFalse(solver.win3())

class SampleCasesTest(unittest.TestCase):
	"""
	def test_solve1(self):

		solver = Solver(Eqn('10*x + 3 + -3=10')) # TODO: change so sp.poly doesn't simplify automatically!
		self.assertEqual(solver.solve(), '10x=10')
	def test_solve2(self):
		solver = Solver(Eqn('10*x + 3 = 10 + 3*x'))
		self.assertEqual(solver.solve(), '7x=7')
		"""
	def test_solve3(self):
		solver = Solver(Eqn('3*x + x**2 = 3 + x**2'))
		expected_steps = 	[	'x**2 + 3*x=x**2 + 3: solve',
								'x**2 + 3*x=x**2 + 3: simp1',
								'x**2 + 3*x + -x**2 - 3=x**2 + 3 + -x**2 - 3: simp4',
								'3*x - 3=x**2 + 3 + -x**2 - 3: simp2',
								'3*x - 3=0: simp2'
							]
		#pdb.set_trace()
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)

	def test_solve4(self):
		solver = Solver(Eqn('(x**2*(x+1)*(x+3))/((x+1)*(x+3)) + 3*x = 3 + x**2'))
		expected_steps =	[	'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3: solve',
								'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3: simp1',
								'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x + -x**2 - 3=x**2 + 3 + -x**2 - 3: simp4',
								'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x + -x**2 - 3=0: simp2',
								'x**2 + 3*x + -x**2 - 3=0: simp5',
								'-3 + 3*x=0: simp2'
							]
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)
		self.assertEqual(str(solver.eqn), '-3 + 3*x=0')
	def test_solve5(self):
		solver = Solver(Eqn('(x-1)/(x**2-2*x-3) + (x+2)/(x**2-9) = (2*x+5)/(x**2 + 4*x + 3)'))
		expected_steps =	[	'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9)=(2*x + 5)/(x**2 + 4*x + 3): solve',
								'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9)=(2*x + 5)/(x**2 + 4*x + 3): simp1',
								'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=(2*x + 5)/(x**2 + 4*x + 3) + (-2*x - 5)/(x**2 + 4*x + 3): simp4',
								'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=(0)/(x**2 + 4*x + 3): simp2',
								'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=0: simp6',
								'(x - 1)/((x + 1) * (x - 3)) + (x + 2)/((x - 3) * (x + 3)) + (-2*x - 5)/((x + 1) * (x + 3))=0: heur1',
								'((x + 1) * (x - 3) * (x + 3) * (x - 1))/((x + 1) * (x - 3)) + ((x + 1) * (x - 3) * (x + 3) * (x + 2))/((x - 3) * (x + 3)) + ((x + 1) * (x - 3) * (x + 3) * (-2*x - 5))/((x + 1) * (x + 3))=0: mult2',
								'(x + 3) * (x - 1) + (x + 1) * (x + 2) + (x - 3) * (-2*x - 5)=0: simp5',
								'x**2 + 2*x - 3 + x**2 + 3*x + 2 + -2*x**2 + x + 15=0: mult5',
								'6*x + 14=0: simp2'
							]
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)
		self.assertEqual(str(solver.eqn), '6*x + 14=0')


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
		self.solver = Solver(Eqn('x+1 = x+3'))

	def test_simp0(self):
		self.solver.eqn = Eqn('((x + 1) + 0 )/(x+3)= x+3')
		# TODO: revised test until string parser works completely
		self.solver.eqn.left = RatPoly(SumPoly([Eqn.strToPolyTree('x+1'),Eqn.strToPolyTree('0')]), Eqn.strToPolyTree('x+3'))
		self.assertTrue(self.solver.simp0())
		self.assertEqual(self.solver.eqn.left, RatPoly(self.sp1, self.sp2))

	def test_simp1(self):
		self.solver.eqn = Eqn('(x+1)*(x+3) = x + 3')
		self.assertTrue(self.solver.simp1())
		self.assertTrue(self.solver.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO))
	def test_simp2(self):
		self.solver.eqn = Eqn('(x+1) + (x+3) = x + 1')
		# TODO: revised test until string parser works completely
		self.solver.eqn.left = SumPoly([Eqn.strToPolyTree('x+1'),Eqn.strToPolyTree('x+3')])
		self.assertTrue(self.solver.simp2())
		self.assertEqual(self.solver.eqn.left, StdPoly(2*x_symb + 4) )
	def test_simp3(self):
		self.solver.eqn = Eqn('3*x + 3 = x +1')
		self.assertTrue(self.solver.simp3())
		left = SumPoly([self.sp3, StdPoly(-1*x_symb-3)])
		right = SumPoly([self.sp1, StdPoly(-3-1*x_symb)])
		self.assertEqual(str(self.solver.eqn.left), str(left)) # (3x + 2 - x -2 )
		self.assertEqual(self.solver.eqn.right, right) # (x + 1 -2 -x)
	def test_simp4(self):
		self.solver.working_mem = WorkingMem() # ensure goal is in wm
		self.solver.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)

		self.solver.eqn = Eqn('(x+1)*(x+3) = (3*x + 3)')
		prod_poly = Eqn.strToPolyTree('(x+1)*(x+3)') 
		self.assertTrue( self.solver.simp4())
		self.assertEqual(self.solver.eqn.left, prod_poly.subtract(self.sp3))
		self.assertTrue(self.sp3.subtract(self.sp3.copy()))
	def test_simp5(self):
		self.solver.eqn = Eqn('(3*x+3)/((x+1)*(3*x+3)) = (x+3)')
		self.assertTrue(self.solver.simp5())
		self.assertEqual(self.solver.eqn.left, RatPoly(StdPoly(1, x_symb), self.sp1) )

	def test_mult1(self):
		self.solver.eqn = Eqn('(x+3) / ((3*x+3)/((x+1)*(3*x+3))) = (x+3)')
		self.assertTrue(self.solver.mult1())
		rp = Eqn.strToPolyTree('(3*x+3)/((x+1)*(3*x+3))')
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp2, rp.reciprocal()]))

	def test_mult2(self):
		solver = Solver(Eqn('1/(x+1) = 1/(x+3)'))
		solver.mult2()
		# result should be...
		left, right = Eqn.strToPolyTree('1/(x+1)'), Eqn.strToPolyTree('1/(x+3)')
		left = ProdPoly([left, ProdPoly([self.sp2, self.sp1]) ])
		right = ProdPoly([right, ProdPoly([self.sp2, self.sp1]) ])

		self.assertEqual(solver.eqn.left, left)
		self.assertEqual(solver.eqn.right, right)

	# TODO: test for mult4
	def test_mult5(self):
		self.solver.eqn = Eqn('(x+1)*(x+3) = 3*x+3')
		self.assertTrue(self.solver.mult5())
		self.assertEqual(self.solver.eqn.left, StdPoly(1*x_symb**2 + 3*x_symb + 1*x_symb + 3) ) # x^2 + 3x + x + 3

	def test_heur1(self):
		self.solver.eqn = Eqn('(x**2 + 4*x + 3) = 3*x + 3')
		self.assertTrue(self.solver.heur1())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp1, self.sp2]))

	def test_heur2(self):
		self.solver.eqn = Eqn('x**2 + 2*x + 3 = 3*x + 3')
		self.assertTrue(self.solver.heur2())
		self.assertEqual(str(self.solver.eqn.left), str(SumPoly([ProdPoly([self.sp1, self.sp1]), StdPoly(2, x_symb)])))

	def test_heur3(self):
		# factor down to linear terms
		self.solver.eqn = Eqn('x**3 + x**2 - 9*x - 9 = 3*x + 3')
		self.assertTrue(self.solver.heur3())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp1, self.sp2, self.sp7]))

		# factor when can't reduce a quadratic term
		self.solver.eqn = Eqn('x**3 + 5*x**2 + x  + 5 = 3*x+3')
		self.assertTrue(self.solver.heur3())
		self.assertEqual(self.solver.eqn.left, ProdPoly([self.sp9, self.sp10]))

class SolverUtilTest(unittest.TestCase):
	def test_computeLCM(self):
		sp1, sp2 = StdPoly(x_symb + 3), StdPoly(x_symb + 2)
		solver = Solver(Eqn('x+3 = x + 2'))
		self.assertEqual( ProdPoly( solver.computeLCM([sp1, sp2]) ), ProdPoly([sp1, sp2]))
		self.assertEqual(solver.computeLCM([sp1, sp1.copy()]), [sp1])

#pdb.set_trace()

# TODO: untested code
# mult1
# sample test cases!
