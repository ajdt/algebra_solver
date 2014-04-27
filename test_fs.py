from forward_solver import *
import unittest
import pdb

#class PolyTest(unittest.TestCase):
    ## TODO: revise these tests to test other polynomial types
	#def setUp(self):
		## monomials
		#self.x_term = StdPoly(x_symb)
		#self.const = StdPoly(3, x_symb)
		#self.x2_term = StdPoly(5*x_symb**2)
		#self.x_term2 = StdPoly(2*x_symb)

		#self.prod_poly = ProdPoly([self.x_term, self.const, self.x2_term])
		#self.new_term = self.prod_poly.divide(ProdPoly([self.x_term2, self.const]) )

	#def test_eq(self):
		#self.assertTrue( self.x_term == self.x_term ) 				# equal to self
		#self.assertTrue( self.x_term  == StdPoly(x_symb) )	# equal to unique copy of self
		#self.assertFalse(self.x_term == self.const)

		#self.assertFalse(self.x_term == self.prod_poly)
		#self.assertFalse(self.new_term == self.prod_poly)
	#def test_ne(self):
		#self.assertTrue( self.x_term != self.const)
		#self.assertTrue( self.x_term != self.prod_poly)
		#self.assertTrue( self.new_term != self.prod_poly)
	#def test_commonTerms(self):
		#sum_poly = SumPoly([self.x_term, self.const, self.x_term2])
		#self.assertTrue(sum_poly.hasCommonTerms())

		## now add the common terms, and test for equality
		#result = SumPoly([self.const, StdPoly(3*x_symb)])
		#self.assertTrue(sum_poly.sumCommonTerms() == result)
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
	def test_solve1(self):

		solver = Solver(Eqn('10*x + 3 + -3=10')) 
		expected_steps = 	[	'10*x - 3 + 3=10: solve',
								'10*x=10: simp2',
								'x=1: simp8',
								'x=1: win2'
							]
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)
		self.assertEqual(str(solver.eqn), 'x=1')
	def test_solve2(self):
		solver = Solver(Eqn('10*x + 3 = 10 + 3*x'))
		expected_steps = 	[	'10*x + 3=3*x + 10: solve', 
								'-3*x + 10*x - 3 + 3=-3*x + 3*x - 3 + 10: simp3', 
								'7*x=-3*x + 3*x - 3 + 10: simp2', 
								'7*x=7: simp2',
								'x=1: simp8', 
								'x=1: win2'
							]
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)
	def test_solve3(self):
		solver = Solver(Eqn('3*x + x**2 = 3 + x**2'))
		expected_steps = 	[	'x**2 + 3*x=x**2 + 3: solve',
								'x**2 + 3*x=x**2 + 3: simp1',
								'-x**2 + x**2 + 3*x - 3=-x**2 + x**2 - 3 + 3: simp4',
								'3*x - 3=-x**2 + x**2 - 3 + 3: simp2',
								'3*x - 3=0: simp2',
								'3*x - 3 + 3=3: simp3',
								'3*x=3: simp2',
								'x=1: simp8',
								'x=1: win2'
							]
		steps = solver.solve()
		self.assertEqual(steps, expected_steps)

	#def test_solve4(self):
		#solver = Solver(Eqn('(x**2*(x+1)*(x+3))/((x+1)*(x+3)) + 3*x = 3 + x**2'))
		#expected_steps =	[	'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3: solve',
								#'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x=x**2 + 3: simp1',
								#'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x + -x**2 - 3=x**2 + 3 + -x**2 - 3: simp4',
								#'((x**2) * (x + 1) * (x + 3))/((x + 1) * (x + 3)) + 3*x + -x**2 - 3=0: simp2',
								#'x**2 + 3*x + -x**2 - 3=0: simp5',
								#'-3 + 3*x=0: simp2',
								#'-3 + 3*x + 3=0 + 3: simp3',
								#'-3 + 3*x + 3=3: simp0',
								#'0 + 3*x=3: simp2',
								#'x=(3)/(3): simp8',
								#'x=(3)/(3): win2'
							#]
		#print "##############################"
		#for s in expected_steps:
			#print s
		#print "------------------------------"
		#for s in steps:
			#print s
		#print steps
		#print "##############################"
		#steps = solver.solve()
		#self.assertEqual(steps, expected_steps)
		#self.assertEqual(str(solver.eqn), 'x=(3)/(3)')
	#def test_solve5(self):
		#solver = Solver(Eqn('(x-1)/(x**2-2*x-3) + (x+2)/(x**2-9) = (2*x+5)/(x**2 + 4*x + 3)'))
		#expected_steps =	[	'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9)=(2*x + 5)/(x**2 + 4*x + 3): solve',
								#'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9)=(2*x + 5)/(x**2 + 4*x + 3): simp1',
								#'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=(2*x + 5)/(x**2 + 4*x + 3) + (-2*x - 5)/(x**2 + 4*x + 3): simp4',
								#'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=(0)/(x**2 + 4*x + 3): simp2',
								#'(x - 1)/(x**2 - 2*x - 3) + (x + 2)/(x**2 - 9) + (-2*x - 5)/(x**2 + 4*x + 3)=0: simp6',
								#'(x - 1)/((x + 1) * (x - 3)) + (x + 2)/((x - 3) * (x + 3)) + (-2*x - 5)/((x + 1) * (x + 3))=0: heur1',
								#'((x + 1) * (x - 3) * (x + 3) * (x - 1))/((x + 1) * (x - 3)) + ((x + 1) * (x - 3) * (x + 3) * (x + 2))/((x - 3) * (x + 3)) + ((x + 1) * (x - 3) * (x + 3) * (-2*x - 5))/((x + 1) * (x + 3))=0: mult2',
								#'(x + 3) * (x - 1) + (x + 1) * (x + 2) + (x - 3) * (-2*x - 5)=0: simp5',
								#'x**2 + 2*x - 3 + x**2 + 3*x + 2 + -2*x**2 + x + 15=0: mult5',
								#'6*x + 14=0: simp2',
								#'6*x + 14 + -14=0 + -14: simp3',
								#'6*x + 14 + -14=-14: simp0',
								#'6*x=-14: simp2',
								#'x=(-14)/(6): simp8',
								#'x=(-14)/(6): win2'
							#]
		#steps = solver.solve()
		#self.assertEqual(steps, expected_steps)
		#self.assertEqual(str(solver.eqn), 'x=(-14)/(6)')


class SolverRuleTest(unittest.TestCase):
	def setUp(self):
		self.sp1 = x_symb + 1 # (x + 1)
		self.sp2 = x_symb + 3 # (x + 3)
		self.sp3 = 3*x_symb + 3 # (3x + 3)
		self.sp4 = x_symb**2 + 4*x_symb + 3 # (x^2 + 4x + 3)
		self.sp5 = x_symb**2 + 2*x_symb + 3 # (x^2 + 2x + 3) , for completing the square test
		self.sp6 = x_symb**3 + x_symb**2 + -9*x_symb + -9 # (x^3 + x^2 -9x -9 )
		self.sp7 = 1*x_symb + -3 # (x - 3)
		self.sp8 = 1*x_symb**3 + 5*x_symb**2 + x_symb + 5 # (x^3 + 5x^2 +x +5 )
		self.sp9 = x_symb + 5 # (x + 3)
		self.sp10 = x_symb**2 + 1 # (x^2 + 1)
		self.solver = Solver(Eqn('x+1 = x+3'))

	def test_simp0(self):
		self.solver.eqn = Eqn('((x + 1) + 0 )/(x+3)= x+3')
		self.solver.eqn.left = sp.sympify('(x+1+0)/(x+3)', evaluate=False)
		#self.assertTrue(SIMP0.checkCondition(self.solver.eqn, self.solver.working_mem)) # TODO: sympy autosimplifies zeroes, causes test to fail
		SIMP0.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(self.solver.eqn.left, self.sp1/self.sp2)

	def test_simp1(self):
		self.solver.eqn = Eqn('(x+1)*(x+3) = x + 3')
		self.assertTrue(SIMP1.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP1.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertTrue(self.solver.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO))
	def test_simp2(self):
		self.solver.eqn = Eqn('x + 1 + x + 3 = x + 1')
		self.assertTrue(SIMP2.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP2.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(self.solver.eqn.left, 2*x_symb + 4 )
	def test_simp3(self):
		self.solver.eqn = Eqn('3*x + 3 = x +1')
		self.assertTrue(SIMP3.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP3.applyAction(self.solver.eqn, self.solver.working_mem)

		# expected results
		left	= sp.sympify('3*x + 3 -x -3 ', evaluate=False)
		right	= sp.sympify('x + 1 -x -3 ', evaluate=False)
		self.assertEqual(str(self.solver.eqn.left), str(left)) # (3x + 2 - x -2 )
		self.assertEqual(self.solver.eqn.right, right) # (x + 1 -2 -x)
		solver = Solver(Eqn('3*x**2 = 0'))
		self.assertFalse( SIMP3.checkCondition(solver.eqn, solver.working_mem))
	def test_simp4(self):
		self.solver.working_mem = WorkingMem() # ensure goal is in wm
		self.solver.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)

		self.solver.eqn = Eqn('(x+1)*(x+3) = (3*x + 3)')
		prod_poly = sp.sympify('(x+1)*(x+3)', evaluate=False) 
		to_subtract = 3*x_symb + 3
		self.assertTrue( SIMP4.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP4.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left),  '-3*x + (x + 1)*(x + 3) - 3') # TODO: equations compared as strings due to sympy peculiarity
		self.assertEqual(str(self.solver.eqn.right), '-3*x + 3*x - 3 + 3')
	def test_simp5(self):
		self.solver.eqn = Eqn('(3*x+3)/((x+1)*(3*x+3)) = (x+3)')
		self.assertTrue(SIMP5.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP5.applyAction(self.solver.eqn, self.solver.working_mem)
		comp_to = sp.sympify('1/(x+1)', evaluate=False)
		self.solver.eqn.left == comp_to
		self.assertEqual(str(self.solver.eqn.left),str(comp_to)) # TODO: have to compare as strings due to hashable_content() issue, tuples have same content but reordered
	def test_simp6(self):
		self.solver.eqn = Eqn('0/(x**2+3*x+3) = (x+3)')
		self.assertTrue(SIMP6.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP6.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertTrue(self.solver.eqn.left.is_zero) 
	def test_simp7(self):
		self.solver.eqn = Eqn('(3*x**2 + 5*x)/(x**2+3*x+3) = 0')
		self.assertTrue(SIMP7.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP7.applyAction(self.solver.eqn, self.solver.working_mem)
		new_lhs = 3*x_symb**2 + 5* x_symb
		self.assertTrue(self.solver.eqn.left, new_lhs) 
	def test_simp8(self):
		self.solver.eqn = Eqn('3*x = 5')
		self.assertTrue(SIMP8.checkCondition(self.solver.eqn, self.solver.working_mem))
		SIMP8.applyAction(self.solver.eqn, self.solver.working_mem)
		new_lhs = x_symb 
		self.assertTrue(self.solver.eqn.left, new_lhs) 
	def test_simp9(self):
		solver = Solver(Eqn('3*x = 5'))
		solver.working_mem.addGoal(WorkingMem.SET_RHS_ZERO)
		self.assertTrue(SIMP9.checkCondition(solver.eqn, solver.working_mem))
		SIMP9.applyAction(solver.eqn, solver.working_mem)
		self.assertFalse(solver.working_mem.hasGoal(WorkingMem.SET_RHS_ZERO))

	def test_mult1(self):
		self.solver.eqn = Eqn('(x+3) / ((3*x+3)/((x+1)*(3*x+3))) = (x+3)')
		self.assertTrue(MULT1.checkCondition(self.solver.eqn, self.solver.working_mem))
		MULT1.applyAction(self.solver.eqn, self.solver.working_mem)
		new_lhs = (x_symb+3)*(x_symb+1)*(3*x_symb+3)/(3*x_symb+3)
		self.assertEqual(self.solver.eqn.left, new_lhs)

	def test_mult2(self):
		solver = Solver(Eqn('1/(x+1) = 1/(x+3)'))
		self.assertTrue(MULT2.checkCondition(solver.eqn, solver.working_mem))
		MULT2.applyAction(solver.eqn, solver.working_mem)

		self.assertEqual(str(solver.eqn.left), '(x + 1)*(x + 3)/(x + 1)') 
		self.assertEqual(str(solver.eqn.right), '(x + 1)*(x + 3)/(x + 3)')

	def test_mult4(self):
		solver = Solver(Eqn('1/(x+1) + 1/(x+3) = 0 '))
		self.assertTrue(MULT4.checkCondition(solver.eqn, solver.working_mem))
		MULT4.applyAction(solver.eqn, solver.working_mem)

		self.assertEqual(str(solver.eqn.left), '(x + 3)/((x + 1)*(x + 3)) + (x + 1)/((x + 1)*(x + 3))') 

	def test_mult5(self):
		self.solver.eqn = Eqn('(x+1)*(x+3) = 3*x+3')
		self.assertTrue(MULT5.checkCondition(self.solver.eqn, self.solver.working_mem))
		MULT5.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), 'x**2 + x + 3*x + 3') 

	def test_heur1(self):
		self.solver.eqn = Eqn('(x**2 + 4*x + 3) = 3*x + 3')
		self.assertTrue(HEUR1.checkCondition(self.solver.eqn, self.solver.working_mem))
		HEUR1.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), '(x + 1)*(x + 3)')

	def test_heur2(self):
		self.solver.eqn = Eqn('x**2 + 2*x + 3 = 3*x + 3')
		self.assertTrue(HEUR2.checkCondition(self.solver.eqn, self.solver.working_mem))
		HEUR2.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), '(x + 1)*(x + 1) + 2') 

	def test_heur3(self):
		# factor down to linear terms
		self.solver.eqn = Eqn('x**3 + x**2 - 9*x - 9 = 3*x + 3')
		self.assertTrue(HEUR3.checkCondition(self.solver.eqn, self.solver.working_mem))
		HEUR3.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), '(x - 3)*(x + 1)*(x + 3)')

		# factor when can't reduce a quadratic term
		self.solver.eqn = Eqn('x**3 + 5*x**2 + x  + 5 = 3*x+3')
		self.assertTrue(HEUR3.checkCondition(self.solver.eqn, self.solver.working_mem))
		HEUR3.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), '(x + 5)*(x**2 + 1)')
	
	def test_heur4(self):
		self.solver.eqn = Eqn('x**2 - 25 = 0')
		self.assertTrue(HEUR4.checkCondition(self.solver.eqn, self.solver.working_mem))
		HEUR4.applyAction(self.solver.eqn, self.solver.working_mem)
		self.assertEqual(str(self.solver.eqn.left), '(x - 5)*(x + 5)')

class RuleHelperTest(unittest.TestCase):
	def test_computeLCM(self):
		sp1, sp2 = x_symb + 3, x_symb + 2
		self.assertEqual( RuleHelper.computeLCM([sp1, sp2]), sp1*sp2)
		self.assertEqual(RuleHelper.computeLCM([sp1, sp1.copy()]), sp1)

#pdb.set_trace()

# TODO: untested code
# mult1
# sample test cases!
