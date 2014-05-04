% -- Armando Diaz Tolentino <ajdt@uw.edu> 
% 
%

% components of a monomial
coefficient(-10..10).
exponent(0..3).
symbol(x).
side(left;right).

% eqn definitions
side(left;right).
operator(add;mul;div).

% parameters for subpolynomials
%depth(1..3).
depth(1).
num_terms(2..3).

eqn(0). % equation(round_number)

% depth, id number, side

% for each side create at least one poly tree

1 {_poly_tree(Op, Children, 0, 1, 0, 0, Side) : operator(Op) : num_terms(Children) } 1 :- side(Side).

% make nodes with a given id and a parent_id
%{make_node(id(D+1, Num-1), Side, Parent) : num_terms(Num) }  :- poly_tree(Op, Child, id(D,Num), Side), Parent := id(D,Num), D < 3.

% why doesn't poly_tree(Op, Child, id(X,Y), Side) work? Nesting predicates
%1 {blah(Num) : num_terms(Num) } 1 :- poly_tree(Op, Child, X, Side), num_terms(N), depth(D), X == id(D,N). %, Parent := id(D,Num), D < 3.

% make_node is a placeholder for a subtree to be generated.
% it's like a directive that says: make a subtree with this id_num, at this depth, for this parent node, on this side.
% NOTE: a node's id is given by both it's depth and id number
NumChild {make_node(Num, Depth, ID, Side) : num_terms(Num) }  NumChild :- _poly_tree(Op, NumChild, Depth, ID, PDepth, PID, Side).


% generating subtrees
1	{	_poly_tree(Op, Child, PDepth + 1, NodeID, PDepth, PID, Side) 
			: operator(Op) 
			: num_terms(Child),
		_monomial(C,B,E, PDepth+1, NodeID, PDepth, PID, Side) 
			: coefficient(C) 
			: symbol(B) 
			: exponent(E)
	} 1 :- make_node(NodeID, PDepth, PID, Side), PDepth < 1. % should be 2

% generating leaf nodes
1	{	_monomial(C,B,E, PDepth+1, NodeID, PDepth, PID, Side) 
			: coefficient(C) 
			: symbol(B) 
			: exponent(E)
	} 1 :- make_node(NodeID, PDepth, PID, Side), PDepth == 1.

% constraints
:- 0{ _monomial(_, _, 1, _, _, _, _, _) } 0.
%:- { make_node(A,B,C) } 0 .

% REMOVE THIS. Just used for testing --project option
optype(add) :- _poly_tree(add,_, _, _, _, _, _).
%optype(div) :- _poly_tree(div,_, _, _, _, _, _).
optype(mul) :- _poly_tree(mul,_, _, _, _, _, _).
optype(blah) :- _monomial(2,x,0,1,_,0,_,left), _monomial(2,x,0,1,_,0,_,right), _poly_tree(add,_,0,_,0,0,left) .

#hide.
#show _poly_tree/7.
#show optype/1.
%#show make_node/4.
#show _monomial/8.
