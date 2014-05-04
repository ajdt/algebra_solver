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

1 {_poly_tree(Op, NumChild, id(0, 1), parent(0, 0), Side) : operator(Op) : num_terms(NumChild) } 1 :- side(Side).

% make nodes with a given id and a parent_id
%{make_node(id(D+1, Num-1), Side, Parent) : num_terms(Num) }  :- poly_tree(Op, Child, id(D,Num), Side), Parent := id(D,Num), D < 3.

% why doesn't poly_tree(Op, Child, id(X,Y), Side) work? Nesting predicates
%1 {blah(Num) : num_terms(Num) } 1 :- poly_tree(Op, Child, X, Side), num_terms(N), depth(D), X == id(D,N). %, Parent := id(D,Num), D < 3.

% make_node is a placeholder for a subtree to be generated.
% it's like a directive that says: make a subtree with this id_num, at this depth, for this parent node, on this side.
% NOTE: a node's id is given by both it's depth and id number
NumChild {make_node(Num, Depth, ID, Side) : num_terms(Num) }  NumChild :- _poly_tree(Op, NumChild, id(Depth, ID), parent(PDepth, PID), Side).


% generating subtrees
1	{	_poly_tree(Op, Child, id(PDepth + 1, NodeID), parent(PDepth, PID), Side) 
			: operator(Op) 
			: num_terms(Child),
		_monomial(monom(C,B,E), id(PDepth+1, NodeID), parent(PDepth, PID), Side) 
			: coefficient(C) 
			: symbol(B) 
			: exponent(E)
	} 1 :- make_node(NodeID, PDepth, PID, Side), PDepth < 1. % should be 2

% generating leaf nodes
1	{	_monomial(monom(C,B,E), id(PDepth+1, NodeID), parent(PDepth, PID), Side) 
			: coefficient(C) 
			: symbol(B) 
			: exponent(E)
	} 1 :- make_node(NodeID, PDepth, PID, Side), PDepth == 1.

% constraints
:- 0{ _monomial(monom(_, _, 1),id( _, _), parent(_, _), _) } 0.
%:- { make_node(A,B,C) } 0 .

% REMOVE THIS. Just used for testing --project option
optype(add) :- _poly_tree(add,_,id( _, _),parent( _, _), _).
%optype(div) :- _poly_tree(div,_, _, _, _, _, _).
optype(mul) :- _poly_tree(mul,_,id( _, _), parent(_, _), _).
optype(blah) :- _monomial(monom(2,x,0),id(1,_),parent(0,_),left), _monomial(monom(2,x,0),id(1,_),parent(0,_),right), _poly_tree(add,_,id(0,_),parent(0,0),left) .


% TODO: REMOVE. These rules exist just so I don't have to modify my python visualization tool
_poly_tree(Op, Num, D, IDNum, PD, PIDNum, Side) :- _poly_tree(Op,Num,id(D, IDNum),parent( PD, PIDNum), Side).
_monomial(C,B,E,D,IDNum,PD,PIDNum,Side) :- _monomial(monom(C,B,E),id(D,IDNum),parent(PD,PIDNum),Side).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% win conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#hide.
%#show _poly_tree/5.
%#show _monomial/4.
#show optype/1.
%#show make_node/4.
#show _monomial/8.
#show _poly_tree/7.
