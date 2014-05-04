% -- Armando Diaz Tolentino <ajdt@uw.edu> 
% 
%

% monomial components and value ranges
coefficient(-10..10).
exponent(0..3).
symbol(x). % include other symbols for multi-variate equations
side(left;right).

% eqn definitions
side(left;right).
operator(add;mul;div).

num_terms(2..3). % number of terms per poly_tree

% NOTES:
% make_node --	is a placeholder for a subtree to be generated.
% 				it's like a directive that says: make a subtree with 
%				this id_num, at this depth, for this parent node, on this side.
%
% NOTE: a node's id is the predicate id(depth, id_number)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EQN GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each side create at least one poly tree
1 {_poly_tree(Op, NumChild, id(0, 1), parent(0, 0), Side) : operator(Op) : num_terms(NumChild) } 1 :- side(Side).

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

% generation constraints
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
% add round number to generated polytree
_poly_tree(Op,Num,id(D, IDNum),parent( PD, PIDNum), Side, 0) :- _poly_tree(Op,Num,id(D, IDNum),parent( PD, PIDNum), Side).
_monomial(monom(C,B,E),id(D,IDNum),parent(PD,PIDNum),Side, 0) :- _monomial(monom(C,B,E),id(D,IDNum),parent(PD,PIDNum),Side).

#const max_steps = 10. % max number of steps to try in solving an equation. REMOVE LATER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RULES AND TRACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove(_monomial(monom(0,X,Y), ID2, PID, S, R)) :- _monomial(X, ID1, PID, S, R), _monomial(monom(0,X,Y), ID2, PID, S, R). % remove zeroes

% remove poly_tree that has only one child
% TODO: remove the poly_tree
%_monomial(Params, ID, PID, Side, R+1) :- _polytree(Op NumC, ID, PID, Side, R),  [ _monomial(_, _, ID, Side, R) ] 1, _monomial(Params, MonID, ID, Side, R).

% rules for removing a monomial

% moving on to the next round: 
_monomial(Params, ID, PID, Side, R+1) :- _monomial(Params, ID, PID, Side, R), R < max_steps, not remove( _monomial(Params, ID, PID, Side, R) ).
_poly_tree(Op, NumC, ID, PID, Side, R+1) :- _poly_tree(Op, NumC, ID, PID, Side, R), R < max_steps, not remove( _poly_tree(Op, NumC, ID, PID, Side, R) ).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% win conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eqn has form x = b
solved :- _monomial(monom(1,x,1), id(0,_), parent(0,0), left, _), _monomial(monom(_,x,0), id(0,_), parent(0,0), right, _).

#hide.
%#show _poly_tree/5.
%#show _monomial/4.
#show optype/1.
%#show make_node/4.
#show _monomial/8.
#show _poly_tree/7.
