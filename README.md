# Neron-Severi lattice 


## Introduction

NS-Lattice is a Python library for doing calculations in Neron-Severi lattices of real weak del Pezzo surfaces.

This library depends on [SageMath](https://SageMath.org) libraries.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install the `ns_lattice` package: 
```    
sage -pip install ns_lattice
```    
If you do not have root access use the following command instead:
```    
sage -pip install --user ns_lattice
```    

* We advice to upgrade the `ns_lattice` package regularly:
```
sage -pip install --upgrade ns_lattice
```
 
* To execute some [usecases](https://github.com/niels-lubbes/ns_lattice/blob/master/ns_lattice/src/ns_lattice/__main__.py) type:
```    
sage -python -m ns_lattice
```

* For showing which files were installed 
or for uninstalling the `ns_lattice` package, 
use one of the following commands:
```
sage -pip show --files ns_lattice
sage -pip uninstall ns_lattice
```


## Examples

See also [this file](https://github.com/niels-lubbes/linear_series/blob/master/ns_lattice/src/ns_lattice/__main__.py) 
for example usecases. 
See the [source code](https://github.com/niels-lubbes/linear_series/blob/master/ns_lattice/src/ns_lattice)
the io-specification of each function.
The [test functions](https://github.com/niels-lubbes/linear_series/blob/master/ns_lattice/src/tests)
might be informative for how to call function.

For running the examples below, either copy paste the code into the Sage interface or run them as a Python module:

    sage -python -m my_module_name.py


__Example 1: Finding classes in Neron-Severi lattice with prescribed intersection products__
```python
# import required libraries
from ns_lattice.class_div import Div
from ns_lattice.div_in_lattice import get_divs
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_div', 'Div', Div )

# find (-1)-classes in Neron-Severi lattice of cubic del Pezzo surface
h = Div.new('3e0-e1-e2-e3-e4-e5-e6',7)
print( get_divs( h, 1, -1, False ) )
print( get_divs( h, 1, -1, True ) )
print( len(get_divs( h, 1, -1, True ))==27 )
```
Output:

    [e1, e0-e1-e2, 2e0-e1-e2-e3-e4-e5]
    [e1, e2, e3, e4, e5, e6, e0-e1-e2, e0-e1-e3, e0-e2-e3, e0-e1-e4, e0-e2-e4, e0-e3-e4, e0-e1-e5, e0-e2-e5, e0-e3-e5, e0-e4-e5, e0-e1-e6, e0-e2-e6, e0-e3-e6, e0-e4-e6, e0-e5-e6, 2e0-e1-e2-e3-e4-e5, 2e0-e1-e2-e3-e4-e6, 2e0-e1-e2-e3-e5-e6, 2e0-e1-e2-e4-e5-e6, 2e0-e1-e3-e4-e5-e6, 2e0-e2-e3-e4-e5-e6]    
    True
    
    
__Example 2: Equivalence classes of Neron-Severi lattice of real sextic weak del Pezzo surface__
```python
# import required libraries
from ns_lattice.class_dp_lattice import DPLattice
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_dp_lattice', 'DPLattice', DPLattice )

# classification of rank 4 lattices
for dpl in DPLattice.get_reduced_cls( 4 ): print dpl.get_marked_Mtype(),'\t', dpl.type
```
Output:    

    A0      A0
    A0      A1
    A0      A1
    A0      2A1
    A0      A2
    A0      A1+A2
    A1'     A0
    A1      A0
    A1'     A1
    A1      A1
    A1'     A2
    2A1     A0


 
__Example 3: Create Neron-Severi lattice of weak del Pezzo surface of degree 4__
 ```python
# import required libraries
from ns_lattice.class_div import Div
from ns_lattice.dp_involutions import basis_to_involution
from ns_lattice.class_dp_lattice import DPLattice
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_div', 'Div', Div )
sage_register_unpickle_override( 'class_dp_lattice', 'DPLattice', DPLattice )

# construct lattice
a1=Div.new('e2-e4',6)
a2=Div.new('e3-e5',6)
a3=Div.new('e0-e1-e2-e4',6)
a4=Div.new('e0-e1-e3-e5',6)

b1=Div.new('e4-e5',6)
b2=Div.new('e0-e1-e2-e3',6)
M = basis_to_involution( [b1,b2], 6 )

dpl = DPLattice( [a1,a2,a3,a4], [b1,b2], M )
print( dpl ) 
```
Output:

    ==================================================
    Degree          = 4
    Rank            = 6
    Intersection    = [(1, 0, 0, 0, 0, 0), (0, -1, 0, 0, 0, 0), (0, 0, -1, 0, 0, 0), (0, 0, 0, -1, 0, 0), (0, 0, 0, 0, -1, 0), (0, 0, 0, 0, 0, -1)]
    Real structure  = 2A1
    Singularities   = 4A1
    Cardinalities   = (0, 12)
    Real involution:
            e0  --->  2e0-e1-e2-e3
            e1  --->  e0-e2-e3
            e2  --->  e0-e1-e3
            e3  --->  e0-e1-e2
            e4  --->  e5
            e5  --->  e4
    Indecomposable (-2)-classes:
            e2-e4  --->  e0-e1-e3-e5
            e3-e5  --->  e0-e1-e2-e4
            e0-e1-e2-e4  --->  e3-e5
            e0-e1-e3-e5  --->  e2-e4
            #real = 0
    Indecomposable (-1)-classes:
            e1  --->  e0-e2-e3
            e4  --->  e5
            e5  --->  e4
            e0-e2-e3  --->  e1
            #real = 0
    Classes of conical families:
            e0-e1  --->  e0-e1
            e0-e2  --->  e0-e2
            e0-e3  --->  e0-e3
            2e0-e2-e3-e4-e5  --->  2e0-e2-e3-e4-e5
            #real = 4
    ==================================================


__Example 4: Determine isomorphism of Neron-Severi lattices of real weak del Pezzo surfaces__
 ```python
# import required libraries
from ns_lattice.class_div import Div
from ns_lattice.dp_involutions import basis_to_involution
from ns_lattice.class_dp_lattice import DPLattice
from ns_lattice.sage_interface import sage_register_unpickle_override
from ns_lattice.class_ns_tools import NSTools
sage_register_unpickle_override( 'class_div', 'Div', Div )
sage_register_unpickle_override( 'class_dp_lattice', 'DPLattice', DPLattice )
NSTools.filter( 'no-output' ) # disable debug info

# construct first lattice
a1=Div.new('e2-e4',6)
a2=Div.new('e3-e5',6)
a3=Div.new('e0-e1-e2-e4',6)
a4=Div.new('e0-e1-e3-e5',6)

b1=Div.new('e4-e5',6)
b2=Div.new('e0-e1-e2-e3',6)
M = basis_to_involution( [b1,b2], 6 )

dpl1 = DPLattice( [a1,a2,a3,a4], [b1,b2], M )

# construct second lattice
c1=Div.new('e1-e2',6)
c2=Div.new('e3-e4',6)
c3=Div.new('e0-e1-e2-e3',6)
c4=Div.new('e0-e3-e4-e5',6)

d1=Div.new('e1-e2',6)
d2=Div.new('e3-e4',6)
M = basis_to_involution( [b1,b2], 6 )

dpl2 = DPLattice( [a1,a2,a3,a4], [b1,b2], M )

# check equivalence
print( dpl1 == dpl2 )

```
Output:

    True
     
    
__Example 5: Change basis of Neron-Severi lattice__
```python
# import required libraries
from ns_lattice.class_div import Div
from ns_lattice.dp_involutions import basis_to_involution
from ns_lattice.class_dp_lattice import DPLattice
from ns_lattice.ns_basis import get_bases_lst
from ns_lattice.class_ns_tools import NSTools
from ns_lattice.sage_interface import sage_matrix
from ns_lattice.sage_interface import sage_ZZ
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_div', 'Div', Div )
sage_register_unpickle_override( 'class_dp_lattice', 'DPLattice', DPLattice )
NSTools.filter( 'no-output' ) # disable debug info

# construct lattice
a = Div.new('e0-e1-e2-e3',4)
M = basis_to_involution( [a], 4 )
dpl1 = DPLattice( [], [a], M )
dpl1.set_attributes()
print 'dpl1 =', dpl1

# obtain base changes
g1=Div.new('e0-e1',4)
g2=Div.new('e0-e2',4)
bas_lst = get_bases_lst( [g1,g2], dpl1.M, dpl1.d_lst, dpl1.m1_lst )
print 'bas_lst =', bas_lst 

# apply base change to dpl1
B = sage_matrix( sage_ZZ, [ d.e_lst for d in bas_lst[0] ] )
dpl2 = dpl1.get_basis_change( B )
print 'B =\n', B
print 'dpl2 =', dpl2

```
Output:

    dpl1 =
    ==================================================
    Degree          = 6
    Rank            = 4
    Intersection    = [(1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1)]
    Real structure  = A1
    Singularities   = A0
    Cardinalities   = (4, 0)
    Real involution:
            e0  --->  2e0-e1-e2-e3
            e1  --->  e0-e2-e3
            e2  --->  e0-e1-e3
            e3  --->  e0-e1-e2
    Indecomposable (-2)-classes:
            #real = 0
    Indecomposable (-1)-classes:
            e1  --->  e0-e2-e3
            e2  --->  e0-e1-e3
            e3  --->  e0-e1-e2
            e0-e1-e2  --->  e3
            e0-e1-e3  --->  e2
            e0-e2-e3  --->  e1
            #real = 0
    Classes of conical families:
            e0-e1  --->  e0-e1
            e0-e2  --->  e0-e2
            e0-e3  --->  e0-e3
            #real = 3
    ==================================================

    bas_lst = [(e0-e1, e0-e2, e3, e0-e1-e2)]

    B =
    [ 1 -1  0  0]
    [ 1  0 -1  0]
    [ 0  0  0  1]
    [ 1 -1 -1  0]

    dpl2 =
    ==================================================
    Degree          = 6
    Rank            = 4
    Intersection    = [(0, 1, 0, 0), (1, 0, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1)]
    Real structure  = A1
    Singularities   = A0
    Cardinalities   = (4, 0)
    Real involution:
            e0  --->  e0
            e1  --->  e1
            e2  --->  e3
            e3  --->  e2
    Indecomposable (-2)-classes:
            #real = 0
    Indecomposable (-1)-classes:
            e1-e3  --->  e1-e2
            e0-e3  --->  e0-e2
            e2  --->  e3
            e3  --->  e2
            e0-e2  --->  e0-e3
            e1-e2  --->  e1-e3
            #real = 0
    Classes of conical families:
            e0  --->  e0
            e1  --->  e1
            e0+e1-e2-e3  --->  e0+e1-e2-e3
            #real = 3
    ==================================================

