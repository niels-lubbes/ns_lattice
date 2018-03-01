# Neron-Severi lattice 


## Introduction

NS-Lattice is a Python library for doing calculations in Neron-Severi lattices of real weak del Pezzo surfaces.   
 
This library depends on [SageMath](https://SageMath.org) libraries.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install the `NS-lattice` package: 
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

For running the examples below, either copy paste the code into the Sage interface or run them as a Python module:

    sage -python -m my_module_name.py


__Example 1: Divisor arithmetic__
```python
from ns_lattice.class_div import Div
from ns_lattice.div_in_lattice import get_divs
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_div', 'Div', Div )

h = Div.new('3e0-e1-e2-e3-e4-e5-e6',7)
a = Div.new('e0-e1-e2',7)

print( (h*h,h*a,a*a,a+a) )
print( get_divs( h, 1, -1, False ) )
print( len(get_divs( h, 1, -1, True ))==27 )
```
    (3, 1, -1, 2e0-2e1-2e2)
    [e1, e0-e1-e2, 2e0-e1-e2-e3-e4-e5]
    True
    
__Example 2: __

    
__Example 3: Create Neron-Severi lattice of weak del Pezzo surface of degree 4__
 ```python
from ns_lattice.class_div import Div
from ns_lattice.dp_involutions import basis_to_involution
from ns_lattice.class_dp_lattice import DPLattice
from ns_lattice.sage_interface import sage_register_unpickle_override
sage_register_unpickle_override( 'class_div', 'Div', Div )
sage_register_unpickle_override( 'class_dp_lattice', 'DPLattice', DPLattice )

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


    



     
    