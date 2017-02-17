'''
Created on Aug 11, 2016
@author: Niels Lubbes

Classification of involutions of Neron-Severi group of 
weak del Pezzo surfaces.

See [http://arxiv.org/abs/1302.6678] for more info.
'''
from sage.all import *

from div_set import *
from dp_root_bases import *

from class_ns_tools import NSTools
from class_div import Div


nt = NSTools()


def complete_basis( d_lst ):
    '''
    INPUT:
        - "d_lst" -- A list of "Div" objects.
    OUTPUT:
        - Returns a square matrix over QQ of full rank. The first columns
          correspond to the elements in d_lst (where d_lst is sorted). 
          The appended columns are orthogonal to the first "len(d_lst)" columns.
         
    EXAMPLES:
        - We explain with 3 examples where the dotted (:) vertical lines 
          denote which columns are appended. 
          
          | 0  0 : 1  0  0  0 | | 0  0  0 : 1  0  0 | | 0  0  0  1 : -3  0 |  
          | 0  0 : 0 -1  0  0 | | 0  0  0 : 0 -1  0 | | 1  0  0 -1 :  1  0 |
          | 0  0 : 0  0 -1  0 | | 1  0  0 : 0  0 -1 | |-1  1  0 -1 :  1  0 |
          | 1  0 : 0  0  0 -1 | |-1  1  0 : 0  0 -1 | | 0 -1  0 -1 :  1  0 |
          |-1  1 : 0  0  0 -1 | | 0 -1  1 : 0  0 -1 | | 0  0  1  0 :  0  1 |
          | 0 -1 : 0  0  0 -1 | | 0  0 -1 : 0  0 -1 | | 0  0 -1  0 :  0  1 |
             
    '''
    # sort
    d_lst = copy( d_lst )
    d_lst.sort()

    # extend with orthogonal vectors
    row_lst = [ d.e_lst for d in d_lst]
    ext_lst = []
    for v_lst in  matrix( ZZ, row_lst ).right_kernel().basis():
        ext_lst += [ [-v_lst[0]] + list( v_lst[1:] ) ]  # accounts for signature (1,rank-1)
    mat = matrix( QQ, row_lst + ext_lst ).transpose()

    # verify output
    if mat.rank() < d_lst[0].rank():
        raise Error( 'Matrix expected to have full rank: ', d_lst, '\n' + str( mat ) )
    de_lst = [ Div( ext ) for ext in ext_lst ]
    for de in de_lst:
        for d in d_lst:
            if d * de != 0:
                raise Error( 'Extended columns are expected to be orthogonal: ', de, d, de_lst, d_lst, list( mat ) )

    return mat


def is_integral_involution( M ):
    '''
    INPUT:
        - A matrix M.
    OUTPUT:
        - Returns True if the matrix is an involution, 
          preserves inner product with signature (1,r) 
          and has integral coefficients.
    '''

    nrows, ncols = M.dimensions()

    # check whether involution
    if M * M != identity_matrix( nrows ):
        return False

    # check whether inner product is preserved
    S = diagonal_matrix( [1] + ( ncols - 1 ) * [-1] )
    if M.transpose() * S * M != S:
        return False

    # check whether coefficients are integral
    for r in range( nrows ):
        for c in range( ncols ):
            if M[r][c] not in ZZ:
                return False

    return True


def get_cls_involutions( max_rank = 9 ):
    '''
    See [http://arxiv.org/abs/1302.6678] for more info.
    
    INPUT:
        - "max_rank" -- An integer in [3,...,9].
    
    OUTPUT:
        - Returns a dictionary "invo_cls_dct" with keys in [3,...,max_rank+1].
          The values for each key consist of a list of matrices over
          ZZ that correspond to (non-equivalent) involutions of 
              ZZ<h,e1,...,er>
          where r=rank-1. Up to equivalence the list contains all possible 
          involutions. Thus a key for "invo_dct" denotes the rank.          
    '''

    # classification of involutions in cache?
    for rank in range( max_rank, 9 + 1 ):
        key = 'get_cls_involutions_' + str( rank )
        if key in nt.get_tool_dct():
            return nt.get_tool_dct()[key]
    key = 'get_cls_involutions_' + str( max_rank )

    # compute classification of involutions
    invo_cls_dct = {}
    for rank in range( 3, max_rank + 1 ):

        nt.p( rank )

        M_lst = [( identity_matrix( ZZ, rank ), [] )]  # include identity involution
        for d_lst in get_cls_root_bases( max_rank )[rank]:

            if d_lst == []:
                continue

            l = len( d_lst )
            V = complete_basis( d_lst )
            D = diagonal_matrix( l * [-1] + ( rank - l ) * [1] )
            M = V * D * V.inverse()  # MV=VD

            if is_integral_involution( M ):
                M_lst += [( M, d_lst )]

        invo_cls_dct[rank] = M_lst

    # store classification
    nt.get_tool_dct()[key] = invo_cls_dct
    nt.save_tool_dct()

    return invo_cls_dct


def get_involutions( rank ):
    '''
    INPUT:
        - rank  -- An integer denoting the rank of lattice:
                   2<="rank"<=9.".  
    OUTPUT:
        - A list [(<M>,<b_lst>),...]
          where <M> is a matrix and <b_lst> a list of Div objects. 
          Each matrix <M> corresponds to a unimodular involution 
          of ZZ<h,e1,...,er> where r=rank-1.
          Up to isomorphism of the eigenspace of 1, 
          each such involution is contained.   
          In order to include all involutions the method 
          "complete_basis()" should be replaced by a methods
          which computes all possible completions of a basis.
    '''
    key = 'get_involutions_' + str( rank )

    # involutions in cache?
    if key in nt.get_tool_dct():
        return nt.get_tool_dct()[key]

    # construct list of all involutions
    MB_lst = [( identity_matrix( ZZ, rank ), [] )]  # include identity involution
    for d_lst in get_root_bases( rank, False ):

        # identity matrix already included
        if d_lst == []:
            continue

        nt.p( 'rank =', rank, 'd_lst = ', d_lst )

        # construct involution matrix corresponding to root basis
        l = len( d_lst )
        V = complete_basis( d_lst )
        D = diagonal_matrix( l * [-1] + ( rank - l ) * [1] )
        M = V * D * V.inverse()  # M*V==V*D

        # if unimodular then store involution
        if is_integral_involution( M ):
            MB_lst += [( M, d_lst )]


    # store involutions
    nt.get_tool_dct()[key] = MB_lst
    nt.save_tool_dct()

    return MB_lst






