#include <iostream>
#include <memory>
#include <string.h>

#include "kmBlock_types.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

Array funByBlocks( const Array & M, const IVector & clu, int dimensions, DMatrix & p_pSepare, Diagonale sDiagonal = Diagonale::Same );
double meanMatrix( const DMatrix & p_matrix );

// [[Rcpp::export]]
Rcpp::List kmBlock( const Array & M, const IVector & clu, Array & weights, const IVector & n, const IVector & nClu )
{
    DMatrix pSeparate;
    Array aRes = funByBlocks( M, clu, Rcpp::sum( nClu ), pSeparate, Diagonale::Ignore );
    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }

    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

Array funByBlocks( const Array & M, const IVector & clu, int dimensions, DMatrix & p_pSepare, Diagonale sDiagonal )
{

//    Rcpp::Rcout << "funByBlocks: begin" << std::endl;

    Array aRes( dimensions, dimensions, M.n_slices, arma::fill::zeros );
    Array S = aRes;
    Array N = aRes;
    DMatrix mDiagonalRes;
    DMatrix mSseprateDiagonal;
    DMatrix mNseprateDiagonal;
    if( sDiagonal == Diagonale::Seperate ) {
        mDiagonalRes = DMatrix( dimensions, M.n_slices, arma::fill::zeros );
        mSseprateDiagonal = mDiagonalRes;
        mNseprateDiagonal = mDiagonalRes;
    }
//    Rcpp::Rcout << "funByBlocks: dimensions= " << M.n_rows << " " << M.n_cols << " " << M.n_slices << std::endl;
//    Rcpp::Rcout << "funByBlocks: Res dimensions= " << aRes.n_rows << " " << aRes.n_cols << " " << aRes.n_slices << std::endl;

    for( size_t r = 0; r < M.n_slices; ++r ) {
        for( size_t i = 0; i < M.n_rows; ++i ) {
            for( size_t j = 0; j < M.n_cols; ++j ) {
                if( sDiagonal == Diagonale::Ignore && i == j ) { // Ignore diagonal
                    continue;
                }
                else if( sDiagonal == Diagonale::Seperate && i == j ) { // Calculate diagonal seperately and ignore it in return Array
                    mSseprateDiagonal( clu.at( i ), r ) += M( i, j, r );
                    mNseprateDiagonal( clu.at( i ), r ) += 1;
                    continue;
                }
//                Rcpp::Rcout << std::endl << "clu(i)=" << clu.at( i ) + 1 << std::endl << "clu(j)=" << clu.at( j ) + 1 << std::endl << "r=" << r << std::endl;
//                Rcpp::Rcout << std::endl << "i=" << i << std::endl << "j=" << j << std::endl << "r=" << r << std::endl;
//                Rcpp::Rcout << "------------------------" << std::endl;
                S( clu.at( i ), clu.at( j ), r ) += M( i, j, r );
                N( clu.at( i ), clu.at( j ), r ) += 1;
            }
        }
    }

    for( size_t r = 0; r < aRes.n_slices; ++r ) {
        double diagMean( sDiagonal == Diagonale::Ignore ? meanMatrix( M.slice( r ) ) : 0 ); // calculate it only once for each iteration
        for( size_t i = 0; i < aRes.n_rows; ++i ) {
            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
                mDiagonalRes( i, r ) = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
            }
            for( size_t j = 0; j < aRes.n_cols; ++j ) {
                double dVal( S( i, j, r ) );
                if( !dVal && sDiagonal == Diagonale::Ignore ) { // If value of the block is and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
                    aRes(i, j, r ) = diagMean;
                }
                else {
                    aRes( i, j, r ) = dVal / N( i, j, r );
                }
            }
        }
    }

    if( sDiagonal == Diagonale::Seperate ) { // Save seperate digaonal values to input parameter
        p_pSepare = std::move( mDiagonalRes );
    }


//    Rcpp::Rcout << "funByBlocks: N= " << N << std::endl << std::endl;
//    Rcpp::Rcout << "funByBlocks: S= " << S << std::endl << std::endl;
//    Rcpp::Rcout << "funByBlocks: Diagonal= " << std::endl << mDiagonalRes << std::endl << std::endl;
//    Rcpp::Rcout << "funByBlocks: slice( 1 )= " << std::endl << M.slice( 0 ) << std::endl << std::endl;
//    Rcpp::Rcout << "funByBlocks: mean( slice( ( 1 ) )= " << meanMatrix( M.slice( 0 ) ) << std::endl;
//    Rcpp::Rcout << "funByBlocks: end" << std::endl;
    return aRes;
}


double meanMatrix( const DMatrix & p_matrix )
{
    size_t sElementsN = 0;
    size_t sum = 0;
    for( size_t i = 0; i < static_cast<size_t>( p_matrix.n_rows ); ++i ) {
        for( size_t j = 0; j < static_cast<size_t>( p_matrix.n_cols ); ++j ) {
            sum += p_matrix.at( i, j );
            ++sElementsN;
        }
    }
    return static_cast<double>( sum ) / sElementsN;
}
