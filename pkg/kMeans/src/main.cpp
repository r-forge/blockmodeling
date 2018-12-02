#include <iostream>
#include <string.h>

#include "kmBlock_types.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

Array funByBlocks( const Array & M, const IVector & clu, int dimensions, std::string sDiagonal = "default" );

// [[Rcpp::export]]
Array kmBlock( const Array & M, const IVector & clu, Array & weights, const IVector & n, const IVector & nClu )
{
    return funByBlocks( M, clu, Rcpp::sum( nClu ) );
}

Array funByBlocks( const Array & M, const IVector & clu, int dimensions, std::string sDiagonal )
{

    Rcpp::Rcout << "funByBlocks: begin" << std::endl;

    Array aRes( dimensions, dimensions, M.n_slices, arma::fill::zeros );
    Array S = aRes;
    Array N = aRes;
    Rcpp::Rcout << "funByBlocks: dimensions= " << M.n_rows << " " << M.n_cols << " " << M.n_slices << std::endl;
    Rcpp::Rcout << "funByBlocks: Res dimensions= " << aRes.n_rows << " " << aRes.n_cols << " " << aRes.n_slices << std::endl;

    for( size_t i = 0; i < M.n_rows; ++i ) {
        for( size_t j = 0; j < M.n_cols; ++j ) {
//            if( sDiagonal == "default" && i == j )
//                continue;
            for( size_t r = 0; r < M.n_slices; ++r ) {
                S( clu.at( i ) - 1, clu.at( j ) - 1, r ) += M( i, j, r );
                N( clu.at( i ) - 1, clu.at( j ) - 1, r ) += 1;
            }
        }
    }

    for( size_t i = 0; i < aRes.n_rows; ++i ) {
        for( size_t j = 0; j < aRes.n_cols; ++j ) {
            for( size_t r = 0; r < aRes.n_slices; ++r ) {
                aRes( i, j, r ) = double( S( i, j, r ) ) / N( i, j, r );
            }
        }
    }


//    Rcpp::Rcout << "funByBlocks: N= " << N << std::endl << std::endl;
//    Rcpp::Rcout << "funByBlocks: S= " << S << std::endl << std::endl;
    Rcpp::Rcout << "funByBlocks: end" << std::endl;
    return aRes;
}
