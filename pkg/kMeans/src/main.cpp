#include <iostream>
#include <memory>
#include <string.h>

#include "kmBlock_types.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Exposed functions
// [[Rcpp::export]]
Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const int dimensions, const std::string diagonal = "same" );
// [[Rcpp::export]]
Rcpp::List kmBlock( const Array & M, const IVector & clu, Array & weights, const IVector & n, const IVector & nClu );
// [[Rcpp::export]]
double criterialFunction( const Array & M, const IVector & clu, const Array & weights, int dimensions );

// Function declarations
Array meansByBlocks( const Array & M, const IVector & clu, int dimensions, DMatrix & p_pSepare, Diagonale sDiagonal = Diagonale::Same );
double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat );
DVector setGroups( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat, const size_t K );
double meanMatrix( const DMatrix & p_matrix );



Rcpp::List kmBlock( const Array & M, const IVector & clu, Array & weights, const IVector & n, const IVector & nClu )
{
    DMatrix pSeparate;
    Array aRes = meansByBlocks( M, clu, Rcpp::sum( nClu ), pSeparate, Diagonale::Ignore );
    double cf = criterialFunction( M, clu, weights, aRes );
    Rcpp::Rcout << "Criterial function value: " << cf << std::endl;
    // Rcpp::Rcout << "Set Groups: " << setGroups( M, clu, weights, aRes, Rcpp::sum( nClu ) ) << std::endl;
    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }

    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const int dimensions, const std::string diagonal )
{
    Diagonale dDiag( Diagonale::Same );
    if( diagonal == "same" );
    else if( diagonal == "seperate" ) dDiag = Diagonale::Seperate;
    else if( diagonal == "ignore"   ) dDiag = Diagonale::Ignore;
    else Rcpp::stop( "Unknow diagonal parameter\nOptions are: [ same, ignore, seperate ]\n" );

    DMatrix pSeparate;
    Array aRes = meansByBlocks( M, clu, dimensions, pSeparate, dDiag );
    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }
    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

double criterialFunction( const Array & M, const IVector & clu, const Array & weights, int dimensions )
{
    DMatrix pSeparate;
    Array aRes = meansByBlocks( M, clu, dimensions, pSeparate, Diagonale::Ignore );
    return criterialFunction( M, clu, weights, aRes );
}

Array meansByBlocks( const Array & M, const IVector & clu, int dimensions, DMatrix & p_pSepare, Diagonale sDiagonal )
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
        double diagMean( sDiagonal == Diagonale::Ignore ? meanMatrix( M.slice( r ) ) : 0 ); // calculate it only once per each iteration
        for( size_t i = 0; i < aRes.n_rows; ++i ) {
            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
                mDiagonalRes( i, r ) = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
            }
            for( size_t j = 0; j < aRes.n_cols; ++j ) {
                double dVal( S( i, j, r ) );
                if( !dVal && sDiagonal == Diagonale::Ignore ) { // If value of the block is 0 and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
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


double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat )
{
    double dRet = 0;

    for( size_t i = 0; i < M.n_rows; ++i ) {
        for( size_t j = 0; j < M.n_cols; ++j ) {
            for( size_t r = 0; r < M.n_slices; ++r ) {
                dRet += weights( i, j, r ) * std::pow( M( i, j, r ) - meansMat( clu.at( i ), clu.at( j ), r ), 2 );
            }
        }
    }

    return  dRet;
}


/*  tole vse sem zakomentiral, ker sem nekaj spreminjal in popravlja, a so zagotovo napake, ker se nisem pretirano
ukvarjal s tem, da bi bila pravilna c++ koda. To poskrbite vi. Zaradi tega sem tudi zakomentiral klic te funkciej v kmBlock

DVector setGroups( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat, const size_t K )
{
	
//  tale vektor e se računa za vsako enoto posebej. Torej, ko računate e, je i fiksen. Je pa potrebno potem to ponoviti za vse enote
//  nato je potrebno za ta i izračunati, pri katerem indeksu (k) je e najmanjši. Pravzaprav, zdaj ko razmišljam, bi bilo celo bolj učinkovito tako, kot sem naredil jaz spodaj. Prosim preveriti, če je prav, ker je dodan komentar //###, pomeni, da sem to dodal ali spreminjal

	double eMin = Inf; //###
	double eTmp //###
	int kMin;//###
	DVector eVec(clu.size()); //###
    IVector vRet;
//  DVector e( K );
    for( size_t i = 0; i < clu.size(); ++i ) {//###
		eTmp = 0;//###
		// predlagam, da zaradi  večje učinkovitosti, clu.at( i ) tu shranite v eno spremenljivko in jo potem v sledečih zankah uporabljate leto (razen, če menite, da se pri izvajajnju to skoraj ne pozna
		for( size_t k = 0; k < K; ++k ) {
			for( size_t j = 0; j < clu.size(); ++j ) { // tu sem i spremenil v j
				for( size_t r = 0; r < M.n_slices; ++r ) {
					eTmp += weights( i, j, r ) * std::pow( M( i, j, r ) - meansMat( clu.at( i ), clu.at( j ), r ), 2 );
					eTmp += weights( j, i, r ) * std::pow( M( j, i, r ) - meansMat( clu.at( j ), clu.at( i ), r ), 2 );
					Rcpp::Rcout << "k = " << k << ", eTmp " << eTmp << ", eMin " << eMin << ", kMin " << kMin << std::endl;
				}
			}
			if (eTmp < eMin){//###
				kMin = k;//###
				eMin = eTmp		//###	
			}//###
		}
		clu[i] = kMin //###	 Tole in spodnje še posebej ne vem, če je prav	
		eVec[i] = eMin //###		
	}

	//    return vRet;
	// Kar tukaj sedaj še manjka je, da se prepričamo,da nobena skupina ni prazna. 
	// Za to morate nekako šteti, koliko enot je v vsaki skupini (od 0 do K)
	// če je kakšna prazna, date notri enoto, ki ima največjo vrednost v eVec
	// to ponavljate, dokler ni nobene prazne skupine.
	
    return clu; //### oziroma, bistvo je, da se posodobi clu - lahko tudi nič ne vrača, le pregleda se
}
*/

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
