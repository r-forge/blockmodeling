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
Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const IVector & n, const IVector & nClu );
// [[Rcpp::export]]
double critFunction( const Array & M, const IVector & clu, const Array & weights, const int dimensions );

// Function declarations
void meansByBlocks( const Array & M, Array & res, const IVector & clu, const int dimensions, DMatrix & p_pSepare, const Diagonale sDiagonal = Diagonale::Same );
double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat );
void setGroups( const Array & M, IVector & clu, const Array & weights, const Array & meansMat, const IVector & nClu, const IVector & n );
unsigned int belongsTo( const int & group, const IVector & borders );
double meanMatrix( const DMatrix & p_matrix );


Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const IVector & n, const IVector & nClu )
{
    const int K = Rcpp::sum( nClu );

    DMatrix pSeparate;
    Array meanBlocks;
    meansByBlocks( M, meanBlocks, clu, K, pSeparate, Diagonale::Ignore );


//    IVector newClu = setGroups( M, clu, weights, meanBlocks, K );
    IVector newClu = Rcpp::clone( clu );
    setGroups( M, newClu, weights, meanBlocks, nClu, n );
    IVector bestClu;

    double newCf = criterialFunction( M, clu, weights, meanBlocks );
    double bestCf = DBL_MAX;

    while( newCf < bestCf ) {
        bestClu = newClu;
        bestCf = newCf;
        meansByBlocks( M, meanBlocks, newClu, K, pSeparate, Diagonale::Ignore );
        setGroups( M, newClu, weights, meanBlocks, nClu, n );
        newCf = criterialFunction( M, newClu, weights, meanBlocks );
    Rcpp::Rcout << "newCf: " << newCf << std::endl;
    Rcpp::Rcout << "newClu: " << newClu<< std::endl;
		
    }

    Rcpp::Rcout << "BestCf: " << bestCf << std::endl;
    Rcpp::Rcout << "BestClu: " << bestClu << std::endl;

//    if( !pSeparate.is_empty() ) {
//        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = meanBlocks, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
//    }

    return Rcpp::List::create( Rcpp::Named( "bestCf" ) = bestCf, Rcpp::Named( "bestClu" ) = bestClu );

}

Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const int dimensions, const std::string diagonal )
{
    Diagonale dDiag( Diagonale::Same );
    if( diagonal == "same" );
    else if( diagonal == "seperate" ) dDiag = Diagonale::Seperate;
    else if( diagonal == "ignore"   ) dDiag = Diagonale::Ignore;
    else Rcpp::stop( "Unknow diagonal parameter\nOptions are: [ same, ignore, seperate ]\n" );

    DMatrix pSeparate;
    Array aRes;
    meansByBlocks( M, aRes, clu, dimensions, pSeparate, dDiag );
    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }
    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

double critFunction( const Array & M, const IVector & clu, const Array & weights, const int dimensions )
{
    DMatrix pSeparate;
    Array aRes;
    meansByBlocks( M, aRes, clu, dimensions, pSeparate, Diagonale::Ignore );
    return criterialFunction( M, clu, weights, aRes );
}

void meansByBlocks( const Array & M, Array & res, const IVector & clu, const int dimensions, DMatrix & p_pSepare, const Diagonale sDiagonal )
{
//    Rcpp::Rcout << "meansByBlocks: begin" << std::endl;
    if( res.is_empty() ) {
        res = Array( dimensions, dimensions, M.n_slices, arma::fill::zeros );
    }
    else {
        res.fill( 0 );
    }
    Array S = res;
    Array N = res;
    DMatrix mDiagonalRes;
    DMatrix mSseprateDiagonal;
    DMatrix mNseprateDiagonal;

    if( sDiagonal == Diagonale::Seperate ) {
        mDiagonalRes = DMatrix( dimensions, M.n_slices, arma::fill::zeros );
        mSseprateDiagonal = mDiagonalRes;
        mNseprateDiagonal = mDiagonalRes;
    }

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
                S( clu.at( i ), clu.at( j ), r ) += M( i, j, r );
                N( clu.at( i ), clu.at( j ), r ) += 1;
            }
        }
    }

    for( size_t r = 0; r < res.n_slices; ++r ) {
        double diagMean( sDiagonal == Diagonale::Ignore ? meanMatrix( M.slice( r ) ) : 0 ); // calculate it only once per each iteration
        for( size_t i = 0; i < res.n_rows; ++i ) {
            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
                mDiagonalRes( i, r ) = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
            }
            for( size_t j = 0; j < res.n_cols; ++j ) {
                double dVal( S( i, j, r ) );
                if( !dVal && sDiagonal == Diagonale::Ignore ) { // If value of the block is 0 and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
                    res(i, j, r ) = diagMean;
                }
                else {
                    res( i, j, r ) = dVal / N( i, j, r );
                }
            }
        }
    }

    if( sDiagonal == Diagonale::Seperate ) { // Save seperate digaonal values to input parameter
        p_pSepare = std::move( mDiagonalRes );
    }
//    Rcpp::Rcout << "meansByBlocks: end" << std::endl;
}


double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat )
{
    double dRet = 0;
    for( size_t i = 0; i < M.n_rows; ++i ) {
        for( size_t j = 0; j < M.n_cols; ++j ) {
            for( size_t r = 0; r < M.n_slices; ++r ) {
                if( i == j ) { // ignore diagonal
                    continue;
                }
                dRet += weights( i, j, r ) * std::pow( ( M( i, j, r ) - meansMat( clu.at( i ), clu.at( j ), r ) ), 2 );
            }
        }
    }

    return  dRet;
}

void setGroups( const Array & M, IVector & clu, const Array & weights, const Array & meansMat, const IVector & nClu, const IVector & n )
{
    IVector borders = Rcpp::cumsum( nClu );
//    Rcpp::Rcout << "Cumsum: " << borders << std::endl;
    double eMin;
    double eTmp;
    int K = Rcpp::sum( nClu );
    int kMin( 0 );
    DVector eVec( clu.size() );
//    IVector vRet = Rcpp::clone( clu );
    DVector e( K );
    IVector countGroups( K );
	
    for( unsigned int i = 0; i < static_cast<unsigned int>( clu.size() ); ++i ) {
        eMin = DBL_MAX;
        int group( clu.at( i ) );
        unsigned int iBelongsTo( belongsTo( group, borders ) );
//        Rcpp::Rcout << "iBelongsTo=" << iBelongsTo << std::endl;
        unsigned int k = 0;
        if( iBelongsTo ) {
            k = borders.at( iBelongsTo - 1 );
        }
//        Rcpp::Rcout << "for k in " << k << ":" << borders.at( iBelongsTo ) << std::endl;
        for( ; k < static_cast<unsigned int>( borders.at( iBelongsTo ) ); ++k ) {
            eTmp = 0;
            for( unsigned int j = 0; j < static_cast<unsigned int>( clu.size() ); ++j ) {
                unsigned int cluJ = clu.at( j );
                if( i == j ) { // ignore diagonal
                    continue;
                }
                for( unsigned int r = 0; r < M.n_slices; ++r ) {
                    eTmp += weights( i, j, r ) * std::pow( M( i, j, r ) - meansMat( k, cluJ, r ), 2 );
                    eTmp += weights( j, i, r ) * std::pow( M( j, i, r ) - meansMat( cluJ, k, r ), 2 );
				}
			}
            if (eTmp < eMin){
                kMin = k;
                eMin = eTmp;
            }

        }
        clu.at( i ) = kMin;
        eVec.at( i ) = eMin;
		++countGroups.at(kMin);
	}
//    IVector countGroups( K );
//    for( int i = 0; i < countGroups.size(); ++i ) {
//        countGroups.at( i ) = std::count( clu.begin(), clu.end(), i );
//    }
//
    IVector nBorders( Rcpp::cumsum( n ) );

	int iBegin, iEnd, kBegin, kEnd, iMax;
	double eMax;
    for( int s = 0; s < nBorders.size(); ++s ) {
        if( !s ) {
            iBegin = 0;
			kBegin = 0;
        } else {
            iBegin = nBorders.at( s - 1 );
            kBegin = borders.at( s - 1 );
        }
		
		iEnd = nBorders.at( s);
		kEnd = borders.at( s);
		for(int k = kBegin; k < kEnd; k++){
			if(countGroups.at(k)==0){
				eMax=-1;
				iMax=iBegin;
				for(int i = iBegin; i< iEnd; ++i){
//					Rcpp::Rcout << " clu.at(i) " << clu.at(i) << std::endl;				
					if(eVec.at(i)>eMax && countGroups.at(clu.at(i))>1){
//						Rcpp::Rcout << "i " << i << std::endl;
						eMax=eVec.at(i);
						iMax=i;
					}					
				}
//				Rcpp::Rcout << " old k " << clu.at(iMax) << " new k " << k << std::endl;			
				-- countGroups.at(clu.at(iMax));
				clu.at(iMax)=k;
				countGroups.at(k)=1;				
				eVec.at(iMax)=0; //this is not really needed
			}
		}
    }
	
	
	
//    Rcpp::Rcout << "nBorders: " << nBorders << std::endl;
//    Rcpp::Rcout << "borders: " << borders << std::endl;
//    Rcpp::Rcout << "nClu: " << nClu << std::endl;
//    Rcpp::Rcout << "clu: " << clu << std::endl;
//    Rcpp::Rcout << "eVec: " << eVec << std::endl;
//    Rcpp::Rcout << "-----------------------" << std::endl;

//    int iBegin, iEnd, k;
//    for( unsigned int i = 0; i < nBorders.size(); ++i ){
//        if( !i ) {
//            iBegin = 0;
//            k = 0;
//        }
//        else {
//            iBegin = nBorders.at( i - 1 );
//            k = borders.at( i - 1 );
//        }
//        iEnd = nBorders.at( i );
//        K = borders.at( i );
////        Rcpp::Rcout << "iBegin: " << iBegin << ", iEnd: " << iEnd << std::endl;
////        Rcpp::Rcout << "k: " << k << ", K: " << K << std::endl;
////        Rcpp::Rcout << "CLU: " << clu << std::endl;
////        Rcpp::Rcout << "Begin element: " << *( clu.begin() + iBegin ) << ", End element: " << *( clu.begin() + iEnd - 1 ) << std::endl;
//        for( ; k < K; ++k ) {
////            Rcpp::Rcout << "nBorders: " << nBorders << std::endl;
////            Rcpp::Rcout << "borders: " << borders << std::endl;
////            Rcpp::Rcout << "nClu: " << nClu << std::endl;
////            Rcpp::Rcout << "clu: " << clu << std::endl;
////            Rcpp::Rcout << "countGroups: " << countGroups << std::endl;
////            Rcpp::Rcout << "eVec: " << eVec << std::endl;
////            Rcpp::Rcout << "-----------------------" << std::endl;
//            if( !( std::find( clu.begin() + iBegin, clu.begin() + iEnd, k ) != ( clu.begin() + iEnd ) ) ) {
//                size_t i = std::distance( eVec.begin() + iBegin, std::max_element( eVec.begin() + iBegin, eVec.begin() + iEnd ) );
////                Rcpp::Rcout << "INSIDE IF i="<< i << ", countGroupst(clu(i))=" << countGroups.at( clu.at( i ) ) << std::endl;
////                ce ma countgroups[i] samo 1 skupino, zberem naslednji max iz drugih skupin - torej ce je k 1 - 3 in ima skupina 1 samo 1 skupino izberem max med skupinama 2 in 3
//                if( countGroups.at( clu.at( i ) ) < 2 ) {
//                    k--;
//                    eVec.at( i ) = -1;
//                    continue;
//                }
//                clu.at( i ) = k;
//                eVec.at( i ) = 0;
//                k = i == 0 ? 0 : borders.size() > 1 ? borders.at( i - 1) : 0;
//            }
//        }
//    }
//
//    K = 0;
//    int k = 0;
//    for( int i = 0; i < borders.size(); ++i ) {
//        if( i ) {
//            k = borders[ i - 1 ];
//            K = borders[ i ];
//        }
//        else {
//            k = 0;
//            K = borders[ i ];
//        }
//        Rcpp::Rcout << "for k in " << k << ":" << K << std::endl;
//        for( ; k < K; ++k ) {
//            if( !( std::find( clu.begin(), clu.end(), k ) != clu.end() ) && countGroups.at( k ) > 1 ) { // if k is not in vRet
//                size_t j = std::distance( eVec.begin(), std::max_element( eVec.begin(), eVec.end() ) );
//                clu.at( j ) = k;
//                eVec.at( j ) = 0;
//                k = i == 0 ? 0 : borders[ i - 1]; // if i == 0 then k = 0 else k = borders[ i - 1 ]
//            }
//        }
//    }

//    for( int k = 0; k < K; ++k ) {
//        if( !( std::find( clu.begin(), clu.end(), k ) != clu.end() ) ) { // if k is not in vRet
//            size_t i = std::distance( eVec.begin(), std::max_element( eVec.begin(), eVec.end() ) );
//            clu.at( i ) = k;
//            eVec.at( i ) = 0;
//            k = 0;
//        }
//    }
}

unsigned int belongsTo( const int & group, const IVector & borders )
{
    for( int i = 0; i < borders.size(); ++i ) {
        if( group < borders.at( i ) ) {
            return i;
        }
    }
    return Rcpp::sum( borders );
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
