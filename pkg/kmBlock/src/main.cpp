#include <iostream>
#include <memory>
#include <string.h>

#include "kmBlock_types.h"


template < typename T >
class Borders {

public:

    Borders(){}
    Borders( const size_t rows, const size_t cols, const size_t slices = -1 );
    Borders( const T & p_lower, const T & p_upper ) : m_lower( p_lower ), m_upper( p_upper ) {}
    virtual ~Borders() {}

    T getLower() const { return m_lower; }
    T getUpper() const { return m_upper; }

    double getLowerAt( const size_t i, const size_t j, const size_t r = -1 ) const { return m_lower.at( i, j ); }
    double getUpperAt( const size_t i, const size_t j, const size_t r = -1 ) const { return m_upper.at( i, j ); }

protected:

    T m_lower;
    T m_upper;

};

template < typename T >
Borders<T>::Borders( const size_t rows, const size_t cols, const size_t /* slices */ )
{
    m_lower = T( rows, cols );
    m_upper = T( rows, cols );

    m_lower.fill( R_NegInf );
    m_upper.fill( R_PosInf );
}

template <>
Borders<Array>::Borders( const size_t rows, const size_t cols, const size_t slices )
{
    m_lower = Array( rows, cols, slices );
    m_upper = Array( rows, cols, slices );

    m_lower.fill( R_NegInf );
    m_upper.fill( R_PosInf );
}

template <>
double Borders<Array>::getLowerAt( const size_t i, const size_t j, const size_t r ) const
{
    return m_lower.at( i, j, r );
}

template <>
double Borders<Array>::getUpperAt( const size_t i, const size_t j, const size_t r ) const
{
    return m_upper.at( i, j, r );
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Exposed functions
// [[Rcpp::export]]
Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const int dimensions, const IVector & n, const std::string diagonal = "ignore",
                         const bool useBorders = false, const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                         const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue );
// [[Rcpp::export]]
Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const IVector & n, const IVector & nClu, const std::string & diagonal = "ignore",
                    const bool useBorders = false, const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue );
// [[Rcpp::export]]
double critFunction( const Array & M, const IVector & clu, const Array & weights, const int dimensions, const IVector & n, const std::string & diagonal = "ignore",
                     const bool useBorders = false, const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                     const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue );

// Functions forward declarations
void meansByBlocks( const Array & M, Array & res, const IVector & clu, const int dimensions, DMatrix & p_pSepare, const DMatrix & p_mMeans, const IVector & n,
                    const Borders<Array> & p_btBorders, const Borders<DMatrix> & p_btBordersSeperate, const Diagonale sDiagonal = Diagonale::Ignore );
double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat, const DMatrix & p_mSeparate, const Diagonale p_diagonale );
void setGroups( const Array & M, IVector & clu, const Array & weights, const Array & meansMat, const IVector & nClu, const IVector & n, const DMatrix & p_mSeparate, const Diagonale p_diagonale );
unsigned int belongsTo( const int & group, const IVector & borders );
DMatrix relationsMeans( const Array & M, const IVector & n );
double meanMatrix( const DMatrix & p_matrix );
Diagonale getDiagonale( const std::string & p_sDiagonal );
void checkInputBorders( const Diagonale & p_diagonale,
                    const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper );

std::ostream & operator << ( std::ostream & p_stream, const Diagonale p_diag );

Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const IVector & n, const IVector & nClu, const std::string & diagonal,
                    const bool useBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper )
{
    const Diagonale dDiag = getDiagonale( diagonal );
    const int K = Rcpp::sum( nClu );
    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;

    if( useBorders ) {
        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( K, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( K, K, M.n_slices );
    }


    DMatrix pSeparate;
//    Array pSeparate;
    Array meanBlocks, mSeparate;
//    Rcpp::Rcout << "DEBUG1" << std::endl;

    const DMatrix MEANS = relationsMeans( M, n );
    meansByBlocks( M, meanBlocks, clu, K, pSeparate, MEANS, n, bordersMeanstMat, bordersSeperate, dDiag );
//    Rcpp::Rcout << "DEBUG2" << std::endl << pSeparate << std::endl;

//    Rcpp::Rcout << bordersMeanstMat.getLower() << std::endl << bordersMeanstMat.getUpper() << std::endl;

//    IVector newClu = setGroups( M, clu, weights, meanBlocks, K );
    IVector newClu = Rcpp::clone( clu );
    setGroups( M, newClu, weights, meanBlocks, nClu, n, pSeparate, dDiag );
    IVector bestClu;

//    Rcpp::stop( "DEBUG" );

    meansByBlocks( M, meanBlocks, newClu, K, pSeparate, MEANS, n, bordersMeanstMat, bordersSeperate, dDiag );
    double newCf = criterialFunction( M, clu, weights, meanBlocks, pSeparate, dDiag );
    double bestCf = DBL_MAX;

    while( newCf < bestCf ) {
        bestClu = newClu;
        bestCf = newCf;
        setGroups( M, newClu, weights, meanBlocks, nClu, n, pSeparate, dDiag );
        meansByBlocks( M, meanBlocks, newClu, K, pSeparate, MEANS, n, bordersMeanstMat, bordersSeperate, dDiag );
        newCf = criterialFunction( M, newClu, weights, meanBlocks, pSeparate, dDiag );
    }

//    Rcpp::Rcout << "BestCf: " << bestCf << std::endl;
//    Rcpp::Rcout << "BestClu: " << bestClu << std::endl;

//    if( !pSeparate.is_empty() ) {
//        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = meanBlocks, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
//    }

//    DMatrix d = relationsMeans( M, n );
//    Rcpp::Rcout << "Relations means matrix:" << std::endl << "-----------------" << std::endl << std::endl;
//    Rcpp::Rcout << d << std::endl << "Dimensions: " << d.n_rows << " " << d.n_cols << std::endl;
//    Rcpp::Rcout << "-----------------" << std::endl << std::endl;

    return Rcpp::List::create( Rcpp::Named( "bestCf" ) = bestCf, Rcpp::Named( "bestClu" ) = bestClu, Rcpp::Named( "IM" ) = meanBlocks );

}

Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const int dimensions, const IVector & n, const std::string diagonal,
                         const bool useBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                         const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper )
{
    Diagonale dDiag = getDiagonale( diagonal );

    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;

    if( useBorders ) {
        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
    }

    DMatrix pSeparate;
    Array aRes;
    const DMatrix MEANS = relationsMeans( M, n );
    meansByBlocks( M, aRes, clu, dimensions, pSeparate, MEANS, n, bordersMeanstMat, bordersSeperate, dDiag );
    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }
    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

double critFunction( const Array & M, const IVector & clu, const Array & weights, const int dimensions, const IVector & n, const std::string & diagonal,
                     const bool useBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                     const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper )
{
    Diagonale dDiagonale = getDiagonale( diagonal );

    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;

    if( useBorders ) {
        checkInputBorders( dDiagonale, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiagonale == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiagonale == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
    }

    DMatrix pSeparate;
    Array aRes;
    const DMatrix MEANS = relationsMeans( M, n );
//    Borders<Array> bordersMeanstMat( dimensions, dimensions, M.n_slices );
//    Borders<DMatrix> bordersSeperate( dimensions, M.n_slices );
    meansByBlocks( M, aRes, clu, dimensions, pSeparate, MEANS, n, bordersMeanstMat, bordersSeperate, dDiagonale );


    return criterialFunction( M, clu, weights, aRes, pSeparate, dDiagonale );
}

void meansByBlocks( const Array & M, Array & res, const IVector & clu, const int dimensions, DMatrix & p_pSepare, const DMatrix & p_mMeans, const IVector & n, const Borders<Array> & p_btBorders, const Borders<DMatrix> & p_btBordersSeperate, const Diagonale sDiagonal )
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
//        mDiagonalRes = DMatrix( dimensions, M.n_slices, arma::fill::zeros );
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
        int s = 0;
        int sCount = 0;
        for( size_t i = 0; i < res.n_rows; ++i ) {
            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
                double dVal = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
                if( dVal < p_btBordersSeperate.getLowerAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getLowerAt( i, r );
                if( dVal > p_btBordersSeperate.getUpperAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getUpperAt( i, r );
                else mDiagonalRes( i, r ) = dVal;
            }
            if( sCount >= n.at( s ) ) {
                ++s;
                sCount = 0;
            }
            ++sCount;
            for( size_t j = 0; j < res.n_cols; ++j ) {
//                if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
//                    mDiagonalRes( i, r ) = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
//                }
                double dVal( S( i, j, r ) );
                if( i == j && N( i, j, r ) == 0 && ((sDiagonal == Diagonale::Ignore)||(sDiagonal == Diagonale::Seperate)) ) { // If value of the block is 0 and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
//                    diagMean = p_mMeans.at( s, r ); // p_mMeans.at( getS( i, n ), r );
                    res(i, j, r ) = p_mMeans.at( s, r );
                }
                else {
                    res( i, j, r ) = dVal / N( i, j, r );
                }

                if( res.at( i, j, r ) < p_btBorders.getLowerAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getLowerAt( i, j, r );
                else if( res.at( i, j, r ) > p_btBorders.getUpperAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getUpperAt( i, j, r );
            }
        }
    }

    if( sDiagonal == Diagonale::Seperate ) { // Save seperate digaonal values to input parameter
        p_pSepare = std::move( mDiagonalRes );
    }
//    Rcpp::Rcout << "meansByBlocks: end" << std::endl;
}

double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const Array & meansMat, const DMatrix & p_mSeparate, const Diagonale p_diagonale )
{
//    Rcpp::Rcout << "blala " << std::endl << p_btBordersSeperate.getLower() << std::endl << p_btBordersSeperate.getUpper() << std::endl;
    double dRet = 0;
    for( size_t i = 0; i < M.n_rows; ++i ) {
        for( size_t j = 0; j < M.n_cols; ++j ) {
            for( size_t r = 0; r < M.n_slices; ++r ) {
                if( p_diagonale == Diagonale::Ignore && i == j ) { // ignore diagonal
                    continue;
                }
                else if( p_diagonale == Diagonale::Seperate && i == j ) {
                    double dAvg = p_mSeparate( clu.at( i ), r );
//                    if( dAvg < p_btBordersSeperate.getLowerAt( clu.at( i ), r ) ) dAvg = p_btBordersSeperate.getLowerAt( clu.at( i ), r );
//                    else if( dAvg > p_btBordersSeperate.getUpperAt( clu.at( i ), r ) ) dAvg = p_btBordersSeperate.getUpperAt( clu.at( i ), r );
                    dRet += weights( i, j, r ) * std::pow( M( i, j, r ) - dAvg, 2 );
                }
                else {
                    double dAvg = meansMat( clu.at( i ), clu.at( j ), r );
//                    if( dAvg < p_btBordersMeansMat.getLowerAt( clu.at( i ), clu.at( j ), r ) ) dAvg = p_btBordersMeansMat.getLowerAt( clu.at( i ), clu.at( j ), r );
//                    else if( dAvg > p_btBordersMeansMat.getUpperAt( clu.at( i ), clu.at( j ), r ) ) dAvg = p_btBordersMeansMat.getUpperAt( clu.at( i ), clu.at( j ), r );
                    dRet += weights( i, j, r ) * std::pow( M( i, j, r ) - dAvg, 2 );
                }
            }
        }
    }

    return  dRet;
}

void setGroups( const Array & M, IVector & clu, const Array & weights, const Array & meansMat, const IVector & nClu, const IVector & n, const DMatrix & p_mSeparate, const Diagonale p_diagonale )
{
    IVector borders = Rcpp::cumsum( nClu );
//    Rcpp::Rcout << "Cumsum: " << borders << std::endl;
    double eMin;
    double eTmp;
    int K = Rcpp::sum( nClu );
    int kMin( 0 );
    DVector eVec( clu.size() );
    IVector countGroups( K );
//    IVector vRet = Rcpp::clone( clu );
    DVector e( K );
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
                if( i == j ){
					if( p_diagonale == Diagonale::Ignore) { // ignore diagonal AZ- pogoj je potrebno spremeniti tako, da se se ne "izvede", če je Diagonale:Same
						// AZ Če je Diagonale:Seperate, potem je tu potrebno narediti pravzaprav spodnjo for zanko, le da namesto meansMat za primerjavo uporabite tisto, kar ste v meansByBlocks izračunali kot mDiagonalRes
						continue; 
					}
                    else if( p_diagonale == Diagonale::Seperate ) {
						for( unsigned int r = 0; r < M.n_slices; ++r ) {
	//                        Rcpp::Rcout << "k=" << k << ", cluJ=" << cluJ << ", r=" << r << std::endl;
                            double dAvg = p_mSeparate.at( k, r );
//                            if( dAvg < p_btBordersSeperate.getLowerAt( k, r ) ) dAvg = p_btBordersSeperate.getLowerAt( k, r );
//                            else if ( dAvg > p_btBordersSeperate.getUpperAt( k, r ) ) dAvg = p_btBordersSeperate.getUpperAt( k, r );
                            eTmp += weights( i, j, r ) * std::pow( M( i, j, r ) - dAvg, 2 );
							// eTmp += weights( j, i, r ) * std::pow( M( j, i, r ) - p_mSeparate.at( k, r ), 2 );
							//Pri vrenosti na diagonali je i==j in to vrednost tako kot  vse ostale celice gledamo le enkrat
						}
					}
                    else if( p_diagonale == Diagonale::Same ) {
						for( unsigned int r = 0; r < M.n_slices; ++r ) {
	//                        Rcpp::Rcout << "k=" << k << ", cluJ=" << cluJ << ", r=" << r << std::endl;
                            double dAvg = meansMat( k, k, r );
//                            if( dAvg < p_btBordersMeansMat.getLowerAt( k, k, r ) ) dAvg = p_btBordersMeansMat.getLowerAt( k, k, r );
//                            else if( dAvg > p_btBordersMeansMat.getUpperAt( k, k, r ) ) dAvg = p_btBordersMeansMat.getUpperAt( k, k, r );
                            eTmp += weights( i, j, r ) * std::pow( M( i, j, r ) - dAvg, 2 );
							// eTmp += weights( j, i, r ) * std::pow( M( j, i, r ) - p_mSeparate.at( k, r ), 2 );
							//Pri vrenosti na diagonali je i==j in to vrednost tako kot  vse ostale celice gledamo le enkrat
						}
					}
				}
                else {
                    for( unsigned int r = 0; r < M.n_slices; ++r ) {
                        double dAvg = meansMat( k, cluJ, r );
//                        if( dAvg < p_btBordersMeansMat.getLowerAt( k, cluJ, r ) ) dAvg = p_btBordersMeansMat.getLowerAt( k, cluJ, r );
//                        else if( dAvg > p_btBordersMeansMat.getUpperAt( k, cluJ, r ) ) dAvg = p_btBordersMeansMat.getUpperAt( k, cluJ, r );
                        eTmp += weights( i, j, r ) * std::pow( M( i, j, r ) - dAvg, 2 );

                        dAvg = meansMat( cluJ, k, r );
//                        if( dAvg < p_btBordersMeansMat.getLowerAt( cluJ, k, r ) ) dAvg = p_btBordersMeansMat.getLowerAt(cluJ, k, r );
//                        else if( dAvg > p_btBordersMeansMat.getUpperAt( cluJ, k, r ) ) dAvg = p_btBordersMeansMat.getUpperAt( cluJ, k, r );
                        eTmp += weights( j, i, r ) * std::pow( M( j, i, r ) - dAvg, 2 );
                    }
                }
            }
            if ( eTmp < eMin ){
                kMin = k;
                eMin = eTmp;
            }

        }
        clu.at( i ) = kMin;
        eVec.at( i ) = eMin;
        ++countGroups.at( kMin );
    }

    IVector nBorders( Rcpp::cumsum( n ) );

//    Rcpp::Rcout << "nBorders: " << nBorders << std::endl;
//    Rcpp::Rcout << "borders: " << borders << std::endl;
//    Rcpp::Rcout << "nClu: " << nClu << std::endl;
//    Rcpp::Rcout << "clu: " << clu << std::endl;
//    Rcpp::Rcout << "eVec: " << eVec << std::endl;
//    Rcpp::Rcout << "-----------------------" << std::endl;



    int iBegin, iEnd, k;
    for(int i = 0; i < nBorders.size(); ++i ){
        if( !i ) {
            iBegin = 0;
            k = 0;
        }
        else {
            iBegin = nBorders.at( i - 1 );
            k = borders.at( i - 1 );
        }
        iEnd = nBorders.at( i );
        K = borders.at( i );
//        Rcpp::Rcout << "iBegin: " << iBegin << ", iEnd: " << iEnd << std::endl;
//        Rcpp::Rcout << "k: " << k << ", K: " << K << std::endl;
//        Rcpp::Rcout << "CLU: " << clu << std::endl;
//        Rcpp::Rcout << "Begin element: " << *( clu.begin() + iBegin ) << ", End element: " << *( clu.begin() + iEnd - 1 ) << std::endl;
        for( ; k < K; ++k ) {
//            Rcpp::Rcout << "begin k=" << k << std::endl;
//            Rcpp::Rcout << "nBorders: " << nBorders << std::endl;
//            Rcpp::Rcout << "borders: " << borders << std::endl;
//            Rcpp::Rcout << "nClu: " << nClu << std::endl;
//            Rcpp::Rcout << "clu: " << clu << std::endl;
//            Rcpp::Rcout << "countGroups: " << countGroups << std::endl;
//            Rcpp::Rcout << "eVec: " << eVec << std::endl;
            if( !( std::find( clu.begin() + iBegin, clu.begin() + iEnd, k ) != ( clu.begin() + iEnd ) ) ) {
//                Rcpp::Rcout << "Emtpy Group: " << k << std::endl;
                size_t g = std::distance( eVec.begin(), std::max_element( eVec.begin() + iBegin, eVec.begin() + iEnd ) );
//                Rcpp::Rcout << "INSIDE IF g="<< g << ", countGroupst(clu(g))=" << countGroups.at( clu.at( g ) ) << std::endl;
//                ce ma countgroups[i] samo 1 skupino, zberem naslednji max iz drugih skupin - torej ce je k 1 - 3 in ima skupina 1 samo 1 skupino izberem max med skupinama 2 in 3
                if( countGroups.at( clu.at(g ) ) < 2 ) {
                    k--;
                    eVec.at( g ) = -1;
                    continue;
                }
                clu.at( g ) = k;
                eVec.at( g ) = 0;
                k = i == 0 ? 0 : borders.size() > 1 ? borders.at( i - 1 ) : 0;
                --k;
//                Rcpp::Rcout << "i=" << i << ", k=" << k << std::endl;
            }
//            Rcpp::Rcout << "-----------------------" << std::endl;
        }
    }


//    int iBegin, iEnd, kBegin, kEnd, iMax;
//	double eMax;
//    for( int s = 0; s < nBorders.size(); ++s ) {
//        if( !s ) {
//            iBegin = 0;
//			kBegin = 0;
//        } else {
//            iBegin = nBorders.at( s - 1 );
//            kBegin = borders.at( s - 1 );
//        }

//		iEnd = nBorders.at( s);
//		kEnd = borders.at( s);
//		for(int k = kBegin; k < kEnd; k++){
//			if(countGroups.at(k)==0){
//				eMax=-1;
//				iMax=iBegin;
//				for(int i = iBegin; i< iEnd; ++i){
////					Rcpp::Rcout << " clu.at(i) " << clu.at(i) << std::endl;
//					if(eVec.at(i)>eMax && countGroups.at(clu.at(i))>1){
////						Rcpp::Rcout << "i " << i << std::endl;
//						eMax=eVec.at(i);
//						iMax=i;
//					}
//				}
////				Rcpp::Rcout << " old k " << clu.at(iMax) << " new k " << k << std::endl;
//				-- countGroups.at(clu.at(iMax));
//				clu.at(iMax)=k;
//				countGroups.at(k)=1;
//				eVec.at(iMax)=0; //this is not really needed
//			}
//		}
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


DMatrix relationsMeans( const Array & M, const IVector & n )
{
    const int S = n.size();
    DMatrix mMeans( S, M.n_slices, arma::fill::zeros );
    IVector cumN = Rcpp::cumsum( n );
    cumN.push_front( 0 );
    for( int s = 0; s < S; ++s ) {
        for( unsigned int r = 0; r < M.n_slices; ++r ) {
            Array subArray = M.subcube( cumN.at( s ), cumN.at( s ), r, cumN.at( s + 1 ) - 1, cumN.at( s + 1 ) - 1, r );
            mMeans.at( s, r ) = meanMatrix( subArray.slice( 0 ) );
        }
    }

    return mMeans;
}


double meanMatrix( const DMatrix & p_matrix )
{
    size_t sElementsN = 0;
    size_t sum = 0;
    for( size_t i = 0; i < static_cast<size_t>( p_matrix.n_rows ); ++i ) {
        for( size_t j = 0; j < static_cast<size_t>( p_matrix.n_cols ); ++j ) {
            if( i == j ) continue;
            sum += p_matrix.at( i, j );
            ++sElementsN;
        }
    }
    return static_cast<double>( sum ) / sElementsN;
}

Diagonale getDiagonale( const std::string & p_sDiagonal )
{
    Diagonale dRet = Diagonale::Ignore;
    std::string sTmp = p_sDiagonal;
    std::transform( sTmp.begin(), sTmp.end(), sTmp.begin(), ::tolower );
    if( sTmp == "ignore" );
    else if( sTmp == "same" )
    {
        dRet = Diagonale::Same;
    }
    else if( sTmp == "seperate" )
    {
        dRet = Diagonale::Seperate;
    }
    else Rcpp::stop( "Unknow diagonal parameter\nOptions are: [ same, ignore, seperate ]\n" );

    return dRet;
}

void checkInputBorders( const Diagonale & p_diagonale, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper )
{
    if( bordersMatLower.isNull() ) Rcpp::stop( "Invalid argument: bordersMatLower is null" );
    if( bordersMatUpper.isNull() ) Rcpp::stop( "Invalid argument: bordersMatUpper is null" );
    if( p_diagonale == Diagonale::Seperate ) {
        if( bordersSeperateLower.isNull() ) Rcpp::stop( "Invalid argument: bordersSeperateLower is null" );
        if( bordersSeperateUpper.isNull() ) Rcpp::stop( "Invalid argument: bordersSeperateUpper is null" );
    }
}

std::ostream & operator << ( std::ostream & p_stream, const Diagonale p_diag )
{
    switch ( p_diag ) {
        case Diagonale::Ignore:
            p_stream << "Ignore";
        break;
        case Diagonale::Same:
            p_stream << "Same";
        break;
        case Diagonale::Seperate:
            p_stream << "Seperate";
        break;
    }

    return p_stream;
}
