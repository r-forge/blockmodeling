#include <memory>
#include <float.h>
#include <algorithm>
#include <RcppArmadillo.h>

using DVector = Rcpp::NumericVector;
using IVector = Rcpp::IntegerVector;
using Array   = arma::Cube<double>; // [ Rows, Cols, Slices ]
using DMatrix = arma::Mat<double>;

enum Diagonale
{
    Same     = 0,
    Ignore   = 1,
    Seperate = 2
};

class MeanObject
{
public:

    enum Type
    {
        Vector = 0,
        Matrix = 1
    };


    MeanObject( const Array & p_array, const IVector & p_n );
    MeanObject( const Array & p_array );

    Type type() const;
    double mean( const int r ) const;
    double mean( const int i, const int r ) const;

    IVector n() const;


protected:
    Type m_type;
    DVector m_vMeans;

    DMatrix m_mMeans;
    IVector m_vN;

};

std::ostream & operator << ( std::ostream & stream, const MeanObject::Type & type );
