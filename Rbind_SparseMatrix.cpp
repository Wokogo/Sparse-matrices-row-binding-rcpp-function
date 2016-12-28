
    // [[Rcpp::plugins(cpp11)]]
    #include <math.h>
    #include <RcppArmadillo.h>
    using namespace Rcpp;
    
    // [[Rcpp::depends(RcppArmadillo)]]
    // [[Rcpp::export]]
    arma::sp_mat rbindSpCpp(arma::sp_mat &a, arma::sp_mat &b)
    { // Function for row-binding two sparse matrices

    // Make const iterator for both a and b matrices
    arma::sp_mat::const_iterator starta = a.begin();
    arma::sp_mat::const_iterator enda = a.end();
    arma::sp_mat::const_iterator startb = b.begin();
    arma::sp_mat::const_iterator endb = b.end();
    
    // Counting the number of points
    int Na = std::distance(starta, enda); // Number of non zero elements in a matrix
    int Nb = std::distance(startb, endb); // Number of non zero elements in b matrix
    
    // Make sure that constraints are respected
    if(Na <= 0) stop("the first matrix is a full sparse one. No values found! Fill it please");
    if(Nb <= 0) stop("the second matrix is a full sparse one. No values found! Fill it please");
    if(a.n_cols != b.n_cols) stop("Matrices mismatch their number of column !");
    
    // Construct locations matrix and values vector for the resulting (final) matrix:
    arma::umat locs(Na+Nb, 2); // Build a location storage matrix
    arma::urowvec tempa(2); arma::urowvec tempb(2); // Create vectors to store each row information in. (Row, Col)
    arma::vec rbindValues(Na+Nb); // Build a values storage vector for the resulting matrix.
    
    // Start collecting information: values and its locations 
    arma::sp_mat::const_iterator ita = starta; // Iterator for collecting information from a matrix.
    arma::sp_mat::const_iterator itb = startb; // Iterator for collecting information from b matrix.
    
    // Two needed tests
    int test = 1; // test = 1 when collecting information from a matrix, and test = -1 when collecting from b.
    int endaTest = 0; // endaTest = 1 means that the latest non zero element is reached for a matrix.

    arma::urowvec tempOld(2); // Will be used for storing the column of the previous non zero element of the current matrix.

    //Collecting locations and values
    for(int i = 0; i < Na+Nb; ++i)
    {
        if(test == 1)
        { // Focus on the a matrix.
            tempa(0) = ita.row(); tempa(1) = ita.col(); // Current non zero element coordinates.
            locs.row(i) = tempa; rbindValues(i) = *ita; // *ita is the non zero value of the current point.
            if(*ita == nonzeros(a)[Na-1])
            { // Check if the current element is the last one.
                endaTest = endaTest + 1; test = -1*test;
            }
            else
            {
                tempOld(1) = tempa(1); ita++; tempa(1) = ita.col(); 
                if(tempa(1) != tempOld(1)){test = -1*test;} // If the column changes, jumps in b matrix.
            }
        }
        else
        {
            if(test == -1)
            { // focus on the b matrix
                tempb(0) = itb.row() + a.n_rows; tempb(1) = itb.col(); // Current non zero element coordinates in the final matrix.
                locs.row(i) = tempb; rbindValues(i) = *itb; // *itb is the non zero value of the current point.
                tempOld(1) = tempb(1);
                itb++; tempb(1) = itb.col();
                if(tempb(1) != tempOld(1))
                {
                    if(endaTest == 1){test = 1*test;} else {test = -1*test;} // If endaTest = 1 stays in b matrix, else jumps in a matrix.
                }
            }
        }
    }
    // Final sparse matrix construction
    arma::sp_mat M(trans(locs), rbindValues, a.n_rows + b.n_rows, a.n_cols);
    return M;
    }





