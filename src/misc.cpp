#include <Rcpp.h>
using namespace Rcpp ;
using namespace std ;

// [[Rcpp::export]]
NumericMatrix collapseData(List mat_list, StringVector colnames){
	//create the mat
	unsigned nCol = colnames.size();
	unsigned nRow =0;
	for(unsigned i = 0; i < mat_list.size(); i++){

		NumericMatrix mat = mat_list(i);
		nRow += mat.nrow();
	}
    NumericMatrix out(nRow, nCol);

    out.attr("dimnames") = List(2);
    SET_VECTOR_ELT(out.attr("dimnames"), 1, colnames);
    //concatenate the mats
    int offset = 0;
    for(unsigned ind = 0; ind < mat_list.size(); ind++){
//    	Rcout << "mat.no:" << ind << endl;
    	//get current mat
    	NumericMatrix mat = mat_list(ind);
    	unsigned nrow = mat.nrow();
//    	Rcout << "mat rows:" << nrow << endl;
    	//update the current block
    	for(unsigned i = 0; i < nrow; i++)
    		for(unsigned j = 0; j < nCol; j++){
    			unsigned rowInd = offset + i;
//    			Rcout << "update i,j:" << rowInd << "," << j << endl;
    			out(rowInd, j) =  mat(i, j);
    		}

    	//update the offset
    	offset += nrow;
    }
    return(out);
}

