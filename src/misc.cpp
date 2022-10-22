#include <cpp11.hpp>
#include <vector>
using namespace std ;

[[cpp11::register]]
cpp11::doubles_matrix<> collapseData(cpp11::list mat_list, cpp11::strings colnames){
	//create the mat
	int nCol = colnames.size();
	int nRow =0;
	for(int i = 0; i < mat_list.size(); i++){

		cpp11::doubles_matrix<> mat(mat_list[i]);
		nRow += mat.nrow();
	}
    cpp11::writable::doubles_matrix<> out(nRow, nCol);

    cpp11::writable::list mydims(2);
	mydims[1] = colnames;
	Rf_setAttrib(cpp11::as_sexp(out), cpp11::as_sexp({"dimnames"}), cpp11::as_sexp(mydims));

    //concatenate the mats
    int offset = 0;
    for(int ind = 0; ind < mat_list.size(); ind++){
//    	Rcout << "mat.no:" << ind << endl;
    	//get current mat
    	cpp11::doubles_matrix<> mat(mat_list[ind]);
    	int nrow = mat.nrow();
//    	Rcout << "mat rows:" << nrow << endl;
    	//update the current block
    	for(int i = 0; i < nrow; i++)
    		for(int j = 0; j < nCol; j++){
    			int rowInd = offset + i;
//    			Rcout << "update i,j:" << rowInd << "," << j << endl;
    			out(rowInd, j) =  mat(i, j);
    		}

    	//update the offset
    	offset += nrow;
    }
    return(out);
}

