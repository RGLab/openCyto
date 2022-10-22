#include "Hungarian.h"
#include <vector>
#include <cpp11.hpp>

[[cpp11::register]] std::vector<int> solve_LSAP_cpp(cpp11::doubles_matrix<> mat)
{
    int r = mat.nrow();
    int c = mat.ncol();
    HungarianAlgorithm HungAlgo;
    std::vector<int> assignment;
    //from r's col-major to c's row-major
    std::vector<std::vector<double>> costMatrix(r);
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            costMatrix[i].push_back(mat(i, j));
        }
    }

    HungAlgo.Solve(costMatrix, assignment);
    return assignment;
}