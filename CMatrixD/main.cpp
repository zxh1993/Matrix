#include <iostream>
#include "CMatrixD.h"
using namespace std;

void print(const CMatrixD &mat)
{
    for (uint i = 1; i <= mat.rowSize(); i++)
    {
        for (uint j = 1; j <= mat.colSize(); j++)
        {
            if (mat.colSize() == j)
            {
                std::cout << "|  ";
            }
            double num = ((int)(mat.at(i,j)*1000))/1000.0;
            std::cout << std::right << std::setw(7) << num <<" ";
        }

        std::cout << std::endl;
    }
}

void equation(const CMatrixD &matx)
{
    std::cout << "初始方程" << std::endl;
    CMatrixD mat(matx);
    print(mat);
    uint xNum = mat.colSize() - 1;
    uint rank = xNum < mat.rowSize() ? xNum : mat.rowSize();

    for (uint i = 1; i <= rank; i++)
    {
        uint j = i;
        for(j=i; j <= mat.rowSize(); j++)
        {
            if (fabs(mat.at(j,i))>=MAT_MIN_ELE)
                break;
        }
        if (j > mat.rowSize())
        {
            rank = i - 1;
            break;
        }

        if (i != j)
        {
            mat.changeRow(i,j);
            std::cout << "交换" << i << "行和" << j << "行" << std::endl;
            print(mat);
        }
        if (fabs(mat.at(i,i) - 1) > MAT_MIN_ELE)
        {
            std::cout << "第" << i <<"行除以" << mat.at(i,i) <<std::endl;
            mat.multRow(i,1/mat.at(i,i));
            print(mat);
        }
        for(j = i + 1; j <= mat.rowSize(); j++)
        {
            double multNum = -mat.at(j,i)/mat.at(i,i);
            if (fabs(multNum) < MAT_MIN_ELE)
                continue;
            mat.multPlusRow(j,i,multNum);
            std::cout<< "第" << j <<"行减去第" << i << "行乘以" << -multNum << std::endl;
            print(mat);
        }

    }

    for (uint i = rank + 1; i <= mat.rowSize(); i++)
    {
        if (fabs(mat.at(i, mat.colSize())) > MAT_MIN_ELE)
        {
            std::cout << "方程无解。" << std::endl;
            return;
        }
    }



    for (uint i = rank; i >= 1; i--)
    {
        for (uint j = i -1; j >= 1; j--)
        {
            double multNum = -mat.at(j,i);
            if (fabs(multNum) < MAT_MIN_ELE)
                continue;
            mat.multPlusRow(j,i,multNum);
            std::cout<< "第" << j <<"行减去第" << i << "行乘以" << -multNum << std::endl;
            print(mat);
        }
    }

    if (rank == xNum)
    {
        std::cout << "方程有唯一解:" << std::endl;
        for (uint i = 1; i <= rank; i++)
        {
            std::cout<<"X"<<i<<"="<<mat.at(i,mat.colSize()) << std::endl;
        }
        return;
    }

    std::cout << "方程有无穷解:" << std::endl;

    std::cout << "方程的一个特解:" << std::endl;

    for (uint i = 1; i <= rank; i++)
    {
        std::cout<<"X"<<i<<"="<<mat.at(i,mat.colSize()) << std::endl;
    }

    for (uint i = rank + 1; i <= xNum; i++)
    {
        std::cout << "X" << i <<"=" << 0 << std::endl;
    }

    std::cout << "齐次方程的基础解系有" << xNum - rank << "组:" << std::endl;

    double xx[xNum - rank][xNum+1];
    for (uint i = 0; i < xNum - rank ; i++)
        for (uint j = 1; j <= xNum; j++)
            xx[i][j] = 0;
    for (uint i = rank + 1; i <= xNum; i++)
    {
        uint zu = i - rank -1;
        for (uint j = rank + 1; j <= xNum; j++)
        {
            xx[zu][j] = (i == j) ? 1 : 0;
        }
        for (uint j = rank; j >= 1; j--)
        {
            xx[zu][j] = 0;
            for (uint k = j + 1; k <= xNum; k++)
            {
                xx[zu][j] += mat.at(j,k)*xx[zu][k];
            }
            xx[zu][j] = -xx[zu][j];
        }
    }

    for (uint i = 1; i <= xNum; i++)
    {

        std::cout<<"X"<<i<<" ";
        for (uint j = 0; j < xNum -rank; j++)
            std::cout << std::right << std::setw(7) <<( (fabs(xx[j][i]) < MAT_MIN_ELE) ? 0 : (xx[j][i]) ) << " ";
        std::cout<<std::endl;

    }

}


int main()
{


    CMatrixD mat1{{1,2,3,10},
                  {2,1,7,8},
                  {2,4,5,2}};//唯一解

    CMatrixD mat2{{1,2,3,10},
                  {2,1,7,8},
                  {3,3,10,2}};//无解

    CMatrixD mat3{{1,2,3,1,1,4},
                 {4,5,4,2,4,3},
                 {7,9,9,3,5,10}};//无穷解

    equation(mat3);

    return 0;
}

