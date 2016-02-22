#include <iostream>
#include <tuple>
#include "CMatrixD.h"
#define DEBUG_PRO 1
using Ans = std::tuple<CMatrixD, CMatrixD, uint>;

auto equation(const CMatrixD &matx)
{
    CMatrixD mat(matx);                                    ///<增广矩阵
    uint equationNum = mat.rowSize();                      ///<方程的个数
    uint xNum = mat.colSize() - 1;                         ///<未知数个数
    uint rank = xNum < equationNum ? xNum : equationNum;   ///<方程的秩

    CMatrixD coefficientMatrix(equationNum, xNum);
    for (uint i = 1; i <= equationNum; i++)
        for (uint j = 1;j <= xNum; j++)
            coefficientMatrix.at(i,j) = matx.at(i,j);

    auto print = [&mat]
    {
        if (! (DEBUG_PRO))
            return;
        for (uint i = 1; i <= mat.rowSize(); i++)
        {
            for (uint j = 1; j <= mat.colSize(); j++)
            {
                if (mat.colSize() == j)
                {
                    std::cout << "|  ";
                }
                double num = ((int)(mat.at(i, j) * 1000)) / 1000.0;
                std::cout << std::right << std::setw(7) << num << " ";
            }
            std::cout << std::endl;
        }
    };

    auto y = [&mat](int i)
    {
        return mat.at(i, mat.colSize());
    };

    std::cout << "初始方程" << std::endl;
    print();

    for (uint i = 1; i <= rank; i++)
    {
        uint j = i;
        for(j=i; j <= equationNum; j++)
        {
            if (fabs(mat.at(i, j)) > MAT_MIN_ELE)
                break;
        }

        if (j > equationNum)
        {
            //得到真正的秩
            rank = i - 1;
            break;
        }

        if (i != j)
        {
            mat.changeRow(i, j);
            std::cout << "交换" << i << "行和" << j << "行" << std::endl;
            print();
        }

        if (fabs(mat.at(i, i) - 1) > MAT_MIN_ELE)
        {
            std::cout << "第" << i <<"行除以" << mat.at(i, i) <<std::endl;
            mat.multRow(i, 1 / mat.at(i, i));
            print();
        }

        for(j = i + 1; j <= equationNum; j++)
        {
            double multNum = -mat.at(j, i) / mat.at(i, i);
            if (fabs(multNum) < MAT_MIN_ELE)
                continue;
            mat.multPlusRow(j, i, multNum);
            std::cout<< "第" << j <<"行减去第" << i << "行乘以" << -multNum << std::endl;
            print();
        }
    }

    for (uint i = rank + 1; i <= equationNum; i++)
    {
        //判断增广矩阵的秩是否大于系数矩阵的秩。
        if (fabs(y(i)) > MAT_MIN_ELE)
        {
            std::cout << "方程无解。" << std::endl;
            return Ans(coefficientMatrix, CMatrixD(1,1,0), 0);
        }
    }

    for (uint i = rank; i >= 1; i--)
    {
        for (uint j = i -1; j >= 1; j--)
        {
            double multNum = -mat.at(j, i);
            if (fabs(multNum) < MAT_MIN_ELE)
                continue;
            mat.multPlusRow(j, i, multNum);
            std::cout<< "第" << j <<"行减去第" << i << "行乘以" << -multNum << std::endl;
            print();
        }
    }

    if (rank == xNum)
    {
        std::cout << "方程有唯一解:" << std::endl;
        CMatrixD solMat(xNum,1,0);
        for (uint i = 1; i <= rank; i++)
        {
            std::cout << "X" << i << "=" << y(i) << std::endl;
            solMat.at(i,1) = y(i);
        }
        return Ans(coefficientMatrix, solMat, 1);
    }

    std::cout << "方程有无穷解:" << std::endl;

    uint solNum = xNum - rank;                   ///<方程解的维数
    CMatrixD solMat(xNum,solNum + 1,0);

    std::cout << "方程的一个特解:" << std::endl;

    for (uint i = 1; i <= rank; i++)
    {
        std::cout << "X" << i << "=" << (fabs(y(i)) < MAT_MIN_ELE ? 0 : y(i)) << std::endl;
        solMat.at(i,1) = y(i);
    }

    for (uint i = rank + 1; i <= xNum; i++)
    {
        std::cout << "X" << i <<"=" << 0 << std::endl;
    }


    std::cout << "齐次方程的基础解系有" << solNum << "组:" << std::endl;


    for (uint i = rank + 1; i <= xNum; i++)
    {
        uint zu = i - rank + 1;                      ///<第几组解
        for (uint j = rank + 1; j <= xNum; j++)
        {
             solMat.at(j,zu) = (i == j) ? 1 : 0;
        }
        for (uint j = rank; j >= 1; j--)
        {
            for (uint k = j + 1; k <= xNum; k++)
            {
                solMat.at(j,zu) += mat.at(j,k) * solMat.at(k,zu);
            }
            solMat.at(j,zu) = -solMat.at(j,zu);
        }
    }

    for (uint i = 1; i <= xNum; i++)
    {
        std::cout << "X" << i << " ";
        for (uint j = 2; j < solNum + 2; j++)
            std::cout << std::right << std::setw(7) << ((fabs(solMat.at(i,j)) < MAT_MIN_ELE) ? 0 : (solMat.at(i,j))) << " ";
        std::cout<<std::endl;
    }

    return Ans(coefficientMatrix, solMat, 2);
}

void checkout(const Ans &ans)
{

    if (0 == std::get<2>(ans))
        return;
    CUTTING_LINE;
    CUTTING_LINE;
    CMatrixD coefficientMatrix(std::get<0>(ans));
    CMatrixD solMat(std::get<1>(ans));
    coefficientMatrix.print("系数矩阵");
    CUTTING_LINE;
    solMat.print("解矩阵");
    CUTTING_LINE;
    CMatrixD y =coefficientMatrix * solMat;
    y.print("乘积");
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

    //CMatrixD mat4(100,101,0);

   // mat4.rand(0,20);

    auto ans = equation(mat2);

    checkout(ans);


    return 0;
}

