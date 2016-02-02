/**
 * @file	         CMatrix.h
 * @brief          实现了矩阵的基本功能
 * @author: 	  zhangxinhao
 * @date: 		  2016-01-27
 */

#ifndef CMATRIXD_H
#define CMATRIXD_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <math.h>



#define GET_MIN(x,y)  ((x)<(y)?(x):(y))
#define GET_MAX(x,y)  ((x)>(y)?(x):(y))
#define MAT_MIN_ELE  0.000001

#define CUTTING_LINE std::cout<<"------------------------------------------------------------------------"<<std::endl

#define MAT_DEBUG(x) std::cout<< x <<std::endl
using std::string;

class MatException :public std::exception{
public:
    MatException( const string& rhs ): _err( rhs ){ }
    virtual const char* what( void ) const throw(){
        return _err.c_str();
    }
    ~MatException() throw(){ }
private:
    std::string _err;
};


class CMatrixD
{
public:
    CMatrixD(uint rowSize, uint colSize, double num = 0);

    CMatrixD(const CMatrixD& mat) = default;

    CMatrixD(CMatrixD&& mat);

    CMatrixD(std::initializer_list<std::initializer_list<double>> list2);

    void load(const std::string& path);

    //保存到文件
    void save(const std::string& path);

    //完全重置
    void reset(uint rowSize, uint colSize, double num);

    //不能改变总大小和元素
    void resize(uint rowSize, uint colSize);

    //生成随机数种子
    static void srand();

    void rand(int min, int max);

    //生成数量矩阵，默认单位矩阵
    void unit(uint size,double num = 1);

    void print(const std::string& info="matrix",uint eachSize = 5);

    CMatrixD tra();

    CMatrixD inv();

    double det();

    //初等变换===========================================
    //交换两行
    bool changeRow(uint rowNum1,uint rowNum2);
    //交换两列
    bool changeCol(uint colNum1,uint colNum2);

    //某一行乘以一个数
    bool multRow(uint rowNum,double num);
    //某一列乘以一个数
    bool multCol(uint colNum,double num);

    //某一行乘以一个数加到另一行上
    bool multPlusRow(uint changedRow, uint row, double num);
    //某一列乘以一个数加到另一列上
    bool multPlusCol(uint changedCol, uint col, double num);
    //初等变换===========================================


    double & at(uint rowNum,uint colNum);

    double  at(uint rowNum,uint colNum) const;

    bool isZero(uint rowNum,uint colNum);

    uint rowSize() const;

    uint colSize() const;

    CMatrixD & operator=(const CMatrixD& mat) = default;

    CMatrixD & operator=(std::initializer_list<std::initializer_list<double>> list2);

    CMatrixD & operator=(CMatrixD&& mat);

    CMatrixD  operator+(const CMatrixD& mat);

    CMatrixD  operator-(const CMatrixD& mat);

    CMatrixD  operator*(const CMatrixD& mat);

    CMatrixD  operator*(double num);


private:
    //移动赋值
    void __move(CMatrixD &mat);
    //判断输入是否在矩阵的行范围
    bool __insideRow(uint rowNum) const;
    //判断输入是否在矩阵的列范围
    bool __insideCol(uint colNum) const;

    void __initList(std::initializer_list<std::initializer_list<double>> list2);

    std::vector<double> _arrEle;
    std::vector<uint> _arrRowNum;
    std::vector<uint> _arrColNum;
    uint _rowSize;
    uint _colSize;
    uint _size;
};



#endif // CMATRIXD_H
