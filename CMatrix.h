/**
 * @file	         CMatrix.h
 * @brief          实现了矩阵的基本功能
 * @author: 	  zhangxinhao
 * @date: 		  2015-12-16
 */

#ifndef CMATRIX_H
#define CMATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <math.h>
#include <typeinfo>

namespace Mat {

#define GET_MIN(x,y)  ((x)<(y)?(x):(y))
#define GET_MAX(x,y)  ((x)>(y)?(x):(y))
#define MAT_MIN_ELE  0.00001
#define CUTTING_LINE std::cout<<"------------------------------------------------------------------------"<<std::endl
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

template <class Type>
class CMatrix
{
public:
    static Type out;

    typedef enum{
        REPLACE,
        ADD,
        DEL,
        MIN,
        MAX,
        MULT,
    }EMergeType;

    CMatrix(){}

    //初始化矩阵大小以及矩阵的值
    CMatrix(uint size,Type value);

    //输出成double的矩阵
    void toDouble(CMatrix<double> & mat);

    //从文件导入数组，默认行存储
    void load(const std::string& path,uint num = 0);

    //生成随机数种子
    static void srand();

    //生成随机数矩阵
    void rand(Type min,Type max);

    //生成数量矩阵，默认单位矩阵
    void unit(uint size,Type num = 1);

    //保存到文件
    void save(const std::string& path);

    //打印该矩阵
    void print(const std::string& info="matrix",uint size = 5,bool withNum = false);

    //设置总共行数，如果rowSize不是size的约数，异常。
    void setRowSize(uint rowSize);

    //设置总共列数，如果colSize不是size的约数，异常。
    void setColSize(uint colSize);

    //设置矩阵大小
    void setSize(uint size,Type value);

    Type value(uint place)const{ return _array[place-1]; }

    Type value(uint rowNum,uint colNum)const{ return _array[_colSize*(rowNum-1)+colNum-1 ];}

    //访问一维地址(从1开始)，返回引用,超出范围返回out；
    Type& at(uint place);

    //访问二维地址(从(1,1)开始)，返回引用,超出范围返回out；
    Type& at(uint rowNum,uint colNum);

    //计算矩阵的和
    Type sum();

    //获取某一行的和，超出范围返回out；
    Type sumRow(uint rowNum);

    //获取某一列的和，超出范围返回out；
    Type sumCol(uint colNum);

    //取出一个子矩阵,不合理的位置返回异常
    void getSubMat(uint startRow,uint startCol,uint endRow,uint endCol,CMatrix<Type> & mat);

    //获取矩阵某一行,不合理的位置返回异常
    void getRow(uint rowNum,CMatrix<Type>& mat);

    //获取矩阵某一列,不合理的位置返回异常
    void getCol(uint colNum,CMatrix<Type>& mat);

    //矩阵合并
    void merge(uint startRow,uint startCol,const CMatrix<Type> & mat,EMergeType mergeType);

    //矩阵合并
    void merge(const CMatrix<Type> & mat1,const CMatrix<Type> & mat2,EMergeType mergeType);

    //矩阵合并
    void merge(const CMatrix<Type> & mat,Type num,EMergeType mergeType);

    //初等变换===========================================
    //交换两行
    void changeRow(uint rowNum1,uint rowNum2);
    //交换两列
    void changeCol(uint colNum1,uint colNum2);

    //某一行乘以一个数
    void multRow(uint rowNum,Type num);
    //某一列乘以一个数
    void multCol(uint colNum,Type num);

    //某一行乘以一个数加到另一行上
    void multPlusRow(uint changedRow, uint row, Type num);
    //某一列乘以一个数加到另一列上
    void multPlusCol(uint changedCol, uint col, Type num);
    //初等变换===========================================

    //矩阵转置
    void  getTranspose(const CMatrix<Type> & mat);

    //矩阵行列式
    Type det();

    //矩阵的逆
    bool inverse(CMatrix<double>& mat);

    //获取矩阵大小
    uint getSize() const { return _size; }

    //获取矩阵行数
    uint getRowSize() const { return _rowSize; }

    //获取矩阵列数
    uint getColSize() const { return _colSize; }

    //获取该矩阵的vector数组
    std::vector <Type> &getArray() { return _array; }

    //判断输入是否在矩阵的行范围
    bool insideRow(uint rowNum) { return (1<=rowNum && rowNum<=_rowSize); }

    //判断输入是否在矩阵的列范围
    bool insideCol(uint colNum) { return (1<=colNum && colNum<=_colSize); }

    Type &operator () (uint rowNum,uint colNum ) { return at(rowNum,colNum); }

    CMatrix<Type>  operator+(const CMatrix<Type>& mat);

    CMatrix<Type>  operator+(Type num);

    CMatrix<Type>  operator-(const CMatrix<Type>& mat);

    CMatrix<Type>  operator-(Type num);

    CMatrix<Type>  operator*(const CMatrix<Type>& mat);

    CMatrix<Type>  operator*(Type num);

    //转为double类型矩阵
    //CMatrix<double> dou(CMatrix<Type> &mat);

    //计算行列式
    //Type det(CMatrix<Type> &mat)

    //矩阵求逆
    //CMatrix<double> inv(CMatrix<Type> &mat)

    //矩阵转置
    //CMatrix<Type> tra(CMatrix<Type> &mat);

private:
    std::vector <Type> _array;
    uint _rowSize;
    uint _colSize;
    uint _size;
};

template <typename Type>
CMatrix<Type>::CMatrix(uint size, Type value){
    _array.resize(size,value);
    setRowSize(1);
}

template <typename Type>
void CMatrix<Type>::toDouble(CMatrix<double> &mat){
    mat.setSize(_size,0);
    mat.setColSize(_colSize);
    std::vector<double> & vecMat = mat.getArray();
    for(uint i = 0;i<_size;i++)
        vecMat[i] =(double) _array[i];
}

template <typename Type>
void CMatrix<Type>::load(const std::string &path, uint num)
{
    _array.clear();
    std::ifstream in(path.c_str());
    if(!in.is_open())
        throw  MatException("---load---");

    //读取所有的
    Type t;
    if(num == 0){
        while (in>>t){
            _array.push_back(t);
        }
        in.close();
        setRowSize(1);
        return;
    }

    //读取指定个数，个数不能超过文件的
    _array.resize(num);
    for(uint i=0;i<num;i++){
        in>>t;
        _array[i] = t;
    }
    in.close();
    setRowSize(1);
    return;
}

template <typename Type>
void CMatrix<Type>::srand(){
    ::srand((unsigned)time(NULL));
}

template <typename Type>
void CMatrix<Type>::rand(Type min, Type max){
    if (min>max)
        throw MatException("---rand---");
    if (typeid(Type)==typeid(double) || typeid(Type)==typeid(float)) {
        for(uint i = 0; i < _size;i++ ){
            _array[i] = min+(max-min)*::rand() / double(RAND_MAX);
        }
    }
    else{
        for(uint i = 0; i < _size;i++ ){
            _array[i] = (Type)(::rand() %(uint) (max-min+1))+ min;
        }
    }
}

template <typename Type>
void CMatrix<Type>::unit(uint size,Type num){
    setSize(size*size,0);
    setRowSize(size);
    for (uint i=1;i<=_rowSize;i++)
        at(i,i) = num;
}

template <typename Type>
void CMatrix<Type>::save(const std::string &path){
    std::ofstream out(path.c_str());
    for (uint i=0;i<_size;i++){
        out<<_array[i]<<" ";
        if ((i+1) % _colSize == 0)
            out<<std::endl;
    }
    out.close();
}

template <typename Type>
void CMatrix<Type>::print(const std::string& info,uint size,bool withNum){
    uint i;
    std::cout<<info<<":"<<_rowSize<<"*"<<_colSize<<std::endl;
    if (withNum){
        for (i=0;i<=_colSize;i++)
            std::cout<<std::left<<std::setw(size)<<i<<" ";
        std::cout<<std::endl;
        for (i=0;i<_size;i++){
            if ( i % _colSize == 0)
                std::cout<<std::left<<std::setw(size)<<(i / _colSize +1)<<" ";
            std::cout<<std::left<<std::setw(size)<<_array[i]<<" ";
            if ((i+1) % _colSize == 0)
                std::cout<<std::endl;
        }
    }
    else {
        for (i=0;i<_size;i++){
            std::cout<<std::left<<std::setw(size)<<_array[i]<<" ";
            if ((i+1) % _colSize == 0)
                std::cout<<std::endl;
        }
    }
}

template <typename Type>
void CMatrix<Type>::setRowSize(uint rowSize){
    _size = _array.size();
    _rowSize = rowSize;
    _colSize = _size / _rowSize;
    if (_size % _rowSize != 0)
        throw MatException("---setRowSize---");
}

template <typename Type>
void CMatrix<Type>::setColSize(uint colSize){
    _size = _array.size();
    _colSize = colSize;
    _rowSize = _size / _colSize;
    if (_size % _colSize != 0)
        throw MatException("---setColSize---");
}

template <typename Type>
void CMatrix<Type>::setSize(uint size, Type value){
    _array.resize(size,value);
    setRowSize(1);
}


template <typename Type>
Type &CMatrix<Type>::at(uint place){
    if(place<1 || place >_size)
        return out;
    return _array[place-1];
}

template <typename Type>
Type &CMatrix<Type>::at(uint rowNum, uint colNum){
    if (insideCol(colNum) && insideRow(rowNum))
        return _array[_colSize*(rowNum-1)+colNum-1 ];
    else
        return out;
}

template <typename Type>
Type CMatrix<Type>::sum(){
    Type sum = 0;
    for(uint i=0;i<_size;i++){
        sum += _array[i];
    }
    return sum;
}

template <typename Type>
Type CMatrix<Type>::sumRow(uint rowNum){
    if (!insideRow(rowNum))
        return out;
    uint start = _colSize*(rowNum-1);
    Type sum = 0;
    for(uint i=0;i<_colSize;i++)
        sum += _array[i+start ];
    return sum;
}

template <typename Type>
Type CMatrix<Type>::sumCol(uint colNum){
    if (!insideCol(colNum))
        return out;
    Type sum = 0;
    for(uint i = 0;i<_rowSize;i++)
        sum += _array[i*_colSize+colNum-1];
    return sum;
}

template <typename Type>
void CMatrix<Type>::getSubMat(uint startRow, uint startCol,
                              uint endRow, uint endCol,CMatrix<Type> &mat){
    if (! (1<=startCol && startCol<=endCol && endCol<=_colSize
           && 1<=startRow && startRow<=endRow && endCol<=_rowSize))
        throw MatException("---getSubMat---");
    mat.setSize((endCol-startCol+1)*(endRow-startRow+1),0);
    uint i,j,k=0;
    for (i=startRow;i<=endRow;i++)
        for(j=startCol;j<=endCol;j++)
            mat.at(++k) = _array[(i-1)*_colSize+j-1];
    mat.setColSize(endCol-startCol+1);
}

template <typename Type>
void CMatrix<Type>::getRow(uint rowNum,CMatrix<Type> &mat){
    if (!insideRow(rowNum))
        throw MatException("---getRow---");
    mat.setSize(_colSize,0);
    uint start = _colSize*(rowNum-1);
    for(uint i=0;i<_colSize;i++)
        mat.at(i+1) = _array[i+start ];
    mat.setRowSize(1);
}

template <typename Type>
void CMatrix<Type>::getCol(uint colNum,CMatrix<Type> &mat){
    if (!insideCol(colNum))
        throw MatException("---getCol---");
    mat.setSize(_rowSize,0);
    for(uint i = 0;i<_rowSize;i++)
        mat.at(i+1) = _array[i*_colSize+colNum-1];
    mat.setColSize(1);
}

template <typename Type>
void CMatrix<Type>::merge(uint startRow, uint startCol,const CMatrix<Type> &mat,
                          CMatrix<Type>::EMergeType mergeType){
    if (!( insideCol(startCol) && insideRow(startRow) ))
        throw MatException("---merge---");
    uint endRow =startRow + mat.getRowSize() -1;
    uint endCol =startCol + mat.getColSize() -1;
    if (endRow>_rowSize || endCol>_colSize)
        throw MatException("---merge---");
    uint i,j;
    switch (mergeType) {
    case CMatrix<Type>::REPLACE:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                at(i,j) = mat.value(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::ADD:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = this->at(i,j) + mat.value(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::DEL:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = this->at(i,j) - mat.value(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::MIN:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = GET_MIN( this->at(i,j) , mat.value(i-startRow+1,j-startCol+1));
        break;
    case CMatrix<Type>::MAX:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = GET_MAX( this->at(i,j) , mat.value(i-startRow+1,j-startCol+1));
        break;
    default:
        throw MatException("---merge---");
        break;
    }
}

template <typename Type>
void CMatrix<Type>::merge(const CMatrix<Type> &mat1, const CMatrix<Type> &mat2, CMatrix::EMergeType mergeType)
{
    if ( (mat1.getColSize() != mat2.getColSize()) || (mat1.getRowSize() != mat2.getRowSize()) )
        throw MatException("---merge2---");
    uint i,size = mat1.getSize();
    _array.resize(size,0);
    setColSize(mat1.getColSize());
    switch (mergeType) {
    case CMatrix<Type>::ADD:
        for (i=1;i<=size;i++)
            at(i) = mat1.value(i) + mat2.value(i);
        break;
    case CMatrix<Type>::MULT:
        for (i=1;i<=size;i++)
            at(i) = mat1.value(i) * mat2.value(i);
        break;
    case CMatrix<Type>::DEL:
        for (i=1;i<=size;i++)
            _array[i] = mat1.value(i) - mat2.value(i);
        break;
    case CMatrix<Type>::MIN:
        for (i=1;i<=size;i++)
            _array[i] = GET_MIN(mat1.value(i) , mat2.value(i));
        break;
    case CMatrix<Type>::MAX:
        for (i=1;i<=size;i++)
            _array[i] = GET_MAX(mat1.value(i), mat2.value(i));
        break;
    default:
        throw MatException("---merge2---");
    }
}

template <typename Type>
void CMatrix<Type>::merge(const CMatrix<Type> &mat,Type num, CMatrix::EMergeType mergeType){
    *this = mat;
    switch (mergeType) {
    case CMatrix<Type>::ADD:
        for(uint i=0;i<_size; i++)
            _array[i] = _array[i] + num;
        break;
    case CMatrix<Type>::DEL:
        for(uint i=0;i<_size; i++)
            _array[i] = _array[i] - num;
        break;
    case CMatrix<Type>::MULT:
        for(uint i=0;i<_size; i++)
            _array[i] = _array[i] * num;
        break;
    default:
        throw MatException("---merge3---");
    }
}

template <typename Type>
void CMatrix<Type>::changeRow(uint rowNum1, uint rowNum2){
    if ( !( (rowNum1!=rowNum2) && insideRow(rowNum1) && insideRow(rowNum2) ) )
        return;
    Type t;
    for(uint i=1;i<=_colSize;i++){
        t = at(rowNum1,i);
        at(rowNum1,i) = at(rowNum2,i);
        at(rowNum2,i) = t;
    }
}

template <typename Type>
void CMatrix<Type>::changeCol(uint colNum1, uint colNum2){
    if ( !( (colNum1!=colNum2) && insideCol(colNum1)  && insideCol(colNum2)) )
        return;
    Type t;
    for(uint i=1;i<=_rowSize;i++){
        t = at(i,colNum1);
        at(i,colNum1) = at(i,colNum2);
        at(i,colNum2) = t;
    }
}

template <typename Type>
void CMatrix<Type>::multRow(uint rowNum, Type num){
    if (!insideRow(rowNum))
        return;
    for(uint i=1;i<=_colSize;i++)
        at(rowNum,i) = at(rowNum,i) * num;
}

template <typename Type>
void CMatrix<Type>::multCol(uint colNum, Type num){
    if(!insideCol(colNum))
        return;
    for(uint i=1;i<=_rowSize;i++)
        at(i,colNum) = at(i,colNum) * num;
}

template <typename Type>
void CMatrix<Type>::multPlusRow(uint changedRow, uint row, Type num){
    if (!(insideRow(changedRow)&&insideRow(row)))
        return;
    for(uint i=1;i<=_colSize;i++)
        at(changedRow,i) = at(changedRow,i) + num*at(row,i);
}

template <typename Type>
void CMatrix<Type>::multPlusCol(uint changedCol, uint col, Type num){
    if (!(insideCol(changedCol)&&insideCol(col)))
        return;
    for(uint i=1;i<=_rowSize;i++)
        at(i,changedCol) = at(i,changedCol) + num*at(i,col);
}

template <typename Type>
void CMatrix<Type>::getTranspose(const CMatrix<Type> &mat){
    setSize(mat.getSize(),0);
    setColSize(mat.getRowSize());
    for (uint i=1;i<=_rowSize;i++)
        for(uint j=1;j<=_colSize;j++){
            at(j,i) = mat.value(i,j);
        }
}

template <typename Type>
Type CMatrix<Type>::det(){
    if(_rowSize != _colSize)
        throw MatException("---det---");
    CMatrix<double> mat;
    toDouble(mat);
    //最后得出的积存在sum
    double sum=1;
    for(uint i=1;i<=_colSize;i++){
        uint j;
        for(j=i; j<=_colSize;j++){
            if (fabs(mat.at(i,j))>=MAT_MIN_ELE)
                break;
        }
        if (j>_colSize)
            return (Type)0;
        if (i!=j){
            sum = -sum;
            mat.changeCol(i,j);
        }

        for (j=i+1;j<=_rowSize;j++)
            mat.multPlusRow(j,i,-mat.at(j,i)/mat.at(i,i));
    }
    for (uint i=1;i<=_colSize;i++){
        sum = sum*mat.at(i,i);
    }
    return (Type)sum;
}

template <typename Type>
bool CMatrix<Type>::inverse(CMatrix<double> &invMat){
    if(_rowSize != _colSize)
        throw MatException("---inverse---");
    if (fabs(det())<=MAT_MIN_ELE)
        return false;
    invMat.unit(_rowSize);

    CMatrix<double> mat;
    toDouble(mat);

    //把右上角的变成0，并且把对角线变成1；
    for(uint i=1;i <= _colSize;i++){
        uint j;
        for(j=i; j <=_colSize;j++){
            if (fabs(mat.at(i,j))>=MAT_MIN_ELE)
                break;
        }
        if (j>_colSize)
            return (Type)0;
        mat.changeCol(i,j);
        invMat.changeCol(i,j);
        double ii=1.0/mat.at(i,i);
        mat.multCol(i,ii);
        invMat.multCol(i,ii);
        for(j=i+1;j<=_colSize;j++){
            double multNum = -mat.at(i,j)/mat.at(i,i);
            mat.multPlusCol(j,i,multNum);
            invMat.multPlusCol(j,i,multNum);
        }
    }
    //把左下角的变成0
    for(uint i=_colSize;i>=1;i--){
        for(uint j=i-1;j>=1;j--){
            double multNum = -mat.at(i,j);
            mat.multPlusCol(j,i,multNum);
            invMat.multPlusCol(j,i,multNum);
        }
    }
    return true;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator+(const CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    tmp.merge(*this,mat,CMatrix<Type>::ADD);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator+(Type num){
    CMatrix<Type> tmp;
    tmp.merge(*this,num,CMatrix<Type>::ADD);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator-(const CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    tmp.merge(*this,mat,CMatrix<Type>::DEL);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator-(Type num){
    CMatrix<Type> tmp;
    tmp.merge(*this,num,CMatrix<Type>::DEL);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator*(const CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    if (getColSize()!=mat.getRowSize())
        throw MatException("---operator*---");
    uint plusNum = getColSize();
    tmp.setSize(getRowSize()*mat.getColSize(),0);
    tmp.setRowSize(getRowSize());
    for (uint i=1;i<=_rowSize;i++)
        for(uint j=1;j<=_colSize;j++){
            for (uint k=1;k<=plusNum;k++){
                tmp.at(i,j) = tmp.at(i,j) + value(i,k)*mat.value(k,j);
                if (fabs(tmp.at(i,j))<=MAT_MIN_ELE)
                    tmp.at(i,j)=0;
            }
        }

    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator*(Type num){
    CMatrix<Type> tmp;
    tmp.merge(*this,num,CMatrix<Type>::MULT);
    return tmp;
}

template <typename Type>
Type CMatrix<Type>::out=0;

template <typename Type>
Type det(CMatrix<Type> &mat){
    return mat.det();
}

template <typename Type>
CMatrix<double> inv(CMatrix<Type> &mat){
    CMatrix<double> tmp;
    mat.inverse(tmp);
    return tmp;
}

template <typename Type>
CMatrix<Type> tra(const CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    tmp.getTranspose(mat);
    return tmp;
}

template <typename Type>
CMatrix<double> dou(CMatrix<Type> &mat){
    CMatrix<double> tmp;
    mat.toDouble(tmp);
    return tmp;
}
}
#endif // CMATRIX_H
