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

typedef size_t site;
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

    //拷贝构造函数
    CMatrix(CMatrix<Type> & mat);

    //初始化矩阵大小以及矩阵的值
    CMatrix(site size,Type value);

    //复制一个矩阵
    void copy(CMatrix<Type> & mat);

    //输出成double的矩阵
    void toDouble(CMatrix<double> & mat);

    //从文件导入数组，默认行存储
    void load(const std::string& path,site num = 0);

    //生成随机数种子
    static void srand();

    //生成随机数矩阵
    void rand(Type min,Type max);

    //生成数量矩阵，默认单位矩阵
    void unit(site size,Type num = 1);

    //保存到文件
    void save(const std::string& path);

    //打印该矩阵
    void print(const std::string& info="matrix",site size = 5,bool withNum = false);

    //设置总共行数，如果rowSize不是size的约数，异常。
    void setRowSize(site rowSize);

    //设置总共列数，如果colSize不是size的约数，异常。
    void setColSize(site colSize);

    //设置矩阵大小
    void setSize(site size,Type value);

    //访问一维地址(从1开始)，返回引用,超出范围返回out；
    Type& at(site place);

    //访问二维地址(从(1,1)开始)，返回引用,超出范围返回out；
    Type& at(site rowNum,site colNum);

    //计算矩阵的和
    Type sum();

    //获取某一行的和，超出范围返回out；
    Type sumRow(site rowNum);

    //获取某一列的和，超出范围返回out；
    Type sumCol(site colNum);

    //取出一个子矩阵,不合理的位置返回异常
    void getSubMat(site startRow,site startCol,site endRow,site endCol,CMatrix<Type> & mat);

    //获取矩阵某一行,不合理的位置返回异常
    void getRow(site rowNum,CMatrix<Type>& mat);

    //获取矩阵某一列,不合理的位置返回异常
    void getCol(site colNum,CMatrix<Type>& mat);

    //矩阵合并
    void merge(site startRow,site startCol,CMatrix<Type> & mat,EMergeType mergeType);

    //矩阵合并
    void merge(CMatrix<Type> & mat1,CMatrix<Type> & mat2,EMergeType mergeType);

    //矩阵合并
    void merge(CMatrix<Type> & mat,Type num,EMergeType mergeType);

    //初等变换===========================================
    //交换两行
    void changeRow(site rowNum1,site rowNum2);
    //交换两列
    void changeCol(site colNum1,site colNum2);

    //某一行乘以一个数
    void multRow(site rowNum,Type num);
    //某一列乘以一个数
    void multCol(site colNum,Type num);

    //某一行乘以一个数加到另一行上
    void multPlusRow(site changedRow, site row, Type num);
    //某一列乘以一个数加到另一列上
    void multPlusCol(site changedCol, site col, Type num);
    //初等变换===========================================

    //矩阵转置
    void  transpose(CMatrix<Type> & mat);

    //矩阵行列式
    Type det();

    //矩阵的逆
    bool inverse(CMatrix<double>& mat);

    //获取矩阵大小
    inline site getSize() const;

    //获取矩阵行数
    inline site getRowSize() const;

    //获取矩阵列数
    inline site getColSize() const;

    //获取该矩阵的vector数组
    inline std::vector <Type> &getMat();

    //判断输入是否在矩阵的行范围
    inline bool insideRow(site rowNum);

    //判断输入是否在矩阵的列范围
    inline bool insideCol(site colNum);

    inline Type &operator ()(site rowNum,site colNum);

    inline Type &operator ()(site place);

    inline CMatrix<Type>  operator+(CMatrix<Type>& mat);

    inline CMatrix<Type>  operator+(Type num);

    inline CMatrix<Type>  operator-(CMatrix<Type>& mat);

    inline CMatrix<Type>  operator-(Type num);

    inline CMatrix<Type>  operator*(CMatrix<Type>& mat);

    inline CMatrix<Type>  operator*(Type num);

    //转为double类型矩阵
    //CMatrix<double> dou(CMatrix<Type> &mat);

    //计算行列式
    //Type det(CMatrix<Type> &mat)

    //矩阵求逆
    //CMatrix<double> inv(CMatrix<Type> &mat)

    //矩阵转置
    //CMatrix<Type> tra(CMatrix<Type> &mat);
private:
    std::vector <Type> _mat;
    site _rowSize;
    site _colSize;
    site _size;
};


template <typename Type>
CMatrix<Type>::CMatrix(CMatrix<Type> &mat){
    copy(mat);
}

template <typename Type>
CMatrix<Type>::CMatrix(site size, Type value){
    _mat.resize(size,value);
    setRowSize(1);
}

template <typename Type>
void CMatrix<Type>::copy(CMatrix<Type> &mat){
    _size = mat.getSize();
    _rowSize = mat.getRowSize();
    _colSize = mat.getColSize();
    _mat.resize(_size);
    std::vector<Type> & vecMat = mat.getMat();
    for (site i=0;i<_size;i++)
        _mat[i] = vecMat[i];
}

template <typename Type>
void CMatrix<Type>::toDouble(CMatrix<double> &mat){
    mat.setSize(_size,0);
    mat.setColSize(_colSize);
    std::vector<double> & vecMat = mat.getMat();
    for(site i = 0;i<_size;i++)
        vecMat[i] =(double) _mat[i];
}

template <typename Type>
void CMatrix<Type>::load(const std::string &path, site num)
{
    _mat.clear();
    std::ifstream in(path.c_str());
    if(!in.is_open())
        throw  MatException("---load---");

    //读取所有的
    Type t;
    if(num == 0){
        while (in>>t){
            _mat.push_back(t);
        }
        in.close();
        setRowSize(1);
        return;
    }

    //读取指定个数，个数不能超过文件的
    _mat.resize(num);
    for(site i=0;i<num;i++){
        in>>t;
        _mat[i] = t;
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
        for(site i = 0; i < _size;i++ ){
            _mat[i] = min+(max-min)*::rand() / double(RAND_MAX);
        }
    }
    else{
        for(site i = 0; i < _size;i++ ){
            _mat[i] = (Type)(::rand() %(uint) (max-min+1))+ min;
        }
    }
}

template <typename Type>
void CMatrix<Type>::unit(site size,Type num){
    setSize(size*size,0);
    setRowSize(size);
    for (site i=1;i<=_rowSize;i++)
        at(i,i) = num;
}

template <typename Type>
void CMatrix<Type>::save(const std::string &path){
    std::ofstream out(path.c_str());
    for (site i=0;i<_size;i++){
        out<<_mat[i]<<" ";
        if ((i+1) % _colSize == 0)
            out<<std::endl;
    }
    out.close();
}

template <typename Type>
void CMatrix<Type>::print(const std::string& info,site size,bool withNum){
    site i;
    std::cout<<info<<":"<<_rowSize<<"*"<<_colSize<<std::endl;
    if (withNum){
        for (i=0;i<=_colSize;i++)
            std::cout<<std::left<<std::setw(size)<<i<<" ";
        std::cout<<std::endl;
        for (i=0;i<_size;i++){
            if ( i % _colSize == 0)
                std::cout<<std::left<<std::setw(size)<<(i / _colSize +1)<<" ";
            std::cout<<std::left<<std::setw(size)<<_mat[i]<<" ";
            if ((i+1) % _colSize == 0)
                std::cout<<std::endl;
        }
    }
    else {
        for (i=0;i<_size;i++){
            std::cout<<std::left<<std::setw(size)<<_mat[i]<<" ";
            if ((i+1) % _colSize == 0)
                std::cout<<std::endl;
        }
    }
}

template <typename Type>
void CMatrix<Type>::setRowSize(site rowSize){
    _size = _mat.size();
    _rowSize = rowSize;
    _colSize = _size / _rowSize;
    if (_size % _rowSize != 0)
        throw MatException("---setRowSize---");
}

template <typename Type>
void CMatrix<Type>::setColSize(site colSize){
    _size = _mat.size();
    _colSize = colSize;
    _rowSize = _size / _colSize;
    if (_size % _colSize != 0)
        throw MatException("---setColSize---");
}

template <typename Type>
void CMatrix<Type>::setSize(site size, Type value){
    _mat.resize(size,value);
    setRowSize(1);
}

template <typename Type>
Type &CMatrix<Type>::at(site place){
    if(place<1 || place >_size)
        return out;
    return _mat[place-1];
}

template <typename Type>
Type &CMatrix<Type>::at(site rowNum, site colNum){
    if (insideCol(colNum) && insideRow(rowNum))
        return _mat[_colSize*(rowNum-1)+colNum-1 ];
    else
        return out;
}

template <typename Type>
Type CMatrix<Type>::sum(){
    Type sum = 0;
    for(site i=0;i<_size;i++){
        sum += _mat[i];
    }
    return sum;
}

template <typename Type>
Type CMatrix<Type>::sumRow(site rowNum){
    if (!insideRow(rowNum))
        return out;
    site start = _colSize*(rowNum-1);
    Type sum = 0;
    for(site i=0;i<_colSize;i++)
        sum += _mat[i+start ];
    return sum;
}

template <typename Type>
Type CMatrix<Type>::sumCol(site colNum){
    if (!insideCol(colNum))
        return out;
    Type sum = 0;
    for(site i = 0;i<_rowSize;i++)
        sum += _mat[i*_colSize+colNum-1];
    return sum;
}

template <typename Type>
void CMatrix<Type>::getSubMat(site startRow, site startCol,
                              site endRow, site endCol,CMatrix<Type> &mat){
    if (! (1<=startCol && startCol<=endCol && endCol<=_colSize
           && 1<=startRow && startRow<=endRow && endCol<=_rowSize))
        throw MatException("---getSubMat---");
    mat.setSize((endCol-startCol+1)*(endRow-startRow+1),0);
    site i,j,k=0;
    for (i=startRow;i<=endRow;i++)
        for(j=startCol;j<=endCol;j++)
            mat.at(++k) = _mat[(i-1)*_colSize+j-1];
    mat.setColSize(endCol-startCol+1);
}

template <typename Type>
void CMatrix<Type>::getRow(site rowNum,CMatrix<Type> &mat){
    if (!insideRow(rowNum))
        throw MatException("---getRow---");
    mat.setSize(_colSize,0);
    site start = _colSize*(rowNum-1);
    for(site i=0;i<_colSize;i++)
        mat.at(i+1) = _mat[i+start ];
    mat.setRowSize(1);
}

template <typename Type>
void CMatrix<Type>::getCol(site colNum,CMatrix<Type> &mat){
    if (!insideCol(colNum))
        throw MatException("---getCol---");
    mat.setSize(_rowSize,0);
    for(site i = 0;i<_rowSize;i++)
        mat.at(i+1) = _mat[i*_colSize+colNum-1];
    mat.setColSize(1);
}

template <typename Type>
void CMatrix<Type>::merge(site startRow, site startCol, CMatrix<Type> &mat,
                          CMatrix<Type>::EMergeType mergeType){
    if (!( insideCol(startCol) && insideRow(startRow) ))
        throw MatException("---merge---");
    site endRow =startRow + mat.getRowSize() -1;
    site endCol =startCol + mat.getColSize() -1;
    if (endRow>_rowSize || endCol>_colSize)
        throw MatException("---merge---");
    site i,j;
    switch (mergeType) {
    case CMatrix<Type>::REPLACE:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = mat.at(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::ADD:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = this->at(i,j) + mat.at(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::DEL:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = this->at(i,j) - mat.at(i-startRow+1,j-startCol+1);
        break;
    case CMatrix<Type>::MIN:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = GET_MIN( this->at(i,j) , mat.at(i-startRow+1,j-startCol+1));
        break;
    case CMatrix<Type>::MAX:
        for(i=startRow;i<=endRow;i++)
            for(j=startCol;j<=endCol;j++)
                this->at(i,j) = GET_MAX( this->at(i,j) , mat.at(i-startRow+1,j-startCol+1));
        break;
    default:
        throw MatException("---merge---");
        break;
    }
}

template <typename Type>
void CMatrix<Type>::merge(CMatrix<Type> &mat1, CMatrix<Type> &mat2, CMatrix::EMergeType mergeType)
{
    if(CMatrix<Type>::MULT == mergeType) {
        if (mat1.getColSize()!=mat2.getRowSize())
            throw MatException("---merge2---");
        site plusNum = mat1.getColSize();
        setSize(mat1.getRowSize()*mat2.getColSize(),0);
        setRowSize(mat1.getRowSize());
        for (site i=1;i<=_rowSize;i++)
            for(site j=1;j<=_colSize;j++){
                for (site k=1;k<=plusNum;k++){
                    at(i,j) = at(i,j) + mat1.at(i,k)*mat2.at(k,j);
                    if (fabs(at(i,j))<=MAT_MIN_ELE)
                        at(i,j)=0;
                }
            }
        return;
    }

    if ( (mat1.getColSize() != mat2.getColSize()) || (mat1.getRowSize() != mat2.getRowSize()) )
        throw MatException("---merge2---");
    std::vector <Type>& _vec1 = mat1.getMat();
    std::vector <Type>& _vec2 = mat2.getMat();
    site i,size = _vec1.size();
    _mat.resize(size,0);
    setColSize(mat1.getColSize());
    switch (mergeType) {
    case CMatrix<Type>::ADD:
        for (i=0;i<size;i++)
            _mat[i] = _vec1[i] + _vec2[i];
        break;
    case CMatrix<Type>::DEL:
        for (i=0;i<size;i++)
            _mat[i] = _vec1[i] - _vec2[i];
        break;
    case CMatrix<Type>::MIN:
        for (i=0;i<size;i++)
            _mat[i] = GET_MIN(_vec1[i] , _vec2[i]);
        break;
    case CMatrix<Type>::MAX:
        for (i=0;i<size;i++)
            _mat[i] = GET_MAX(_vec1[i] , _vec2[i]);
        break;
    default:
        throw MatException("---merge2---");
    }
}

template <typename Type>
void CMatrix<Type>::merge(CMatrix<Type> &mat,Type num, CMatrix::EMergeType mergeType){
    copy(mat);
    switch (mergeType) {
    case CMatrix<Type>::ADD:
        for(site i=0;i<_size; i++)
            _mat[i] = _mat[i] + num;
        break;
    case CMatrix<Type>::MULT:
        for(site i=0;i<_size; i++)
            _mat[i] = _mat[i] * num;
        break;
    default:
        throw MatException("---merge3---");
    }
}

template <typename Type>
void CMatrix<Type>::changeRow(site rowNum1, site rowNum2){
    if ( !( (rowNum1!=rowNum2) && insideRow(rowNum1) && insideRow(rowNum2) ) )
        return;
    Type t;
    for(site i=1;i<=_colSize;i++){
        t = at(rowNum1,i);
        at(rowNum1,i) = at(rowNum2,i);
        at(rowNum2,i) = t;
    }
}

template <typename Type>
void CMatrix<Type>::changeCol(site colNum1, site colNum2){
    if ( !( (colNum1!=colNum2) && insideCol(colNum1)  && insideCol(colNum2)) )
        return;
    Type t;
    for(site i=1;i<=_rowSize;i++){
        t = at(i,colNum1);
        at(i,colNum1) = at(i,colNum2);
        at(i,colNum2) = t;
    }
}

template <typename Type>
void CMatrix<Type>::multRow(site rowNum, Type num){
    if (!insideRow(rowNum))
        return;
    for(site i=1;i<=_colSize;i++)
        at(rowNum,i) = at(rowNum,i) * num;
}

template <typename Type>
void CMatrix<Type>::multCol(site colNum, Type num){
    if(!insideCol(colNum))
        return;
    for(site i=1;i<=_rowSize;i++)
        at(i,colNum) = at(i,colNum) * num;
}

template <typename Type>
void CMatrix<Type>::multPlusRow(site changedRow, site row, Type num){
    if (!(insideRow(changedRow)&&insideRow(row)))
        return;
    for(site i=1;i<=_colSize;i++)
        at(changedRow,i) = at(changedRow,i) + num*at(row,i);
}

template <typename Type>
void CMatrix<Type>::multPlusCol(site changedCol, site col, Type num){
    if (!(insideCol(changedCol)&&insideCol(col)))
        return;
    for(site i=1;i<=_rowSize;i++)
        at(i,changedCol) = at(i,changedCol) + num*at(i,col);
}

template <typename Type>
void CMatrix<Type>::transpose(CMatrix<Type> &mat){
    mat.setSize(_size,0);
    mat.setColSize(_rowSize);
    for (site i=1;i<=_rowSize;i++)
        for(site j=1;j<=_colSize;j++){
            mat.at(j,i) = this->at(i,j);
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
    for(site i=1;i<=_colSize;i++){
        site j;
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
    for (site i=1;i<=_colSize;i++){
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
    for(site i=1;i<=_colSize;i++){
        site j;
        for(j=i; j<=_colSize;j++){
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
    for(site i=_colSize;i>=1;i--){
        for(site j=i-1;j>=1;j--){
            double multNum = -mat.at(i,j);
            mat.multPlusCol(j,i,multNum);
            invMat.multPlusCol(j,i,multNum);
        }
    }
    return true;
}


template <typename Type>
site CMatrix<Type>::getSize() const{
    return _size;
}

template <typename Type>
site CMatrix<Type>::getRowSize() const{
    return _rowSize;
}

template <typename Type>
site CMatrix<Type>::getColSize() const {
    return _colSize;
}

template <typename Type>
std::vector<Type> &CMatrix<Type>::getMat(){
    return _mat;
}

template <typename Type>
bool CMatrix<Type>::insideRow(site rowNum){
    return (1<=rowNum && rowNum<=_rowSize);
}

template <typename Type>
bool CMatrix<Type>::insideCol(site colNum){
    return (1<=colNum && colNum<=_colSize);
}

template <typename Type>
Type &CMatrix<Type>::operator ()(site rowNum, site colNum){
    return at(rowNum,colNum);
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator+(CMatrix<Type> &mat){
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
CMatrix<Type> CMatrix<Type>::operator-(CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    tmp.merge(*this,mat,CMatrix<Type>::DEL);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator-(Type num){
    CMatrix<Type> tmp;
    tmp.merge(*this,-num,CMatrix<Type>::ADD);
    return tmp;
}

template <typename Type>
CMatrix<Type> CMatrix<Type>::operator*(CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    tmp.merge(*this,mat,CMatrix<Type>::MULT);
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
CMatrix<Type> tra(CMatrix<Type> &mat){
    CMatrix<Type> tmp;
    mat.transpose(tmp);
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