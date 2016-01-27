#include "CMatrixD.h"

CMatrixD::CMatrixD(uint rowSize, uint colSize, double num)
{
    reset(rowSize, colSize, num);
}

CMatrixD::CMatrixD(CMatrixD &&mat)
{
    __move(mat);
}

void CMatrixD::load(const std::string &path)
{
    reset(_rowSize,_colSize,0);
    std::ifstream in(path.c_str());
    if(!in.is_open())
        throw  MatException("---load---");

    double t;
    uint num = 0;
    while (in>>t)
    {
        _arrEle[num] = t;
        num++;
        if (_rowSize * _colSize == num)
            break;
    }
    if (_rowSize * _colSize != num)
    {
        in.close();
        throw  MatException("---load---");
    }
    in.close();
}

void CMatrixD::save(const std::string &path)
{
    std::ofstream out(path.c_str());
    for (uint i = 1; i <= _rowSize; i++)
    {
        for (uint j = 1; j <= _colSize; j++)
            out<<at(i,j)<<" ";
        out<<std::endl;
    }
    out.close();
}

void CMatrixD::reset(uint rowSize, uint colSize, double num)
{
    _rowSize = rowSize;
    _colSize = colSize;
    _size = _rowSize * _colSize;
    _arrEle.resize(_size, num);
    _arrRowNum.resize(rowSize, 0);
    _arrColNum.resize(colSize, 0);
    for (uint i=0; i<rowSize; _arrRowNum[i] = i, i++);
    for (uint j=0; j<colSize; _arrColNum[j] = j, j++);
}

void CMatrixD::resize(uint rowSize, uint colSize)
{
    if (rowSize * colSize != _size)
        throw MatException("---resize---");
    CMatrixD mat(rowSize, colSize);
    uint k = 0;
    for (uint i = 1; i <= _rowSize; i++)
        for (uint j = 1; j<= _colSize; j++)
            mat._arrEle[k++] = at(i,j);
    __move(mat);
}

void CMatrixD::srand()
{
    ::srand((unsigned)time(NULL));
}

void CMatrixD::rand(int min, int max)
{
    if (min>max)
        throw MatException("---rand---");
    for(uint i = 0; i < _size;i++ )
    {
        _arrEle[i] = (::rand() %(uint) (max-min+1)+ min);
    }
}

void CMatrixD::unit(uint size, double num)
{
    this->reset(size, size, 0);
    for (uint i = 1; i <= size; at(i,i) =num, i++);
}

void CMatrixD::print(const std::string &info, uint eachSize)
{
    std::cout<<info<<":"<<_rowSize<<"*"<<_colSize<<std::endl;
    for (uint i = 1; i <= _rowSize; i++)
    {
        for (uint j = 1; j <= _colSize; j++)
            std::cout << std::left << std::setw(eachSize) << (at(i,j) > MAT_MIN_ELE ? at(i,j) : 0) <<" ";
        std::cout << std::endl;
    }
}

CMatrixD CMatrixD::tra()
{
    CMatrixD mat(_colSize, _rowSize);
    for (uint i = 1; i <= _rowSize; i++)
        for (uint j = 1; j<= _colSize; j++)
            mat.at(j,i) = at(i,j);
    return mat;
}

CMatrixD CMatrixD::inv()
{
    if(_rowSize != _colSize)
        throw MatException("---inverse---");
    if (fabs(det())<=MAT_MIN_ELE)
        throw MatException("---inverse---");
    CMatrixD invMat(_rowSize, _colSize);
    invMat.unit(_rowSize);

    CMatrixD mat(*this);

    //把右上角的变成0，并且把对角线变成1；
    for(uint i=1;i <= _colSize;i++){
        uint j;
        for(j=i; j <=_colSize;j++){
            if (fabs(mat.at(i,j))>=MAT_MIN_ELE)
                break;
        }
        if (j>_colSize)
            throw MatException("---inverse---");
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
    return invMat;
}

double CMatrixD::det()
{
    if(_rowSize != _colSize)
        throw MatException("---det---");
    CMatrixD mat(*this);
    //最后得出的积存在sum
    double sum=1;
    for(uint i=1; i <= _colSize; i++)
    {
        uint j;
        for(j=i; j <= _colSize;j++)
        {
            if (fabs(mat.at(i, j))>=MAT_MIN_ELE)
                break;
        }
        if (j>_colSize)
            return 0;
        if (i!=j)
        {
            sum = -sum;
            mat.changeCol(i,j);
        }

        for (j=i+1;j<=_rowSize;j++)
            mat.multPlusRow(j,i,-mat.at(j,i)/mat.at(i,i));
    }
    for (uint i=1;i<=_colSize;i++)
    {
        sum = sum*mat.at(i,i);
    }
    return sum;
}

bool CMatrixD::changeRow(uint rowNum1, uint rowNum2)
{
    if ( !( (rowNum1!=rowNum2) && __insideRow(rowNum1) && __insideRow(rowNum2) ) )
        return false;
    uint tmp = _arrRowNum[rowNum1-1];
    _arrRowNum[rowNum1 - 1] = _arrRowNum[rowNum2 - 1];
    _arrRowNum[rowNum2 - 1] = tmp;
    return true;
}


bool CMatrixD::changeCol(uint colNum1, uint colNum2)
{
    if ( !( (colNum1!=colNum2) && __insideCol(colNum1) && __insideCol(colNum2) ) )
        return false;
    uint tmp = _arrColNum[colNum1-1];
    _arrColNum[colNum1 - 1] = _arrColNum[colNum2 - 1];
    _arrColNum[colNum2 - 1] = tmp;
    return true;
}

bool CMatrixD::multRow(uint rowNum, double num)
{
    if (!(__insideRow(rowNum) && (fabs(num) >= MAT_MIN_ELE)))
        return false;
    for (uint i = 1; i <= _colSize; at(rowNum, i) *= num, i++);
    return true;
}

bool CMatrixD::multCol(uint colNum, double num)
{
    if (!(__insideCol(colNum) && (fabs(num) >= MAT_MIN_ELE)))
        return false;
    for (uint i = 1; i <= _rowSize; at(i, colNum) *= num, i++);
    return true;
}

bool CMatrixD::multPlusRow(uint changedRow, uint row, double num)
{
    if (!( __insideRow(changedRow) && __insideRow(row) && (changedRow != row) && (fabs(num) > MAT_MIN_ELE)))
        return false;
    for(uint i=1;i<=_colSize; at(changedRow,i) = at(changedRow,i) + num*at(row,i), i++);
    return true;
}

bool CMatrixD::multPlusCol(uint changedCol, uint col, double num)
{
    if (!(__insideCol(changedCol)&&__insideCol(col) && (changedCol != col) && (fabs(num) > MAT_MIN_ELE)))
        return false;
    for (uint i=1; i<=_rowSize; at(i,changedCol) = at(i,changedCol) + num*at(i,col), i++);
    return true;
}

bool CMatrixD::__insideRow(uint rowNum) const
{
    return (1<=rowNum && rowNum<=_rowSize);
}

bool CMatrixD::__insideCol(uint colNum) const
{
    return (1<=colNum && colNum<=_colSize);
}


double & CMatrixD::at(uint rowNum,uint colNum)
{
    if (__insideCol(colNum) && __insideRow(rowNum))
        return _arrEle[_colSize*(_arrRowNum[rowNum - 1])+_arrColNum[colNum - 1]];
    else
        throw MatException("---at---");
}

double  CMatrixD::at(uint rowNum,uint colNum) const
{
    if (__insideCol(colNum) && __insideRow(rowNum))
        return _arrEle[_colSize*(_arrRowNum[rowNum - 1])+_arrColNum[colNum - 1]];
    else
        throw MatException("---at---");
}

uint CMatrixD::rowSize() const
{
    return _rowSize;
}

uint CMatrixD::colSize() const
{
    return _colSize;
}

void CMatrixD::__move(CMatrixD &mat)
{
    _arrEle = std::move(mat._arrEle);
    _arrRowNum = std::move(mat._arrRowNum);
    _arrColNum = std::move(mat._arrColNum);
    _rowSize = mat.rowSize();
    _colSize = mat.colSize();
    _size = mat._size;
    mat._rowSize = 0;
    mat._colSize = 0;
    mat._size = 0;
}

CMatrixD &CMatrixD::operator=(CMatrixD &&mat)
{
    __move(mat);
    return *this;
}

CMatrixD CMatrixD::operator+(const CMatrixD &mat)
{
    if (!((mat.rowSize() == _rowSize) && (mat.colSize() == _colSize)))
        throw MatException("+ err");
    CMatrixD matAns(_rowSize, _colSize);
    for (uint i = 1; i <= _rowSize; i++)
        for (uint j=1; j <= _colSize; j++)
            matAns.at(i,j) = this->at(i,j) + mat.at(i,j);
    return matAns;
}

CMatrixD CMatrixD::operator-(const CMatrixD &mat)
{
    if (!((mat.rowSize() == _rowSize) && (mat.colSize() == _colSize)))
        throw MatException("- err");
    CMatrixD matAns(_rowSize, _colSize);
    for (uint i = 1; i <= _rowSize; i++)
        for (uint j=1; j <= _colSize; j++)
            matAns.at(i,j) = this->at(i,j) - mat.at(i,j);
    return matAns;
}

CMatrixD CMatrixD::operator*(const CMatrixD &mat)
{
    if (!((mat.rowSize() == _colSize) && (mat.colSize() == _rowSize)))
        throw MatException("* err");
    uint rowSize = _rowSize;
    uint colSize = mat.colSize();
    uint plusTimes = _colSize;
    CMatrixD matAns(rowSize, colSize);
    for (uint i = 1; i <= _rowSize; i++)
        for (uint j=1; j <= _colSize; j++)
            for (uint k=1; k <= plusTimes; k++)
                matAns.at(i,j) += this->at(i,k) * mat.at(k,j);
    return matAns;
}

CMatrixD CMatrixD::operator*(double num)
{
    CMatrixD matAns(*this);
    for (auto& it : matAns._arrEle)
    {
        it = it * num;
    }
    return matAns;
}
