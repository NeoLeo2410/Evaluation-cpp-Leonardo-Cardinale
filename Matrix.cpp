#include "Matrix.hpp"
#include <stdexcept>

// Constructeur

Matrix::Matrix(unsigned n, unsigned p, double value){
    this->rows = n;
    this->columns = p;
    this->matrix.resize(n);
    for (unsigned i = 0; i < this->matrix.size(); i++){
        this->matrix[i].resize(p, value);
    }
}

// Conversion d'un vecteur 1D en matrice 1D

Matrix::Matrix(std::vector<double> vec, unsigned k){
    unsigned n = vec.size();
    this->rows = 1;
    this->columns = n;
    this->matrix.resize(1);
    this->matrix.resize(n);
    this->matrix[0] = vec;
}

// Somme de deux matrices

Matrix Matrix::operator+(Matrix B){
    if ((this->rows != B.getrows()) || (this->columns != B.getcolumns())){
        throw std::runtime_error("Incompatible matrix sizes");
    }
    else{
        Matrix C(this->rows, this->columns, 0.0);
        for (unsigned i = 0; i < this->rows; i++){
            for (unsigned j = 0; j < this->columns; j++){
                C(i,j) = this->matrix[i][j] + B(i,j);
            }
        }
        return C;
    }
}

// Différence de deux matrices

Matrix Matrix::operator-(Matrix B){
    if ((this->rows != B.getrows()) || (this->columns != B.getcolumns())){
        throw std::runtime_error("Incompatible matrix sizes");
    }
    else{
        Matrix C(this->rows, this->columns, 0.0);
        for (unsigned i = 0; i < this->rows; i++){
            for (unsigned j = 0; j < this->columns; j++){
                C(i,j) = this->matrix[i][j] - B(i,j);
            }
        }
        return C;
    }
}

// Produit de deux matrices

Matrix Matrix::operator*(Matrix B){
    if (this->columns != B.getrows()){
        throw std::runtime_error("Incompatible matrix sizes");
    }
    else{
        double sum;
        Matrix C(this->rows, B.getcolumns(), 0.0);
        for (unsigned i = 0; i < this->rows; i++){
            for (unsigned j = 0; j < B.getcolumns(); j++){
                sum = 0.0;
                for (unsigned k = 0; k < this->columns; k++){
                    sum += this->matrix[i][k] * B(k,j);
                }
                C(i,j) = sum;
            }
        }
        return C;
    }
}

// Multiplication d'une matrice par un scalaire

Matrix Matrix::operator*(double scalar){
    Matrix B(this->rows, this->columns, 0.0);
    unsigned i;
    unsigned j;
    for (i = 0; i < this->rows; i++){
        for (j = 0; j < this->columns; j++){
            B(i,j) = scalar * this->matrix[i][j];
        }
    }
    return B;
}

// M(i,j)

double& Matrix::operator()(unsigned row, unsigned column){
    return this->matrix[row][column];
}

// Nombre de lignes

unsigned Matrix::getrows(){
    return this->rows;
}

// Nombre de colonnes

unsigned Matrix::getcolumns(){
    return this->columns;
}

// Transposée

Matrix Matrix::transpose(){
    Matrix T(this->columns, this->rows, 0.0);
    for (unsigned i = 0; i < this->columns; i++){
        for (unsigned j = 0; j < this->rows; j++){
            T(i,j) = this->matrix[j][i];
        }
    }
    return T;
}

// Affichage

void Matrix::print(){
    std::cout << "Matrix: " << std::endl;
    for (unsigned i = 0; i < this->rows; i++){
        for (unsigned j = 0; j < this->columns; j++){
            std::cout << "[" << matrix[i][j] << "]";
        }
        std::cout << std::endl;
    }
}

// Constructeur de copie

Matrix::Matrix(Matrix& B){
    unsigned n = B.getrows();
    unsigned p = B.getcolumns();
    this->rows = n;
    this->columns = p;
    this->matrix.resize(n);
    for (unsigned i = 0; i < this->matrix.size(); i++){
        this->matrix[i].resize(p);
    }
    for (unsigned i = 0; i < n; i++){
        for (unsigned j = 0; j < p; j++){
            this->matrix[i][j] = B(i,j);
        }
    }
}

// Conversion Matrix -> std::vector

std::vector<double> Matrix::tovector(){
    if (this->rows != 1){
        throw std::runtime_error("Matrix needs to be a row");
    }
    else{
        std::vector<double> result = this->matrix[0];
        return result;
    }
}

// Echange de lignes pour l'algorithme de Gauss

void Matrix::swap(unsigned i, unsigned j){
    std::vector<double> vec1 = this->matrix[i];
    std::vector<double> vec2 = this->matrix[j];
    this->matrix[i] = vec2;
    this->matrix[j] = vec1;
}