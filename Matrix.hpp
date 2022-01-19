#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <functional>
#include <cmath>
#include <chrono>
#include <fstream>

class Matrix{
    private: // une matrice est ici caractérisée par sa forme et des données stockées dans un vecteur de vecteurs (les lignes) de doubles
        unsigned rows; 
        unsigned columns;
        std::vector<std::vector<double> > matrix;
    public:
        Matrix(const unsigned, const unsigned, const double); // pour construire une matrice de taille donnée initialisée avec une valeur
        Matrix(Matrix&); // constructeur de copie
        Matrix(std::vector<double>, const unsigned); // convertisseur vecteur -> matrice
        Matrix operator+(Matrix); // on utilise la surcharge d'opérateurs
        Matrix operator-(Matrix);
        Matrix operator*(Matrix);
        Matrix operator*(const double);
        Matrix transpose();
        void print();
        std::vector<double> tovector(); // convertisseur matrice -> vecteur
        double& operator()(unsigned, unsigned); // pour avoir l'écriture M(i,j) des éléments d'une matrice M
        unsigned getrows(); // pour obtenir le nombre de lignes
        unsigned getcolumns(); // pour obtenir le nombre de colonnes
};
