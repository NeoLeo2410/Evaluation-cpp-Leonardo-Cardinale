#include "QBonus2_implicite.cpp"

// A l'exécution, on obtient des temps de calcul pour Euler explicite/implicite avec D = 1, version sparse ou non.
// Cela correspond donc aux questions 2, 4, Bonus 2 explicite/implicite.

int main(){

    unsigned nt = 1001;
    double horiz = 0.5; // Horizon temporelle
    double dt = horiz/(nt-1); // Pas temporel
    
    // Question 2

    auto start_Q2 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution_Q2 = euler_explicite_Q2(dt,horiz);
    auto stop_Q2 = std::chrono::high_resolution_clock::now();
    auto duration_Q2 = std::chrono::duration_cast<std::chrono::microseconds>(stop_Q2 - start_Q2);
    std::cout << "Le calcul pour la question 2 prend " << duration_Q2.count() << " μs." << std::endl;
    exportsolution_Q2(solution_Q2);

    // Question 4

    auto start_Q4 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution_Q4 = euler_implicite_Q4(dt,horiz);
    auto stop_Q4 = std::chrono::high_resolution_clock::now();
    auto duration_Q4 = std::chrono::duration_cast<std::chrono::microseconds>(stop_Q4 - start_Q4);
    std::cout << "Le calcul pour la question 4 prend " << duration_Q4.count() << " μs." << std::endl;
    exportsolution_Q4(solution_Q4);

    // Question Bonus 1 (explicite)

    std::vector<std::vector<double> > solution_QB1e = euler_explicite_QB1e(dt,horiz);
    exportsolution_QB1e(solution_QB1e);

    // Question Bonus 1 (implicite)

    std::vector<std::vector<double> > solution_QB1i = euler_implicite_QB1i(dt,horiz);
    exportsolution_QB1i(solution_QB1i);

    // Question Bonus 2 (explicite)

    auto start_QB2e = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution_QB2e = euler_explicite_QB2e(dt,horiz);
    auto stop_QB2e = std::chrono::high_resolution_clock::now();
    auto duration_QB2e = std::chrono::duration_cast<std::chrono::microseconds>(stop_QB2e - start_QB2e);
    std::cout << "Le calcul pour la question bonus 2 (explicite) prend " << duration_QB2e.count() << " μs." << std::endl;
    exportsolution_QB2e(solution_QB2e);

    // Question Bonus 2 (implicite)

    auto start_QB2i = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double> > solution_QB2i = euler_implicite_QB2i(dt,horiz);
    auto stop_QB2i = std::chrono::high_resolution_clock::now();
    auto duration_QB2i = std::chrono::duration_cast<std::chrono::microseconds>(stop_QB2i - start_QB2i);
    std::cout << "Le calcul pour la question bonus 2 (implicite) prend " << duration_QB2i.count() << " μs." << std::endl;
    exportsolution_QB2i(solution_QB2i);

    return EXIT_SUCCESS;
}