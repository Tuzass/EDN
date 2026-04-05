#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

double f(double x, double y){
    return -y * y;
}

double phi(double x){
    return 1.0 / x;
}

void writeResultPoints(std::vector<double> x, std::vector<double> y, int N){
    std::ofstream results_file("results.txt", std::ios::out);
    if (!results_file)
        throw std::runtime_error("Failed to open results.txt for writing");

    for (int j = 0; j <= N; j++){
        results_file << x[j] << ", " << y[j] << "\n";
    }
}

void writeDifferences(std::vector<double> x, std::vector<double> y, int N){
    std::ofstream differences_file("differences.txt", std::ios::out);
    if (!differences_file)
        throw std::runtime_error("Failed to open differences.txt for writing");

    for (int j = 0; j <= N; j++){
        differences_file << x[j] << ", " << std::abs(y[j] - phi(x[j])) << "\n";
    }
}

void writeSolutionPoints(double x_0, double X, int N){
    std::ofstream differences_file("solution.txt", std::ios::out);
    if (!differences_file)
        throw std::runtime_error("Failed to open solution.txt for writing");

    // how many more points to calculate to display the exact solution in the same interval (default: 5x)
    const int solution_point_density {5};
    double h = (X - x_0) / (solution_point_density * N);

    for (double x {x_0}; x < X; x += h){
        differences_file << x << ", " << phi(x) << "\n";
    }
}

std::string stripString(std::string str){
    if (str.size() == 0) return "";

    std::size_t start {};
    while (str[start] == ' ' && start < str.size()) start++;

    if (start == str.size()) return "";

    std::size_t end {str.size() - 1};
    while (str[end] == ' ' && end > 0) end--;

    return str.substr(start, end - start + 1);
}

std::map<std::string, std::string> loadSettings(){
    std::map<std::string, std::string> settings {};

    std::ifstream settings_file("tp1-settings.ini", std::ios::in);
    if (!settings_file)
        throw std::runtime_error("Failed to open settings.ini");
    
    std::string line {};
    while (std::getline(settings_file, line)){
        if (line[0] == '#') continue;

        std::size_t delim_index {line.find('=')};
        std::string key = stripString(line.substr(0, delim_index));
        std::string value = stripString(line.substr(delim_index + 1));
        settings[key] = value;
    }

    return settings;
}

void eulerMethod(double x_0, double X, double y_0, double h){
    const int N = (X - x_0) / h;
    std::vector<double> x(N + 1);
    std::vector<double> y(N + 1);
    x[0] = x_0;
    y[0] = y_0;

    for (int j = 1; j <= N; j++){
        x[j] = x[0] + j * h;
        y[j] = y[j - 1] + h * f(x[j - 1], y[j - 1]);
    }

    writeResultPoints(x, y, N);
    writeDifferences(x, y, N);
    writeSolutionPoints(x_0, X, N);
}

// Runge-Kutta method with 2 stages
void rk2Method(double x_0, double X, double y_0, double h, double b_2){
    const int N = (X - x_0) / h;
    std::vector<double> x(N + 1);
    std::vector<double> y(N + 1);
    x[0] = x_0;
    y[0] = y_0;

    for (int j = 1; j <= N; j++){
        double k_1 {f(x[j - 1], y[j - 1])};
        double k_2 {f(x[j - 1] + (0.5 * h / b_2), y[j - 1] + 0.5 * h * k_1 / b_2)};
        
        x[j] = x[0] + j * h;
        y[j] = y[j - 1] + h * (1 - b_2) * k_1 + h * b_2 * k_2;
    }

    writeResultPoints(x, y, N);
    writeDifferences(x, y, N);
    writeSolutionPoints(x_0, X, N);
}

// Runge-Kutta method with 4 stages
void rk4Method(double x_0, double X, double y_0, double h){
    const int N = (X - x_0) / h;
    std::vector<double> x(N + 1);
    std::vector<double> y(N + 1);
    x[0] = x_0;
    y[0] = y_0;

    for (int j = 1; j <= N; j++){
        double k_1 {f(x[j - 1], y[j - 1])};
        double k_2 {f(x[j - 1] + 0.5 * h, y[j - 1] + 0.5 * h * k_1)};
        double k_3 {f(x[j - 1] + 0.5 * h, y[j - 1] + 0.5 * h * k_2)};
        double k_4 {f(x[j - 1] + h, y[j - 1] + h * k_3)};
        
        x[j] = x[0] + j * h;
        y[j] = y[j - 1] + h * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
    }

    writeResultPoints(x, y, N);
    writeDifferences(x, y, N);
    writeSolutionPoints(x_0, X, N);
}

// implicit theta method; num_corrections indicates how many corrections in the predictor-corrector scheme
void thetaMethod(double x_0, double X, double y_0, double h, double theta, int num_corrections){
    const int N = (X - x_0) / h;
    std::vector<double> x(N + 1);
    std::vector<double> y(N + 1);
    x[0] = x_0;
    y[0] = y_0;

    for (int j = 1; j <= N; j++){
        double y_approx = y[j - 1] + h * f(x[j - 1], y[j - 1]); // first approximation given by Euler's method
        x[j] = x[0] + j * h;

        for (int k = 1; k <= num_corrections; k++){
            y_approx = y[j - 1] + h * ((1 - theta) * f(x[j - 1], y[j - 1]) + theta * f(x[j], y_approx));
        }

        y[j] = y_approx;
    }

    writeResultPoints(x, y, N);
    writeDifferences(x, y, N);
    writeSolutionPoints(x_0, X, N);
}

int main(){
    std::map<std::string, std::string> settings = loadSettings();

    double h {std::stod(settings["h"])};
    double x_0 {std::stod(settings["x_0"])};
    double X {std::stod(settings["X"])};
    double y_0 {std::stod(settings["y_0"])};
    int num_corrections {std::stoi(settings["num_corrections"])};

    if (settings["method"] == "euler"){
        eulerMethod(x_0, X, y_0, h);
    } else if (settings["method"] == "heun"){
        rk2Method(x_0, X, y_0, h, 0.5); // Heun's method is RK2 with b_2 = 0.5
    } else if (settings["method"] == "midpoint"){
        rk2Method(x_0, X, y_0, h, 1); // Midpoint method is RK2 with b_2 = 1
    } else if (settings["method"] == "rk4"){
        rk4Method(x_0, X, y_0, h);
    } else if (settings["method"] == "implicit_euler"){
        thetaMethod(x_0, X, y_0, h, 1, num_corrections); // Implicit Euler's method is Theta with theta = 1
    } else if (settings["method"] == "trapezoidal"){
        thetaMethod(x_0, X, y_0, h, 0.5, num_corrections); // Trapezoidal method is Theta with theta = 0.5
    }

    return 0;
}
