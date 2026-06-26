#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

double f(double x){
    if (x <= 0.5) return 2 * x;
    return 2 - 2 * x;
}

double solution(double x, double t){
    double result {};

    if (std::abs(x) < 1e-14 || std::abs(x - 1.0) < 1e-14)
        return 0.0;

    for (int i = 0; i < 21; i++){
        int index_squared {(2 * i + 1) * (2 * i + 1)};
        double exp_term = exp(- index_squared * M_PI * M_PI * t);
        double sin_term = sin((2 * i + 1) * M_PI * x);

        if (i % 2 == 0) result += (1 / index_squared) * exp_term * sin_term;
        else result -= (1.0 / index_squared) * exp_term * sin_term;
    }

    return (8.0 / (M_PI * M_PI)) * result;
}

void writeResultPoints(const std::vector<double> &x, const std::vector<double> &t, const std::vector<std::vector<double>> &u, int M, int N){
    std::ofstream results_file("results.txt", std::ios::out);
    if (!results_file)
        throw std::runtime_error("Failed to open results.txt for writing");

    for (int j = 0; j <= N; j++)
        for (int m = 0; m <= M; m++)
            results_file << t[j] << ", " << x[m] << ", " << u[j][m] << "\n";
}

void writeDifferences(const std::vector<double> &x, const std::vector<double> &t, const std::vector<std::vector<double>> &u, int M, int N){
    std::ofstream differences_file("differences.txt", std::ios::out);
    if (!differences_file)
        throw std::runtime_error("Failed to open differences.txt for writing");

    for (int j {0}; j <= N; j++)
        for (int m {0}; m <= M; m++)
            differences_file << t[j] << ", " << x[m] << ", " << std::abs(u[j][m] - solution(x[m], t[j])) << "\n";
}

void writeSolutionPoints(double x_0, double X, double t_0, double T, int M, int N){
    std::ofstream differences_file("solution.txt", std::ios::out);
    if (!differences_file)
        throw std::runtime_error("Failed to open solution.txt for writing");

    // how many more points to calculate to display the exact solution in the same interval (default: 2x)
    const int solution_point_density {2};
    double dx = (X - x_0) / (solution_point_density * M);
    double dt = (T - t_0) / (solution_point_density * N);

    for (int j {}; j <= solution_point_density * N; j++){
        double t = t_0 + j * dt;
        for (int m {}; m <= solution_point_density; m++){
            double x = x_0 + m * dx;
            differences_file << t << ", " << x << ", " << solution(x, t) << "\n";
        }
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

    std::ifstream settings_file("settings.ini", std::ios::in);
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

std::vector<double> solveTridiagonal(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d){
    int n = d.size();

    std::vector<double> c_(n);
    std::vector<double> d_(n);
    std::vector<double> x(n);

    c_[0] = c[0] / b[0];
    d_[0] = d[0] / b[0];

    for (int i = 1; i < n; i++){
        double denom = b[i] - a[i] * c_[i - 1];

        c_[i] = c[i] / denom;
        d_[i] =
            (d[i] - a[i] * d_[i - 1]) / denom;
    }

    x[n - 1] = d_[n - 1];

    for (int i = n - 2; i >= 0; i--){
        x[i] =
            d_[i]
            - c_[i] * x[i + 1];
    }

    return x;
}

void ftcs(double x_0, double X, double t_0, double T, double dx, double dt){
    const int M = (X - x_0) / dx;
    const int N = (T - t_0) / dt;
    const double mu {dt / (dx * dx)};

    if (mu > 0.5)
        throw std::runtime_error("Method is unstable with current settings");

    std::vector<double> x(M + 1);
    std::vector<double> t(N + 1);
    std::vector<std::vector<double>> u(N + 1, std::vector<double>(M + 1));

    for (int m {}; m <= M; m++)
        x[m] = x_0 + m * dx;

    for (int j {}; j <= N; j++)
        t[j] = t_0 + j * dt;

    for (int m {}; m <= M; m++)
        u[0][m] = f(x_0 + m * dx);
    
    for (int j {}; j <= N; j++){
        u[j][0] = 0;
        u[j][M] = 0;
    }

    for (int j {}; j < N; j++)
        for (int m {1}; m <= M - 1; m++)
            u[j + 1][m] = u[j][m] + mu * (u[j][m - 1] - 2 * u[j][m] + u[j][m + 1]);

    writeResultPoints(x, t, u, M, N);
    writeDifferences(x, t, u, M, N);
    writeSolutionPoints(x_0, X, t_0, T, M, N);
}

void btcs(double x_0, double X, double t_0, double T, double dx, double dt){
    const int M = (X - x_0) / dx;
    const int N = (T - t_0) / dt;
    const double mu {dt / (dx * dx)};
   
    std::vector<double> x(M + 1);
    std::vector<double> t(N + 1);
    std::vector<std::vector<double>> u(N + 1, std::vector<double>(M + 1));

    for (int m {}; m <= M; m++)
        x[m] = x_0 + m * dx;

    for (int j {}; j <= N; j++)
        t[j] = t_0 + j * dt;
    
    for (int m {}; m <= M; m++)
        u[0][m] = f(x_0 + m * dx);
    
    for (int j {}; j <= N; j++){
        u[j][0] = 0;
        u[j][M] = 0;
    }

    std::vector<double> a(M - 1, -mu / 2.0);
    std::vector<double> b(M - 1, 1 + mu);
    std::vector<double> c(M - 1, -mu / 2.0);

    for (int j {}; j < N; j++){
        a[0] = 0;
        c[M - 2] = 0;
        
        std::vector<double> d(M - 1);
        for (int m {1}; m < M; m++)
            d[m - 1] = (mu / 2.0) * u[j][m - 1] + (1 - mu) * u[j][m] + (mu / 2.0) * u[j][m + 1];
        
        std::vector<double> solution = solveTridiagonal(a, b, c, d);
        for (int m {1}; m < M; m++)
            u[j + 1][m] = solution[m - 1];
    }

    writeResultPoints(x, t, u, M, N);
    writeDifferences(x, t, u, M, N);
    writeSolutionPoints(x_0, X, t_0, T, M, N);
}

void crankNicolson(double x_0, double X, double t_0, double T, double dx, double dt){
    const int M = (X - x_0) / dx;
    const int N = (T - t_0) / dt;
    const double mu {dt / (dx * dx)};
   
    std::vector<double> x(M + 1);
    std::vector<double> t(N + 1);
    std::vector<std::vector<double>> u(N + 1, std::vector<double>(M + 1));

    for (int m {}; m <= M; m++)
        x[m] = x_0 + m * dx;

    for (int j {}; j <= N; j++)
        t[j] = t_0 + j * dt;
    
    for (int m {}; m <= M; m++)
        u[0][m] = f(x_0 + m * dx);
    
    for (int j {}; j <= N; j++){
        u[j][0] = 0;
        u[j][M] = 0;
    }

    std::vector<double> a(M - 1, -mu);
    std::vector<double> b(M - 1, 1 + 2 * mu);
    std::vector<double> c(M - 1, -mu);

    for (int j {}; j < N; j++){
        a[0] = 0;
        c[M - 2] = 0;
        
        std::vector<double> d(M - 1);
        for (int m {1}; m < M; m++)
            d[m - 1] = u[j][m];
        
        std::vector<double> solution = solveTridiagonal(a, b, c, d);
        for (int m {1}; m < M; m++)
            u[j + 1][m] = solution[m - 1];
    }

    writeResultPoints(x, t, u, M, N);
    writeDifferences(x, t, u, M, N);
    writeSolutionPoints(x_0, X, t_0, T, M, N);
}

int main(){
    std::map<std::string, std::string> settings = loadSettings();

    double x_0 {std::stod(settings["x_0"])};
    double X {std::stod(settings["X"])};
    double t_0 {std::stod(settings["t_0"])};
    double T {std::stod(settings["T"])};
    double dx {std::stod(settings["dx"])};
    double dt {std::stod(settings["dt"])};
    int num_corrections {std::stoi(settings["num_corrections"])};

    if (settings["method"] == "ftcs"){
        ftcs(x_0, X, t_0, T, dx, dt);
    } else if (settings["method"] == "btcs"){
        btcs(x_0, X, t_0, T, dx, dt);
    } else if (settings["method"] == "crank-nicolson"){
        crankNicolson(x_0, X, t_0, T, dx, dt);
    }

    return 0;
}
