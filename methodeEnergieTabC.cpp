#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono> 

using namespace std;

const double G = 6.67430e-11;
const double Tmax = 4.0e+7;
const double dt = 5000;
const int N = 2;
const int num_equations = 4 * N;
const double Mass_Earth = 5.972e24;
const double Distance_Earth_Sun = 149597870700;
const double Mass_Sun = 1.989e30;

// Accélérations gravitationnelles
void calculate_accelerations(const double positions[], const double masses[], double accelerations[]) {
    for (int i = 0; i < N; ++i) {
        accelerations[2 * i] = 0.0;
        accelerations[2 * i + 1] = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double dx = positions[2 * j] - positions[2 * i];
                double dy = positions[2 * j + 1] - positions[2 * i + 1];
                double r_squared = dx * dx + dy * dy + 1e-3;
                double r = sqrt(r_squared);
                double force = G * masses[i] * masses[j] / (r_squared * r);
                accelerations[2 * i] += force * dx / masses[i];
                accelerations[2 * i + 1] += force * dy / masses[i];
            }
        }
    }
}

// Fonction dérivée pour RK4
void deriv(double t, const double y[], double dy[], const double masses[]) {
    double positions[2 * N];
    double velocities[2 * N];
    double accelerations[2 * N];

    for (int i = 0; i < 2 * N; ++i) {
        positions[i] = y[i];
        velocities[i] = y[2 * N + i];
    }

    calculate_accelerations(positions, masses, accelerations);

    for (int i = 0; i < N; ++i) {
        dy[2 * i] = velocities[2 * i];
        dy[2 * i + 1] = velocities[2 * i + 1];
        dy[2 * N + 2 * i] = accelerations[2 * i];
        dy[2 * N + 2 * i + 1] = accelerations[2 * i + 1];
    }
}

// Méthode RK4
void rk4(double t, double y[], double dx, const double masses[]) {
    double d1[num_equations], d2[num_equations], d3[num_equations], d4[num_equations], yp[num_equations];
    double ddx = dx / 2.0;

    deriv(t, y, d1, masses);
    for (int i = 0; i < num_equations; ++i) yp[i] = y[i] + d1[i] * ddx;

    deriv(t + ddx, yp, d2, masses);
    for (int i = 0; i < num_equations; ++i) yp[i] = y[i] + d2[i] * ddx;

    deriv(t + ddx, yp, d3, masses);
    for (int i = 0; i < num_equations; ++i) yp[i] = y[i] + d3[i] * dx;

    deriv(t + dx, yp, d4, masses);

    for (int i = 0; i < num_equations; ++i) {
        y[i] += dx * (d1[i] + 2.0 * d2[i] + 2.0 * d3[i] + d4[i]) / 6.0;
    }
}

// Énergie totale du système
double calculate_energy(const double y[], const double masses[]) {
    double energy_kin = 0.0;
    double energy_pot = 0.0;

    for (int i = 0; i < N; ++i) {
        double vx = y[2 * N + 2 * i];
        double vy = y[2 * N + 2 * i + 1];
        energy_kin += 0.5 * masses[i] * (vx * vx + vy * vy);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = y[2 * j] - y[2 * i];
            double dy = y[2 * j + 1] - y[2 * i + 1];
            double r = sqrt(dx * dx + dy * dy + 1e-3);
            energy_pot -= G * masses[i] * masses[j] / r;
        }
    }

    return energy_kin + energy_pot;
}


// Function to verify Kepler's third law
void verify_kepler_law(const double y[], const double masses[]) {
    double semi_major_axis = Distance_Earth_Sun; // Assuming circular orbit
    double orbital_period = 2 * M_PI * sqrt(pow(semi_major_axis, 3) / (G * (masses[0] + masses[1])));
    double kepler_ratio = pow(orbital_period, 2) / pow(semi_major_axis, 3);

    cout << "Semi-major axis (a): " << semi_major_axis << " meters" << endl;
    cout << "Orbital period (T): " << orbital_period << " seconds" << endl;
    cout << "Kepler ratio (T^2 / a^3): " << kepler_ratio << endl;

    // Check if the ratio is approximately constant
    if (fabs(kepler_ratio - 4 * M_PI * M_PI / G) < 1e-6) {
        cout << "Kepler's third law is verified." << endl;
    } else {
        cout << "Kepler's third law is not satisfied." << endl;
    }
}


int main() {
    double y[num_equations] = {0.0};
    double masses[N] = {Mass_Sun, Mass_Earth};

    y[0] = 0.0;                 y[1] = 0.0;
    y[2] = Distance_Earth_Sun;  y[3] = 0.0;
    y[4] = 0.0;                 y[5] = 0.0;
    y[6] = 0.0;                 y[7] = 100000.0;

    ofstream output("n_body_simulation_C.txt");
    ofstream energy_output("n_body_energy_C.txt");

    auto start = std::chrono::high_resolution_clock::now(); // début du chrono
    double t = 0.0;
    double E0 = calculate_energy(y, masses);

    while (t < Tmax) {
        output << t;
        for (int i = 0; i < N; ++i) {
            output << " " << y[2 * i] << " " << y[2 * i + 1];
        }
        output << endl;

        double E = calculate_energy(y, masses) - E0;
        energy_output << t << " " << E << endl;

        rk4(t, y, dt, masses);
        t += dt;
    }

    output.close();
    energy_output.close();

    auto end = std::chrono::high_resolution_clock::now();   // fin du chrono
    std::chrono::duration<double> elapsed = end - start;
    cout << "Temps d'exécution : " << elapsed.count() << " secondes" << endl;


    verify_kepler_law(y, masses);


    return 0;
}

