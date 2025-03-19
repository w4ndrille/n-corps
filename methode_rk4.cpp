#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

const double G = 6.67430e-11; // Gravitational constant
const double Tmax = 4.0e+7;
const double dt = 10000;
const int N = 2; // Number of particles
const double Mass_Earth = 5.972e24; // Mass of the Earth
const double Distance_Earth_Moon = 384400000; // Distance between the Earth and the Moon
const double Distance_Earth_Sun = 149597870700; // Distance between the Earth and the Sun
const double Mass_Sun = 1.989e30; // Mass of the Sun

// Function to calculate gravitational accelerations
void calculate_accelerations(const vector<double>& positions, const vector<double>& masses, vector<double>& accelerations) {
    int num_particles = masses.size();
    for (int i = 0; i < num_particles; ++i) {
        accelerations[2 * i] = 0.0;     // Reset x-acceleration
        accelerations[2 * i + 1] = 0.0; // Reset y-acceleration
        for (int j = 0; j < num_particles; ++j) {
            if (i != j) {
                double dx = positions[2 * j] - positions[2 * i];
                double dy = positions[2 * j + 1] - positions[2 * i + 1];
                double r_squared = dx * dx + dy * dy + 1e-3; // Add softening factor
                double r = sqrt(r_squared);
                double force = G * masses[i] * masses[j] / (r_squared * r);
                accelerations[2 * i] += force * dx / masses[i];
                accelerations[2 * i + 1] += force * dy / masses[i];
            }
        }
    }
}

// Derivative function for RK4
void deriv(int n, double t, const vector<double>& y, vector<double>& dy, const vector<double>& masses) {
    vector<double> positions(y.begin(), y.begin() + 2 * N);
    vector<double> velocities(y.begin() + 2 * N, y.end());
    vector<double> accelerations(2 * N, 0.0);

    calculate_accelerations(positions, masses, accelerations);

    for (int i = 0; i < N; ++i) { // Loop over particles
        dy[2 * i] = velocities[2 * i];         // dx/dt = vx
        dy[2 * i + 1] = velocities[2 * i + 1]; // dy/dt = vy
        dy[2 * N + 2 * i] = accelerations[2 * i];     // dvx/dt = ax
        dy[2 * N + 2 * i + 1] = accelerations[2 * i + 1]; // dvy/dt = ay
    }
}

// RK4 integration method
void rk4(int n, double t, vector<double>& y, double dx, void deriv(int, double, const vector<double>&, vector<double>&, const vector<double>&), const vector<double>& masses) {
    vector<double> d1(n), d2(n), d3(n), d4(n), yp(n);
    double ddx = dx / 2.0;

    deriv(n, t, y, d1, masses);
    for (int i = 0; i < n; ++i) yp[i] = y[i] + d1[i] * ddx;
    deriv(n, t + ddx, yp, d2, masses);
    for (int i = 0; i < n; ++i) yp[i] = y[i] + d2[i] * ddx;
    deriv(n, t + ddx, yp, d3, masses);
    for (int i = 0; i < n; ++i) yp[i] = y[i] + d3[i] * dx;
    deriv(n, t + dx, yp, d4, masses);

    for (int i = 0; i < n; ++i) {
        y[i] += dx * (d1[i] + 2.0 * d2[i] + 2.0 * d3[i] + d4[i]) / 6.0;
    }
}




int main() {
    int num_equations = 4 * N; // 2 for position and 2 for velocity per particle
    vector<double> y(num_equations, 0.0); // State vector: [x1, y1, x2, y2, ..., vx1, vy1, vx2, vy2, ...]
    // Initialisation des masses pour le système Alpha Centauri
    vector<double> masses = {Mass_Sun, Mass_Sun}; // Masses en multiples de la masse terrestre

    // Initial conditions
    // y[0] = 0.0;                 y[1] = 0.0; // Particle 1 position
    // y[2] = Distance_Earth_Sun;  y[3] = 0.0; // Particle 2 position
    // y[4] = 0.0;                 y[5] = Distance_Earth_Sun; // Particle 3 position
    // y[6] = 0.0;                 y[7] = 0.0; // Particle 1 velocity
    // y[8] = 0.0;                 y[9] = 30000.0; // Particle 2 velocity
    // y[10] = 0.0;                y[11] = 0.0; // Particle 3 velocity

    y[0] = 0.0; y[1] = 0.0; // Particle 1 position
    y[2] = Distance_Earth_Sun/100; y[3] = 0.0; // Particle 2 position
    y[4] = 0.0; y[5] = 0.0; // Particle 1 velocity
    y[6] = 0.0; y[7] = 100000.0; // Particle 2 velocity

    ofstream output("n_body_simulation.txt");

    double t = 0.0;
    while (t < Tmax) {
        output << t;
        for (int i = 0; i < N; ++i) {
            output << " " << y[2 * i] << " " << y[2 * i + 1]; // Output positions
        }
        output << endl;

        rk4(num_equations, t, y, dt, deriv, masses);
        t += dt;
    }

    output.close();
    return 0;
}
